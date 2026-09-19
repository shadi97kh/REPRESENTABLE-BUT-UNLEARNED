"""Shared training loop, so the monotone model and its unconstrained baseline are
trained by identical code and any accuracy difference is attributable to the sign
constraint alone rather than to two separate training scripts."""
import json
import numpy as np


def fit(X, A, y, idx_tr, idx_va, d_hid=32, n_layers=2, epochs=250, lr=0.01,
        batch=128, seed=0, nonneg=True, verbose=False):
    """X (N,L,F) non-negative, A (N,L,L) with self-loops, y (N,). Returns (model, best_val)."""
    import torch
    from .models import TorchMonoMPNN
    torch.manual_seed(seed)
    Xt = torch.tensor(np.asarray(X), dtype=torch.float32)
    At = torch.tensor(np.asarray(A), dtype=torch.float32)
    Y = torch.tensor(np.asarray(y), dtype=torch.float32)
    tr, va = torch.tensor(np.asarray(idx_tr)), torch.tensor(np.asarray(idx_va))
    model = TorchMonoMPNN(Xt.shape[2], d_hid, n_layers, agg="sum", seed=seed, nonneg=nonneg)

    with torch.no_grad():                       # put the readout on the target's scale
        raw = model(Xt[tr], At[tr])
        s = float(Y[tr].std() / raw.std().clamp(min=1e-12))
        model.set_output_affine(s, float(Y[tr].mean() - s * raw.mean()))

    opt = torch.optim.Adam(model.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=20)
    g = torch.Generator().manual_seed(seed)
    best, best_state, best_ep = float("inf"), None, -1
    for ep in range(epochs):
        model.train()
        perm = tr[torch.randperm(len(tr), generator=g)]
        for i in range(0, len(perm), batch):
            b = perm[i:i + batch]
            opt.zero_grad()
            torch.nn.functional.mse_loss(model(Xt[b], At[b]), Y[b]).backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            vl = float(torch.nn.functional.mse_loss(model(Xt[va], At[va]), Y[va]))
        sched.step(vl)
        if vl < best - 1e-7:
            best, best_ep = vl, ep
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        if verbose and ep % 50 == 0:
            print(f"      epoch {ep:4d}  val {vl:.4f}")
    model.load_state_dict(best_state)
    return model, best, best_ep


def predict(model, X, A, idx):
    import torch
    with torch.no_grad():
        return model(torch.tensor(np.asarray(X)[idx], dtype=torch.float32),
                     torch.tensor(np.asarray(A)[idx], dtype=torch.float32)).numpy()


def save_numpy_model(m, path):
    json.dump({"W": [w.tolist() for w in m.W], "B": [b.tolist() for b in m.B],
               "out": m.out.tolist(), "scale": m.scale, "shift": m.shift,
               "agg": m.agg, "nonneg": m.nonneg}, open(path, "w"))


def load_numpy_model(path):
    from .models import MonoMPNN
    d = json.load(open(path))
    m = MonoMPNN(1, 1, 1, agg=d["agg"], nonneg=d["nonneg"], seed=0)
    m.W = [np.array(w) for w in d["W"]]
    m.B = [np.array(b) for b in d["B"]]
    m.out = np.array(d["out"]); m.scale = d["scale"]; m.shift = d["shift"]
    return m

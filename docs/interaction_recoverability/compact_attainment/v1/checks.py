"""Fixed analytic checks only; run under resource_gate.py. Never samples."""
import json
from pathlib import Path
from fractions import Fraction as Q
import mpmath as mp
import estimator as e


def main():
    out = Path(__file__).with_name('checks.json')
    if out.exists():
        raise FileExistsError('Retain all checks: choose a fresh version for another run')
    mp.mp.dps = 48
    tol = mp.mpf('1e-43')
    result = {'status': 'started', 'scope': 'Fixed analytic functions and hand-written polygon; no observations or sampling.',
              'backend': 'mpmath diagnostics, not certified interval arithmetic'}
    try:
        rational = {'endpoint_lower_upper_bound_below_5': 4+Q(153,160) < 5,
                    'endpoint_upper_lower_bound_above_5': 4+Q(184,160) > 5,
                    'hp_upper_squared': 570**2 < 71**2*65,
                    'hp_lower_squared': 89**2 > 81*97,
                    'condition_factor_below_69': Q(6675,98) < 69,
                    'projected_contraction_below_3_4': Q(71,98) < Q(3,4)}
        assert all(rational.values())
        result['rational_certificates'] = rational
        # r(0)=0; e^x sin(x)/1000 obeys the weighted C3 local envelope.
        g1 = lambda x: mp.exp(x)*(1+mp.sin(x)/1000)
        g2 = lambda x: mp.exp(x)*(1+mp.sin(2*x)/2000)
        r1 = lambda x: mp.exp(x)*mp.sin(x)/1000
        r2 = lambda x: mp.exp(x)*mp.sin(2*x)/2000
        aa = (mp.mpf(2)/3, mp.mpf(3)/4)
        bb = (mp.mpf(3)/5, mp.mpf(4)/5)
        shifts = ((0,0),(aa[0],bb[0]),(aa[1],bb[1]),(1,0),(0,1))
        curves = [lambda z,a=a,b=b: g1(z+a)+g2(e.h(z)+b) for a,b in shifts]
        tangents = [lambda z,a=a,b=b: r1(z+a)+r2(e.h(z)+b) for a,b in shifts]
        op = e.CompactOperators(curves, 24, tol)
        tangent = e.CompactOperators(tangents, 24, tol)
        x, y = -mp.mpf(3)/2, -mp.mpf(1)
        pop_error = abs(op.E(x)-(g1(op.C(x))-g1(x)))
        tan_error = abs(tangent.E(x)-(r1(op.C(x))-r1(x)))
        xx, yy = x, y
        for _ in range(op.depth):
            xx, yy = op.C(xx), op.C(yy)
        finite_error = abs(op.difference(x,y)-(g1(x)-g1(y)-g1(xx)+g1(yy)))
        assert max(pop_error, tan_error, finite_error) < mp.mpf('1e-36')
        result['identity'] = {'population_error':str(pop_error),'tangent_error':str(tan_error),
                              'finite_depth':op.depth,'finite_remainder_identity_error':str(finite_error),
                              'normalization_error':str(abs(1-(curves[3](-1)-curves[0](-1))-g1(-1)))}
        # Fixed population calibration fixture, not an empirical/model fit.
        coeff_records = []
        for pair in (1,2):
            M, c0 = e.public_matrix(pair)
            true = (c0[0]+mp.mpf('.004'),c0[1]-mp.mpf('.003'))
            B = lambda c: mp.matrix([g1(z+c[0])+g2(e.h(z)+c[1]) for z in (-1,1)])
            recovered = e.projected_iteration(B, B(true), pair, 80)
            error = max(abs(recovered[j]-true[j]) for j in range(2))
            condition = mp.sqrt(mp.fsum(v*v for v in M))*mp.sqrt(mp.fsum(v*v for v in M**-1))
            q = mp.exp(mp.mpf('.02'))-1+mp.mpf('.01')*mp.exp(mp.mpf('.02'))*condition
            assert q < mp.mpf(3)/4 and error < mp.mpf('1e-35')
            coeff_records.append({'pair':pair,'iterations':80,'q_bound':str(q),'error':str(error)})
        result['population_iteration'] = coeff_records
        # Hand-written quantile polygon, no random sample. Known kink locations.
        poly = e.polygon([0,1,4])
        cuts = [mp.mpf(0),mp.mpf(1)/3,mp.mpf(2)/3,mp.mpf(1)]
        value = bound = mp.mpf(0)
        panels = 0
        for left,right in zip(cuts,cuts[1:]):
            v,b,k = e.gauss_two_piece(lambda p: poly(p)*p**3, left,right,mp.mpf('0.04'),216)
            value += v; bound += b; panels += k
        # On [1/3,2/3], Q=3p-1; on [2/3,1], Q=9p-5.
        ant = lambda p,a,b: a*p**5/5+b*p**4/4
        exact = ant(cuts[2],3,-1)-ant(cuts[1],3,-1)+ant(cuts[3],9,-5)-ant(cuts[2],9,-5)
        assert abs(value-exact) <= bound
        lo,hi = e.monotone_knot_bracket(lambda t:t*t,mp.mpf('.5'),mp.mpf(0),mp.mpf(1),mp.mpf('1e-30'))
        assert lo <= mp.sqrt(mp.mpf('.5')) <= hi
        result['piecewise_quadrature'] = {'absolute_error':str(abs(value-exact)),
            'proved_truncation_bound':str(bound),'panels':panels,'knot_bracket_width':str(hi-lo)}
        def core_bound(left,right,maps,lip,lam):
            B=max(abs(left),abs(right))
            f1=mp.mpf('1.01')*(mp.exp(right+1)+e.hp(B)*mp.exp(e.h(right)+1))
            f2=mp.mpf('1.01')*(mp.exp(right+1)+(mp.mpf('.2')+e.hp(B)**2)*mp.exp(e.h(right)+1))
            invphi=1/e.phi(B)
            b=f1*invphi
            bp=(f2+B*f1)*invphi
            C=mp.sqrt(5)*(mp.sqrt(right-left)*bp/2+b*mp.sqrt(e.phi(0)))
            return maps*C*mp.sqrt(lip)/(1-mp.sqrt(lam))
        old=core_bound(mp.mpf(3),mp.mpf(7),6,mp.mpf('2.3'),mp.mpf(575)/613)
        new=core_bound(-mp.mpf(7)/4,mp.mpf(13)/4,18,mp.mpf(271)/200,mp.mpf(271)/290)
        result['central_norm_envelopes']={'old_bound':str(old),'compact_bound':str(new),
             'ratio_of_bounds':str(new/old),'interpretation':'Comparable conservative bounds, not variances or lower bounds.'}
        # Concrete arithmetic counts, not runs at these sizes or proof-onset gates.
        records=[]
        for N in (1000,1000000):
            s=e.schedule(N)
            records.append({'N':N,'depth':s.depth,'iterations':s.iterations,'lattice':s.lattice,
                'radius':str(s.radius),'central_minimum_probability':str(e.Phi(-mp.mpf(13)/4)),
                'tail_guard_probability':str(e.Phi(-s.radius-1)),
                'central_reference_evaluation_bound':36*s.depth,
                'central_knot_upper_bound':int(mp.ceil(36*(N//5)/(1-mp.mpf(271)/290)+36*s.depth+2*(N//5)+4))})
        result['schedule_counts_only']=records
        try:
            e.estimate_from_sources(None)
            raise AssertionError('missing sampled-data guard')
        except PermissionError:
            result['sample_entry_refused']=True
        result['status']='passed'
    except BaseException as exc:
        result['status']='failed'
        result['failure']=repr(exc)
        raise
    finally:
        out.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()

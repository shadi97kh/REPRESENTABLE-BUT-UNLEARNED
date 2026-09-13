import argparse,subprocess,sys,json
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);root=ap.parse_args().run
code=Path(__file__).resolve().parent
cmd=[sys.executable,'-m','pip','install','--disable-pip-version-check','--target',str(code/'_deps'),'beautifulsoup4==4.14.3','lxml==6.0.2']
subprocess.run(cmd,check=True)
(root/'installed_dependencies.json').write_text(json.dumps({'command':cmd,'scope':'Fresh task-local parser dependencies; no global environment changes.'},indent=2)+'\n')
(root/'dependencies.outputs.json').write_text(json.dumps(['installed_dependencies.json']))

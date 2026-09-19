"""Inline the preserved manuscript into a fresh four-entry source directory."""
from pathlib import Path
import re,shutil,json,hashlib
R=Path(__file__).resolve().parents[2];old=R/'papers/interaction_recoverability_iclr2027/v8';new=R/'papers/interaction_recoverability_iclr2027/v9';source=new/'source';source.mkdir(parents=True,exist_ok=False);typesetting=new/'typesetting';typesetting.mkdir();(new/'artifacts').mkdir()
def inline(text,stack=()):
 def repl(m):
  name=m[1];p=old/(name if name.endswith('.tex') else name+'.tex');assert p not in stack,('cyclic include',p);return inline(p.read_text(),stack+(p,))
 return re.sub(r'\\input\{([^}]+)\}',repl,text)
main=(old/'main.tex').read_text();begin=main.index('\\appendix');end=main.index('\\end{document}');appendix=inline(main[begin:end]);main=inline(main[:begin])+r'\input{appendix}'+'\n'+main[end:]
(source/'main.tex').write_text(main);(source/'appendix.tex').write_text(appendix)
# Parse balanced entries and reject conflicting duplicate citation keys.
def entries(text):
 pattern=re.compile(r'@(\w+)\s*\{\s*([^,\s]+)\s*,');offset=0
 while m:=pattern.search(text,offset):
  i=m.end();depth=1;quoted=False
  while i<len(text) and depth:
   c=text[i]
   if c=='\\':i+=2;continue
   if c=='"':quoted=not quoted
   if not quoted:
    if c=='{':depth+=1
    elif c=='}':depth-=1
   i+=1
  assert depth==0,m[2];yield m[2],text[m.start():i];offset=i
bib={};conflicts=[]
for name in ['references.bib','references_new.bib']:
 for key,value in entries((old/name).read_text()):
  if key in bib:assert re.sub(r'\s+','',bib[key])==re.sub(r'\s+','',value),('conflicting citation key',key)
  else:bib[key]=value
(source/'references.bib').write_text('\n\n'.join(bib.values())+'\n')
main=(source/'main.tex').read_text().replace(r'\bibliography{references,references_new}',r'\bibliography{references}');(source/'main.tex').write_text(main)
figs=set(re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',main+appendix));(source/'figures').mkdir()
for name in figs:
 p=old/name;assert p.is_file() and Path(name).parts[0]=='figures';shutil.copy2(p,source/name)
for name in ['iclr2027_conference.sty','iclr2027_conference.bst']:shutil.copy2(old/name,typesetting/name)
(new/'typesetting/build.sh').write_text('''#!/usr/bin/env bash
set -euo pipefail
TYPESETTING_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="${1:-$TYPESETTING_DIR/../source}"
OUTPUT_DIR="${2:-$TYPESETTING_DIR/../artifacts/build}"
SOURCE_DIR="$(realpath "$SOURCE_DIR")"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(realpath "$OUTPUT_DIR")"
export TEXINPUTS="$TYPESETTING_DIR:$SOURCE_DIR:${TEXINPUTS:-}"
export BSTINPUTS="$TYPESETTING_DIR:${BSTINPUTS:-}"
export BIBINPUTS="$SOURCE_DIR:${BIBINPUTS:-}"
cd "$SOURCE_DIR"
tectonic -Z "search-path=$TYPESETTING_DIR" --outdir "$OUTPUT_DIR" --keep-logs --keep-intermediates main.tex
''')
(new/'typesetting/README.md').write_text('External unmodified official anonymous ICLR 2027 dependencies. Build with `bash typesetting/build.sh SOURCE_DIR OUTPUT_DIR`. Tectonic 0.15.0 is required; installed TeX resources/cache remain external. The source directory has exactly main.tex, appendix.tex, references.bib and figures/. No build output belongs there. The additional Tectonic search-path flag complements TEXINPUTS/BSTINPUTS because Tectonic handles its own resource bundle.\n')
(new/'artifacts/initial_inlining_audit.json').write_text(json.dumps(dict(parent='v8',source_entries=sorted(p.name for p in source.iterdir()),bibliography_entries=len(bib),figures=len(figs),official_style_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in typesetting.glob('*.st*')},scope='Initial source layout only; scientific revision and final layout validation pending'),indent=2)+'\n')
print(source, 'inlined',len(bib),'bibliography entries and',len(figs),'figure assets')

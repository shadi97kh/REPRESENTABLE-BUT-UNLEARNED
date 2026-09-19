"""Validate the exact four-entry manuscript layout, includes, bibliography and assets."""
from pathlib import Path
import argparse,re,json,zipfile,tempfile,stat
ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('--report',type=Path);args=ap.parse_args()
def validate(root):
 expected={'main.tex','appendix.tex','references.bib','figures'};assert {p.name for p in root.iterdir()}==expected,'source must have exactly four entries'
 assert (root/'figures').is_dir()
 for p in root.rglob('*'):assert not p.is_symlink(),p
 for n in ['main.tex','appendix.tex','references.bib']:assert (root/n).is_file()
 raw={n:(root/n).read_text() for n in ['main.tex','appendix.tex']};tex={n:re.sub(r'(?<!\\)%[^\n]*','',s) for n,s in raw.items()};includes={n:re.findall(r'\\(?:input|include)\s*\{([^}]+)\}',s) for n,s in tex.items()};assert includes['main.tex'] in [['appendix'],['appendix.tex']] and not includes['appendix.tex'],includes
 assert not any(re.search(r'\\(?:subfile|import|includeonly|lstinputlisting|VerbatimInput)\b',s) for s in tex.values())
 assert re.findall(r'\\bibliography\s*\{([^}]+)\}',tex['main.tex'])==['references']
 figures=set(re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}','\n'.join(tex.values())));actual={str(p.relative_to(root)) for p in (root/'figures').rglob('*') if p.is_file()};assert figures==actual,(figures-actual,actual-figures)
 for n in actual:assert Path(n).suffix.lower() in ['.pdf','.png','.svg','.jpg','.jpeg','.eps'] and len(Path(n).parts)==2,n
 keys=re.findall(r'@\w+\s*\{\s*([^,\s]+)\s*,',(root/'references.bib').read_text());assert len(keys)==len(set(keys))
 cited={k.strip() for s in tex.values() for v in re.findall(r'\\cite\w*(?:\[[^]]*\])?\{([^}]+)\}',s) for k in v.split(',')};assert cited<=set(keys),cited-set(keys)
 for name,s in tex.items():
  for kind,body in re.findall(r'\\begin\{(figure\*?|table\*?)\}(.*?)\\end\{\1\}',s,re.S):
   if '\\caption' in body:
    caption=body.index('\\caption');objects=[m.end() for m in re.finditer(r'\\includegraphics(?:\[[^]]*\])?\{[^}]+\}|\\end\{(?:tabular|tabularx)\}',body)];assert not objects or caption>max(objects),(name,kind,'caption before object')
 return dict(valid=True,source_entries=sorted(expected),tex_files=2,figures=len(figures),bibliography_entries=len(keys),cited_entries=len(cited),external_typesetting_dependencies_required=True,build_outputs_inside_source=False)
if args.source.suffix=='.zip':
 with zipfile.ZipFile(args.source) as z,tempfile.TemporaryDirectory() as d:
  for i in z.infolist():assert not Path(i.filename).is_absolute() and '..' not in Path(i.filename).parts and not stat.S_ISLNK(i.external_attr>>16)
  z.extractall(d);result=validate(Path(d))
else:result=validate(args.source)
if args.report:args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

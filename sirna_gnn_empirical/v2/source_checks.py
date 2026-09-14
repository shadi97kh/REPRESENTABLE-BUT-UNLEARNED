from core import *
import re,shutil
(RUN/'source_audit').mkdir(parents=True,exist_ok=True)
p=OLD/'raw/davis2025_supplements/Supplemental Tables.xlsx';checks=[];rc=lambda s:s.translate(str.maketrans('ACGUT','UGCAA'))[::-1];km=lambda s:{s[i:i+13] for i in range(len(s)-12)}
for sheet,header,cols in [('Supplemental Table 1',5,(4,5,6)),('Supplemental Table 6',4,(2,3,4)),('Supplemental Table 7',4,(2,3,4))]:
 df=pd.read_excel(p,sheet_name=sheet,header=header)
 for j,x in df.iterrows():
  ga,ps,ta=[str(x.iloc[k]) for k in cols];g=''.join(re.findall(r'\([mfr]([ACGU])\)',ga));s=''.join(re.findall(r'\([mfr]([ACGU])\)',ps));t=ta.upper().replace('T','U');checks.append(dict(sheet=sheet,excel_row=j+header+2,compound=str(x.iloc[0]),raw_guide=ga,raw_passenger=ps,raw_target=ta,guide=g,passenger=s,target=t,guide_target_sense_13mer=bool(km(g)&km(t)),guide_target_reverse_complement_13mer=bool(km(g)&km(rc(t))),guide_passenger_reverse_complement_13mer=bool(km(g)&km(rc(s)))))
pd.DataFrame(checks).to_csv(RUN/'source_audit/davis_orientation_rows.csv',index=False)
df=pd.DataFrame(checks);print(df.groupby('sheet')[['guide_target_sense_13mer','guide_target_reverse_complement_13mer','guide_passenger_reverse_complement_13mer']].sum().to_string())
shutil.copy2(p,RUN/'source_audit/davis_original_supplement.xlsx')
write(RUN/'source_audit/davis_source_identity.json',dict(workbook_sha256=sha(p),primary_html_sha256=sha(OLD/'raw/davis2025.html'),primary_url='https://academic.oup.com/nar/article/53/12/gkaf479/8171869',conversion='Decode literal (m/f/r BASE) annotations in printed order only; T->U only in target matching. No reverse/complement repair of an annotated strand.',S1='All 1248 printed guide strings share target-sense 13mers; quarantine retained. No correction validating a conversion was found in checked publisher/PubMed text.',S7='Audit each printed duplex against target reverse complement and its printed passenger. Only structurally consistent rows may enter separately qualified endpoint prediction. No prediction-accuracy-based orientation choice.'))

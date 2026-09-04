#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, shutil
from pathlib import Path
import pandas as pd


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input',required=True)
    p.add_argument('--catalog',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--expected-count',type=int,default=19)
    a=p.parse_args()

    source=Path(a.input); out=Path(a.output)
    cat=pd.read_csv(a.catalog)
    expected=set(cat.candidate_id.astype(int))
    if out.exists(): shutil.rmtree(out)
    (out/'candidates').mkdir(parents=True,exist_ok=True)

    rows=[]
    for sf in source.rglob('summary.json'):
        d=json.loads(sf.read_text(encoding='utf-8'))
        cid=int(d['candidate_id']); rows.append(d)
        shutil.copytree(sf.parent,out/'candidates'/f'candidate_{cid}',dirs_exist_ok=True)

    df=pd.DataFrame(rows)
    if len(df):
        df['candidate_id']=pd.to_numeric(df['candidate_id'],errors='coerce').astype('Int64')
        df=df.sort_values('candidate_id').reset_index(drop=True)
    df.to_csv(out/'ALL19_PM_AUDIT.csv',index=False)

    got=set(df.candidate_id.dropna().astype(int)) if len(df) else set()
    missing=sorted(expected-got)
    if len(expected)!=a.expected_count or missing:
        raise SystemExit(f'QC failed: expected={len(expected)} missing={missing}')

    pair_found=int(df['pair_status'].astype(str).eq('PAIR_FOUND').sum()) if len(df) else 0
    valid_pm=int((~df['classification'].astype(str).eq('INSUFFICIENT_DATA')).sum()) if len(df) else 0
    counts=df['classification'].astype(str).value_counts() if len(df) else pd.Series(dtype=int)
    reasons=df.get('reason',pd.Series(dtype=object)).fillna('NONE').astype(str).value_counts() if len(df) else pd.Series(dtype=int)
    counts.rename_axis('classification').reset_index(name='n_candidates').to_csv(out/'CLASSIFICATION_COUNTS.csv',index=False)
    reasons.rename_axis('reason').reset_index(name='n_candidates').to_csv(out/'REASON_COUNTS.csv',index=False)

    lines=[
        '# Remaining 19 candidates — independent PM search\n\n',
        f'Candidate summaries: **{len(got)}/{len(expected)}**\n\n',
        f'Archive-selected NIRCam epoch pair found: **{pair_found}**\n\n',
        f'Candidates reaching a non-INSUFFICIENT PM classification: **{valid_pm}**\n\n',
        '## Classification counts\n\n'
    ]
    for k,v in counts.items(): lines.append(f'- `{k}`: {int(v)}\n')
    (out/'SUMMARY.md').write_text(''.join(lines),encoding='utf-8')
    (out/'RUN_QC.md').write_text(f'# Run QC\n\nCandidate summaries present: **{len(got)}/{len(expected)}**\n\nMissing: {missing}\n',encoding='utf-8')
    print(df.to_string(index=False) if len(df) else 'No summaries')

if __name__=='__main__': main()

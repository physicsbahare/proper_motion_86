#!/usr/bin/env python3
"""Small read-only VizieR TAP reconnaissance for COSMOS2025 cross-filter controls."""
import json, os, sys, requests
ENDPOINTS=[
 "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync",
 "https://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync",
 "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync",
]
def query(adql):
    failures=[]
    for endpoint in ENDPOINTS:
        try:
            response=requests.get(endpoint,params={
              "REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":adql
            },timeout=35,headers={"User-Agent":"PM86-public-research/1.0"})
            response.raise_for_status()
            if "<INFO" in response.text[:300] or response.text.lstrip().lower().startswith("<html"):
                raise RuntimeError("TAP response not CSV: "+response.text[:280])
            if "ERROR" in response.text[:280] and "column_name" not in response.text[:280]:
                raise RuntimeError("TAP query error: "+response.text[:280])
            return {"endpoint":endpoint,"csv":response.text}
        except Exception as e: failures.append(str(e)[:250])
    raise RuntimeError("All endpoints failed: "+"; ".join(failures))
def main():
    rows={}
    for short,tab in [("phot","J/A+A/704/A339/phot"),("lephare","J/A+A/704/A339/lephare")]:
        q="SELECT column_name, datatype FROM TAP_SCHEMA.columns WHERE table_name='"+tab+"'"
        try:
            result=query(q)
            import pandas as pd
            from io import StringIO
            df=pd.read_csv(StringIO(result["csv"]))
            rows[short]={"endpoint":result["endpoint"],"n_fields":len(df),"columns":df.to_dict(orient="records")}
            print(short, "N_COLUMNS",len(df), "COLUMNS",df.to_dict(orient="records"),flush=True)
        except Exception as e:
            rows[short]={"error":repr(e)}
            print(short,"FAILED",e,flush=True)
    os.makedirs("results_stationary_galaxy_controls",exist_ok=True)
    with open("results_stationary_galaxy_controls/tap_probe.json","w") as f:json.dump(rows,f,indent=2)
    if any("error" in x for x in rows.values()):sys.exit(1)
if __name__=="__main__":main()

"""Task: merge the source spreadsheets (landmark coordinates, centroid size, capture dates) into one table per individual
(outputs/individuals_full.csv). Author: RNV."""
import os, sys
ROOT = (sys.argv[1].rstrip("/\\") + "/") if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__))).replace("\\", "/") + "/"
RAW = ROOT + "data/raw/"; OUTD = ROOT + "outputs/"; os.makedirs(OUTD + "figures/", exist_ok=True)
D = RAW; OUT = OUTD
import openpyxl, pandas as pd, numpy as np, re, datetime as dt
def sheet(path, name):
    ws=openpyxl.load_workbook(path, data_only=True)[name]
    return [list(r) for r in ws.iter_rows(values_only=True)]
pref=lambda s: re.sub(r'\d+$','',str(s))
c1=[r for r in sheet(D+"DADOS MORFOMETRIA (1).xlsx","coordenadas dos marcos anatomic")[1:] if r[0] is not None]
c2=[r for r in sheet(D+"DADOS MORFOMETRIA_TODOS_MOD.xlsx","coordenadas dos marcos anatomic")[1:] if r[0] is not None]
key=lambda r: tuple(float(v) for v in r[1:25])
m2={}
for r in c2: m2.setdefault(key(r),[]).append(r[0])
dup2=sum(1 for v in m2.values() if len(v)>1); print("duplicate coordinate rows in TODOS_MOD:", dup2)
recs=[]
for r in c1:
    k=key(r); ids=m2.get(k,[])
    recs.append(dict(id_orig=r[0],pref_orig=pref(r[0]),id_AI=ids[0] if ids else None,n_match=len(ids)))
df=pd.DataFrame(recs); print("unmatched:", (df.n_match==0).sum(), "multi:", (df.n_match>1).sum())
df["frag"]=df.id_AI.map(pref)
print(pd.crosstab(df.pref_orig, df.frag))
X=np.array([r[1:25] for r in c1],dtype=float).reshape(-1,12,2)
df["cs"]=np.round(np.sqrt(((X-X.mean(axis=1,keepdims=True))**2).sum(axis=(1,2)))*1e-5,3)
# compare with Centroide columns
rows=sheet(D+"DADOS MORFOMETRIA_TODOS.xlsx","Centroide"); hdr=rows[0]
cent={h:sorted(r[j] for r in rows[1:] if isinstance(r[j],(int,float))) for j,h in enumerate(hdr[:9])}
for f,g in df.groupby("frag"):
    best=max(cent, key=lambda c: len(set(np.round(cent[c],3)) & set(g.cs)))
    print(f, "n",len(g), "-> site code", best, "shared values", len(set(np.round(cent[best],3)) & set(g.cs)), "of", len(cent[best]))
dates=[(r[0],r[1]) for r in sheet(D+"Data_coleta_El_cingulata_ morfometria.xlsx","Plan1")[2:] if r[0] is not None and isinstance(r[1],dt.datetime)]
df=df.merge(pd.DataFrame(dates,columns=["id_orig","date"]).drop_duplicates("id_orig"),on="id_orig",how="left")
print("dates matched:", df.date.notna().sum())
df["year"]=df.date.dt.year; df["month"]=df.date.dt.month
df["season"]=np.where(df.month.isin([10,11,12,1,2,3]),"hot-wet","cold-dry")
coords=pd.DataFrame(X.reshape(len(X),24),columns=[f"{a}{i}" for i in range(1,13) for a in ("x","y")])
full=pd.concat([df,coords],axis=1).sort_values(["frag","id_AI"]).reset_index(drop=True)
full.to_csv(OUT+"individuals_full.csv",index=False)
pd.set_option("display.width",200)
print(full.groupby("frag").agg(n=("cs","size"),mean=("cs","mean"),sd=("cs","std"),med=("cs","median"),mn=("cs","min"),mx=("cs","max"),ndates=("date","nunique")).round(2))
print(full.groupby(["frag","year","season"]).size().unstack([1,2]).fillna(0).astype(int))
print("per-fragment dates:"); print(full.groupby("frag").date.agg(lambda s: ", ".join(sorted(set(s.dt.strftime("%d/%m/%Y"))))).to_string())
print("overall:", round(full.cs.mean(),3), round(full.cs.std(),3), round(full.cs.median(),3), full.cs.min(), full.cs.max())
print("HB sheet 22: mean", round(np.mean(cent['HB']),2),"sd",round(np.std(cent['HB'],ddof=1),2))

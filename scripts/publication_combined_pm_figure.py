#!/usr/bin/env python3
"""Compact two-candidate figure from original JWST CAL stamps + independently registered astrometry.
Source image panels are NOT co-registered; apparent astrometric displacement comes
only from frozen independent PM table, and parallax comes from JPL Horizons.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse

SOURCES={
 387132:dict(name="387132 / COSMOS2025 333719",d=19.1,
             first="jw01727139001_04101_00003_nrcblong.npz",
             last="jw06434296001_06201_00002_nrcalong.npz",
             csv="387132_jwst_l2_conditional_motion.csv"),
 1651050:dict(name="1651050 / COSMOS2025 726310",d=11.2,
             first="jw01727065001_04101_00003_nrcblong.npz",
             last="jw05893016008_05101_00001_nrcblong.npz",
             csv="jwst_l2_parallax_conditional_pm.csv"),
}

def load_source(file):
    with np.load(file,allow_pickle=False) as z:
        d={key:z[key].item() if z[key].ndim==0 else z[key].copy() for key in z.files}
    # The 1651050 CAL-QA archive reports coordinates in parent 401px cutout,
    # whereas the 387132 archive reports coordinates in its exported 73px stamp.
    if float(d["catalog_x_local"])>=d["SCI"].shape[1]:
        for key in ("catalog_x_local","gaussian_x_local","com_x_local"):
            d[key]=float(d[key])-float(d["pixel_window_x0"])
        for key in ("catalog_y_local","gaussian_y_local","com_y_local"):
            d[key]=float(d[key])-float(d["pixel_window_y0"])
    return d

def image(ax,d,label):
    arr=d["SCI"].astype(float);err=d["ERR"].astype(float);dq=d["DQ"].astype(np.uint32)
    x=int(round(float(d["catalog_x_local"])));y=int(round(float(d["catalog_y_local"])))
    h=10
    if x-h<0 or y-h<0 or x+h>=arr.shape[1] or y+h>=arr.shape[0]:
        raise ValueError("Catalog position beyond source stamp")
    yy=slice(y-h,y+h+1);xx=slice(x-h,x+h+1)
    image=arr[yy,xx]-np.nanmedian(arr)
    ok=np.isfinite(image)&np.isfinite(err[yy,xx])&(err[yy,xx]>0)&((dq[yy,xx]&(1|2|4|16))==0)
    lo,hi=np.percentile(image[ok],[10,99])
    ax.imshow(np.where(ok,image,np.nan),origin="lower",cmap="gray",vmin=lo,vmax=max(hi,lo+1e-6),interpolation="nearest")
    ax.scatter(float(d["catalog_x_local"])-(x-h),float(d["catalog_y_local"])-(y-h),
               marker="+",s=95,c="#39c3a3",linewidths=2)
    ax.scatter(float(d["gaussian_x_local"])-(x-h),float(d["gaussian_y_local"])-(y-h),
               marker="o",s=105,facecolors="none",edgecolors="#e67e22",linewidths=1.7)
    ax.set_xticks([]);ax.set_yticks([])
    ax.set_title(f"{label}: {d['filter']} (MJD {float(d['mjd']):.1f})",fontsize=10)
    ax.text(.02,.04,"original detector pixels",transform=ax.transAxes,color="white",fontsize=8,
            bbox=dict(facecolor="black",alpha=.5,edgecolor="none",pad=2))

def arrows(ax,row,q,dt):
    obs=dt*np.array([float(row.mu_alpha_cosdec_masyr),float(row.mu_delta_masyr)])
    pi=np.array([float(q.parallax_displacement_east_mas),float(q.parallax_displacement_north_mas)])
    residual=obs-pi
    expected=dt*np.array([float(q.conditional_pm_east_masyr),float(q.conditional_pm_north_masyr)])
    if not np.allclose(residual,expected,rtol=0,atol=0.08):
        raise ValueError(f"Parallax correction and raw measurement disagree: {residual}, {expected}")
    def draw(a,b,color,width):
        ax.annotate("",xy=b,xytext=a,arrowprops=dict(arrowstyle="-|>",color=color,lw=width,mutation_scale=13,shrinkA=0,shrinkB=0))
    ax.axhline(0,c=".8",lw=.8);ax.axvline(0,c=".8",lw=.8)
    draw([0,0],obs,"#252c37",3)
    draw([0,0],pi,"#c47a24",2.3)
    draw(pi,obs,"#127a86",2.5)
    ax.scatter(*obs,s=32,color="#252c37")
    ax.scatter(*pi,s=29,color="#c47a24")
    ax.add_patch(Ellipse(obs,width=2*float(row.sigma_mu_alpha_cosdec_masyr)*dt,
                            height=2*float(row.sigma_mu_delta_masyr)*dt,
                            facecolor="#252c37",edgecolor="none",alpha=.12))
    anchors=np.vstack(([0,0],obs,pi))
    rx=max(80,(np.ptp(anchors[:,0])+2*row.sigma_mu_alpha_cosdec_masyr*dt)*1.35)
    ry=max(70,(np.ptp(anchors[:,1])+2*row.sigma_mu_delta_masyr*dt)*1.5)
    ax.set_xlim(np.min(anchors[:,0])-.22*rx,np.max(anchors[:,0])+.20*rx)
    ax.set_ylim(np.min(anchors[:,1])-.32*ry,np.max(anchors[:,1])+.27*ry)
    ax.set_aspect("equal",adjustable="box")
    ax.grid(alpha=.15)
    ax.set_xlabel(r"$\Delta\alpha^*$ (mas)")
    ax.set_ylabel(r"$\Delta\delta$ (mas)")
    ax.set_title(f"Registered shift, assumed d = {float(q.assumed_distance_pc):.1f} pc",fontsize=10)
    return dict(observed_east_mas=obs[0],observed_north_mas=obs[1],
        parallax_east_mas=pi[0],parallax_north_mas=pi[1],
        conditional_motion_east_mas=residual[0],conditional_motion_north_mas=residual[1])

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,default=Path("publication_inputs"))
    p.add_argument("--output",type=Path,default=Path("publication_check/combined_pm_figure"))
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    cat=pd.read_csv("data/two_pm_candidates_125.csv").set_index("sample_candidate_id")
    fig,axs=plt.subplots(2,3,figsize=(14.3,7.1),gridspec_kw=dict(width_ratios=[1.02,1.02,1.65]),layout="constrained")
    rows=[]
    for i,cid in enumerate((387132,1651050)):
        m=SOURCES[cid]; c=cat.loc[cid]
        input_dir=a.input/str(cid)
        early=load_source(input_dir/"stamps"/"pixel_stamps"/m["first"])
        late=load_source(input_dir/"stamps"/"pixel_stamps"/m["last"])
        if str(early["filter"])!=str(c.filter_early) or str(late["filter"])!=str(c.filter_late):
            raise RuntimeError("Unexpected filter vs independently measured astrometric pair")
        image(axs[i,0],early,"Early"); image(axs[i,1],late,"Late")
        q=pd.read_csv(input_dir/"parallax"/m["csv"])
        q=q.loc[np.isclose(q.assumed_distance_pc,m["d"])]
        if len(q)!=1:raise RuntimeError("Missing unique nominal distance/observer geometry")
        log=arrows(axs[i,2],c,q.iloc[0],float(c.baseline_days)/365.25)
        rows.append(dict(candidate_id=cid,**log))
        axs[i,0].text(-.13,.5,m["name"],rotation=90,ha="center",va="center",
                      transform=axs[i,0].transAxes,fontsize=11,fontweight="bold")
    fig.suptitle("COSMOS-Web: observed displacement and conditional proper motion",fontsize=14)
    items=[
       Line2D([0],[0],marker="+",linestyle="none",markersize=10,color="#39c3a3",label="Catalog / WCS"),
       Line2D([0],[0],marker="o",linestyle="none",markersize=9,markerfacecolor="none",markeredgecolor="#e67e22",label="Gaussian centroid"),
       Line2D([0],[0],color="#252c37",lw=2.8,label="Observed displacement"),
       Line2D([0],[0],color="#c47a24",lw=2.8,label="Assumed-distance parallax"),
       Line2D([0],[0],color="#127a86",lw=2.8,label="Conditional residual motion")]
    fig.legend(handles=items,loc="lower center",bbox_to_anchor=(.52,-.02),ncol=5,fontsize=8,frameon=False)
    for ext in ("pdf","png"):
        fig.savefig(a.output/f"pm_two_candidates_combined.{ext}",dpi=210,bbox_inches="tight",pad_inches=.28)
    pd.DataFrame(rows).to_csv(a.output/"figure_values.csv",index=False)
    (a.output/"READ_ME.txt").write_text(
       "The early and late original CAL detector images are NOT registered to each other. "
       "The displayed source centroids are per-exposure measurements, not a newly "
       "derived displacement. Motion arrows use the published locally registered "
       "two-epoch results. Orange vector uses true JWST Horizons observer parallax "
       "at an ASSUMED photometric distance. Teal vector is distance-conditional "
       "residual, not a fitted trigonometric parallax or independent proper motion. "
       "Gray ellipse only formal apparent-shift errors, excluding cross-filter PSF "
       "systematics and distance uncertainty; the two source plots have different ranges.\n")
    plt.close(fig)
    print(pd.DataFrame(rows).to_string(index=False))

if __name__=="__main__":
    main()

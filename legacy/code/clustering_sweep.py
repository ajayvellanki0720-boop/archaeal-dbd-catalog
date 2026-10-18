#!/usr/bin/env python
# Clustering-stability sweep: re-cluster 25,472 orphans across id/cov/mode; compare to canonical 250 families.
import subprocess, os, csv, itertools, time, shutil
BASE="/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/workspaces/624da9a9-94cf-4316-b2b2-684c2f60f45d"
os.chdir(BASE)
BIN="/Users/ajayvellanki/.claude-science/conda/envs/archaea-dormancy/bin"
FAA=f"{BASE}/round2/orphan_proteins_refetched.faa"
SW=f"{BASE}/round2/sweep"; os.makedirs(SW,exist_ok=True)
os.environ["PATH"]=BIN+":"+os.environ["PATH"]

# canonical family membership: map member-protein -> family_id (canonical 250)
canon={}
import csv as _csv
with open(f"{BASE}/repro_pkg/archaeal_regulators/data/novel_family_membership.csv") as fh:
    for r in _csv.DictReader(fh):
        canon[r["member"]]=r["family_id"]
print("canonical members mapped:",len(canon),flush=True)

# build mmseqs DB once
db=f"{SW}/orphanDB"
subprocess.run([f"{BIN}/mmseqs","createdb",FAA,db],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)

# canonical setting for reference: id0.30 cov0.6 cov-mode0 cluster-mode1 (greedy set-cover)
grid=[]
for idt in [0.25,0.30,0.35,0.40]:
    for cov in [0.5,0.6,0.8]:
        for cmode in [0,1]:   # 0=set-cover(greedy), 1=connected-component
            grid.append((idt,cov,cmode))

def run_cluster(idt,cov,cmode,tag):
    out=f"{SW}/clu_{tag}"; tmp=f"{SW}/tmp_{tag}"
    os.makedirs(tmp,exist_ok=True)
    cmd=[f"{BIN}/mmseqs","cluster",db,out,tmp,
         "--min-seq-id",str(idt),"-c",str(cov),"--cov-mode","0","--cluster-mode",str(cmode),
         "--threads","6","-v","0"]
    subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
    tsv=f"{SW}/clu_{tag}.tsv"
    subprocess.run([f"{BIN}/mmseqs","createtsv",db,db,out,tsv],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
    # parse rep->members
    m2rep={}
    for ln in open(tsv):
        rep,mem=ln.rstrip("\n").split("\t")
        m2rep[mem]=rep
    shutil.rmtree(tmp,ignore_errors=True)
    return m2rep

# stability metric: for each canonical family, do its members stay together?
# We compute: (a) # new clusters, (b) size>=4 clusters, and per canonical family the
# fraction of its members that fall in a single dominant new cluster (cohesion),
# and whether that dominant cluster is "pure" (mostly this family).
from collections import defaultdict
def compare(m2rep):
    # group canonical members by family; only members present in clustering
    fam_members=defaultdict(list)
    for mem,fam in canon.items():
        if mem in m2rep: fam_members[fam].append(mem)
    # new cluster sizes
    clsize=defaultdict(int)
    for mem,rep in m2rep.items(): clsize[rep]+=1
    n_clusters=len(set(m2rep.values()))
    n_ge4=sum(1 for c,s in clsize.items() if s>=4)
    # per family cohesion + fate
    fates={"stable":0,"split":0,"merged":0,"lost":0}
    rows=[]
    # reverse: new cluster -> canonical families present
    rep2fams=defaultdict(lambda: defaultdict(int))
    for mem,fam in canon.items():
        if mem in m2rep: rep2fams[m2rep[mem]][fam]+=1
    for fam,mems in fam_members.items():
        reps=[m2rep[m] for m in mems]
        from collections import Counter
        rc=Counter(reps)
        dom_rep,dom_n=rc.most_common(1)[0]
        cohesion=dom_n/len(mems)               # frac of family in its biggest new cluster
        # purity of that dominant cluster: is it mostly this family?
        domcl_fams=rep2fams[dom_rep]
        purity=domcl_fams[fam]/sum(domcl_fams.values())
        # classify
        if cohesion>=0.8 and purity>=0.8: fate="stable"
        elif cohesion<0.8 and purity>=0.8: fate="split"
        elif cohesion>=0.8 and purity<0.8: fate="merged"
        else: fate="split_and_merged"
        rows.append((fam,len(mems),round(cohesion,3),round(purity,3),fate))
    return n_clusters,n_ge4,rows

summ=open("round2/clustering_stability_summary.csv","w",newline="")
sw=csv.writer(summ); sw.writerow(["min_seq_id","coverage","cluster_mode","mode_name","n_clusters","n_clusters_ge4","n_families_stable","n_families_split","n_families_merged","n_families_other","frac_stable"])
perfam=open("round2/clustering_stability_perfamily.csv","w",newline="")
pw=csv.writer(perfam); pw.writerow(["min_seq_id","coverage","cluster_mode","family_id","n_members","cohesion","purity","fate"])
t0=time.time()
for idt,cov,cmode in grid:
    tag=f"id{int(idt*100)}_c{int(cov*100)}_m{cmode}"
    try:
        m2rep=run_cluster(idt,cov,cmode,tag)
        nc,nge4,rows=compare(m2rep)
        from collections import Counter
        fc=Counter(r[4] for r in rows)
        stable=fc.get("stable",0); split=fc.get("split",0); merged=fc.get("merged",0)
        other=len(rows)-stable-split-merged
        sw.writerow([idt,cov,cmode,"setcover" if cmode==0 else "conncomp",nc,nge4,stable,split,merged,other,round(stable/max(len(rows),1),3)])
        summ.flush()
        for r in rows: pw.writerow([idt,cov,cmode]+list(r))
        perfam.flush()
        # cleanup big files
        for f in os.listdir(SW):
            if f.startswith(("clu_","tmp_")): 
                pp=os.path.join(SW,f)
                if os.path.isfile(pp): os.remove(pp)
                elif os.path.isdir(pp): shutil.rmtree(pp,ignore_errors=True)
        print(f"{tag}: clusters={nc} ge4={nge4} stable={stable}/{len(rows)} {time.time()-t0:.0f}s",flush=True)
    except Exception as e:
        print(f"{tag} FAILED {e}",flush=True)
summ.close(); perfam.close()
print("SWEEP DONE",round(time.time()-t0),"s",flush=True)

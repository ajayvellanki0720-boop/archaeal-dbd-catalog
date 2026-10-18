#!/usr/bin/env python
import json, os, gzip, urllib.request, time, re
from collections import defaultdict
BASE="/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/workspaces/624da9a9-94cf-4316-b2b2-684c2f60f45d"
os.chdir(BASE)
WORK=f"{BASE}/round2/tgt_work"; os.makedirs(WORK,exist_ok=True)
J=json.load(open("handoff/target_harvest_tasks.json"))
tasks=J["tasks"]; memmap=J["memmap"]; targets=J["targets"]
CATS={"TA_toxin": re.compile(r"toxin|VapC|antitoxin|HicA|HicB|RelE|ParE|MazF|PIN ",re.I),
      "CRISPR_defense": re.compile(r"CRISPR|Cas[0-9]|argonaute|defense",re.I)}
def rd(url,dst):
    try:
        with urllib.request.urlopen(urllib.request.Request(url),timeout=120) as r: open(dst,"wb").write(r.read()); return True
    except Exception: return False
def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTN","TGCAN"))
def parse_ft(path):
    rows=[]
    with gzip.open(path,"rt") as f:
        hd=f.readline().rstrip("\n").split("\t"); ix={c:i for i,c in enumerate(hd)}
        for line in f:
            p=line.rstrip("\n").split("\t")
            if not p or p[0]!="CDS": continue
            try: rows.append((p[ix["genomic_accession"]],int(p[ix["start"]]),int(p[ix["end"]]),p[ix["strand"]],p[ix["product_accession"]],p[ix["name"]]))
            except: continue
    return rows
def load_fna(path):
    seqs={}; name=None; buf=[]
    with gzip.open(path,"rt") as f:
        for line in f:
            if line.startswith(">"):
                if name: seqs[name]="".join(buf)
                name=line[1:].split()[0]; buf=[]
            else: buf.append(line.strip())
        if name: seqs[name]="".join(buf)
    return seqs
def upstream(seq,start,end,strand,n=250):
    # start/end 1-based inclusive; strand + -> upstream is before start; - -> after end (revcomp)
    if strand=="+":
        a=max(0,start-1-n); return seq[a:start-1]
    else:
        return revcomp(seq[end:end+n])

REG=open("round2/regulator_upstreams.tsv","w"); REG.write("family\tgenome\tprotein\tupstream\n")
TGT=open("round2/target_upstreams.tsv","w"); TGT.write("family\ttarget_cat\tgenome\ttarget_protein\tupstream\n")
t0=time.time(); done=0; miss=0; nreg=0; ntgt=0
for i,t in enumerate(tasks):
    g=t["genome"]; fturl=t["url"]
    fnaurl=fturl.replace("_feature_table.txt.gz","_genomic.fna.gz")
    ftgz=f"{WORK}/{g}.ft.gz"; fnagz=f"{WORK}/{g}.fna.gz"
    if not (rd(fturl,ftgz) and rd(fnaurl,fnagz)): miss+=1; 
    else:
        try:
            cds=parse_ft(ftgz); fna=load_fna(fnagz)
        except Exception: cds=None
        for p in (ftgz,fnagz):
            if os.path.exists(p): os.remove(p)
        if cds and fna:
            tgtfams=memmap.get(g,{})
            bycontig=defaultdict(list)
            for r in cds: bycontig[r[0]].append(r)
            for contig,lst in bycontig.items():
                lst.sort(key=lambda x:x[1])
                seq=fna.get(contig)
                if not seq: continue
                for j,r in enumerate(lst):
                    fam=tgtfams.get(r[4])
                    if not fam: continue
                    # regulator upstream
                    us=upstream(seq,r[1],r[2],r[3])
                    if len(us)>=200: REG.write(f"{fam}\t{g}\t{r[4]}\t{us}\n"); nreg+=1
                    # target neighbors within +-5
                    tcat=targets[fam]; pat=CATS[tcat]
                    for off in range(-5,6):
                        k=j+off
                        if k<0 or k>=len(lst) or off==0: continue
                        nb=lst[k]
                        if pat.search(nb[5] or ""):
                            tus=upstream(seq,nb[1],nb[2],nb[3])
                            if len(tus)>=200: TGT.write(f"{fam}\t{tcat}\t{g}\t{nb[4]}\t{tus}\n"); ntgt+=1
            done+=1
    if (i+1)%40==0:
        REG.flush(); TGT.flush(); print(f"{i+1}/{len(tasks)} {time.time()-t0:.0f}s done={done} miss={miss} reg={nreg} tgt={ntgt}",flush=True)
REG.close(); TGT.close()
print(f"TARGET HARVEST DONE done={done} miss={miss} reg={nreg} tgt={ntgt} wall={time.time()-t0:.0f}s",flush=True)

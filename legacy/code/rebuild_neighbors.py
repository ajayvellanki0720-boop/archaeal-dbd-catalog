#!/usr/bin/env python
# Reconstruct +-5 CDS neighborhoods, then dependence-corrected enrichment.
import json, os, gzip, urllib.request, time, re
from collections import defaultdict
BASE="/Users/ajayvellanki/.claude-science/orgs/a85da623-e9bf-4a84-bae1-c872a64348e0/workspaces/624da9a9-94cf-4316-b2b2-684c2f60f45d"
os.chdir(BASE)
WORK=f"{BASE}/round2/ft_work"; os.makedirs(WORK,exist_ok=True)
tasks=json.load(open("handoff/ft_tasks.json"))["tasks"]
memmap={g:v for g,v in json.load(open("handoff/member_family_map.json")).items()}

CATS = {
 "transporter_ABC": r"ABC transporter|ATP-binding|permease",
 "transporter_MFS": r"MFS transporter|major facilitator",
 "transport_other": r"transporter|symporter|antiporter|channel|efflux",
 "TA_toxin": r"toxin|VapC|antitoxin|HicA|HicB|RelE|ParE|MazF|PIN ",
 "stress_USP": r"universal stress",
 "methyltransferase": r"methyltransferase|SAM-dependent",
 "acetyltransferase": r"acetyltransferase|GNAT",
 "oxidoreductase": r"oxidoreductase|dehydrogenase|reductase|FAD|NAD",
 "metallohydrolase": r"metallo-hydrolase|MBL fold|hydrolase",
 "peptidase": r"peptidase|protease",
 "kinase": r"kinase|phosphatase",
 "signal_HalOD1": r"HalOD1|output domain|PAS|GAF|signal",
 "chaperone_proteostasis": r"chaperone|protease|Clp|Hsp|GroEL|DnaK|degradation|proteasome",
 "nucleotide_cyclase": r"cyclase|diguanylate|c-di-GMP|adenylate",
 "DNA_repair_replication": r"repair|helicase|nuclease|recombinase|topoisomerase|replicat",
 "metabolism_biosynth": r"synthase|synthetase|biosynth|aminotransferase|decarboxylase|dehydratase",
 "restriction_modification": r"restriction|modification|methylase",
 "CRISPR_defense": r"CRISPR|Cas[0-9]|argonaute|defense",
 "regulator_HTH": r"helix-turn-helix|winged helix|transcriptional regulator|DNA-binding",
 "hypothetical": r"hypothetical protein|DUF|uncharacterized",
}
CATSc={c:re.compile(p,re.I) for c,p in CATS.items()}
def categorize(name):
    if not isinstance(name,str): return ["other"]
    h=[c for c,p in CATSc.items() if p.search(name)]
    return h or ["other"]

def parse_ft(path):
    rows=[]
    with gzip.open(path,"rt") as f:
        header=f.readline().rstrip("\n").split("\t")
        idx={c:i for i,c in enumerate(header)}
        for line in f:
            p=line.rstrip("\n").split("\t")
            if not p or p[0]!="CDS": continue
            try: rows.append((p[idx["genomic_accession"]],int(p[idx["start"]]),int(p[idx["end"]]),p[idx["strand"]],p[idx["product_accession"]],p[idx["name"]]))
            except: continue
    return rows

def rd(url,dst):
    try:
        with urllib.request.urlopen(urllib.request.Request(url),timeout=90) as r: data=r.read()
        open(dst,"wb").write(data); return True
    except Exception: return False

# neighbor records: family_id, genome, genus, group, offset, category (exploded)
NEI=open("round2/neighbor_records.csv","w"); NEI.write("family_id,genome,genus,group,offset,category\n")
# also per-genome CDS category pool for null
POOL=open("round2/genome_cds_pool.csv","w"); POOL.write("genome,category,count\n")
t0=time.time(); done=0; miss=0
for i,t in enumerate(tasks):
    g=t["genome"]; genus=t.get("genus") or "NA"; grp=t.get("group") or "NA"
    gz=f"{WORK}/{g}.gz"
    if not rd(t["url"],gz): miss+=1; continue
    try: cds=parse_ft(gz)
    except Exception: miss+=1; os.remove(gz) if os.path.exists(gz) else None; continue
    os.remove(gz)
    if not cds: continue
    # per-genome pool
    poolc=defaultdict(int)
    bycontig=defaultdict(list)
    for r in cds:
        for c in categorize(r[5]): poolc[c]+=1
        bycontig[r[0]].append(r)
    for c,n in poolc.items(): POOL.write(f"{g},{c},{n}\n")
    tgt=memmap.get(g,{})
    if not tgt: continue
    for contig,lst in bycontig.items():
        lst.sort(key=lambda x:x[1])
        for j,r in enumerate(lst):
            fam=tgt.get(r[4])
            if not fam: continue
            for off in range(-5,6):
                k=j+off
                if k<0 or k>=len(lst) or off==0: continue
                nb=lst[k]
                for c in categorize(nb[5]):
                    NEI.write(f"{fam},{g},{genus},{grp},{off},{c}\n")
    done+=1
    if (i+1)%50==0:
        NEI.flush(); POOL.flush(); print(f"{i+1}/{len(tasks)} {time.time()-t0:.0f}s done={done} miss={miss}",flush=True)
NEI.close(); POOL.close()
print(f"NEIGHBOR REBUILD DONE done={done} miss={miss} wall={time.time()-t0:.0f}s",flush=True)

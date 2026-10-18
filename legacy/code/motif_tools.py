"""Vectorized numpy motif discovery + null controls for archaeal operator analysis."""
import numpy as np
from collections import defaultdict

B="ACGT"; BI={b:i for i,b in enumerate(B)}

def _prep(seqs, L):
    """Keep seqs with len>=L, take last L bases (nearest CDS), one-hot (n,L,4) float32."""
    keep=[s for s in seqs if len(s)>=L]
    oh=np.zeros((len(keep),L,4),dtype=np.float32)
    for i,s in enumerate(keep):
        t=s[-L:]
        for j,c in enumerate(t):
            k=BI.get(c,-1)
            if k>=0: oh[i,j,k]=1.0
    return oh

def _bg(oh):
    c=oh.sum((0,1)); c=np.clip(c,1.0,None); return c/c.sum()

def find_motif_em(oh, w, n_starts=4, n_iter=25, rng=None):
    """Vectorized ZOOPS EM. oh:(n,L,4). Returns (pwm(w,4), ic, ll)."""
    if rng is None: rng=np.random.default_rng(0)
    n,L,_=oh.shape
    if n<4 or L<w: return None,0.0,None
    npos=L-w+1
    bg=_bg(oh); logbg=np.log(bg)
    # precompute bg score per window position (constant across EM)
    bgsc=np.zeros((n,npos),dtype=np.float32)
    for j in range(w): bgsc+=oh[:,j:j+npos,:]@logbg
    best=None
    for st in range(n_starts):
        i0=rng.integers(n); p0=rng.integers(npos)
        pwm=np.full((w,4),0.1,dtype=np.float32)
        seed=oh[i0,p0:p0+w,:]
        pwm+=0.6*seed
        pwm/=pwm.sum(1,keepdims=True)
        for it in range(n_iter):
            logpwm=np.log(pwm)
            sc=np.zeros((n,npos),dtype=np.float32)
            for j in range(w): sc+=oh[:,j:j+npos,:]@logpwm[j]
            llr=sc-bgsc
            llr-=llr.max(1,keepdims=True)
            wt=np.exp(llr); wt/=wt.sum(1,keepdims=True)   # (n,npos) posterior
            counts=np.full((w,4),0.25,dtype=np.float32)
            for j in range(w):
                counts[j]+=np.einsum('np,npb->b', wt, oh[:,j:j+npos,:])
            newpwm=counts/counts.sum(1,keepdims=True)
            d=np.abs(newpwm-pwm).max(); pwm=newpwm
            if d<1e-4: break
        # final LL (best position per seq)
        logpwm=np.log(pwm)
        sc=np.zeros((n,npos),dtype=np.float32)
        for j in range(w): sc+=oh[:,j:j+npos,:]@logpwm[j]
        ll=float((sc-bgsc).max(1).sum())
        ic=float((pwm*np.log2(np.clip(pwm,1e-9,None)/bg)).sum())
        if best is None or ll>best[2]: best=(pwm,ic,ll)
    return best

def consensus(pwm): return "".join("ACGT"[i] for i in pwm.argmax(1))
def palindrome_mismatches(c):
    rc=c.translate(str.maketrans("ACGT","TGCA"))[::-1]
    return sum(1 for a,b in zip(c,rc) if a!=b)
def gc_content(seqs):
    tot=sum(len(s) for s in seqs); gc=sum(s.count('G')+s.count('C') for s in seqs)
    return gc/max(1,tot)

def mono_shuffle(s,rng):
    a=list(s); rng.shuffle(a); return "".join(a)

def dinuc_shuffle(s,rng):
    s=[c for c in s if c in BI]
    if len(s)<4: return "".join(s)
    edges=defaultdict(list)
    for i in range(len(s)-1): edges[s[i]].append(s[i+1])
    for _ in range(20):
        e={k:v[:] for k,v in edges.items()}
        for k in e: rng.shuffle(e[k])
        seq=[s[0]]; cur=s[0]; ok=True
        for _ in range(len(s)-1):
            if not e[cur]: ok=False; break
            nxt=e[cur].pop(); seq.append(nxt); cur=nxt
        if ok and len(seq)==len(s): return "".join(seq)
    return mono_shuffle("".join(s),rng)

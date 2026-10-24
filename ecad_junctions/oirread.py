import re, numpy as np
PAT=re.compile(rb'((?:REF_LSM\d|z\d{3}_\d+_\d+)_([0-9a-f]{8}-[0-9a-f-]{27})_(\d+))\x00?')
def parse_chunks(f):
    out=[]
    for m in re.finditer(rb'(?:REF_LSM\d|z\d{3}_\d+_\d+)_[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}_\d+',f):
        p=m.start(); nl=int.from_bytes(f[p-4:p],'little')
        if nl!=len(m.group()): continue
        name=m.group().decode()
        L=int.from_bytes(f[p+nl:p+nl+4],'little')
        s=p+nl+8
        parts=name.split('_'); idx=int(parts[-1]); uid=parts[-2]; pref='_'.join(parts[:-2])
        out.append((p,pref,uid,idx,s,L))
    return out
def read_oir(path):
    f=open(path,'rb').read()
    ch=parse_chunks(f)
    # runs: consecutive chunks with same pref+uid and increasing idx
    runs=[]
    for p,pref,uid,idx,s,L in ch:
        if runs and runs[-1]['uid']==uid and runs[-1]['pref']==pref and idx==runs[-1]['last']+1:
            r=runs[-1]; r['last']=idx; r['data'].append(f[s:s+L])
        else: runs.append(dict(pref=pref,uid=uid,last=idx,start=p,data=[f[s:s+L]]))
    return f,runs

# channel uuid-order -> name (by FV3000 emission windows / laser lines)
def load_stack(path, nz=None):
    """returns dict ch-> (Z,H,W) uint16 for channels dapi, actin, dnmt3a, ecad"""
    f,runs=read_oir(path)
    planes={}
    for r in runs:
        if not r['pref'].startswith('z'): continue
        n=sum(map(len,r['data']))
        if n!=800*800*2: continue
        z=int(r['pref'][1:4])
        planes[(z,r['uid'])]=np.frombuffer(b''.join(r['data']),dtype='<u2').reshape(800,800)
    uid2name={'bde2e175':'dapi','dda7677d':'actin','7abca5fe':'dnmt3a','98919cc1':'ecad'}
    zs=sorted({k[0] for k in planes})
    out={}
    for u,nm in uid2name.items():
        out[nm]=np.stack([planes[(z,k)] for z in zs for k in [kk for kk in [x[1] for x in planes if x[0]==z] if kk.startswith(u)]])
    return out

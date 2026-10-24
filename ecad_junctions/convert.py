import sys, re, json, os, glob
import numpy as np
from multiprocessing import Pool
from oirread import load_stack
UID={'bde2e175':'dapi','dda7677d':'actin','7abca5fe':'dnmt3a','98919cc1':'ecad'}
def meta(path):
    f=open(path,'rb').read()
    i=f.find(b'<lsmimage:imageProperties'); j=f.find(b'</lsmimage:imageProperties>',i)
    x=re.sub(r'\s+',' ',f[i:j].decode('ascii','ignore'))
    m={}
    for uid,nm in UID.items():
        r=re.search(r'<lsmparam:pmt channelId="'+uid+r'[^"]*" detectorId="[^"]*"> <lsmparam:voltage>([\d.]+)</lsmparam:voltage> <lsmparam:offset>([\d.-]+)</lsmparam:offset> <lsmparam:gain>([\d.]+)',x)
        if r: m[nm+'_pmt']=[float(v) for v in r.groups()]
    for ld in ['LD405','LD488','LD561','LD640','LD594']:
        r=re.search(r'<commonparam:mainLaser laserDataId="'+ld+r'_\d+Imaging_main"> <commonparam:power>([\d.]+)</commonparam:power> <commonparam:transmissivity>([\d.]+)',x)
        if r: m[ld+'_power_transm']=[float(v) for v in r.groups()]
    for k,pat in [('zoom',r'<lsmparam:zoom>([\d.]+)'),('speed',r'<commonparam:speed>([\d.]+)'),('zstep',r'<commonparam:step>([\d.]+)</commonparam:step> <commonparam:maxSize>(\d+)')]:
        r=re.search(pat,x)
        if r: m[k]=r.groups() if len(r.groups())>1 else r.group(1)
    return m
def run(args):
    p,DST=args
    name=os.path.basename(p)[:-4]; out=f"{DST}/{name}.npz"
    if os.path.exists(out): return name,'cached'
    try:
        s=load_stack(p); md=meta(p)
        np.savez_compressed(out,**s,meta=json.dumps(md))
        return name,{k:v.shape for k,v in s.items()}['ecad']
    except Exception as e: return name,'ERR '+repr(e)
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description='Convert Olympus .oir stacks to cached .npz (uint16 stacks + acquisition metadata)')
    ap.add_argument('--src',required=True,help='folder with .oir files'); ap.add_argument('--dst',default='cache'); ap.add_argument('--procs',type=int,default=2)
    a=ap.parse_args(); os.makedirs(a.dst,exist_ok=True)
    files=sorted(glob.glob(os.path.join(a.src,'*.oir')))
    with Pool(a.procs) as pl:
        for r in pl.imap_unordered(run,[(f,a.dst) for f in files]): print(r,flush=True)

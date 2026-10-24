"""E-cadherin junction quantification for wound-proximal confocal images (Olympus FV3000 .oir -> cache npz).
Channels: dapi(405) actin(488) dnmt3a(561) ecad(647).  Run: python ecad_quant.py [--pixel-um 0.265]
"""
import os, re, sys, json, glob, argparse
import numpy as np, pandas as pd
from scipy import ndimage as ndi
from skimage import filters, morphology, measure, segmentation, feature, io
from skimage.filters import sato

PIXEL_UM = 0.265          # ASSUMPTION (FV3000, 60x, zoom 1, 800 px); not stored in these .oir files
ZONE_UM = 50.0            # analysis zone: tissue within this distance of the wound edge
WIN = 3                   # number of z-slices averaged around the E-cad best-focus plane
MIN_SEG_PX = 30           # min skeleton component length (px) kept
CONT_PX = 50              # "continuous" junction segment length (px)
MODE = os.environ.get('ECAD_MODE', 'hyst')   # 'hyst' (default, used for all reported results) | 'global' (legacy single global threshold) | 'localnorm' (experimental, adds noise)
V2 = MODE in ('hyst', 'localnorm'); LOCALNORM = MODE == 'localnorm'   
LOW_K, HIGH_K = 2.0, 3.0   # hysteresis thresholds in MAD units (v2)
LOW_THR = None
BAND_R = 3                # half-width of scaffold band (px)

def norm99(a):
    hi = np.percentile(a, 99.5); return a / (hi + 1e-9)

def load(name, cache):
    d = np.load(f"{cache}/{name}.npz")
    return {k: d[k].astype(np.float32) for k in ['dapi', 'actin', 'dnmt3a', 'ecad']}, json.loads(str(d['meta']))

def best_window(ecad):
    sc = []
    for z in range(ecad.shape[0]):
        a = ecad[z]; hp = filters.gaussian(a, 1) - filters.gaussian(a, 6)
        sc.append(np.var(hp) / (a.mean() + 1e-6) ** 2)
    sc = np.array(sc); w = min(WIN, len(sc))
    tot = [sc[i:i + w].sum() for i in range(len(sc) - w + 1)]
    i = int(np.argmax(tot)); return i, i + w

def tissue_masks(P):
    S = sum(norm99(filters.gaussian(P[c], 2)) for c in ['actin', 'dnmt3a', 'ecad'])
    S = filters.gaussian(S, 3)
    t = filters.threshold_otsu(S); t = 0.6 * t     # permissive: keep dim cytoplasm
    tissue = S > t
    tissue = morphology.remove_small_objects(tissue, 3000)
    tissue = ndi.binary_fill_holes(ndi.binary_closing(tissue, morphology.disk(5)))
    bg = ~tissue
    lab = measure.label(bg); gap = np.zeros_like(tissue)
    border = set(np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]])) - {0}
    cands = [(int((lab == l).sum()), l) for l in border]
    cands = [c for c in cands if c[0] > 4000]
    if cands: gap = lab == max(cands)[1]
    tissue = ~gap
    return tissue, gap

def nuclei_mask(dapi_stack, tissue):
    D = filters.gaussian(dapi_stack.max(0), 1.5)
    t = filters.threshold_otsu(D[tissue]) if tissue.sum() else filters.threshold_otsu(D)
    m = (D > 0.6 * t) & tissue
    m = morphology.remove_small_objects(ndi.binary_opening(m, morphology.disk(2)), 80)
    m = ndi.binary_fill_holes(m)
    dist = ndi.distance_transform_edt(m)
    pk = feature.peak_local_max(filters.gaussian(dist, 1), min_distance=9, labels=measure.label(m), exclude_border=False)
    mk = np.zeros(m.shape, int); mk[tuple(pk.T)] = np.arange(1, len(pk) + 1)
    lab = segmentation.watershed(-dist, mk, mask=m)
    return lab

def zone_metrics(Es, E, Eb, R, thr, zone, nmask, nlab, satmask, bgv, pixel_um, sfx):
    rb = (filters.apply_hysteresis_threshold(R, LOW_THR, thr) if (V2 and LOW_THR is not None) else (R > thr)) & zone
    rb = morphology.remove_small_objects(rb, 40)
    sk = morphology.skeletonize(rb)
    skl = measure.label(sk, connectivity=2)
    sizes = np.bincount(skl.ravel())[1:] if skl.max() else np.array([])
    if len(sizes):
        sk = np.isin(skl, np.where(sizes >= MIN_SEG_PX)[0] + 1); skl = measure.label(sk, connectivity=2)
        sizes = np.bincount(skl.ravel())[1:] if skl.max() else np.array([])
    nz = int(zone.sum()); sk_len = int(sk.sum()); area_um2 = nz * pixel_um ** 2
    skd = ndi.binary_dilation(sk, morphology.disk(1))
    cyt = zone & ~ndi.binary_dilation(sk, morphology.disk(6)) & ~ndi.binary_dilation(nmask, morphology.disk(3)) & ~satmask
    ann = ndi.binary_dilation(sk, morphology.disk(7)) & ~ndi.binary_dilation(sk, morphology.disk(4)) & zone
    Jv = float(Es[skd].mean() - bgv) if sk_len > 50 else np.nan
    Cv = float(Es[cyt].mean() - bgv) if cyt.sum() > 500 else np.nan
    nzone = len(np.unique(nlab[zone & nmask]))
    o = {
        'zone_area_um2': area_um2, 'n_nuclei': nzone,
        'skel_density': sk_len * pixel_um / area_um2 * 100 if nz else np.nan,            # um of junction per 100 um2
        'skel_len_per_nucleus_um': sk_len * pixel_um / max(nzone, 1),
        'n_segments': int(len(sizes)),
        'segments_per_100um': len(sizes) / (sk_len * pixel_um) * 100 if sk_len else np.nan,
        'mean_segment_um': sizes.mean() * pixel_um if len(sizes) else np.nan,
        'continuous_frac': sizes[sizes >= CONT_PX].sum() / sizes.sum() if len(sizes) else np.nan,
        'JC_ratio': Jv / Cv if (Cv and Cv > 0 and not np.isnan(Jv)) else np.nan,
        'junction_bgsub': Jv, 'cytoplasm_bgsub': Cv,
        'ridge_contrast': float((Es[skd].mean() - bgv) / max(Es[ann].mean() - bgv, 1e-6)) if sk_len > 50 and ann.sum() > 100 else np.nan,
        'ridge_frac': float((R[zone] > thr).mean()) if nz else np.nan,
        'ridge_snr_p95': float((np.percentile(R[zone], 95) - np.median(R[zone])) / (np.median(np.abs(R[zone] - np.median(R[zone]))) * 1.4826 + 1e-9)) if nz else np.nan,
        'total_ecad_bgsub': float(Eb[zone].mean()) if nz else np.nan,
        'nuclear_to_cyt_ecad': float((Es[nmask & zone].mean() - bgv) / max(Cv, 1e-6)) if (nmask & zone).sum() > 200 and Cv and Cv > 0 else np.nan,
    }
    return {f'{k}{sfx}': v for k, v in o.items()}, sk


def cell_table(name, Es, R, med, mad, sk, nlab, tissue, gap, bgv, pixel_um):
    """Per-cell (nucleus-seeded territory) metrics. Territory = tissue pixels nearer to this nucleus than to any other."""
    H, W = Es.shape
    ids = np.unique(nlab[nlab > 0])
    if len(ids) < 3: return []
    idx = ndi.distance_transform_edt(nlab == 0, return_distances=False, return_indices=True)
    terr = np.where(tissue, nlab[idx[0], idx[1]], 0)
    nmask = nlab > 0
    # complete territories only: not touching image border or wound gap
    bad = np.zeros((H, W), bool); bad[0, :] = bad[-1, :] = bad[:, 0] = bad[:, -1] = True
    bad |= ndi.binary_dilation(gap, morphology.disk(4))
    bad_ids = set(np.unique(terr[bad])) - {0}
    skl = measure.label(sk, connectivity=2); sizes = np.r_[0, np.bincount(skl.ravel())[1:]] if skl.max() else np.array([0])
    cont = (sizes[skl] >= CONT_PX) & sk
    skd = ndi.binary_dilation(sk, morphology.disk(1))
    cyt_all = ~ndi.binary_dilation(sk, morphology.disk(6)) & ~ndi.binary_dilation(nmask, morphology.disk(3))
    Rn = (R - med) / mad
    areas = ndi.sum(np.ones_like(terr), terr, ids)
    med_area = np.median(areas)
    rows = []
    sk_px = ndi.sum(sk, terr, ids); ct_px = ndi.sum(cont, terr, ids); Rm = ndi.mean(Rn, terr, ids)
    for k, i in enumerate(ids):
        a = areas[k]
        if i in bad_ids or a < 0.3 * med_area or a > 3 * med_area: continue
        m = terr == i
        j = skd & m; c = cyt_all & m
        Jv = Es[j].mean() - bgv if j.sum() > 10 else np.nan
        Cv = Es[c].mean() - bgv if c.sum() > 50 else np.nan
        rows.append(dict(image=name, cell_id=int(i), territory_area_um2=a * pixel_um ** 2,
            junction_len_um=sk_px[k] * pixel_um, junction_density=sk_px[k] * pixel_um / (a * pixel_um ** 2) * 100,
            has_junction=bool(sk_px[k] >= 20), continuity=(ct_px[k] / sk_px[k]) if sk_px[k] >= 20 else np.nan,
            ridge_snr=float(Rm[k]), JC_ratio=(Jv / Cv) if (Cv is not np.nan and Cv and Cv > 0 and not np.isnan(Jv)) else np.nan))
    return rows

def quantify(name, cache, qc_dir=None, pixel_um=PIXEL_UM, cells=None):
    S, meta = load(name, cache)
    z0, z1 = best_window(S['ecad'])
    P = {c: S[c][z0:z1].mean(0) for c in ['actin', 'dnmt3a', 'ecad']}
    E = P['ecad']; Es = filters.gaussian(E, 1.0)
    tissue, gap = tissue_masks(P)
    has_gap = gap.sum() > 4000
    dist = ndi.distance_transform_edt(~gap) if has_gap else np.full(E.shape, 1e9)
    full = tissue.copy()
    edge = tissue & (dist <= ZONE_UM / pixel_um)
    bgv = float(np.median(E[gap])) if has_gap else float(np.percentile(E[tissue], 1))
    nlab = nuclei_mask(S['dapi'], tissue); nmask = nlab > 0
    satmask = S['ecad'][z0:z1].max(0) >= 4095
    Eb = np.clip(Es - bgv, 0, None)
    global LOW_THR
    if LOCALNORM:   # divide by local mean (sigma 12 px ~ 3 um) so dim junctions in dim regions are judged against their own surroundings
        loc = filters.gaussian(Eb, 12); fl = 0.3 * np.median(loc[full]) + 1e-6
        Ein = filters.gaussian(Eb / np.maximum(loc, fl), 0.5)
        R = sato(Ein, sigmas=[1.5, 2.5], black_ridges=False)
    else:
        R = sato(Eb, sigmas=[1.5, 2.5], black_ridges=False)
    rz = R[full]; med = np.median(rz); mad = np.median(np.abs(rz - med)) * 1.4826 + 1e-9
    thr = med + (HIGH_K if V2 else 3.0) * mad; LOW_THR = med + LOW_K * mad
    m_full, sk_full = zone_metrics(Es, E, Eb, R, thr, full, nmask, nlab, satmask, bgv, pixel_um, '')
    m_edge, _ = zone_metrics(Es, E, Eb, R, thr, edge, nmask, nlab, satmask, bgv, pixel_um, '_edge50um')
    out = dict(image=name, n_z=S['ecad'].shape[0], z_window=f"{z0+1}-{z1}", has_wound_gap=bool(has_gap),
               ecad_pmt_V=meta.get('ecad_pmt', [np.nan])[0], ecad_bg=bgv, ecad_sat_frac=float(satmask[full].mean()),
               ridge_thr=float(thr))
    out.update(m_full); out.update(m_edge)
    if cells is not None: cells.extend(cell_table(name, Es, R, med, mad, sk_full, nlab, tissue, gap, bgv, pixel_um))
    if qc_dir:
        os.makedirs(qc_dir, exist_ok=True); qc(name, E, tissue, gap, edge, nlab, sk_full, qc_dir)
    return out

def qc(name, E, tissue, gap, zone, nlab, sk, qc_dir):
    lo, hi = np.percentile(E[zone] if zone.any() else E, [1, 99.7])
    g = np.clip((E - lo) / (hi - lo + 1e-9), 0, 1)
    rgb = np.dstack([g, g, g])
    a = rgb.copy()
    gp = np.pad(gap, 1, constant_values=False)
    gb = segmentation.find_boundaries(gp, mode='inner')[1:-1, 1:-1]
    a[morphology.dilation(gb)] = [1, 1, 0]
    zp = np.pad(zone, 1, constant_values=False)
    zb = segmentation.find_boundaries(zp, mode='inner')[1:-1, 1:-1] & ~gb
    zb[[0, -1], :] = False; zb[:, [0, -1]] = False
    a[zb] = [0, 1, 1]
    b = rgb.copy(); b[segmentation.find_boundaries(nlab, mode='inner')] = [0, 0.6, 1]
    c = rgb.copy(); c[ndi.binary_dilation(sk)] = [0, 1, 0]
    panel = np.hstack([rgb, a, b, c])[::2, ::2]
    io.imsave(f"{qc_dir}/{name}.png", (panel * 255).astype(np.uint8), check_contrast=False)

def group_of(name):
    """treatment, time (min), batch label from the file name, e.g. 'CalA 10 Min WP E-Cad_0003' or '15 Min WP E-Cad Rep2_0004'"""
    treat = 'CalA' if name.startswith('CalA') else 'Control'
    t = int(re.search(r'(\d+) Min', name).group(1))
    batch = 'Rep2label' if 'Rep2' in name else 'main'
    return treat, t, batch

def _job(args):
    n, cache, out, pu = args
    cells = []
    r = quantify(n, cache, qc_dir=os.path.join(out, 'qc'), pixel_um=pu, cells=cells)
    r['treatment'], r['time_min'], r['file_batch'] = group_of(n); r['condition'] = f"{r['treatment']}_{r['time_min']}min"
    return r, cells

if __name__ == '__main__':
    from multiprocessing import Pool
    ap = argparse.ArgumentParser(); ap.add_argument('--cache', default='cache'); ap.add_argument('--out', default='results')
    ap.add_argument('--pixel-um', type=float, default=PIXEL_UM); ap.add_argument('--only', default=None); ap.add_argument('--procs', type=int, default=2)
    ap.add_argument('--cells', action='store_true', help='also write per_cell_metrics.csv')
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    names = sorted(os.path.basename(p)[:-4] for p in glob.glob(f"{a.cache}/*.npz"))
    if a.only: names = [n for n in names if re.search(a.only, n)]
    rows, allcells = [], []
    with Pool(a.procs) as pl:
        for (r, c), n in zip(pl.imap(_job, [(n, a.cache, a.out, a.pixel_um) for n in names]), names):
            rows.append(r); allcells += c; print('ok', n, len(c), 'cells', flush=True)
    pd.DataFrame(rows).to_csv(f"{a.out}/per_image_metrics.csv", index=False)
    if a.cells:
        cd = pd.DataFrame(allcells); cd['treatment'], cd['time_min'], _ = zip(*[group_of(i) for i in cd.image]); cd.to_csv(f"{a.out}/per_cell_metrics.csv", index=False)

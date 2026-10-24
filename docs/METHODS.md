# Methods

## Input
Olympus FV3000 `.oir` stacks, 800x800 px, 12-bit, 4-8 z-slices (0.4 um). Chunks of the OIR container are parsed directly (`oirread.py`): channels are identified by UUID (DAPI, actin, DNMT3a, E-cadherin), and PMT voltage/offset/gain, laser settings, zoom and z-step are read from the embedded XML.

## Per field
1. **Z handling.** The 3 consecutive slices with the sharpest E-cad (high-pass variance) are averaged, so every field has the same optical thickness.
2. **Tissue and wound gap.** Actin, DNMT3a and E-cad channels are normalised, summed and blurred, then Otsu-thresholded at 0.6x. The largest border-touching dark component is the wound gap. Background = median E-cad inside the gap (or the 1st percentile of tissue if no gap).
3. **Nuclei.** DAPI max projection, Otsu, opening, watershed.
4. **Junction tracing.** E-cad is background-subtracted and smoothed (sigma 1). A Sato ridge filter (sigma 1.5 and 2.5) gives a ridge image R. Noise is estimated per field as median and MAD of R over tissue. Junctions are traced by **hysteresis thresholding** (high = median + 3 MAD, low = median + 2 MAD), objects < 40 px removed, skeletonised, and skeleton segments < 30 px dropped. Hysteresis lets a dim, dotted stretch of junction be traced when it connects to a clearly detected one. (The older single-threshold version, `ECAD_MODE=global`, missed such stretches.)
5. **Metrics.** See README. Continuity uses segments >= 50 px (~13 um). `ridge_snr_p95` does not depend on the threshold. Each metric also has an `_edge50um` variant restricted to tissue within 50 um of the wound edge.

## Per cell
Tissue is divided into territories by nearest nucleus (not true cell boundaries). Only complete territories are kept: not touching the image border or the wound gap (dilated 4 px), area 0.3-3x the field median. A cell has a junction if >= 20 skeleton px fall in its territory.

## Statistics
* **Field level (`timecourse_stats.py`):** fields are the unit. Each treatment is compared with its own 0 min (Mann-Whitney, Cliff's delta, bootstrap CI of the median fold, Holm over time points); CalA vs control at each time; a within-time-stratified permutation test for the treatment x time interaction (mean fold at 5/10/15 min, CalA minus control); Spearman trend; equivalence margin +-25%; voltage sensitivity (log metric ~ time + voltage) and a matched-voltage subset (540-575 V).
* **Cell level (`cell_stats.py`):** the naive per-cell Mann-Whitney is shown only to illustrate pseudoreplication. The reported inference swaps whole fields between groups (5,000 permutations) and uses a cluster bootstrap for CIs; Holm correction within each comparison family.

## Known issues
* Detection thresholds (2 / 3 MAD) were chosen by visual QC on CalA 5 min fields and applied to all fields. Results at other settings were not systematically explored.
* In the control 5 min set, the first ~6 fields acquired (`5 Min WP E-Cad` to `_0005`) show nuclear-shaped E-cad signal and speckle with almost no traceable junctions, and junction density rises almost monotonically with acquisition order (Spearman 0.96). This looks like a focus/sample problem early in the session, not biology. Decide whether to exclude them before interpreting that time point.
* The Fiji macro is not validated and is not equivalent to the current Python pipeline.

// Fiji/ImageJ macro: simplified version of the Python E-cadherin junction pipeline.
// NOT validated in Fiji (written without Fiji access). Use to sanity-check the Python numbers on a few images.
// Requires Fiji (Bio-Formats and "Analyze Skeleton (2D/3D)" are included by default).
// CHECK the channel index: E-cad (647) is expected to be the LAST channel (4); the file order is DAPI, DNMT3a(561), actin(488), E-cad(647) by acquisition order.
ECAD_C = 4;  PIXEL_UM = 0.265;  NZ = 3;  K = 3.0;
dir = getDirectory("Choose folder with .oir files");
list = getFileList(dir);
run("Clear Results"); setBatchMode(true);
out = "Name,JunctionLen_um,Area_um2,Density_um_per_100um2,Segments,MeanSegment_um\n";
for (i = 0; i < list.length; i++) {
  if (!endsWith(list[i], ".oir")) continue;
  run("Bio-Formats", "open=[" + dir + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
  id = getImageID(); getDimensions(w, h, c, z, t);
  run("Duplicate...", "title=ec duplicate channels=" + ECAD_C);
  // best-focus 3-slice window (highest high-pass variance)
  best = 1; bs = -1;
  for (s = 1; s <= z - NZ + 1; s++) {
    run("Duplicate...", "title=tmp duplicate range=" + s + "-" + (s + NZ - 1));
    run("Z Project...", "projection=[Average Intensity]"); rename("p");
    run("Duplicate...", "title=a"); run("Gaussian Blur...", "sigma=1");
    run("Duplicate...", "title=b"); run("Gaussian Blur...", "sigma=6");
    imageCalculator("Subtract create 32-bit", "a", "b"); getStatistics(ar, m, mn, mx, sd);
    v = sd * sd / (m * m + 1e-6); close("Result of a"); close("a"); close("b"); close("p"); close("tmp");
    if (v > bs) { bs = v; best = s; }
  }
  selectImage("ec"); run("Z Project...", "start=" + best + " stop=" + (best + NZ - 1) + " projection=[Average Intensity]"); rename("E");
  run("32-bit"); run("Duplicate...", "title=E1"); run("Gaussian Blur...", "sigma=1.5");
  run("Duplicate...", "title=E2"); run("Gaussian Blur...", "sigma=4");
  imageCalculator("Subtract create 32-bit", "E1", "E2"); rename("DoG");   // band-pass ridge proxy
  // tissue mask (exclude wound gap): Otsu on heavily blurred image
  selectImage("E"); run("Duplicate...", "title=T"); run("Gaussian Blur...", "sigma=6"); setAutoThreshold("Otsu dark");
  run("Convert to Mask"); run("Fill Holes"); rename("tissue");
  selectImage("DoG"); run("Duplicate...", "title=D"); run("Select None");
  imageCalculator("Multiply create 32-bit", "D", "tissue"); // tissue*255 scale, only for stats
  selectImage("tissue"); run("Create Selection"); selectImage("DoG"); run("Restore Selection");
  getStatistics(area, mean, mn, mx, sd); thr = mean + K * sd; run("Select None");
  setThreshold(thr, 1e9); run("Convert to Mask"); rename("ridge");
  imageCalculator("AND create", "ridge", "tissue"); rename("rk");
  run("Skeletonize"); run("Analyze Skeleton (2D/3D)", "prune=none");
  nb = nResults; // each row = one skeleton
  tl = 0; ns = 0;
  for (r = 0; r < nb; r++) { tl += getResult("# Branches", r) * getResult("Average Branch Length", r); ns++; }
  selectImage("tissue"); run("Create Selection"); getStatistics(tissueArea);
  aum = tissueArea * PIXEL_UM * PIXEL_UM; lum = tl * PIXEL_UM;
  out += list[i] + "," + lum + "," + aum + "," + (lum / aum * 100) + "," + ns + "," + (lum / maxOf(ns, 1)) + "\n";
  run("Clear Results"); run("Close All");
}
File.saveString(out, dir + "Ecad_junction_quant_fiji.csv");
setBatchMode(false); print("Done -> " + dir + "Ecad_junction_quant_fiji.csv");

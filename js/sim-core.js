/* sim-core.js
 *
 * Statistical core for the Simple Interval Mapping (SIM) Explorer.
 *
 * Implements, for a recombinant inbred line (RIL, by repeated selfing)
 * population on a single chromosome:
 *
 *   1. Population simulation.  The RIL genotype along the chromosome is
 *      modelled as a Markov chain: between two loci at per-meiosis
 *      recombination fraction r (Haldane, no interference) the genotype
 *      switches parental origin with probability R = 2r / (1 + 2r), the
 *      Haldane-Waddington map expansion for RILs by selfing.
 *   2. QTL genotype weights.  At any genome position lambda, the
 *      probability that the (unobserved) QTL carries the Parent A allele,
 *      conditional on the nearest observed flanking markers, from the
 *      Lander-Botstein interval-mapping conditional probabilities.
 *   3. Haley-Knott regression scan.  The unobserved QTL genotype x_i in
 *      {-1, +1} is replaced by its conditional expectation
 *      xhat_i = 2*P(QTL = A | markers) - 1, and y is regressed on xhat
 *      at each grid position.  LOD = (n/2) * log10(RSS0 / RSS1).
 *   4. Churchill-Doerge permutation threshold for genome-wide
 *      significance.
 *
 * Pure functions only; no DOM access.  Usable in the browser as a classic
 * script (functions become globals) and in Node via module.exports.
 *
 * References:
 *   Lander ES, Botstein D (1989) Mapping Mendelian factors underlying
 *     quantitative traits using RFLP linkage maps. Genetics 121:185-199.
 *   Haley CS, Knott SA (1992) A simple regression method for mapping
 *     quantitative trait loci in line crosses using flanking markers.
 *     Heredity 69:315-324.
 *   Churchill GA, Doerge RW (1994) Empirical threshold values for
 *     quantitative trait mapping. Genetics 138:963-971.
 *   Haldane JBS, Waddington CH (1931) Inbreeding and linkage. Genetics
 *     16:357-374.
 */

"use strict";

/* ------------------------------------------------------------------ */
/* Random number generation                                            */
/* ------------------------------------------------------------------ */

/* Deterministic 32-bit PRNG (mulberry32).  Seed input keeps every
 * simulation reproducible. */
function mulberry32(seed) {
  let t = (seed >>> 0) || 1;
  return function () {
    t += 0x6D2B79F5;
    let z = Math.imul(t ^ (t >>> 15), t | 1);
    z ^= z + Math.imul(z ^ (z >>> 7), z | 61);
    return ((z ^ (z >>> 14)) >>> 0) / 4294967296;
  };
}

/* Standard normal draw via Box-Muller. */
function randn(rng) {
  let u = 0, v = 0;
  while (u === 0) u = rng();
  while (v === 0) v = rng();
  return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
}

/* ------------------------------------------------------------------ */
/* Map functions                                                       */
/* ------------------------------------------------------------------ */

/* Haldane map function: per-meiosis recombination fraction r for a
 * distance of d centiMorgans (assumes no interference). */
function haldaneR(dCM) {
  return 0.5 * (1 - Math.exp(-2 * dCM / 100));
}

/* RIL-by-selfing map expansion (Haldane & Waddington 1931): the
 * probability that a RIL switches parental genotype between two loci
 * at per-meiosis recombination fraction r. */
function rilSwitchR(r) {
  return (2 * r) / (1 + 2 * r);
}

/* ------------------------------------------------------------------ */
/* Population simulation                                               */
/* ------------------------------------------------------------------ */

/* Simulate n RILs on one chromosome.
 *
 * params: { n, chromLen, markerSpacing, qtlPos, effect, residSD,
 *           missingRate }
 *   effect    - additive QTL effect a: y = a * x_qtl + e, x in {-1,+1},
 *               so a is half the homozygote difference.
 *   residSD   - residual standard deviation sigma.
 *   missingRate - fraction of marker genotypes masked as missing.
 *
 * Returns { n, m, markerPos, geno, qtlGeno, y, params } where geno is a
 * Float64Array of length n*m holding -1 / +1 / NaN (missing), qtlGeno
 * holds the true (unobserved) QTL genotypes, and y the phenotypes. */
function simulateRIL(params, rng) {
  const n = params.n;
  const chromLen = params.chromLen;
  const spacing = params.markerSpacing;
  const qtlPos = params.qtlPos;
  const effect = params.effect;
  const residSD = params.residSD;
  const missingRate = params.missingRate || 0;

  /* Marker positions: evenly spaced, pinned to both chromosome ends. */
  const markerPos = [];
  for (let p = 0; p < chromLen - 1e-9; p += spacing) markerPos.push(p);
  markerPos.push(chromLen);
  const m = markerPos.length;

  /* Joint ordered position list (markers + QTL) so the QTL genotype is
   * simulated consistently with the flanking markers. */
  const uniq = [];
  const allPos = markerPos.concat([qtlPos]).sort((a, b) => a - b);
  for (const p of allPos) {
    if (uniq.length === 0 || Math.abs(p - uniq[uniq.length - 1]) > 1e-9) uniq.push(p);
  }
  const qtlIndex = uniq.findIndex((p) => Math.abs(p - qtlPos) < 1e-9);
  const markerIdx = markerPos.map((p) =>
    uniq.findIndex((q) => Math.abs(q - p) < 1e-9)
  );

  const geno = new Float64Array(n * m);
  const qtlGeno = new Int8Array(n);
  const y = new Float64Array(n);
  const chain = new Int8Array(uniq.length);

  for (let i = 0; i < n; i++) {
    /* Markov chain along the chromosome. */
    let g = rng() < 0.5 ? 1 : -1;
    chain[0] = g;
    for (let k = 1; k < uniq.length; k++) {
      const r = haldaneR(uniq[k] - uniq[k - 1]);
      if (rng() < rilSwitchR(r)) g = -g;
      chain[k] = g;
    }
    qtlGeno[i] = chain[qtlIndex];
    for (let j = 0; j < m; j++) {
      let og = chain[markerIdx[j]];
      geno[i * m + j] = missingRate > 0 && rng() < missingRate ? NaN : og;
    }
    y[i] = effect * qtlGeno[i] + residSD * randn(rng);
  }
  return { n, m, markerPos, geno, qtlGeno, y, params };
}

/* ------------------------------------------------------------------ */
/* QTL genotype weights: P(QTL = A | flanking marker data)             */
/* ------------------------------------------------------------------ */

/* Conditional probability that the QTL at position lambda carries the
 * Parent A allele, given the nearest observed flanking markers.
 *
 * markerPos - array of marker positions (cM), ascending.
 * geno      - array-like of length m with -1 / +1 / NaN (missing).
 * lambda    - genome position (cM) at which to evaluate.
 *
 * With flanking markers M1 (genotype g1) at distance d1 and M2 (g2) at
 * d2, and RIL switch probabilities R1, R2 (R12 = R1 + R2 - 2*R1*R2):
 *   P(A | A,A) = (1-R1)(1-R2) / (1-R12)
 *   P(A | B,B) = 1 - (1-R1)(1-R2) / (1-R12)
 *   P(A | A,B) = (1-R1) R2 / R12
 *   P(A | B,A) = R1 (1-R2) / R12
 * At an observed marker the probability is 0/1; outside the observed
 * marker range (or with no observed markers) it is the marginal 0.5. */
function qtlWeight(markerPos, geno, lambda) {
  const m = markerPos.length;
  let li = -1, ri = -1;
  for (let j = 0; j < m; j++) {
    const gj = geno[j];
    if (Number.isNaN(gj)) continue;
    if (Math.abs(markerPos[j] - lambda) < 1e-9) return gj > 0 ? 1 : 0;
    if (markerPos[j] < lambda) li = j;
    else { ri = j; break; }
  }
  if (li < 0 || ri < 0) return 0.5;

  const g1 = geno[li], g2 = geno[ri];
  const R1 = rilSwitchR(haldaneR(lambda - markerPos[li]));
  const R2 = rilSwitchR(haldaneR(markerPos[ri] - lambda));
  const R12 = R1 + R2 - 2 * R1 * R2;
  const a1 = 1 - R1, a2 = 1 - R2;

  if (g1 > 0 && g2 > 0) return (a1 * a2) / (1 - R12);
  if (g1 < 0 && g2 < 0) return 1 - (a1 * a2) / (1 - R12);
  if (g1 > 0) return (a1 * R2) / R12; /* g1 = A, g2 = B */
  return (R1 * a2) / R12;            /* g1 = B, g2 = A */
}

/* Expected QTL genotypes xhat = 2w - 1 for every RIL at every grid
 * position, plus the per-position mean information content
 * mean_i |xhat_i| (1 = fully determined, 0 = uninformative).
 * Grid must be ascending. */
function expectedGenotypes(sim, grid) {
  const n = sim.n, m = sim.m, G = grid.length;
  const X = new Float64Array(n * G);
  const info = new Float64Array(G);
  const obsPos = [], obsGeno = [];
  for (let i = 0; i < n; i++) {
    obsPos.length = 0; obsGeno.length = 0;
    for (let j = 0; j < m; j++) {
      const gj = sim.geno[i * m + j];
      if (!Number.isNaN(gj)) { obsPos.push(sim.markerPos[j]); obsGeno.push(gj); }
    }
    const w = weightTraceSweep(obsPos, obsGeno, grid);
    for (let g = 0; g < G; g++) {
      const xh = 2 * w[g] - 1;
      X[i * G + g] = xh;
      info[g] += Math.abs(xh);
    }
  }
  for (let g = 0; g < G; g++) info[g] /= n;
  return { X, info };
}

/* Weight trace for a single RIL via an ascending sweep over the grid
 * (O(K + G) for K observed markers); grid must be ascending. */
function weightTraceSweep(obsPos, obsGeno, grid) {
  const w = new Float64Array(grid.length);
  const K = obsPos.length;
  if (K === 0) { w.fill(0.5); return w; }
  let k = 0;
  for (let g = 0; g < grid.length; g++) {
    const lam = grid[g];
    while (k < K && obsPos[k] < lam - 1e-9) k++;
    if (k < K && Math.abs(obsPos[k] - lam) < 1e-9) {
      w[g] = obsGeno[k] > 0 ? 1 : 0;
      continue;
    }
    const li = k - 1, ri = k;
    if (li < 0 || ri >= K) { w[g] = 0.5; continue; }
    const g1 = obsGeno[li], g2 = obsGeno[ri];
    const R1 = rilSwitchR(haldaneR(lam - obsPos[li]));
    const R2 = rilSwitchR(haldaneR(obsPos[ri] - lam));
    const R12 = R1 + R2 - 2 * R1 * R2;
    const a1 = 1 - R1, a2 = 1 - R2;
    if (g1 > 0 && g2 > 0) w[g] = (a1 * a2) / (1 - R12);
    else if (g1 < 0 && g2 < 0) w[g] = 1 - (a1 * a2) / (1 - R12);
    else if (g1 > 0) w[g] = (a1 * R2) / R12;
    else w[g] = (R1 * a2) / R12;
  }
  return w;
}

/* Weight trace for a single RIL (used by the weight-curve panel). */
function weightTrace(sim, rilIndex, grid) {
  const obsPos = [], obsGeno = [];
  for (let j = 0; j < sim.m; j++) {
    const gj = sim.geno[rilIndex * sim.m + j];
    if (!Number.isNaN(gj)) { obsPos.push(sim.markerPos[j]); obsGeno.push(gj); }
  }
  return weightTraceSweep(obsPos, obsGeno, grid);
}

/* ------------------------------------------------------------------ */
/* Haley-Knott interval-mapping scan                                   */
/* ------------------------------------------------------------------ */

/* Regress y on the expected QTL genotype at each grid position.
 * Returns per-position { lod, beta, se, r2 }.
 *   LOD = (n/2) * log10(RSS0 / RSS1), RSS0 from the intercept-only model. */
function haleyKnott(y, X, n, G) {
  const lod = new Float64Array(G);
  const beta = new Float64Array(G);
  const se = new Float64Array(G);
  const r2 = new Float64Array(G);

  let ybar = 0;
  for (let i = 0; i < n; i++) ybar += y[i];
  ybar /= n;
  let Syy = 0;
  const yc = new Float64Array(n);
  for (let i = 0; i < n; i++) { yc[i] = y[i] - ybar; Syy += yc[i] * yc[i]; }

  for (let g = 0; g < G; g++) {
    let sx = 0, sxx = 0, sxy = 0;
    for (let i = 0; i < n; i++) {
      const x = X[i * G + g];
      sx += x; sxx += x * x; sxy += x * yc[i];
    }
    const Sxx = sxx - (sx * sx) / n;
    if (Sxx < 1e-12 || Syy < 1e-12) {
      lod[g] = 0; beta[g] = 0; se[g] = NaN; r2[g] = 0;
      continue;
    }
    const b = sxy / Sxx;
    const rss1 = Syy - (sxy * sxy) / Sxx;
    lod[g] = (n / 2) * Math.log10(Syy / Math.max(rss1, 1e-300));
    beta[g] = b;
    se[g] = Math.sqrt(rss1 / (n - 2) / Sxx);
    r2[g] = 1 - rss1 / Syy;
  }
  return { lod, beta, se, r2 };
}

/* Index of the maximum LOD (first occurrence on ties). */
function peakIndex(lod) {
  let best = 0;
  for (let g = 1; g < lod.length; g++) if (lod[g] > lod[best]) best = g;
  return best;
}

/* Churchill-Doerge permutation threshold: permute phenotypes, rescan,
 * record the genome-wide maximum LOD each time, return the 95th
 * percentile.  The expected-genotype matrix X is phenotype-free, so it
 * is computed once and reused across permutations. */
function permutationThreshold(y, X, n, G, nPerm, rng, onProgress) {
  const Sx = new Float64Array(G);
  const Sxx = new Float64Array(G);
  for (let g = 0; g < G; g++) {
    let s = 0, ss = 0;
    for (let i = 0; i < n; i++) { const x = X[i * G + g]; s += x; ss += x * x; }
    Sx[g] = s;
    Sxx[g] = ss - (s * s) / n;
  }
  const idx = new Int32Array(n);
  for (let i = 0; i < n; i++) idx[i] = i;
  const yp = new Float64Array(n);
  const maxLods = new Float64Array(nPerm);

  for (let p = 0; p < nPerm; p++) {
    for (let i = n - 1; i > 0; i--) {
      const k = Math.floor(rng() * (i + 1));
      const t = idx[i]; idx[i] = idx[k]; idx[k] = t;
    }
    let ybar = 0;
    for (let i = 0; i < n; i++) ybar += y[idx[i]];
    ybar /= n;
    let Syy = 0;
    for (let i = 0; i < n; i++) {
      const d = y[idx[i]] - ybar;
      yp[i] = d; Syy += d * d;
    }
    let mx = 0;
    for (let g = 0; g < G; g++) {
      if (Sxx[g] < 1e-12 || Syy < 1e-12) continue;
      let sxy = 0;
      for (let i = 0; i < n; i++) sxy += X[i * G + g] * yp[i];
      const rss1 = Syy - (sxy * sxy) / Sxx[g];
      const l = (n / 2) * Math.log10(Syy / Math.max(rss1, 1e-300));
      if (l > mx) mx = l;
    }
    maxLods[p] = mx;
    if (onProgress && p % 25 === 0) onProgress(p / nPerm);
  }
  const sorted = Array.from(maxLods).sort((a, b) => a - b);
  return {
    threshold: sorted[Math.min(nPerm - 1, Math.floor(0.95 * nPerm))],
    maxLods: sorted,
    nPerm,
  };
}

/* Evenly spaced scan grid (1 cM steps) over [0, chromLen]. */
function makeGrid(chromLen, step) {
  step = step || 1;
  const grid = [];
  for (let p = 0; p <= chromLen + 1e-9; p += step) {
    grid.push(Math.min(p, chromLen));
  }
  return grid;
}

/* Node export shim; in the browser the functions are globals. */
if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    mulberry32, randn, haldaneR, rilSwitchR,
    simulateRIL, qtlWeight, weightTraceSweep, expectedGenotypes, weightTrace,
    haleyKnott, peakIndex, permutationThreshold, makeGrid,
  };
}

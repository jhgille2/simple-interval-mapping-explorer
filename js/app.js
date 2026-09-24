/* app.js — UI wiring and canvas rendering for the Simple Interval Mapping Explorer.
 * Depends on js/sim-core.js (loaded first). No external dependencies. */

"use strict";

/* ------------------------------------------------------------------ */
/* State                                                               */
/* ------------------------------------------------------------------ */

const MAXD = 160; /* max RIL rows / traces drawn */

const state = {
  params: null,
  sim: null,
  grid: null,
  X: null,
  info: null,
  scan: null,
  perm: null,
  selectedRil: 0,
};

const $ = (id) => document.getElementById(id);

/* ------------------------------------------------------------------ */
/* Canvas helpers                                                      */
/* ------------------------------------------------------------------ */

function setupCanvas(id, cssHeight) {
  const cv = $(id);
  const dpr = window.devicePixelRatio || 1;
  const wCss = Math.max(60, cv.clientWidth);
  cv.width = Math.round(wCss * dpr);
  cv.height = Math.round(cssHeight * dpr);
  cv.style.height = cssHeight + "px";
  const ctx = cv.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, wCss, cssHeight);
  return { ctx, w: wCss, h: cssHeight };
}

/* Draw axes; returns plot area and data->pixel mappers. */
function drawAxes(ctx, w, h, m, x0, x1, y0, y1, xticks, yticks, xlabel, ylabel) {
  const pw = w - m.l - m.r, ph = h - m.t - m.b;
  const X = (v) => m.l + ((v - x0) / (x1 - x0)) * pw;
  const Y = (v) => m.t + ph - ((v - y0) / (y1 - y0)) * ph;
  ctx.strokeStyle = "#9aa7b1";
  ctx.fillStyle = "#5a6b76";
  ctx.lineWidth = 1;
  ctx.font = "11px Helvetica, Arial, sans-serif";
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  /* x ticks */
  for (const t of xticks) {
    const x = X(t);
    ctx.beginPath(); ctx.moveTo(x, m.t + ph); ctx.lineTo(x, m.t + ph + 5); ctx.stroke();
    ctx.fillText(String(t), x, m.t + ph + 7);
    ctx.strokeStyle = "#e6ebef";
    ctx.beginPath(); ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + ph); ctx.stroke();
    ctx.strokeStyle = "#9aa7b1";
  }
  /* y ticks */
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  for (const t of yticks) {
    const y = Y(t);
    ctx.beginPath(); ctx.moveTo(m.l - 5, y); ctx.lineTo(m.l, y); ctx.stroke();
    ctx.fillText(String(t), m.l - 8, y);
    ctx.strokeStyle = "#e6ebef";
    ctx.beginPath(); ctx.moveTo(m.l, y); ctx.lineTo(m.l + pw, y); ctx.stroke();
    ctx.strokeStyle = "#9aa7b1";
  }
  /* frame */
  ctx.strokeStyle = "#9aa7b1";
  ctx.strokeRect(m.l, m.t, pw, ph);
  /* labels */
  ctx.fillStyle = "#5a6b76";
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  if (xlabel) ctx.fillText(xlabel, m.l + pw / 2, h - 14);
  ctx.save();
  ctx.translate(12, m.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = "bottom";
  if (ylabel) ctx.fillText(ylabel, 0, 0);
  ctx.restore();
  return { X, Y, pw, ph };
}

function vline(ctx, X, Y, y0, y1, x, style, dash) {
  ctx.save();
  ctx.strokeStyle = style;
  ctx.lineWidth = 1.5;
  if (dash) ctx.setLineDash(dash);
  ctx.beginPath();
  ctx.moveTo(X(x), Y(y0));
  ctx.lineTo(X(x), Y(y1));
  ctx.stroke();
  ctx.restore();
}

function niceTicks(x0, x1, target) {
  const span = x1 - x0;
  const raw = span / target;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  const step = (norm >= 5 ? 10 : norm >= 2 ? 5 : norm >= 1 ? 2 : 1) * mag;
  const ticks = [];
  for (let t = Math.ceil(x0 / step) * step; t <= x1 + 1e-9; t += step) {
    ticks.push(Math.abs(t) < 1e-9 ? 0 : Math.round(t * 100) / 100);
  }
  return ticks;
}

const fmt = (x, d) => (Number.isFinite(x) ? x.toFixed(d) : "n/a");

/* ------------------------------------------------------------------ */
/* Rendering                                                           */
/* ------------------------------------------------------------------ */

const COL_A = "#1f6fb2", COL_B = "#d95f02", COL_MISS = "#cfd6dc";

function renderHap() {
  const sim = state.sim, p = state.params;
  const R = Math.min(MAXD, sim.n);
  const rowH = 4;
  const cssH = R * rowH + 64;
  const { ctx, w, h } = setupCanvas("cv-hap", cssH);
  const m = { l: 46, r: 14, t: 10, b: 40 };
  const { X } = drawAxes(ctx, w, h, m, 0, p.chromLen, 0, 1, niceTicks(0, p.chromLen, 10), [], "Position (cM)", null);

  const hw = p.markerSpacing / 2;
  for (let r = 0; r < R; r++) {
    const yTop = m.t + r * rowH;
    for (let j = 0; j < sim.m; j++) {
      const g = sim.geno[r * sim.m + j];
      const pos = sim.markerPos[j];
      const x0 = X(Math.max(0, pos - hw)), x1 = X(Math.min(p.chromLen, pos + hw));
      ctx.fillStyle = Number.isNaN(g) ? COL_MISS : g > 0 ? COL_A : COL_B;
      ctx.fillRect(x0, yTop, Math.max(1, x1 - x0), rowH - 0.6);
    }
  }
  /* selected row highlight */
  if (state.selectedRil < R) {
    ctx.strokeStyle = "#111";
    ctx.lineWidth = 1.5;
    ctx.strokeRect(m.l, m.t + state.selectedRil * rowH - 1, w - m.l - m.r, rowH + 1);
  }
  vline(ctx, X, (v) => m.t + v * (R * rowH), 0, 1, p.qtlPos, "#333", [5, 4]);
  ctx.fillStyle = "#5a6b76";
  ctx.font = "11px Helvetica, Arial, sans-serif";
  ctx.textAlign = "left";
  ctx.fillText(sim.n > R ? `first ${R} of ${sim.n} RILs` : `${sim.n} RILs`, m.l, h - 26);
}

function renderWeights() {
  const sim = state.sim, p = state.params, grid = state.grid;
  const { ctx, w, h } = setupCanvas("cv-w", 300);
  const m = { l: 46, r: 14, t: 12, b: 40 };
  const { X, Y } = drawAxes(ctx, w, h, m, 0, p.chromLen, 0, 1,
    niceTicks(0, p.chromLen, 10), [0, 0.25, 0.5, 0.75, 1],
    "Position (cM)", "P(QTL = A | markers)");

  const R = Math.min(MAXD, sim.n);
  const showAll = $("chk-traces").checked;

  const traceOf = (ril) => {
    const obsPos = [], obsGeno = [];
    for (let j = 0; j < sim.m; j++) {
      const g = sim.geno[ril * sim.m + j];
      if (!Number.isNaN(g)) { obsPos.push(sim.markerPos[j]); obsGeno.push(g); }
    }
    return weightTraceSweep(obsPos, obsGeno, grid);
  };

  if (showAll) {
    ctx.strokeStyle = "rgba(31,111,178,0.07)";
    ctx.lineWidth = 1;
    for (let r = 0; r < R; r++) {
      if (r === state.selectedRil) continue;
      const wt = traceOf(r);
      ctx.beginPath();
      for (let g = 0; g < grid.length; g++) {
        const x = X(grid[g]), y = Y(wt[g]);
        if (g === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }
  }
  /* mean information content */
  if ($("chk-info").checked) {
    ctx.save();
    ctx.strokeStyle = COL_B;
    ctx.lineWidth = 1.5;
    ctx.setLineDash([6, 4]);
    ctx.beginPath();
    for (let g = 0; g < grid.length; g++) {
      const x = X(grid[g]), y = Y(state.info[g]);
      if (g === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.restore();
  }
  /* selected RIL */
  const wSel = traceOf(state.selectedRil);
  ctx.strokeStyle = COL_A;
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let g = 0; g < grid.length; g++) {
    const x = X(grid[g]), y = Y(wSel[g]);
    if (g === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.stroke();
  /* truth */
  vline(ctx, X, Y, 0, 1, p.qtlPos, "#333", [5, 4]);
  const qg = sim.qtlGeno[state.selectedRil] > 0 ? 1 : 0;
  ctx.fillStyle = "#333";
  ctx.beginPath();
  ctx.arc(X(p.qtlPos), Y(qg), 4, 0, 2 * Math.PI);
  ctx.fill();
  ctx.fillStyle = "#5a6b76";
  ctx.font = "11px Helvetica, Arial, sans-serif";
  ctx.textAlign = "left";
  ctx.fillText(`RIL ${state.selectedRil} (true QTL genotype: ${qg > 0 ? "AA" : "BB"})`, m.l, h - 26);
}

function renderLOD() {
  const sim = state.sim, p = state.params, grid = state.grid, scan = state.scan;
  const { ctx, w, h } = setupCanvas("cv-lod", 300);
  const pk = peakIndex(scan.lod);
  const thr = state.perm ? state.perm.threshold : 3;
  const yMax = Math.max(scan.lod[pk], thr, 3.2) * 1.15;
  const m = { l: 46, r: 14, t: 12, b: 40 };
  const { X, Y } = drawAxes(ctx, w, h, m, 0, p.chromLen, 0, yMax,
    niceTicks(0, p.chromLen, 10), niceTicks(0, yMax, 5),
    "Position (cM)", "LOD");

  /* marker rug */
  ctx.strokeStyle = "#7d8b96";
  ctx.lineWidth = 1;
  for (const pos of sim.markerPos) {
    const x = X(pos);
    ctx.beginPath(); ctx.moveTo(x, Y(0)); ctx.lineTo(x, Y(0) - 6); ctx.stroke();
  }
  /* LOD curve */
  ctx.strokeStyle = "#1b7a3d";
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let g = 0; g < grid.length; g++) {
    const x = X(grid[g]), y = Y(scan.lod[g]);
    if (g === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.stroke();
  /* truth + threshold + peak */
  vline(ctx, X, Y, 0, yMax, p.qtlPos, "#333", [5, 4]);
  ctx.save();
  ctx.strokeStyle = state.perm ? "#c0392b" : "#999999";
  ctx.setLineDash(state.perm ? [6, 4] : [2, 3]);
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(m.l, Y(thr)); ctx.lineTo(w - m.r, Y(thr));
  ctx.stroke();
  ctx.restore();
  ctx.fillStyle = "#c0392b";
  ctx.beginPath();
  ctx.arc(X(grid[pk]), Y(scan.lod[pk]), 4, 0, 2 * Math.PI);
  ctx.fill();
  ctx.fillStyle = "#1e2a32";
  ctx.font = "11px Helvetica, Arial, sans-serif";
  ctx.textAlign = "left";
  const label = `peak ${fmt(grid[pk], 0)} cM, LOD ${fmt(scan.lod[pk], 2)}`;
  let lx = X(grid[pk]) + 8;
  if (lx + ctx.measureText(label).width > w - m.r) lx = X(grid[pk]) - 8 - ctx.measureText(label).width;
  ctx.fillText(label, lx, Math.max(m.t + 10, Y(scan.lod[pk]) - 10));
  ctx.fillStyle = "#5a6b76";
  ctx.fillText(state.perm
    ? `red dashed: 5% permutation threshold (${state.perm.nPerm} perms)`
    : "grey dotted: LOD = 3 reference (no permutation threshold computed)",
    m.l, h - 26);
}

function renderSummary() {
  const p = state.params, scan = state.scan, grid = state.grid;
  const pk = peakIndex(scan.lod);
  const rows = [
    ["RILs (n)", p.n],
    ["Markers", `${state.sim.m} (spacing ${p.markerSpacing} cM)`],
    ["True QTL position", `${p.qtlPos} cM`],
    ["True additive effect a", fmt(p.effect, 3)],
    ["Peak position", `${fmt(grid[pk], 0)} cM`],
    ["Peak LOD", fmt(scan.lod[pk], 2)],
    ["Significance threshold", state.perm ? `${fmt(state.perm.threshold, 2)} (${state.perm.nPerm} permutations)` : "not computed"],
    ["Estimated effect at peak", `${fmt(scan.beta[pk], 3)} ± ${fmt(scan.se[pk], 3)}`],
    ["R² at peak", fmt(scan.r2[pk], 3)],
    ["Mean information at peak", fmt(state.info[pk], 3)],
  ];
  $("sum-body").innerHTML = rows
    .map(([k, v]) => `<tr><td>${k}</td><td>${v}</td></tr>`)
    .join("");
  $("marker-count-hint").textContent =
    `${state.sim.m} markers at ${p.markerSpacing} cM spacing on a ${p.chromLen} cM chromosome.`;
}

function renderAll() {
  renderHap();
  renderWeights();
  renderLOD();
  renderSummary();
}

/* ------------------------------------------------------------------ */
/* Simulation run                                                      */
/* ------------------------------------------------------------------ */

function readParams() {
  return {
    n: parseInt($("ctl-n").value, 10),
    chromLen: parseInt($("ctl-len").value, 10),
    markerSpacing: parseInt($("ctl-spacing").value, 10),
    qtlPos: parseInt($("ctl-qtlpos").value, 10),
    effect: parseFloat($("ctl-effect").value),
    residSD: parseFloat($("ctl-sigma").value),
    missingRate: parseFloat($("ctl-missing").value),
    seed: parseInt($("ctl-seed").value, 10) || 1,
  };
}

function run() {
  const p = readParams();
  const rng = mulberry32(p.seed);
  state.params = p;
  state.sim = simulateRIL(p, rng);
  state.grid = makeGrid(p.chromLen, 1);
  const eg = expectedGenotypes(state.sim, state.grid);
  state.X = eg.X;
  state.info = eg.info;
  state.scan = haleyKnott(state.sim.y, state.X, p.n, state.grid.length);
  state.perm = null;
  setPermStatus("No permutation threshold computed; dashed line shows LOD = 3.");
  state.selectedRil = Math.max(0, Math.min(p.n - 1, state.selectedRil));
  $("ctl-ril").max = p.n - 1;
  $("ctl-ril").value = state.selectedRil;
  renderAll();
}

/* ------------------------------------------------------------------ */
/* Permutation threshold                                               */
/* ------------------------------------------------------------------ */

function setPermStatus(text) {
  $("perm-status").textContent = text;
}

function runPermutation() {
  const btn = $("btn-perm");
  btn.disabled = true;
  setPermStatus("Computing permutations…");
  /* let the UI update before the blocking computation */
  setTimeout(() => {
    try {
      const p = state.params, G = state.grid.length;
      const requested = parseInt($("sel-perms").value, 10);
      /* bound total work to ~1.5e8 inner iterations */
      const eff = Math.max(50, Math.min(requested, Math.floor(1.5e8 / (p.n * G))));
      const res = permutationThreshold(state.sim.y, state.X, p.n, G, eff,
        mulberry32(p.seed * 100003 + 17));
      state.perm = res;
      setPermStatus(`95th percentile of genome-wide max LOD over ${eff} permutations: ` +
        `LOD = ${fmt(res.threshold, 2)}.` +
        (eff < requested ? ` (Reduced from ${requested} to bound runtime.)` : ""));
      renderLOD();
      renderSummary();
    } finally {
      btn.disabled = false;
    }
  }, 40);
}

/* ------------------------------------------------------------------ */
/* Control wiring                                                      */
/* ------------------------------------------------------------------ */

function bindRange(id, outId, format, onChange) {
  const el = $(id), out = $(outId);
  const update = () => { out.textContent = format(parseFloat(el.value)); };
  el.addEventListener("input", () => { update(); if (onChange) onChange(); scheduleRun(); });
  update();
}

let debounceTimer = null;
function scheduleRun() {
  clearTimeout(debounceTimer);
  debounceTimer = setTimeout(run, 160);
}

function updateVarExp() {
  const a = parseFloat($("ctl-effect").value);
  const s = parseFloat($("ctl-sigma").value);
  $("val-varexp").textContent = (100 * a * a / (a * a + s * s)).toFixed(1) + "%";
}

function init() {
  bindRange("ctl-n", "val-n", (v) => v.toFixed(0));
  bindRange("ctl-len", "val-len", (v) => v.toFixed(0), () => {
    const L = parseInt($("ctl-len").value, 10);
    const qp = $("ctl-qtlpos");
    qp.max = L;
    if (parseInt(qp.value, 10) > L) qp.value = L;
    $("val-qtlpos").textContent = qp.value;
  });
  bindRange("ctl-spacing", "val-spacing", (v) => v.toFixed(0));
  bindRange("ctl-missing", "val-missing", (v) => v.toFixed(2));
  bindRange("ctl-qtlpos", "val-qtlpos", (v) => v.toFixed(0));
  bindRange("ctl-effect", "val-effect", (v) => v.toFixed(2), updateVarExp);
  bindRange("ctl-sigma", "val-sigma", (v) => v.toFixed(2), updateVarExp);
  updateVarExp();

  $("ctl-seed").addEventListener("change", () => {
    $("val-seed").textContent = $("ctl-seed").value;
    run();
  });
  $("ctl-ril").addEventListener("change", () => {
    const p = state.params;
    state.selectedRil = Math.max(0, Math.min(p.n - 1, parseInt($("ctl-ril").value, 10) || 0));
    $("ctl-ril").value = state.selectedRil;
    renderHap();
    renderWeights();
  });
  $("chk-traces").addEventListener("change", renderWeights);
  $("chk-info").addEventListener("change", renderWeights);
  $("btn-run").addEventListener("click", run);
  $("btn-perm").addEventListener("click", runPermutation);

  /* click a haplotype row to select that RIL */
  $("cv-hap").addEventListener("click", (ev) => {
    const cv = $("cv-hap");
    const rect = cv.getBoundingClientRect();
    const y = ev.clientY - rect.top;
    const row = Math.floor((y - 10) / 4);
    const R = Math.min(MAXD, state.sim.n);
    if (row >= 0 && row < R) {
      state.selectedRil = row;
      $("ctl-ril").value = row;
      renderHap();
      renderWeights();
    }
  });

  let resizeTimer = null;
  window.addEventListener("resize", () => {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(() => { if (state.sim) renderAll(); }, 200);
  });

  run();
}

init();

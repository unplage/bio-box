"use strict";
const AA = "ARNDCQEGHILKMFPSTWYVBZX*";
const AA_IDX = {};
for (let i = 0; i < AA.length; i++) AA_IDX[AA[i]] = i;

const BLOSUM62 = [
[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4],
[-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4],
[-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4],
[-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4],
[0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
[-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4],
[-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
[0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4],
[-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4],
[-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4],
[-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4],
[-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4],
[-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4],
[-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4],
[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
[1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4],
[0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4],
[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4],
[-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4],
[0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4],
[-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4],
[-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
[0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4],
[-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1]
];

const PAM250 = [
[2,-2,0,0,-3,1,0,0,-2,-1,-1,-1,-2,-1,1,1,1,-6,-3,0,0,0,0,-8],
[-2,6,0,-1,-4,3,1,-2,0,-3,-2,2,-1,-3,0,0,1,-2,-4,-1,0,1,0,-8],
[0,0,2,2,-3,0,1,1,2,-2,-1,2,0,-2,1,0,1,-2,-1,1,1,1,0,-8],
[0,-1,2,5,-4,2,3,0,1,-3,-2,3,1,-2,0,0,1,-5,-3,0,3,3,0,-8],
[-3,-4,-3,-4,10,-4,-5,-4,-3,-3,-2,-4,-2,-3,0,0,-6,-4,-3,-3,-4,-4,0,-8],
[1,3,0,2,-4,4,2,-1,1,-2,-1,2,0,-2,0,0,0,-3,-3,0,2,3,0,-8],
[0,1,1,3,-5,2,4,-1,1,-3,-2,3,0,-2,0,0,0,-4,-3,1,3,4,0,-8],
[0,-2,1,0,-4,-1,-1,5,-2,-3,-2,1,-1,-1,1,0,0,-4,-4,-1,1,1,0,-8],
[-2,0,2,1,-3,1,1,-2,6,-3,-2,1,-2,-2,0,0,1,-4,-3,1,1,1,0,-8],
[-1,-3,-2,-3,-3,-2,-3,-3,-3,4,2,-3,2,1,-3,-2,-1,-5,-3,3,-2,-3,0,-8],
[-1,-2,-1,-2,-2,-1,-2,-2,-3,2,3,-2,3,2,-3,-1,-1,-4,-3,2,-1,-2,0,-8],
[-1,2,2,3,-4,2,3,1,1,-3,-2,4,1,-2,0,0,0,-4,-3,1,2,3,0,-8],
[-2,-1,0,1,-2,-2,0,-1,-2,2,3,-1,3,2,-2,-1,-1,-3,-3,2,0,-1,0,-8],
[-1,-3,-2,-2,-3,-2,-2,-1,-2,1,2,-2,2,5,-3,-2,-1,-4,-3,1,-2,-2,0,-8],
[1,0,1,0,0,0,0,1,0,-3,-3,0,-2,-3,6,1,0,-7,-4,0,0,0,0,-8],
[1,0,0,0,0,0,0,0,0,-2,-1,0,-1,-2,1,2,1,-3,-3,0,0,0,0,-8],
[1,1,1,1,-1,0,0,0,1,-1,-1,0,-1,-1,0,1,4,-5,-3,0,1,1,0,-8],
[-6,-2,-2,-5,-4,-3,-4,-4,-4,-5,-4,-4,-3,-4,-7,-3,-5,14,3,-5,-4,-4,0,-8],
[-3,-4,-1,-3,-3,-3,-3,-4,-3,-3,-3,-3,-3,3,-4,-3,-3,3,8,-3,-2,-3,0,-8],
[0,-1,1,0,-3,0,1,-1,1,3,2,1,2,1,0,0,0,-5,-3,3,1,1,0,-8],
[0,0,1,3,-4,0,3,1,1,-2,-1,2,0,-2,0,0,1,-4,-2,1,3,3,0,-8],
[0,1,0,3,-4,3,4,1,1,-3,-2,3,-1,-2,0,0,1,-4,-3,1,4,4,0,-8],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8],
[-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,1]
];

function makeScorer(params) {
  const mismatch = Number(params.mismatch);
  if (params.type === "aa") {
    const M = params.matrix === "pam250" ? PAM250 : BLOSUM62;
    return (x, y) => {
      const i = AA_IDX[x], j = AA_IDX[y];
      if (i === undefined || j === undefined) return mismatch;
      return M[i][j];
    };
  }
  const match = Number(params.match);
  return (x, y) => (x === y ? match : mismatch);
}

const NEG = -1e18;

function nwAlign(a, b, sf, go, ge, freeEnds) {
  const n = a.length, m = b.length, W = m + 1;
  const size = (n + 1) * W;
  const M = new Float64Array(size);
  const X = new Float64Array(size);
  const Y = new Float64Array(size);
  const pM = new Int8Array(size);
  const pX = new Int8Array(size);
  const pY = new Int8Array(size);

  M[0] = 0;
  for (let i = 1; i <= n; i++) {
    const k = i * W;
    if (freeEnds) {
      M[k] = 0; X[k] = 0; Y[k] = NEG;
      pM[k] = 0; pX[k] = 0;
    } else {
      M[k] = NEG;
      X[k] = -(go + (i - 1) * ge);
      pX[k] = 1;
      Y[k] = NEG;
    }
  }
  for (let j = 1; j <= m; j++) {
    if (freeEnds) {
      M[j] = 0; Y[j] = 0; X[j] = NEG;
      pM[j] = 0; pY[j] = 0;
    } else {
      M[j] = NEG;
      Y[j] = -(go + (j - 1) * ge);
      pY[j] = 2;
      X[j] = NEG;
    }
  }
  X[0] = NEG; Y[0] = NEG;

  for (let i = 1; i <= n; i++) {
    const ai = a[i - 1];
    const row = i * W, prow = (i - 1) * W;
    for (let j = 1; j <= m; j++) {
      const idx = row + j;
      const prev = prow + j - 1;
      let bs = M[prev], st = 0;
      if (X[prev] > bs) { bs = X[prev]; st = 1; }
      if (Y[prev] > bs) { bs = Y[prev]; st = 2; }
      M[idx] = bs + sf(ai, b[j - 1]);
      pM[idx] = st;

      const a1 = M[prow + j] - go, a2 = X[prow + j] - ge, a3 = Y[prow + j] - go;
      if (a1 >= a2 && a1 >= a3) { X[idx] = a1; pX[idx] = 0; }
      else if (a2 >= a3) { X[idx] = a2; pX[idx] = 1; }
      else { X[idx] = a3; pX[idx] = 2; }

      const b1 = M[row + j - 1] - go, b2 = Y[row + j - 1] - ge, b3 = X[row + j - 1] - go;
      if (b1 >= b2 && b1 >= b3) { Y[idx] = b1; pY[idx] = 0; }
      else if (b2 >= b3) { Y[idx] = b2; pY[idx] = 1; }
      else { Y[idx] = b3; pY[idx] = 2; }
    }
  }

  let i = n, j = m, state = 0, finalScore;
  let endI = n, endJ = m;
  if (freeEnds) {
    finalScore = NEG;
    for (let jj = 0; jj <= m; jj++) {
      const k = n * W + jj;
      if (M[k] > finalScore) { finalScore = M[k]; i = n; j = jj; state = 0; }
      if (X[k] > finalScore) { finalScore = X[k]; i = n; j = jj; state = 1; }
    }
    for (let ii = 0; ii <= n; ii++) {
      const k = ii * W + m;
      if (M[k] > finalScore) { finalScore = M[k]; i = ii; j = m; state = 0; }
      if (Y[k] > finalScore) { finalScore = Y[k]; i = ii; j = m; state = 2; }
    }
    endI = i; endJ = j;
  } else {
    const k = n * W + m;
    finalScore = M[k]; state = 0;
    if (X[k] > finalScore) { finalScore = X[k]; state = 1; }
    if (Y[k] > finalScore) { finalScore = Y[k]; state = 2; }
  }

  const ra = [], rb = [];
  let matches = 0, both = 0, guard = 0;
  while ((i > 0 || j > 0) && guard++ < n + m + 5) {
    if (i === 0) {
      ra.push("-"); rb.push(b[j - 1]); j--;
      continue;
    }
    if (j === 0) {
      ra.push(a[i - 1]); rb.push("-"); i--;
      continue;
    }
    const idx = i * W + j;
    if (state === 0) {
      ra.push(a[i - 1]); rb.push(b[j - 1]);
      both++;
      if (a[i - 1] === b[j - 1]) matches++;
      const st = pM[idx];
      i--; j--; state = st;
    } else if (state === 1) {
      ra.push(a[i - 1]); rb.push("-");
      const st = pX[idx];
      i--; state = st;
    } else {
      ra.push("-"); rb.push(b[j - 1]);
      const st = pY[idx];
      j--; state = st;
    }
  }
  ra.reverse(); rb.reverse();
  if (freeEnds && endI < n && endJ >= m) {
    for (let t = endI; t < n; t++) { ra.push(a[t]); rb.push("-"); }
  } else if (freeEnds && endJ < m && endI >= n) {
    for (let t = endJ; t < m; t++) { ra.push("-"); rb.push(b[t]); }
  }
  return { a: ra.join(""), b: rb.join(""), score: finalScore, matches, both };
}

function colScore(ca, cb, sf) {
  let sum = 0, cnt = 0;
  for (let i = 0; i < ca.length; i++) {
    const x = ca[i];
    if (x === "-") continue;
    for (let j = 0; j < cb.length; j++) {
      const y = cb[j];
      if (y === "-") continue;
      sum += sf(x, y);
      cnt++;
    }
  }
  return cnt ? sum / cnt : 0;
}

function alignProfiles(rowsA, rowsB, sf, go, ge, freeEnds) {
  const n = rowsA[0].length, m = rowsB[0].length, W = m + 1;
  const size = (n + 1) * W;
  const M = new Float64Array(size);
  const X = new Float64Array(size);
  const Y = new Float64Array(size);
  const pM = new Int8Array(size);
  const pX = new Int8Array(size);
  const pY = new Int8Array(size);

  const colsA = [], colsB = [];
  for (let i = 0; i < n; i++) {
    let s = "";
    for (const r of rowsA) s += r[i];
    colsA.push(s);
  }
  for (let j = 0; j < m; j++) {
    let s = "";
    for (const r of rowsB) s += r[j];
    colsB.push(s);
  }

  M[0] = 0;
  for (let i = 1; i <= n; i++) {
    const k = i * W;
    if (freeEnds) { M[k] = 0; X[k] = 0; Y[k] = NEG; pM[k] = 0; pX[k] = 0; }
    else { M[k] = NEG; X[k] = -(go + (i - 1) * ge); pX[k] = 1; Y[k] = NEG; }
  }
  for (let j = 1; j <= m; j++) {
    if (freeEnds) { M[j] = 0; Y[j] = 0; X[j] = NEG; pM[j] = 0; pY[j] = 0; }
    else { M[j] = NEG; Y[j] = -(go + (j - 1) * ge); pY[j] = 2; X[j] = NEG; }
  }
  X[0] = NEG; Y[0] = NEG;

  for (let i = 1; i <= n; i++) {
    const row = i * W, prow = (i - 1) * W;
    const ca = colsA[i - 1];
    for (let j = 1; j <= m; j++) {
      const idx = row + j;
      const prev = prow + j - 1;
      let bs = M[prev], st = 0;
      if (X[prev] > bs) { bs = X[prev]; st = 1; }
      if (Y[prev] > bs) { bs = Y[prev]; st = 2; }
      M[idx] = bs + colScore(ca, colsB[j - 1], sf);
      pM[idx] = st;

      const a1 = M[prow + j] - go, a2 = X[prow + j] - ge, a3 = Y[prow + j] - go;
      if (a1 >= a2 && a1 >= a3) { X[idx] = a1; pX[idx] = 0; }
      else if (a2 >= a3) { X[idx] = a2; pX[idx] = 1; }
      else { X[idx] = a3; pX[idx] = 2; }

      const b1 = M[row + j - 1] - go, b2 = Y[row + j - 1] - ge, b3 = X[row + j - 1] - go;
      if (b1 >= b2 && b1 >= b3) { Y[idx] = b1; pY[idx] = 0; }
      else if (b2 >= b3) { Y[idx] = b2; pY[idx] = 1; }
      else { Y[idx] = b3; pY[idx] = 2; }
    }
  }

  let i = n, j = m, state = 0;
  let endI = n, endJ = m;
  if (freeEnds) {
    let best = NEG;
    for (let jj = 0; jj <= m; jj++) {
      const k = n * W + jj;
      if (M[k] > best) { best = M[k]; i = n; j = jj; state = 0; }
      if (X[k] > best) { best = X[k]; i = n; j = jj; state = 1; }
    }
    for (let ii = 0; ii <= n; ii++) {
      const k = ii * W + m;
      if (M[k] > best) { best = M[k]; i = ii; j = m; state = 0; }
      if (Y[k] > best) { best = Y[k]; i = ii; j = m; state = 2; }
    }
    endI = i; endJ = j;
  } else {
    const k = n * W + m;
    let best = M[k]; state = 0;
    if (X[k] > best) { best = X[k]; state = 1; }
    if (Y[k] > best) { state = 2; }
  }

  const ops = [];
  let guard = 0;
  while ((i > 0 || j > 0) && guard++ < n + m + 5) {
    if (i === 0) { ops.push({ t: "Y", j: j - 1 }); j--; continue; }
    if (j === 0) { ops.push({ t: "X", i: i - 1 }); i--; continue; }
    const idx = i * W + j;
    if (state === 0) {
      ops.push({ t: "M", i: i - 1, j: j - 1 });
      const st = pM[idx];
      i--; j--; state = st;
    } else if (state === 1) {
      ops.push({ t: "X", i: i - 1 });
      const st = pX[idx];
      i--; state = st;
    } else {
      ops.push({ t: "Y", j: j - 1 });
      const st = pY[idx];
      j--; state = st;
    }
  }
  ops.reverse();
  if (freeEnds && endI < n && endJ >= m) {
    for (let t = endI; t < n; t++) ops.push({ t: "X", i: t });
  } else if (freeEnds && endJ < m && endI >= n) {
    for (let t = endJ; t < m; t++) ops.push({ t: "Y", j: t });
  }

  const outA = rowsA.map(() => []);
  const outB = rowsB.map(() => []);
  for (const op of ops) {
    if (op.t === "M") {
      for (let k = 0; k < rowsA.length; k++) outA[k].push(rowsA[k][op.i]);
      for (let k = 0; k < rowsB.length; k++) outB[k].push(rowsB[k][op.j]);
    } else if (op.t === "X") {
      for (let k = 0; k < rowsA.length; k++) outA[k].push(rowsA[k][op.i]);
      for (let k = 0; k < rowsB.length; k++) outB[k].push("-");
    } else {
      for (let k = 0; k < rowsA.length; k++) outA[k].push("-");
      for (let k = 0; k < rowsB.length; k++) outB[k].push(rowsB[k][op.j]);
    }
  }
  const resA = outA.map(arr => arr.join(""));
  const resB = outB.map(arr => arr.join(""));
  return resA.concat(resB);
}

function upgma(dist, n) {
  const clusters = [];
  for (let i = 0; i < n; i++) clusters.push({ size: 1, left: i, right: -1 });
  const d = [];
  for (let i = 0; i < n; i++) d.push(Array.from(dist[i]));
  let active = [];
  for (let i = 0; i < n; i++) active.push(i);
  const merges = [];

  while (active.length > 1) {
    let ia = active[0], ib = active[1], best = Infinity;
    for (let a = 0; a < active.length; a++) {
      for (let b = a + 1; b < active.length; b++) {
        const v = d[active[a]][active[b]];
        if (v < best) { best = v; ia = active[a]; ib = active[b]; }
      }
    }
    const ca = clusters[ia], cb = clusters[ib];
    const newIdx = clusters.length;
    clusters.push({ size: ca.size + cb.size, left: ia, right: ib });
    d.push(new Array(newIdx + 1).fill(0));
    const rest = active.filter(x => x !== ia && x !== ib);
    for (const k of rest) {
      const nd = (ca.size * d[ia][k] + cb.size * d[ib][k]) / (ca.size + cb.size);
      while (d[k].length <= newIdx) d[k].push(0);
      d[k][newIdx] = nd;
      d[newIdx][k] = nd;
    }
    merges.push({ a: ia, b: ib, node: newIdx });
    rest.push(newIdx);
    active = rest;
  }
  return { clusters, merges, root: active[0] };
}

function newickOf(clusters, root, n, nameOf) {
  const build = (idx) => {
    if (idx < n) return nameOf(idx).replace(/[(),:;\s]/g, "_");
    const c = clusters[idx];
    return "(" + build(c.left) + "," + build(c.right) + ")";
  };
  return build(root) + ";";
}

function detectType(seqs) {
  let total = 0, ntOk = 0;
  const ntSet = new Set("ACGTUNRYSWKMBDHVacgtunryswkmbdhv");
  const all = new Set();
  for (const s of seqs) {
    for (const ch of s) {
      if (ch === "-") continue;
      total++;
      all.add(ch.toUpperCase());
      if (ntSet.has(ch)) ntOk++;
    }
  }
  if (total === 0) return "nt";
  const onlyNt = [...all].every(c => ntSet.has(c));
  return onlyNt && ntOk / total >= 0.95 ? "nt" : "aa";
}

function runMSA(seqs, params, post) {
  const n = seqs.length;
  if (n < 2) throw new Error("至少需要 2 条序列");
  for (let i = 0; i < n; i++) {
    if (!seqs[i] || !seqs[i].length) throw new Error("第 " + (i + 1) + " 条序列为空");
  }

  let type = params.type;
  if (type === "auto") type = detectType(seqs);
  let matrixName = params.matrix;
  if (matrixName === "auto") matrixName = type === "aa" ? "blosum62" : "custom";
  if (type === "nt" && (matrixName === "blosum62" || matrixName === "pam250")) matrixName = "custom";

  const eff = {
    ...params, type, matrix: matrixName,
    match: Number(params.match), mismatch: Number(params.mismatch),
    gapOpen: Number(params.gapOpen), gapExtend: Number(params.gapExtend),
    freeEnds: !params.endGap
  };
  const sf = makeScorer(eff);
  const go = eff.gapOpen, ge = eff.gapExtend, freeEnds = eff.freeEnds;

  post({ type: "progress", pct: 1, stage: `准备：${n} 条序列 · ${type === "nt" ? "核酸" : "蛋白质"} · ${matrixName}` });

  const dist = Array.from({ length: n }, () => new Float64Array(n));
  const pairs = [];
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) pairs.push([i, j]);
  const totalPairs = pairs.length;
  const step = Math.max(1, Math.floor(totalPairs / 25));
  for (let t = 0; t < totalPairs; t++) {
    const [i, j] = pairs[t];
    const r = nwAlign(seqs[i], seqs[j], sf, go, ge, freeEnds);
    const id = r.both > 0 ? r.matches / r.both : 0;
    dist[i][j] = dist[j][i] = 1 - id;
    if ((t + 1) % step === 0 || t === totalPairs - 1) {
      post({ type: "progress", pct: 2 + ((t + 1) / totalPairs) * 56, stage: `两两比对 ${t + 1}/${totalPairs}` });
    }
  }

  post({ type: "progress", pct: 60, stage: "构建 UPGMA 引导树…" });
  const tree = upgma(dist, n);

  post({ type: "progress", pct: 64, stage: "渐进比对…" });
  const prof = new Array(tree.clusters.length);
  const profIds = new Array(tree.clusters.length);
  for (let i = 0; i < n; i++) { prof[i] = [seqs[i]]; profIds[i] = [i]; }

  for (let mi = 0; mi < tree.merges.length; mi++) {
    const { a, b, node } = tree.merges[mi];
    prof[node] = alignProfiles(prof[a], prof[b], sf, go, ge, freeEnds);
    profIds[node] = profIds[a].concat(profIds[b]);
    if (mi % 2 === 0 || mi === tree.merges.length - 1) {
      post({ type: "progress", pct: 64 + ((mi + 1) / tree.merges.length) * 32, stage: `渐进比对 ${mi + 1}/${tree.merges.length}` });
    }
  }

  const rootRows = prof[tree.root];
  const rootIds = profIds[tree.root];
  const aligned = new Array(n);
  for (let k = 0; k < n; k++) aligned[rootIds[k]] = rootRows[k];
  if (aligned.some(r => r === undefined)) throw new Error("渐进比对失败：行序恢复出错");
  if (new Set(aligned.map(r => r.length)).size !== 1) throw new Error("渐进比对失败：列长不一致");

  post({ type: "progress", pct: 97, stage: "计算 consensus…" });
  const L = aligned[0].length;
  const cons = [], consFrac = [];
  for (let c = 0; c < L; c++) {
    const count = Object.create(null);
    let nonGap = 0;
    for (const row of aligned) {
      const ch = row[c];
      if (ch === "-" || ch === ".") continue;
      nonGap++;
      count[ch] = (count[ch] || 0) + 1;
    }
    let bestCh = "-", bestN = 0;
    for (const k in count) if (count[k] > bestN) { bestN = count[k]; bestCh = k; }
    cons.push(nonGap === 0 ? "-" : bestCh);
    consFrac.push(bestN / n);
  }

  const nk = newickOf(tree.clusters, tree.root, n, (i) => "seq" + (i + 1));
  post({ type: "progress", pct: 100, stage: "完成" });
  post({
    type: "done", aligned, consensus: cons.join(""), consFrac,
    newick: nk, seqType: type, matrix: matrixName, params: eff
  });
}

if (typeof self !== "undefined" && typeof importScripts !== "undefined") {
  self.onmessage = (e) => {
    try {
      runMSA(e.data.seqs, e.data.params, (msg) => self.postMessage(msg));
    } catch (err) {
      self.postMessage({ type: "error", message: err.message || String(err) });
    }
  };
}

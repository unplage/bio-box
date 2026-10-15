"use strict";
const fs = require("fs");
const path = require("path");

const WORKER_PATH = process.env.MSA_WORKER || path.join(__dirname, "..", "msa-worker.js");

function loadWorker() {
  const src = fs.readFileSync(WORKER_PATH, "utf8").replace(/^\s*"use strict";\s*/m, "");
  const names = ["makeScorer", "nwAlign", "alignProfiles", "neighborJoining", "newickOf", "detectType", "runMSA"];
  const fn = new Function("self", src + "\nreturn {" + names.join(",") + "};");
  return fn({});
}

const W = loadWorker();

let passed = 0, failed = 0;
const failures = [];
const filter = process.argv[2] || "";

function assert(cond, msg) {
  if (!cond) throw new Error(msg || "assertion failed");
}
function assertEq(got, want, msg) {
  if (got !== want) throw new Error((msg || "") + " expected " + JSON.stringify(want) + ", got " + JSON.stringify(got));
}
function assertClose(got, want, tol, msg) {
  if (!(Math.abs(got - want) <= tol)) throw new Error((msg || "") + " expected ≈" + want + ", got " + got);
}

function runTest(name, fn) {
  if (filter && !name.includes(filter)) return;
  const t0 = Date.now();
  try {
    fn();
    passed++;
    console.log("  ✓ " + name + " (" + (Date.now() - t0) + "ms)");
  } catch (e) {
    failed++;
    failures.push({ name, err: e });
    console.log("  ✗ " + name);
    console.log("    " + (e && e.message));
  }
}
function suite(name) {
  console.log("\n[" + name + "]");
}

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function randSeq(len, alphabet, rng) {
  let s = "";
  for (let i = 0; i < len; i++) s += alphabet[Math.floor(rng() * alphabet.length)];
  return s;
}

function gapRunList(row) {
  const runs = [];
  let run = 0, start = 0;
  for (let k = 0; k <= row.length; k++) {
    if (k < row.length && row[k] === "-") {
      if (run === 0) start = k;
      run++;
    } else if (run) {
      runs.push({ start, len: run });
      run = 0;
    }
  }
  return runs;
}
function scoreAlignment(ar, br, sf, go, ge, freeEnds) {
  let s = 0;
  for (let k = 0; k < ar.length; k++) {
    if (ar[k] !== "-" && br[k] !== "-") s += sf(ar[k], br[k]);
  }
  const cost = (row) => {
    let c = 0;
    for (const r of gapRunList(row)) {
      const leading = r.start === 0;
      const trailing = r.start + r.len === row.length;
      if (freeEnds && (leading || trailing)) continue;
      c += go + (r.len - 1) * ge;
    }
    return c;
  };
  return s - cost(ar) - cost(br);
}
function bruteMaxScore(a, b, sf, go, ge, freeEnds) {
  let best = -Infinity;
  const walk = (i, j, pa, pb) => {
    if (i === a.length && j === b.length) {
      const sc = scoreAlignment(pa, pb, sf, go, ge, freeEnds);
      if (sc > best) best = sc;
      return;
    }
    if (i < a.length && j < b.length) walk(i + 1, j + 1, pa + a[i], pb + b[j]);
    if (i < a.length) walk(i + 1, j, pa + a[i], pb + "-");
    if (j < b.length) walk(i, j + 1, pa + "-", pb + b[j]);
  };
  walk(0, 0, "", "");
  return best;
}
function sisterPairs(newick) {
  const pairs = [];
  const re = /\(seq(\d+):[^,()]+,seq(\d+):[^,()]+\)/g;
  let m;
  while ((m = re.exec(newick))) {
    pairs.push([Number(m[1]), Number(m[2])].sort((x, y) => x - y).join("-"));
  }
  return pairs;
}
function parseNewickTree(newick) {
  const s = newick.trim();
  let pos = 0;
  let nextId = 0;
  const nodes = [];
  function parseNode() {
    const id = nextId++;
    nodes.push({ id, name: null, children: [] });
    if (s[pos] === "(") {
      pos++;
      for (;;) {
        const child = parseNode();
        let len = 0;
        if (s[pos] === ":") {
          pos++;
          const st = pos;
          while (pos < s.length && /[-\d.eE+]/.test(s[pos])) pos++;
          len = Number(s.slice(st, pos));
        }
        nodes[id].children.push({ id: child, len });
        if (s[pos] === ",") { pos++; continue; }
        break;
      }
      if (s[pos] === ")") pos++;
    } else {
      const st = pos;
      while (pos < s.length && s[pos] !== ":" && s[pos] !== "," && s[pos] !== ")" && s[pos] !== ";") pos++;
      nodes[id].name = s.slice(st, pos);
    }
    return id;
  }
  const root = parseNode();
  const parent = new Map();
  function collect(id) {
    for (const c of nodes[id].children) {
      parent.set(c.id, { from: id, len: c.len });
      collect(c.id);
    }
  }
  collect(root);
  const leafId = {};
  for (const nd of nodes) if (nd.name) leafId[nd.name] = nd.id;
  function pathLen(a, b) {
    const chains = [];
    for (const start of [a, b]) {
      const chain = [];
      let cur = start;
      while (cur !== root) {
        const p = parent.get(cur);
        chain.push({ node: cur, len: p.len });
        cur = p.from;
      }
      chains.push(chain);
    }
    const [ca, cb] = chains;
    let sum = 0;
    const inA = new Set(ca.map(x => x.node));
    let shared = null;
    for (const x of cb) if (inA.has(x.node)) { shared = x.node; break; }
    for (const x of ca) { if (x.node === shared) break; sum += x.len; }
    for (const x of cb) { if (x.node === shared) break; sum += x.len; }
    return sum;
  }
  return { leafId, pathLen };
}

const AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY";
const NT_ALPHABET = "ACGT";
const sfB62 = W.makeScorer({ type: "aa", matrix: "blosum62", match: 2, mismatch: -1 });
const sfPam = W.makeScorer({ type: "aa", matrix: "pam250", match: 2, mismatch: -1 });
const sfCustom = W.makeScorer({ type: "aa", matrix: "custom", match: 2, mismatch: -1 });
const sfNt = W.makeScorer({ type: "nt", matrix: "custom", match: 1, mismatch: -4 });

suite("makeScorer");
runTest("BLOSUM62 对角线为标准值", () => {
  const diag = { A: 4, R: 5, N: 6, D: 6, C: 9, Q: 5, E: 5, G: 6, H: 8, I: 4, L: 4, K: 5, M: 5, F: 6, P: 7, S: 4, T: 5, W: 11, Y: 7, V: 4, B: 4, Z: 4 };
  for (const [ch, v] of Object.entries(diag)) assertEq(sfB62(ch, ch), v, "B62 " + ch);
});
runTest("PAM250 抽样锁定现表", () => {
  assertEq(sfPam("A", "A"), 2);
  assertEq(sfPam("W", "W"), 14);
  assertEq(sfPam("M", "M"), 3);
  assertEq(sfPam("Y", "Y"), 8);
});
runTest("自定义 aa 回退 BLOSUM62，nt 用 match/mismatch", () => {
  assertEq(sfCustom("A", "A"), 4, "aa+custom → BLOSUM62");
  assertEq(sfCustom("?", "?"), -1, "未知残基按 mismatch");
  assertEq(sfNt("A", "A"), 1);
  assertEq(sfNt("A", "T"), -4);
});

suite("nwAlign 基础");
runTest("相同序列全匹配", () => {
  const r = W.nwAlign("MKTAYIAKQR", "MKTAYIAKQR", sfB62, 10, 0.5, false);
  assertEq(r.both, 10);
  assertEq(r.matches, 10);
  assertEq(r.score, 49);
  assertEq(r.a, "MKTAYIAKQR");
});
runTest("单残基", () => {
  const r = W.nwAlign("W", "W", sfB62, 10, 0.5, false);
  assertEq(r.score, 11);
  assertEq(r.both, 1);
});
runTest("AAAA/AA 全局 = -2.5", () => {
  const r = W.nwAlign("AAAA", "AA", sfCustom, 10, 0.5, false);
  assertClose(r.score, -2.5, 1e-9);
  assertEq(scoreAlignment(r.a, r.b, sfCustom, 10, 0.5, false), r.score, "回溯自洽");
});
runTest("自由末端：前导缺口免费（AAAM/M = 5，全局 = -6）", () => {
  const free = W.nwAlign("AAAM", "M", sfB62, 10, 0.5, true);
  assertClose(free.score, 5, 1e-9);
  assertEq(scoreAlignment(free.a, free.b, sfB62, 10, 0.5, true), free.score);
  const glob = W.nwAlign("AAAM", "M", sfB62, 10, 0.5, false);
  assertClose(glob.score, -6, 1e-9);
  assertEq(scoreAlignment(glob.a, glob.b, sfB62, 10, 0.5, false), glob.score);
});

suite("穷举 oracle（捕获回溯指针类 bug）");
runTest("ACG 字母表 len≤3 全部序列对 × 全局/自由末端", () => {
  const seqs = [];
  for (const len of [1, 2, 3]) {
    const count = Math.pow(3, len);
    for (let x = 0; x < count; x++) {
      let s = "", t = x;
      for (let i = 0; i < len; i++) { s += "ACG"[t % 3]; t = Math.floor(t / 3); }
      seqs.push(s);
    }
  }
  let checked = 0;
  for (const a of seqs) {
    for (const b of seqs) {
      for (const free of [false, true]) {
        const r = W.nwAlign(a, b, sfCustom, 10, 0.5, free);
        const oracle = bruteMaxScore(a, b, sfCustom, 10, 0.5, free);
        assertClose(r.score, oracle, 1e-9, `score ${a} vs ${b} freeEnds=${free}`);
        const recon = scoreAlignment(r.a, r.b, sfCustom, 10, 0.5, free);
        assertClose(recon, r.score, 1e-9, `reconstruction ${a} vs ${b} freeEnds=${free}`);
        checked++;
      }
    }
  }
  assert(checked > 1000, "pairs=" + checked);
});
runTest("随机 len4-5 × BLOSUM62 vs 穷举", () => {
  const rng = mulberry32(42);
  for (let t = 0; t < 30; t++) {
    const a = randSeq(4 + Math.floor(rng() * 2), "ACG", rng);
    const b = randSeq(4 + Math.floor(rng() * 2), "ACG", rng);
    for (const free of [false, true]) {
      const r = W.nwAlign(a, b, sfCustom, 10, 0.5, free);
      assertClose(r.score, bruteMaxScore(a, b, sfCustom, 10, 0.5, free), 1e-9, `${a} vs ${b}`);
    }
  }
});
runTest("随机蛋白对：回溯重建 = DP 分数（30 对 × 双模式）", () => {
  const rng = mulberry32(7);
  for (let t = 0; t < 30; t++) {
    const a = randSeq(10 + Math.floor(rng() * 50), AA_ALPHABET, rng);
    const b = randSeq(10 + Math.floor(rng() * 50), AA_ALPHABET, rng);
    for (const free of [false, true]) {
      const r = W.nwAlign(a, b, sfB62, 10, 0.5, free);
      const recon = scoreAlignment(r.a, r.b, sfB62, 10, 0.5, free);
      assertClose(recon, r.score, 1e-9, `recon ${t} freeEnds=${free}`);
    }
  }
});

suite("历史 bug 回归（pY 指针互换）");
const E51H2 = "QVQLVQSGAEVKKPGSSVKVSCKDSDGTVFPIAYMSWVRQAPGQGLEWMGGIFPSIGRTIYGEKFEDRVTITADESTSTAYMELSSLRSEDTAVYYCARGRTYWEYYHAMDNWGQGTTVTVSS";
const h82H = "QVQLVQSGAEVKKPGASVKVSCKASGFYTTTYWMHWVRQAPGQGLEWMGMIHPSDGETRLNQNFKDRVTMTRDTSTSTVYMELSSLRSEDTAVYYCAREDMGGYPWFAYWGQGTLVTVSS";
runTest("E51H2 vs h82H：全局 400/117/87", () => {
  const r = W.nwAlign(E51H2, h82H, sfB62, 10, 0.5, false);
  assertEq(r.score, 400);
  assertEq(r.both, 117);
  assertEq(r.matches, 87);
  assertClose(scoreAlignment(r.a, r.b, sfB62, 10, 0.5, false), 400, 1e-9, "重建");
});
runTest("E51H2 vs h82H：自由末端同样 400/117/87", () => {
  const r = W.nwAlign(E51H2, h82H, sfB62, 10, 0.5, true);
  assertEq(r.score, 400);
  assertEq(r.both, 117);
  assertEq(r.matches, 87);
  assertClose(scoreAlignment(r.a, r.b, sfB62, 10, 0.5, true), 400, 1e-9, "重建");
});

suite("neighborJoining");
runTest("4 分类可加矩阵：拓扑 + 路径距离还原", () => {
  const d = [
    [0, 0.3, 0.5, 0.35],
    [0.3, 0, 0.6, 0.45],
    [0.5, 0.6, 0, 0.45],
    [0.35, 0.45, 0.45, 0]
  ];
  const tree = W.neighborJoining(d, 4);
  const nk = W.newickOf(tree.clusters, tree.root, 4, (i) => "seq" + (i + 1));
  const pairs = sisterPairs(nk);
  assert(pairs.includes("1-2"), "sisters=" + pairs.join(" ") + " nk=" + nk);
  assert(pairs.includes("3-4"), "sisters=" + pairs.join(" "));
  const { leafId, pathLen } = parseNewickTree(nk);
  for (let i = 0; i < 4; i++) {
    for (let j = i + 1; j < 4; j++) {
      const got = pathLen(leafId["seq" + (i + 1)], leafId["seq" + (j + 1)]);
      assertClose(got, d[i][j], 1e-4, `path seq${i + 1}-seq${j + 1}`);
    }
  }
});

suite("runMSA 不变量");
function msaParams(extra) {
  return Object.assign({ type: "aa", matrix: "blosum62", match: 2, mismatch: -1, gapOpen: 10, gapExtend: 0.5, endGap: false }, extra || {});
}
function collect(seqs, params) {
  let done = null;
  W.runMSA(seqs, params, (m) => { if (m.type === "done") done = m; if (m.type === "error") throw new Error(m.message); });
  assert(done, "未收到 done");
  return done;
}
function checkInvariants(d, n, alphabet) {
  const L = d.aligned[0].length;
  assertEq(d.aligned.length, n, "行数");
  for (const row of d.aligned) {
    assertEq(row.length, L, "列长一致");
    for (const ch of row) assert(ch === "-" || alphabet.includes(ch), "非法字符 " + ch);
  }
  assertEq(d.consensus.length, L, "consensus 长度");
  assertEq(d.consFrac.length, L, "consFrac 长度");
  for (const f of d.consFrac) assert(f >= 0 && f <= 1, "consFrac 范围");
  let allGap = 0;
  for (let c = 0; c < L; c++) {
    if (d.aligned.every(row => row[c] === "-")) allGap++;
  }
  assertEq(allGap, 0, "全 gap 列数");
  const m = d.newick.match(/seq(\d+)/g) || [];
  const ids = m.map(x => Number(x.slice(3))).sort((a, b) => a - b);
  assertEq(ids.join(","), Array.from({ length: n }, (_, i) => i + 1).join(","), "newick 叶子齐全");
}
runTest("相同序列 ×3：L=长度、无全 gap、consensus=序列", () => {
  const s = "MKTAYIAKQRQISFVKSHFSRQ";
  const d = collect([s, s, s], msaParams());
  assertEq(d.aligned[0].length, s.length);
  assertEq(d.consensus, s);
  checkInvariants(d, 3, AA_ALPHABET);
});
runTest("随机 nt 集合不变量 + 确定性", () => {
  const rng = mulberry32(99);
  const base = randSeq(30, NT_ALPHABET, rng);
  const seqs = [base];
  for (let i = 0; i < 4; i++) {
    const mutated = base.split("");
    for (let k = 0; k < 3; k++) {
      const p = Math.floor(rng() * mutated.length);
      mutated[p] = NT_ALPHABET[Math.floor(rng() * 4)];
    }
    seqs.push(mutated.join(""));
  }
  const params = msaParams({ type: "nt", matrix: "custom", match: 1, mismatch: -4 });
  const d1 = collect(seqs, params);
  checkInvariants(d1, 5, NT_ALPHABET);
  const d2 = collect(seqs, params);
  assertEq(JSON.stringify({ a: d1.aligned, c: d1.consensus, n: d1.newick }),
           JSON.stringify({ a: d2.aligned, c: d2.consensus, n: d2.newick }), "两次运行完全一致");
});
runTest("抗 6 序列夹具：L=137、无全 gap、正确姐妹对", () => {
  const seqs = [
    "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYWMHWVRQAPGQGLEWIGQIDPSDSYTYYNEDFKDRATLTVDQSTSTAYMELSSLRSEDTAVYYCAKGYYDYDWGYAMDYWGQGTTVTVSS",
    E51H2,
    "DIQMTQSPSSLSASVGDRVTITCRASQGISSYLAWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKLEIK",
    "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYWMHWVRQAPGQGLEWMGQIDPSDSYTYYNEDFKDRATMTVDTSTSTAYMELSSLRSEDTAVYYCAKGYYDYDWGYAMDYWGQGTTVTVSS",
    h82H,
    "DIVMTQSPDSLAVSLGERATISCKASKSVSSSAYSYVHWYQQKPGQPPKLLIYLASNLGSGVPARFSGSGSGTDFTLTISSVQAEDVAVYYCQHSWALPYTFGGGTKLEIK"
  ];
  const d = collect(seqs, msaParams());
  assertEq(d.aligned[0].length, 137, "L");
  assertEq(d.seqType, "aa");
  checkInvariants(d, 6, AA_ALPHABET);
  const pairs = sisterPairs(d.newick);
  assert(pairs.includes("1-4"), "C45H4-C45H3 sisters: " + pairs.join(" "));
  assert(pairs.includes("2-5"), "E51H2-h82H sisters: " + pairs.join(" "));
  assert(pairs.includes("3-6"), "VL-h82L sisters: " + pairs.join(" "));
});
runTest("输入校验：1 条抛错、空序列抛错", () => {
  let e1 = null, e2 = null;
  try { W.runMSA(["ACGT"], msaParams(), () => {}); } catch (e) { e1 = e; }
  assert(e1 && e1.message.includes("至少"), "n<2 应抛错");
  try { W.runMSA(["ACGT", ""], msaParams(), () => {}); } catch (e) { e2 = e; }
  assert(e2 && e2.message.includes("为空"), "空序列应抛错");
});

suite("detectType");
runTest("纯 nt / 蛋白 / 含歧义", () => {
  assertEq(W.detectType(["ACGTACGTNN", "ACGTACGTAA"]), "nt");
  assertEq(W.detectType(["MKTAYIAKQR", "MKTAYIAKQRWY"]), "aa");
  assertEq(W.detectType(["ACGTACGT??", "ACGTACGTTT"]), "aa", "? 不属于核酸字母表");
});

suite("alignProfiles");
runTest("相同 profile：无全 gap 列、列长一致", () => {
  const rows = ["ACDE", "ACDE", "ACDE"];
  const out = W.alignProfiles(rows, ["ACDE"], sfCustom, 10, 0.5, false);
  assertEq(out.length, 4);
  const L = out[0].length;
  for (const row of out) assertEq(row.length, L);
  for (let c = 0; c < L; c++) assert(out.some(row => row[c] !== "-"), "全 gap 列 at " + c);
});

console.log("\n" + "=".repeat(50));
if (filter) console.log("过滤: " + filter);
console.log("通过 " + passed + " / 失败 " + failed + " / 共 " + (passed + failed));
if (failed) {
  console.log("\n失败详情:");
  for (const f of failures) {
    console.log("- " + f.name);
    console.log("  " + f.err.message);
    if (f.err.stack) console.log("  " + f.err.stack.split("\n").slice(1, 3).join("\n  "));
  }
  process.exit(1);
}

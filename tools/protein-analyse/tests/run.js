"use strict";
/*
 * protein-analyse 零依赖回归测试（Node，无 npm 依赖）
 *   node tools/protein-analyse/tests/run.js          跑全部用例
 *   node tools/protein-analyse/tests/run.js fasta    只跑名字含 "fasta" 的用例
 *   PROTEIN_ANALYSE_HTML=<path> 指定被测 index.html
 *
 * 覆盖重点：批量 FASTA 解析（按 > 标题行区分序列）与批量 pI / MW 计算，
 * 并断言批量结果与单条计算路径逐位一致。
 */
const fs = require("fs");
const path = require("path");

const HTML_PATH = process.env.PROTEIN_ANALYSE_HTML || path.join(__dirname, "..", "index.html");

/* ── 最小 DOM 桩：让 index.html 内联脚本能在 Node 中加载 ── */
function makeStubEl(id) {
    return {
        id: id,
        value: "",
        checked: false,
        textContent: "",
        innerHTML: "",
        className: "",
        style: {},
        dataset: {},
        files: [],
        classList: { add() {}, remove() {}, toggle() {}, contains() { return false; } },
        addEventListener() {}, removeEventListener() {},
        appendChild() {}, removeChild() {}, insertBefore() {}, replaceChild() {},
        click() {}, focus() {}, blur() {}, scrollIntoView() {},
        querySelector() { return null; }, querySelectorAll() { return []; },
        getAttribute() { return null; }, setAttribute() {}, removeAttribute() {}
    };
}

function loadApi() {
    const html = fs.readFileSync(HTML_PATH, "utf8");
    const m = html.match(/<script>([\s\S]*?)<\/script>/);
    if (!m) throw new Error("未在 " + HTML_PATH + " 中找到内联 <script>");
    const cache = new Map();
    const document = {
        getElementById(id) {
            if (!cache.has(id)) cache.set(id, makeStubEl(id));
            return cache.get(id);
        },
        querySelector(sel) { return makeStubEl(sel); },
        querySelectorAll() { return []; },
        createElement() { return makeStubEl("created"); },
        createElementNS() { return makeStubEl("created"); },
        body: makeStubEl("body"),
        documentElement: makeStubEl("html"),
        addEventListener() {}
    };
    const navigator = {};
    const location = { protocol: "file:" };
    const names = [
        "AA_MASS", "WATER", "PKA_SETS", "PKA_SET_ORDER",
        "parseSeq", "countResidues", "molecularWeight", "effectiveCysCount",
        "netCharge", "isoelectricPoint",
        "piCleanLine", "piParseFasta", "piComputeBatch", "calcPiBatch"
    ];
    const factory = new Function(
        "document", "window", "navigator", "location", "console", "alert",
        m[1] + "\nreturn {" + names.join(",") + ", get batchLastRun() { return piBatchLastRun; }};"
    );
    const api = factory(document, globalThis, navigator, location, console, alert);
    api.document = document;
    api.alerts = alertStub.messages;
    return api;
}

const alertStub = { messages: [], fn: (msg) => { alertStub.messages.push(String(msg)); } };
function alert(msg) { alertStub.fn(msg); }

const A = loadApi();

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
    try {
        fn();
        passed++;
        console.log("  ✓ " + name);
    } catch (e) {
        failed++;
        failures.push({ name, err: e });
        console.log("  ✗ " + name + " — " + e.message);
    }
}

const FASTA2 = ">seq1 first\nMALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT\n>seq2 second\nEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK";

/* ── 解析 ── */
runTest("fasta: 标题行区分多条序列", () => {
    const entries = A.piParseFasta(FASTA2, "fasta");
    assertEq(entries.length, 2, "序列条数");
    assertEq(entries[0].id, "seq1 first", "第 1 条 ID");
    assertEq(entries[1].id, "seq2 second", "第 2 条 ID");
    assertEq(entries[0].seq.length, 54, "第 1 条长度");
    assert(entries[0].seq.startsWith("MALWMRLLPLL"), "第 1 条序列开头");
});

runTest("fasta: 多行序列拼接为一条", () => {
    const entries = A.piParseFasta(">a\nMALWM\nRLLPL\nLALLAL\n>b\nPEPTIDE\nPEPTIDE", "fasta");
    assertEq(entries.length, 2, "序列条数");
    assertEq(entries[0].seq, "MALWMRLLPLLALLAL", "多行拼接");
    assertEq(entries[1].seq, "PEPTIDEPEPTIDE", "多行拼接（第 2 条）");
});

runTest("fasta: 无标题行时整段视为一条序列", () => {
    const entries = A.piParseFasta("MALWMRLLPLLALLAL\nPEPTIDE", "fasta");
    assertEq(entries.length, 1, "序列条数");
    assertEq(entries[0].id, "unnamed", "ID");
    assertEq(entries[0].seq, "MALWMRLLPLLALLALPEPTIDE", "整段拼接");
});

runTest("fasta: 每行一条模式", () => {
    const entries = A.piParseFasta("PEPTIDE\nPEPTIDEK\nPEPTIDEKEE", "lines");
    assertEq(entries.length, 3, "序列条数");
    assertEq(entries[2].id, "seq3", "ID 编号");
    assertEq(entries[2].seq, "PEPTIDEKEE", "序列内容");
});

runTest("fasta: 非标准字符被过滤并计数", () => {
    const entries = A.piParseFasta(">x\nPEPT*IDE X", "fasta");
    assertEq(entries.length, 1, "序列条数");
    assertEq(entries[0].seq, "PEPTIDE", "过滤后序列");
    assertEq(entries[0].invalid, "X", "非标准字符记录");
});

runTest("fasta: 空输入返回空数组", () => {
    assertEq(A.piParseFasta("", "fasta").length, 0, "空输入");
    assertEq(A.piParseFasta("   \n  ", "fasta").length, 0, "空白输入");
});

runTest("fasta: 与单条 parseSeq 的清洗规则一致", () => {
    const single = A.parseSeq(">x\nPEPT*IDE X").seq;
    const batch = A.piParseFasta(">x\nPEPT*IDE X", "fasta")[0].seq;
    assertEq(batch, single, "批量/单条清洗结果");
});

/* ── 计算 ── */
runTest("calc: MW 与残基质量求和一致（PEPTIDE ≈ 799.83 Da）", () => {
    const mw = A.molecularWeight("PEPTIDE");
    assertClose(mw, 799.8328, 0.001, "PEPTIDE MW");
});

runTest("calc: 批量 MW/pI 与单条路径逐位一致", () => {
    const entries = A.piParseFasta(FASTA2, "fasta");
    const run = A.piComputeBatch(entries, "bjellqvist", "free");
    assertEq(run.successful, 2, "成功条数");
    for (const r of run.results) {
        assertEq(r.mw, A.molecularWeight(r.seq), "MW 一致性");
        assertEq(r.pI, A.isoelectricPoint(r.seq, "bjellqvist", "free"), "pI 一致性");
        for (const k of A.PKA_SET_ORDER) {
            assertEq(r.pkaPi[k], A.isoelectricPoint(r.seq, k, "free"), "pKa 参数集一致性：" + k);
        }
    }
});

runTest("calc: 预计算计数不改变 netCharge 结果", () => {
    const seq = A.piParseFasta(FASTA2, "fasta")[0].seq;
    const cnt = A.countResidues(seq);
    for (const pH of [0, 3.5, 7, 10.5, 14]) {
        assertEq(A.netCharge(seq, pH, "bjellqvist", "free", cnt),
            A.netCharge(seq, pH, "bjellqvist", "free"), "netCharge 一致性 pH=" + pH);
    }
});

runTest("calc: Cys 处理方式影响 pI", () => {
    const entries = A.piParseFasta(">c\nPEPTIDECPEPTIDEC", "fasta");
    const free = A.piComputeBatch(entries, "bjellqvist", "free").results[0].pI;
    const pair = A.piComputeBatch(entries, "bjellqvist", "pair").results[0].pI;
    const none = A.piComputeBatch(entries, "bjellqvist", "none").results[0].pI;
    assert(!(free === pair && pair === none), "三种 Cys 模式应给出不同 pI");
});

runTest("calc: 批量保留各参数集 pI 且顺序稳定", () => {
    const entries = A.piParseFasta(FASTA2, "fasta");
    const run = A.piComputeBatch(entries, "ipc_protein", "free");
    for (const r of run.results) {
        for (const k of A.PKA_SET_ORDER) {
            assert(typeof r.pkaPi[k] === "number" && isFinite(r.pkaPi[k]), "pKa pI 有效：" + k);
        }
        assertEq(r.pI, r.pkaPi.ipc_protein, "所选参数集即对应列");
    }
});

runTest("calc: 空序列被计入失败而非中断整批", () => {
    const run = A.piComputeBatch([{ id: "ok", seq: "PEPTIDE" }, { id: "bad", seq: "" }], "bjellqvist", "free");
    assertEq(run.total, 2, "总数");
    assertEq(run.successful, 1, "成功数");
    assertEq(run.failed, 1, "失败数");
    assertEq(run.errors[0].id, "bad", "失败 ID");
});

runTest("calc: 批量规模性能（200 条 × 约 300 aa）", () => {
    const base = A.piParseFasta(FASTA2, "fasta")[0].seq.repeat(6).slice(0, 300);
    const entries = [];
    for (let i = 0; i < 200; i++) entries.push({ id: "s" + i, seq: base });
    const t0 = Date.now();
    const run = A.piComputeBatch(entries, "bjellqvist", "free");
    const dt = Date.now() - t0;
    assertEq(run.successful, 200, "全部成功");
    assert(dt < 5000, "耗时 " + dt + "ms 应 < 5000ms");
});

runTest("ui: 批量计算渲染表格与汇总", () => {
    const doc = A.document;
    doc.getElementById("piBatchInput").value = FASTA2;
    doc.getElementById("piBatchMode").value = "fasta";
    doc.getElementById("pkaSet").value = "bjellqvist";
    doc.getElementById("pkaCysMode").value = "free";
    A.calcPiBatch();
    const table = doc.getElementById("piBatchTable").innerHTML;
    assert(table.includes("seq1 first") && table.includes("seq2 second"), "表格含两条序列 ID");
    assert(table.includes("54"), "表格含第 1 条长度");
    assert(doc.getElementById("piBatchSummary").innerHTML.includes("共 2 条序列"), "汇总文案");
    assertEq(A.batchLastRun.run.successful, 2, "批量结果条数");
});

runTest("ui: 空输入触发提示且不产生结果", () => {
    const doc = A.document;
    const before = A.alerts.length;
    doc.getElementById("piBatchInput").value = "   ";
    A.calcPiBatch();
    assert(A.alerts.length > before, "应弹出提示");
    assertEq(A.batchLastRun, null, "无结果");
});

console.log("");
console.log(failed === 0
    ? "全部通过：" + passed + " 个用例"
    : "通过 " + passed + " 个，失败 " + failed + " 个");
if (failures.length) {
    for (const f of failures) console.log("  - " + f.name + ": " + f.err.stack);
    process.exit(1);
}

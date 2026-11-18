"use strict";

const fs = require("fs");
const path = require("path");
const assert = require("assert");
const vm = require("vm");

const root = path.resolve(__dirname, "..");
const dataDir = path.join(root, "data");
const read = name => fs.readFileSync(path.join(root, name), "utf8");
const readJson = name => JSON.parse(read(path.join("data", name + ".json")));
const tests = [];

function test(name, fn) {
  tests.push({ name, fn });
}

const meta = readJson("meta");
const refs = readJson("refs");
const nodes = readJson("nodes");
const edges = readJson("edges");
const drugs = readJson("drugs");

const moduleIds = new Set(meta.modules.map(m => m.id));
const nodeIds = new Set(nodes.map(n => n.id));
const refIds = new Set(Object.keys(refs));
const nodeKinds = new Set(["receptor", "ligand", "cytokine", "complex", "cell", "drug"]);
const relTypes = new Set(["bind", "activate", "inhibit", "compete", "block", "agonist", "regulate", "recruit"]);
const mechanisms = new Set(["block", "agonist"]);

function checkRefs(list, owner) {
  assert(Array.isArray(list), `${owner} 的 refs 必须是数组`);
  for (const id of list) assert(refIds.has(id), `${owner} 引用了不存在的文献 ${id}`);
}

function checkModules(list, owner) {
  assert(Array.isArray(list) && list.length, `${owner} 的 modules 必须是非空数组`);
  for (const id of list) assert(moduleIds.has(id), `${owner} 引用了未知模块 ${id}`);
}

test("数据文件存在且为合法 JSON", () => {
  for (const name of ["meta", "refs", "nodes", "edges", "drugs"]) {
    assert(fs.existsSync(path.join(dataDir, name + ".json")), `缺少 data/${name}.json`);
  }
});

test("meta 字段与模块定义完整", () => {
  for (const key of ["title", "snapshot", "disclaimer", "modules"]) assert(meta[key], `meta.${key} 不能为空`);
  assert(/^\d{4}-\d{2}$/.test(meta.snapshot), "meta.snapshot 应为 YYYY-MM");
  assert.strictEqual(moduleIds.size, meta.modules.length, "meta.modules 的 id 不能重复");
  for (const mod of meta.modules) {
    assert(mod.id && mod.name && mod.desc, `模块 ${mod.id || "?"} 字段不完整`);
    assert(Array.isArray(mod.groups) && mod.groups.length, `模块 ${mod.id} 缺少 groups`);
  }
});

test("节点 ID 唯一且字段有效", () => {
  const seen = new Set();
  for (const node of nodes) {
    assert(node.id && !seen.has(node.id), `节点 ID 重复: ${node.id}`);
    seen.add(node.id);
    assert(node.label && node.kind && node.group && node.note, `节点 ${node.id} 字段不完整`);
    assert(nodeKinds.has(node.kind), `节点 ${node.id} 的 kind 无效: ${node.kind}`);
    checkModules(node.modules, `节点 ${node.id}`);
    checkRefs(node.refs || [], `节点 ${node.id}`);
    const groups = meta.modules.filter(m => node.modules.includes(m.id)).flatMap(m => m.groups);
    assert(groups.includes(node.group), `节点 ${node.id} 的 group ${node.group} 不在所属模块分组中`);
  }
});

test("边端点存在、关系类型有效且 ID 唯一", () => {
  const seen = new Set();
  for (const edge of edges) {
    const key = `${edge.source}|${edge.target}|${edge.rel}`;
    assert(!seen.has(key), `重复关系: ${key}`);
    seen.add(key);
    assert(nodeIds.has(edge.source), `关系源节点不存在: ${edge.source}`);
    assert(nodeIds.has(edge.target), `关系目标节点不存在: ${edge.target}`);
    assert(relTypes.has(edge.rel), `关系类型无效: ${edge.rel}`);
    assert(edge.label && edge.note, `关系 ${key} 缺少 label 或 note`);
    assert(["A", "B", "C"].includes(edge.evidence), `关系 ${key} 的 evidence 无效: ${edge.evidence}`);
    checkModules(edge.modules, `关系 ${key}`);
    checkRefs(edge.refs || [], `关系 ${key}`);
  }
});

test("药物靶点存在、机制有效且 ID 唯一", () => {
  for (const drug of drugs) {
    assert(drug.id && drug.label && drug.inn, `药物 ${drug.id || "?"} 字段不完整`);
    assert(mechanisms.has(drug.mechanism), `药物 ${drug.id} 的 mechanism 无效: ${drug.mechanism}`);
    assert(drug.modality && drug.stage, `药物 ${drug.id} 缺少 modality 或 stage`);
    assert(Array.isArray(drug.targets) && drug.targets.length, `药物 ${drug.id} 缺少 targets`);
    for (const target of drug.targets) assert(nodeIds.has(target), `药物 ${drug.id} 的靶点不存在: ${target}`);
    checkModules(drug.modules, `药物 ${drug.id}`);
    checkRefs(drug.refs || [], `药物 ${drug.id}`);
    if (drug.graph !== undefined) assert.strictEqual(typeof drug.graph, "boolean", `药物 ${drug.id} 的 graph 应为布尔值`);
  }
  const seen = new Set();
  for (const drug of drugs) {
    assert(!seen.has(drug.id), `药物 ID 重复: ${drug.id}`);
    seen.add(drug.id);
  }
});

test("每个模块均具备节点、关系与药物", () => {
  for (const id of moduleIds) {
    const modNodes = nodes.filter(n => n.modules.includes(id));
    const modEdges = edges.filter(e => e.modules.includes(id));
    const modDrugs = drugs.filter(d => d.modules.includes(id));
    assert(modNodes.length >= 15, `模块 ${id} 节点过少: ${modNodes.length}`);
    assert(modEdges.length >= 10, `模块 ${id} 关系过少: ${modEdges.length}`);
    assert(modDrugs.length >= 10, `模块 ${id} 药物过少: ${modDrugs.length}`);
    assert(modDrugs.some(d => mechanisms.has(d.mechanism)), `模块 ${id} 没有有效机制药物`);
  }
  for (const mech of mechanisms) assert(drugs.some(d => d.mechanism === mech), `全局缺少 ${mech} 药物`);
});

test("没有孤立节点", () => {
  const used = new Set();
  for (const edge of edges) {
    used.add(edge.source);
    used.add(edge.target);
  }
  const orphans = nodes.filter(node => !used.has(node.id)).map(node => node.id);
  assert.deepStrictEqual(orphans, [], `存在孤立节点: ${orphans.join(", ")}`);
});

test("已生成 index.html 且无模板占位符", () => {
  const out = path.join(root, "index.html");
  assert(fs.existsSync(out), "缺少 index.html，请先运行 node build.js");
  const html = fs.readFileSync(out, "utf8");
  assert(!html.includes("__IMMUNE_DATA__"), "index.html 仍包含 __IMMUNE_DATA__ 占位符");
  for (const marker of [
    'id="graph"', 'id="detail"', 'id="drugBody"', 'id="tableCount"', 'id="sourceList"',
    'id="snapshot"', 'id="snapshot2"', 'id="tabs"', 'id="search"', 'id="stats"',
    'id="relChips"', 'id="legend"', 'id="immune-data"',
    'rel="manifest" href="manifest.json"', 'rel="apple-touch-icon" href="icon-180.png"',
    'navigator.serviceWorker.register(swUrl, { scope })'
  ]) assert(html.includes(marker), `index.html 缺少关键标记 ${marker}`);
});

test("PWA 清单与图标完整且作用域受限", () => {
  const manifest = JSON.parse(read("manifest.json"));
  assert.strictEqual(manifest.start_url, "./index.html", "manifest.start_url 必须限制在本目录");
  assert.strictEqual(manifest.scope, "./", "manifest.scope 必须限制在本目录");
  assert.strictEqual(manifest.display, "standalone", "PWA 应使用 standalone 显示模式");
  assert(Array.isArray(manifest.icons) && manifest.icons.length, "manifest.icons 不能为空");
  assert(manifest.icons.some(icon => icon.sizes === "192x192" && icon.type === "image/png"), "缺少 192x192 PNG 图标");
  assert(manifest.icons.some(icon => icon.sizes === "512x512" && icon.type === "image/png"), "缺少 512x512 PNG 图标");
  assert(manifest.icons.some(icon => icon.sizes === "512x512" && String(icon.purpose).includes("maskable")), "缺少 maskable 图标");
  for (const icon of manifest.icons) assert(fs.existsSync(path.join(root, icon.src)), `图标文件不存在: ${icon.src}`);
  assert(fs.existsSync(path.join(root, "sw.js")), "缺少 sw.js");
});

test("SW 只拦截本目录请求", () => {
  const listeners = {};
  const scope = "https://example.com/tools/immune-target-map/";
  const cache = { match: async () => null, put: async () => {} };
  const context = {
    URL,
    Response,
    self: {
      registration: { scope },
      addEventListener: (type, fn) => { listeners[type] = fn; },
      skipWaiting: async () => {},
      clients: { claim: async () => {} }
    },
    caches: { open: async () => cache, keys: async () => [], delete: async () => true },
    fetch: async () => new Response("ok")
  };
  vm.runInNewContext(read("sw.js"), context);
  assert.strictEqual(typeof listeners.fetch, "function", "sw.js 未注册 fetch 事件");
  const dispatch = (url, method = "GET") => {
    let intercepted = false;
    listeners.fetch({ request: { url, method }, respondWith: () => { intercepted = true; }, waitUntil: () => {} });
    return intercepted;
  };
  assert.strictEqual(dispatch(scope + "index.html"), true, "SW 应处理本目录页面");
  assert.strictEqual(dispatch(scope + "index.html", "POST"), false, "SW 不应处理非 GET 请求");
  assert.strictEqual(dispatch("https://example.com/tools/protein-analyse/index.html"), false, "SW 不应拦截其他 PWA 目录");
  assert.strictEqual(dispatch("https://example.com/tools/immune-target-map-old/index.html"), false, "SW 不应拦截相似前缀目录");
  assert.strictEqual(dispatch("https://other.example/tools/immune-target-map/index.html"), false, "SW 不应拦截跨域请求");
});

test("构建产物内嵌数据与源文件一致", () => {
  const html = read("index.html");
  const match = html.match(/<script id="immune-data" type="application\/json">([\s\S]*?)<\/script>/);
  assert(match, "index.html 中找不到内嵌数据脚本");
  const embedded = JSON.parse(match[1]);
  assert.strictEqual(embedded.nodes.length, nodes.length, "内嵌节点数量不一致");
  assert.strictEqual(embedded.edges.length, edges.length, "内嵌关系数量不一致");
  assert.strictEqual(embedded.drugs.length, drugs.length, "内嵌药物数量不一致");
  assert.strictEqual(Object.keys(embedded.refs).length, refIds.size, "内嵌文献数量不一致");
  assert.strictEqual(embedded.meta.snapshot, meta.snapshot, "内嵌快照日期不一致");
  assert.deepStrictEqual(embedded.modules.map(m => m.id), meta.modules.map(m => m.id), "内嵌模块顺序不一致");
});

test("页面脚本语法有效", () => {
  const html = read("index.html");
  const start = html.lastIndexOf("<script>", html.indexOf("const DATA = JSON.parse"));
  const end = html.indexOf("</script>", start);
  assert(start >= 0 && end > start, "找不到页面主脚本");
  assert.doesNotThrow(() => new Function(html.slice(start + 8, end)), "页面主脚本存在语法错误");
});

test("页面无外部运行时依赖", () => {
  const html = read("index.html");
  assert(!/<script[^>]+src=/i.test(html), "页面不应引用外部脚本");
  assert(!/<link[^>]+rel=["']stylesheet["'][^>]+href=["']https?:/i.test(html), "页面不应引用外部样式表");
  assert(!/@import\s+url\(/i.test(html), "页面不应使用 CSS @import");
});

let failed = 0;
for (const { name, fn } of tests) {
  try {
    fn();
    console.log(`[ok] ${name}`);
  } catch (error) {
    failed += 1;
    console.error(`[fail] ${name}`);
    console.error(`       ${error.message}`);
  }
}

console.log(`\n${tests.length - failed}/${tests.length} 项通过`);
if (failed) process.exit(1);

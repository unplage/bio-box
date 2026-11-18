"use strict";

const fs = require("fs");
const path = require("path");

const root = __dirname;
const dataDir = path.join(root, "data");
const readJson = name => JSON.parse(fs.readFileSync(path.join(dataDir, name + ".json"), "utf8"));

const meta = readJson("meta");
const refs = readJson("refs");
const nodes = readJson("nodes");
const edges = readJson("edges");
const drugs = readJson("drugs");
const templatePath = path.join(root, "template.html");
const outputPath = path.join(root, "index.html");
const placeholder = "__IMMUNE_DATA__";

const missing = ["title", "snapshot", "disclaimer", "modules"].filter(key => meta[key] === undefined);
if (missing.length) throw new Error(`meta.json 缺少字段: ${missing.join(", ")}`);
if (!Array.isArray(meta.modules) || !meta.modules.length) throw new Error("meta.modules 必须是非空数组");
if (!templatePath || !fs.existsSync(templatePath)) throw new Error("找不到 template.html");

const template = fs.readFileSync(templatePath, "utf8");
const count = template.split(placeholder).length - 1;
if (count !== 1) throw new Error(`template.html 中占位符 ${placeholder} 应出现 1 次，实际 ${count} 次`);

const data = { meta, modules: meta.modules, refs, nodes, edges, drugs };
const json = JSON.stringify(data, null, 2).replace(/<\/script/gi, "<\\/script");
const html = template.replace(placeholder, json);

fs.writeFileSync(outputPath, html);
console.log(`已生成 ${path.relative(process.cwd(), outputPath)}`);
console.log(`节点 ${nodes.length} · 关系 ${edges.length} · 药物 ${drugs.length} · 文献 ${Object.keys(refs).length}`);

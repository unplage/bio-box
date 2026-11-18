# 免疫靶点关系图谱

单文件、离线可用的免疫靶点机制图谱，覆盖检查点/共刺激与细胞因子/趋化因子两个模块。

`index.html` 可直接用浏览器打开，不需要服务器、npm、CDN 或网络连接。

## 文件

| 文件 | 作用 |
|:---|:---|
| `index.html` | 构建产物，最终交付的单文件页面，请勿手工修改 |
| `template.html` | 页面模板，包含样式、交互和网络图逻辑 |
| `data/meta.json` | 标题、数据快照日期、免责声明、模块与分组定义 |
| `data/refs.json` | 文献与监管数据源，键名为引用 ID |
| `data/nodes.json` | 受体、配体、细胞因子、受体复合物与细胞效应节点 |
| `data/edges.json` | 结合、激活、抑制、竞争、表达调控等关系及亲和力/证据 |
| `data/drugs.json` | 药物、靶点、机制、临床阶段、关键试验与备注 |
| `manifest.json` | PWA 清单，作用域限制在本目录 |
| `sw.js` | 本目录专用 Service Worker，只缓存本目录资源 |
| `icon.svg` / `icon-*.png` | 页面图标与 PWA 安装图标 |
| `build.js` | 零依赖构建脚本，将数据和模板合并为 `index.html` |
| `tests/run.js` | 零依赖数据、构建产物与 PWA 作用域校验 |

## 使用方法

```bash
cd tools/immune-target-map
node build.js
node tests/run.js
```

然后直接打开 `index.html`：

```bash
xdg-open index.html
```

修改 `data/*.json` 或 `template.html` 后，必须重新运行 `node build.js`。浏览器只读取 `index.html` 中的内嵌数据副本。

通过 HTTP/HTTPS 部署后，首次在线访问会注册 `sw.js` 并缓存本页资源；之后断网仍可打开该 URL，浏览器也可按 PWA 方式安装。直接以 `file://` 打开时不会注册 Service Worker，但单文件本身仍可离线使用。

## PWA 作用域

- `manifest.json` 使用 `start_url: "./index.html"`、`scope: "./"`。
- `sw.js` 注册时显式限定到自身所在目录，并在 `fetch` 中再次检查路径前缀与同源。
- 本目录之外的请求会直接放行，不会拦截 `tools/protein-analyse/`、`tools/target-info/` 或其他目录。
- 修改缓存策略时，请同步更新 `sw.js` 中的 `CACHE` 版本号。

## 页面功能

- 两个模块切换：检查点/共刺激、细胞因子/趋化因子。
- 按关系类型、证据等级筛选，支持显示/隐藏药物节点与边标签。
- 搜索靶点、配体、药物或机制关键词。
- 点击节点查看表达、机制说明、关系、相关药物和文献。
- 点击关系查看亲和力、测定体系、证据等级与解读。
- 药物表可按状态、机制筛选，点击行定位到网络或查看详情。
- 导出当前网络 PNG、证据 CSV、药物 CSV。
- 网络支持缩放、拖拽和重置视图。

## 数据口径

- 关系类型：`bind`（结合）、`activate`（激活/共刺激）、`inhibit`（抑制）、`compete`（竞争）、`regulate`（表达调控）。
- 药物机制：`block`（阻断）或 `agonist`（激动）。
- 证据等级：
  - `A`：结构与生化、功能实验及人体/临床药物数据相互支持。
  - `B`：体外或动物功能证据较强，人体因果性有限。
  - `C`：相关性、单一体系或机制仍有争议。
- 亲和力必须结合 `assay` 字段阅读。二价融合蛋白、单体 SPR、细胞结合实验的结果不可直接横向比较。
- 药物监管状态与临床阶段是快照信息，实际以最新说明书、NMPA/FDA 公告和注册临床试验为准。

## 扩展数据

新增节点时：

1. 在 `data/nodes.json` 添加唯一 `id`，并填写 `label`、`kind`、`group`、`modules`、`expr`、`note`、`refs`。
2. 在 `data/edges.json` 添加关系，确保 `source`、`target` 均存在。
3. 在 `data/refs.json` 补充被引用的文献 ID。
4. 运行 `node build.js && node tests/run.js`。

新增药物时，`targets` 必须使用已存在的节点 ID；`modules` 必须属于 `meta.json` 中定义的模块。若想让药物出现在网络图中，将 `graph` 设为 `true` 并在 `data/edges.json` 中补充药物干预边。

## 免责声明

本工具用于机制梳理与靶点组合假设生成，不构成临床决策或用药建议。

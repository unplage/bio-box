# NVIDIA Bio Tools

Unified CLI and Python API for NVIDIA BioNeMo NIM models:

- **ProteinMPNN** (ipd/proteinmpnn): Predict amino acid sequences from protein backbones
- **Boltz-2** (mit/boltz2): Predict protein/protein-ligand/protein-nucleic acid structures

## ⚠️ Hosted API Limitations (Important)

> **关键提示**：NVIDIA 托管 API 为精简版，生产环境建议自托管。

### ProteinMPNN (托管版 vs 本地版)

| 功能 | 托管 API | 本地开源版 |
|------|----------|------------|
| 序列输出 (mfasta) | ✅ | ✅ |
| Per-position log-probs | ❌ | ✅ `.npz` |
| Per-position 概率 (20aa) | ❌ | ✅ |
| Per-position entropy | ❌ | ✅ |
| Sequence recovery | ❌ | ✅ |
| Tied positions (对称设计) | ❌ | ✅ |
| Fixed positions (固定位点) | ❌ | ✅ |
| 随机种子控制 | ❌ | ✅ `--seed` |
| Chain-specific design | ❌ | ✅ |
| Soluble model | ✅ | ✅ |

**托管 API 仅返回 mfasta 字符串，无任何置信度/概率/熵指标。**

### Boltz-2 (托管版 vs 自托管)

| 功能 | 托管 API | 自托管 NIM |
|------|----------|------------|
| 结构预测 | ✅ | ✅ |
| pLDDT 输出 | ❌ (仅 B-factor 归一化) | ✅ 完整 pLDDT |
| PAE/PDE 矩阵 | ❌ | ✅ 完整矩阵 |
| B-factor | 归一化 (0-1) | 实验级 (Å²) |
| 随机种子控制 | ❌ | ✅ `--seed` |
| 完整置信度分解 | ❌ | ✅ |
| 自定义模型权重 | ❌ | ✅ |

**托管 B-factor 已归一化 (0.7-1.0)，非实验 B-factor (Å²)，无法用于 pLDDT 阈值筛选。**

---

## Installation

```bash
cd tools/nvidia-bio
pip install -r requirements.txt
```

## Requirements

- Python 3.8+
- NVIDIA API key (from https://build.nvidia.com)
- For self-hosted: Docker with NVIDIA Container Toolkit

## Quick Start

### ProteinMPNN: Sequence design from backbone

```bash
# Using hosted NIM (requires API key)
python main.py proteinmpnn input.pdb --api-key $NVIDIA_API_KEY -o output.fa

# Self-hosted (no API key needed)
python main.py proteinmpnn input.pdb --self-hosted -o output.fa

# Dry run to see request
python main.py proteinmpnn input.pdb --dry-run
```

### Boltz-2: Structure prediction

```bash
# Single protein sequence
python main.py boltz2 --sequence "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN" --api-key $NVIDIA_API_KEY -o output/

# Multi-chain / ligand (via JSON)
python main.py boltz2 --chains chains.json --api-key $NVIDIA_API_KEY -o output/

# Self-hosted
python main.py boltz2 --sequence "..." --self-hosted -o output/
```

### chains.json format

```json
[
  {"id": "A", "molecule_type": "protein", "sequence": "MALWMRLLPLL..."},
  {"id": "B", "molecule_type": "protein", "sequence": "MKTAYIAKQRQ..."},
  {"id": "C", "molecule_type": "ligand", "ccd": "ATP", "smiles": "CC1=CC=CC=C1"}
]
```

## Python API

```python
from proteinmpnn import ProteinMPNNClient, predict_sequences
from boltz2 import Boltz2Client, predict_structure

# ProteinMPNN
client = ProteinMPNNClient(api_key="your-key")
result = client.predict_from_file("input.pdb")
print(result["mfasta"])

# Boltz-2
client = Boltz2Client(api_key="your-key")
result = client.predict_protein("MALWMRLLPLL...")
client.save_all_structures(result, "output/")

# Convenience functions
result = predict_sequences("input.pdb", api_key="key")
result = predict_structure("MALWMRLLPLL...", api_key="key")

# Multi-chain / ligand
chains = [
    {"id": "A", "molecule_type": "protein", "sequence": "..."},
    {"id": "B", "molecule_type": "ligand", "ccd": "ATP", "smiles": "..."}
]
result = predict_structure(chains=chains)
```

## API Key

Get your API key from https://build.nvidia.com

- Never hardcode keys in code
- Use `--api-key` parameter or `NVIDIA_API_KEY` environment variable
- Keys are never logged or saved to files

## Self-Hosted NIM

For production use, deploy self-hosted NIMs:

```bash
# ProteinMPNN
docker run --rm --gpus all -p 8000:8000 \
  -e NGC_API_KEY -v ~/.cache/nim:/opt/nim/.cache \
  nvcr.io/nim/ipd/proteinmpnn:1.0.0

# Boltz-2
docker run --rm --gpus all -p 8000:8000 \
  -e NGC_API_KEY -v ~/.cache/nim:/opt/nim/.cache --shm-size=16g \
  nvcr.io/nim/mit/boltz2:1.10.0
```

Then use `--self-hosted` flag (no API key needed):
```bash
python main.py proteinmpnn input.pdb --self-hosted
python main.py boltz2 --sequence "..." --self-hosted
```

## Endpoints

| Model | Hosted Endpoint | Self-Hosted Path |
|-------|----------------|------------------|
| ProteinMPNN | `https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict` | `/biology/ipd/proteinmpnn/predict` |
| Boltz-2 | `https://health.api.nvidia.com/v1/biology/mit/boltz2/predict` | `/biology/mit/boltz2/predict` |

## Dry Run

Use `--dry-run` to print the request payload without sending:

```bash
python main.py proteinmpnn input.pdb --dry-run
python main.py boltz2 --sequence "..." --dry-run
```

## Requirements

```
requests>=2.31.0
```

Install with:
```bash
pip install -r requirements.txt
```
"""Read confidence metrics from AF3-like predictors (Boltz-2, AlphaFold3, ...).

Different tools name their confidence JSON differently:
- Boltz-2: ``confidence_<job>_model_0.json`` with ``ptm`` / ``iptm`` /
  ``complex_plddt`` (0-1 scale).
- AlphaFold3: ``<job>/summary_confidences.json`` with ``ptm`` / ``iptm`` /
  ``ranking_score``.
Only pTM/ipTM/ranking score are taken here; pLDDT is read from the structure
B-factor so it works even without a JSON sidecar.
"""
import json
import re
from pathlib import Path

_ALIASES = {
    'iptm': ('iptm', 'interface_ptm', 'interface_iptm'),
    'ptm': ('ptm',),
    'ranking_score': ('ranking_score', 'confidence_score'),
}


def _first(mapping, keys):
    for key in keys:
        value = mapping.get(key)
        if value is None:
            continue
        try:
            return float(value)
        except (TypeError, ValueError):
            continue
    return None


def _job_base(stem):
    return re.split(r'[._-](?:unrelaxed|relaxed|scores|model|rank|pred)',
                    stem)[0]


def _json_candidates(model_path):
    model_path = Path(model_path)
    stem = model_path.stem
    base = _job_base(stem)
    patterns = ('confidence_*.json', 'summary_confidences*.json',
                '*confidences*.json', '*scores*.json')
    found = []
    for folder in (model_path.parent, model_path.parent.parent):
        if not folder.exists():
            continue
        for pattern in patterns:
            for candidate in folder.glob(pattern):
                name = candidate.name
                if base and base in name:
                    rank = 0
                elif stem in name:
                    rank = 1
                else:
                    rank = 2
                found.append((rank, candidate))
    found.sort(key=lambda item: (item[0], item[1].name))
    # Prefer JSONs that share the model's job base; never grab an unrelated
    # design's scores from a shared directory.
    exact = [c for r, c in found if r <= 1]
    if exact:
        return exact
    best = found[0][0] if found else 0
    return [c for r, c in found if r == best]


def load_confidence(model_path):
    """Return {'ptm','iptm','ranking_score'} (may be None) for a model file."""
    result = {'ptm': None, 'iptm': None, 'ranking_score': None,
              'source': None}  # type: dict
    for candidate in _json_candidates(model_path):
        try:
            payload = json.loads(candidate.read_text())
        except (OSError, ValueError):
            continue
        if not isinstance(payload, dict):
            continue
        for key, aliases in _ALIASES.items():
            value = _first(payload, aliases)
            if value is not None:
                result[key] = value
        result['source'] = str(candidate)
        if any(result[k] is not None for k in _ALIASES):
            break
    return result

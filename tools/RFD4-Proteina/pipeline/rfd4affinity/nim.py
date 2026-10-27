# -*- coding: utf-8 -*-
"""RFD4-Proteina NIM client."""

from __future__ import annotations

import json
from pathlib import Path

from .net import post_json
from .util import PipelineError, log, read_json, require

RFD4_PATH = "/biology/ipd/rfd4-proteina/generative"


def endpoint(base_url: str) -> str:
    return base_url.rstrip("/") + RFD4_PATH


def run_design(
    request: dict,
    api_key: str,
    base_url: str = "https://health.api.nvidia.com/v1",
    timeout: float = 1800.0,
) -> dict:
    require(bool(api_key), "NVIDIA_API_KEY is required (trial keys start with 'nvapi-')")
    log(
        f"Calling RFD4-Proteina: {request.get('sampling', {}).get('num_replicates', 1)} "
        f"design(s), strategy={request.get('sampling', {}).get('recycle_strategy')}"
    )
    response = post_json(endpoint(base_url), request, api_key=api_key, timeout=timeout)
    if "designs" not in response:
        raise PipelineError(
            "Unexpected response (no 'designs'): "
            + json.dumps(response, ensure_ascii=False)[:1500]
        )
    designs = response.get("designs") or []
    log(f"Received {len(designs)} design(s) in {response.get('duration_ms', 0) / 1000:.1f}s")
    for warning in response.get("warnings") or []:
        log(f"service warning: {warning}", "WARN")
    return response


def load_replay(path: Path) -> dict:
    require(path.exists(), f"--replay file not found: {path}")
    data = read_json(path)
    require("designs" in data, f"{path} does not look like a NIM response (no 'designs')")
    return data

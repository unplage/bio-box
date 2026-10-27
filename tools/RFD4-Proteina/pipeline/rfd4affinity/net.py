# -*- coding: utf-8 -*-
"""Tiny HTTP helper with retries.  Standard library only (no requests)."""

from __future__ import annotations

import json
import random
import re
import time
import urllib.error
import urllib.request

from .util import PipelineError, warn

USER_AGENT = "rfd4affinity/0.1 (+https://github.com/unplage/bio-box)"


def _retry_after(value: str | None) -> float | None:
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _sleep(attempt: int, retry_after: float | None = None) -> None:
    if retry_after is not None:
        delay = min(max(retry_after, 1.0), 120.0)
    else:
        delay = min(2.0 ** attempt + random.random(), 60.0)
    time.sleep(delay)


def get_text(url: str, timeout: float = 60.0, retries: int = 3) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code in (429, 500, 502, 503, 504) and attempt + 1 < retries:
                warn(f"HTTP {exc.code} from {url}; retrying")
                retry_after = _retry_after(exc.headers.get("Retry-After"))
                _sleep(attempt, retry_after)
                continue
            raise PipelineError(f"HTTP {exc.code} fetching {url}: {exc.reason}") from exc
        except urllib.error.URLError as exc:
            last_error = exc
            if attempt + 1 < retries:
                warn(f"Network error fetching {url}: {exc.reason}; retrying")
                _sleep(attempt)
                continue
            raise PipelineError(f"Network error fetching {url}: {exc.reason}") from exc
        except OSError as exc:
            last_error = exc
            if attempt + 1 < retries:
                warn(f"Network error fetching {url}: {exc}; retrying")
                _sleep(attempt)
                continue
            raise PipelineError(f"Network error fetching {url}: {exc}") from exc
    raise PipelineError(f"Failed to fetch {url}: {last_error}")


def _is_fatal_server_error(detail: str) -> bool:
    """Detect deterministic backend crashes that retrying cannot fix."""
    lowered = detail.lower()
    return any(
        marker in lowered
        for marker in (
            "cuda error",
            "acceleratorerror",
            "illegal memory access",
            "out of memory",
            "outofmemory",
        )
    )


def _metric_tag_hint(detail: str) -> str:
    """Extract the offending metric tags from a 422 payload."""
    tags = sorted(set(re.findall(r"'([a-z0-9_]+)'", detail)))
    if not tags:
        return ""
    return (
        " The server rejected metric tag(s) "
        + ", ".join(tags)
        + "; drop them with --metrics-tags (for example --metrics-tags design)."
    )


def post_json(
    url: str,
    payload: dict,
    api_key: str | None = None,
    timeout: float = 900.0,
    retries: int = 4,
) -> dict:
    body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": USER_AGENT,
    }
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    last_error: Exception | None = None
    for attempt in range(retries):
        request = urllib.request.Request(url, data=body, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                raw = response.read().decode("utf-8", errors="replace")
            try:
                return json.loads(raw)
            except json.JSONDecodeError as exc:
                raise PipelineError(
                    f"Server returned non-JSON ({len(raw)} bytes): {raw[:400]}"
                ) from exc
        except urllib.error.HTTPError as exc:
            last_error = exc
            try:
                detail = exc.read().decode("utf-8", errors="replace")[:1200]
            except Exception:
                detail = ""
            if exc.code == 422 and "metric tag" in detail:
                raise PipelineError(
                    f"HTTP 422 from {url}: the deployment rejected the requested "
                    f"metric tag(s).{_metric_tag_hint(detail)}"
                ) from exc
            if exc.code >= 500 and _is_fatal_server_error(detail):
                raise PipelineError(
                    f"HTTP {exc.code} from {url}: the NIM backend crashed while "
                    "running the model "
                    f"({detail[:300]}). This is a server-side failure; retrying "
                    "will not help. Check the service status or use a self-hosted "
                    "NIM (nvcr.io/nim/ipd/rfd4-proteina)."
                ) from exc
            if exc.code in (429, 500, 502, 503, 504) and attempt + 1 < retries:
                retry_after = _retry_after(exc.headers.get("Retry-After"))
                warn(
                    f"HTTP {exc.code} (attempt {attempt + 1}/{retries}) from "
                    f"{url.rsplit('/', 1)[-1]}; backing off"
                )
                _sleep(attempt, retry_after)
                continue
            hint = ""
            if exc.code in (401, 403):
                hint = " (check NVIDIA_API_KEY; trial keys start with 'nvapi-')"
            raise PipelineError(
                f"HTTP {exc.code} from {url}{hint}: {detail or exc.reason}"
            ) from exc
        except urllib.error.URLError as exc:
            last_error = exc
            if attempt + 1 < retries:
                warn(f"Network error (attempt {attempt + 1}/{retries}): {exc.reason}")
                _sleep(attempt)
                continue
            raise PipelineError(f"Network error calling {url}: {exc.reason}") from exc
        except OSError as exc:
            last_error = exc
            if attempt + 1 < retries:
                warn(f"Network error (attempt {attempt + 1}/{retries}): {exc}")
                _sleep(attempt)
                continue
            raise PipelineError(f"Network error calling {url}: {exc}") from exc
    raise PipelineError(f"Request failed after {retries} attempts: {last_error}")

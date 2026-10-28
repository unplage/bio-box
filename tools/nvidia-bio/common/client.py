"""Common NVIDIA NIM API client with authentication and retry logic."""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any

import requests


class NVIDIAAPIError(Exception):
    """Exception raised for NVIDIA API errors."""

    def __init__(self, message: str, status_code: int | None = None, response: Any = None):
        super().__init__(message)
        self.status_code = status_code
        self.response = response


class NVIDIAClient:
    """Base client for NVIDIA NIM API calls."""

    def __init__(
        self,
        base_url: str,
        api_key: str | None = None,
        timeout: float = 300.0,
        max_retries: int = 3,
        retry_backoff: float = 2.0,
    ):
        """
        Initialize the NVIDIA API client.

        Args:
            base_url: Base URL for the NIM endpoint (e.g., "https://health.api.nvidia.com/v1")
            api_key: NVIDIA API key (can also be set via NVIDIA_API_KEY env var)
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
            retry_backoff: Exponential backoff factor for retries
        """
        self.base_url = base_url.rstrip("/")
        self.api_key = api_key or os.getenv("NVIDIA_API_KEY")
        self.timeout = timeout
        self.max_retries = max_retries
        self.retry_backoff = retry_backoff

        if not self.api_key:
            raise ValueError(
                "NVIDIA API key is required. Pass api_key parameter or set NVIDIA_API_KEY environment variable."
            )

        self.session = requests.Session()
        self.session.headers.update(
            {
                "Authorization": f"Bearer {self.api_key}",
                "Content-Type": "application/json",
                "Accept": "application/json",
            }
        )

    def _request(
        self,
        method: str,
        endpoint: str,
        json_data: dict | None = None,
        params: dict | None = None,
        headers: dict | None = None,
        stream: bool = False,
    ) -> requests.Response:
        """Make an HTTP request with retry logic."""
        url = f"{self.base_url}{endpoint}"
        request_headers = dict(self.session.headers)
        if headers:
            request_headers.update(headers)

        last_exception = None
        for attempt in range(self.max_retries + 1):
            try:
                response = self.session.request(
                    method=method,
                    url=url,
                    json=json_data,
                    params=params,
                    headers=request_headers,
                    timeout=self.timeout,
                    stream=stream,
                )

                # Handle different status codes
                if response.status_code == 200:
                    return response

                if response.status_code == 401:
                    raise NVIDIAAPIError(
                        "Authentication failed: Invalid or expired API key",
                        status_code=401,
                        response=response,
                    )
                if response.status_code == 403:
                    raise NVIDIAAPIError(
                        "Access forbidden: API key may not have required permissions",
                        status_code=403,
                        response=response,
                    )
                if response.status_code == 422:
                    try:
                        error_detail = response.json()
                    except Exception:
                        error_detail = response.text
                    raise NVIDIAAPIError(
                        f"Validation error: {error_detail}",
                        status_code=422,
                        response=response,
                    )
                if response.status_code == 429:
                    # Rate limited - retry with backoff
                    retry_after = response.headers.get("Retry-After")
                    wait_time = float(retry_after) if retry_after else self.retry_backoff ** attempt
                    if attempt < self.max_retries:
                        time.sleep(wait_time)
                        continue
                    raise NVIDIAAPIError(
                        f"Rate limited after {self.max_retries} retries",
                        status_code=429,
                        response=response,
                    )
                if 500 <= response.status_code < 600:
                    # Server error - retry with backoff
                    if attempt < self.max_retries:
                        wait_time = self.retry_backoff ** attempt
                        time.sleep(wait_time)
                        continue
                    raise NVIDIAAPIError(
                        f"Server error: {response.status_code} - {response.text[:500]}",
                        status_code=response.status_code,
                        response=response,
                    )

                # Other client errors - don't retry
                raise NVIDIAAPIError(
                    f"Request failed: {response.status_code} - {response.text[:500]}",
                    status_code=response.status_code,
                    response=response,
                )

            except requests.exceptions.Timeout:
                last_exception = NVIDIAAPIError(f"Request timeout after {self.timeout}s")
            except requests.exceptions.ConnectionError as e:
                last_exception = NVIDIAAPIError(f"Connection error: {e}")
            except requests.exceptions.RequestException as e:
                last_exception = NVIDIAAPIError(f"Request failed: {e}")

            if attempt < self.max_retries:
                wait_time = self.retry_backoff ** attempt
                time.sleep(wait_time)
            else:
                raise last_exception or NVIDIAAPIError("Max retries exceeded")

        raise last_exception or NVIDIAAPIError("Max retries exceeded")

    def get(self, endpoint: str, params: dict | None = None) -> requests.Response:
        return self._request("GET", endpoint, params=params)

    def post(self, endpoint: str, json_data: dict | None = None) -> requests.Response:
        return self._request("POST", endpoint, json_data=json_data)

    def health_check(self, endpoint: str = "/v1/health/ready") -> bool:
        """Check if the NIM endpoint is ready."""
        try:
            response = self.get(endpoint)
            return response.status_code == 200 and response.json() is True
        except Exception:
            return False

    def wait_for_ready(self, endpoint: str = "/v1/health/ready", timeout: float = 300.0, poll_interval: float = 5.0) -> bool:
        """Wait for the NIM endpoint to become ready."""
        start_time = time.time()
        while time.time() - start_time < timeout:
            if self.health_check(endpoint):
                return True
            time.sleep(poll_interval)
        return False


def get_api_key(api_key: str | None = None) -> str:
    """Get API key from parameter or environment variable."""
    key = api_key or os.getenv("NVIDIA_API_KEY")
    if not key:
        raise ValueError(
            "NVIDIA API key is required. Pass api_key parameter or set NVIDIA_API_KEY environment variable."
        )
    return key
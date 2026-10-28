"""ProteinMPNN client for NVIDIA NIM API."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

# Add parent directory to path for common module
sys.path.insert(0, str(Path(__file__).parent.parent))

from common.client import NVIDIAClient, NVIDIAAPIError, get_api_key


class ProteinMPNNClient(NVIDIAClient):
    """Client for ProteinMPNN NIM API."""

    DEFAULT_HOSTED_URL = "https://health.api.nvidia.com/v1"
    DEFAULT_SELF_HOSTED_URL = "http://localhost:8000"
    PREDICT_ENDPOINT = "/biology/ipd/proteinmpnn/predict"
    HEALTH_ENDPOINT = "/v1/health/ready"

    def __init__(
        self,
        api_key: str | None = None,
        base_url: str | None = None,
        use_self_hosted: bool = False,
        timeout: float = 300.0,
        max_retries: int = 3,
    ):
        """
        Initialize ProteinMPNN client.

        Args:
            api_key: NVIDIA API key (or set NVIDIA_API_KEY env var)
            base_url: Custom base URL (overrides use_self_hosted)
            use_self_hosted: Use self-hosted endpoint (http://localhost:8000)
            timeout: Request timeout in seconds
            max_retries: Maximum retry attempts
            retry_backoff: Exponential backoff factor
        """
        if base_url is None:
            base_url = self.DEFAULT_SELF_HOSTED_URL if use_self_hosted else self.DEFAULT_HOSTED_URL

        super().__init__(
            base_url=base_url,
            api_key=api_key,
            timeout=timeout,
            max_retries=max_retries,
            retry_backoff=2.0,
        )

    def predict(
        self,
        pdb_content: str,
        ca_only: bool = False,
        use_soluble_model: bool = False,
        sampling_temp: list[float] | None = None,
        num_sequences: int = 1,
        chain_ids: list[str] | None = None,
    ) -> dict:
        """
        Predict amino acid sequences for a protein backbone.

        Args:
            pdb_content: PDB file content (ATOM records only, max ~200 lines recommended)
            ca_only: Use only CA atoms (faster, less accurate)
            use_soluble_model: Use soluble protein model variant
            sampling_temp: List of sampling temperatures (default: [0.1])
            num_sequences: Number of sequences to generate per temperature
            chain_ids: Specific chain IDs to design (default: all)

        Returns:
            Dictionary with 'mfasta' (multi-FASTA) and metadata
        """
        if sampling_temp is None:
            sampling_temp = [0.1]

        # Reduce PDB to first ~200 ATOM lines if too large
        lines = pdb_content.strip().split("\n")
        atom_lines = [line for line in lines if line.startswith("ATOM")]
        if len(atom_lines) > 200:
            atom_lines = atom_lines[:200]
        reduced_pdb = "\n".join(atom_lines)

        payload = {
            "input_pdb": reduced_pdb,
            "ca_only": ca_only,
            "use_soluble_model": use_soluble_model,
            "sampling_temp": sampling_temp,
        }

        response = self.post(self.PREDICT_ENDPOINT, json_data=payload)
        return response.json()

    def predict_from_file(
        self,
        pdb_path: str | Path,
        ca_only: bool = False,
        use_soluble_model: bool = False,
        sampling_temp: list[float] | None = None,
        num_sequences: int = 1,
        max_atom_lines: int = 200,
    ) -> dict:
        """Predict sequences from a PDB file."""
        path = Path(pdb_path)
        if not path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")

        content = path.read_text(encoding="utf-8")
        lines = content.strip().split("\n")
        atom_lines = [line for line in lines if line.startswith("ATOM")]
        if len(atom_lines) > max_atom_lines:
            # Keep first max_atom_lines ATOM records
            kept = []
            count = 0
            for line in lines:
                if line.startswith("ATOM"):
                    if count >= max_atom_lines:
                        continue
                    count += 1
                kept.append(line)
            pdb_content = "\n".join(kept)
        else:
            pdb_content = "\n".join(atom_lines)

        return self.predict(
            pdb_content=pdb_content,
            ca_only=ca_only,
            use_soluble_model=use_soluble_model,
            sampling_temp=sampling_temp,
        )

    

    def save_mfasta(self, result: dict, output_path: str | Path) -> Path:
        """Save multi-FASTA output from result."""
        output_path = Path(output_path)
        mfasta = result.get("mfasta", "")
        if not mfasta:
            raise ValueError("No 'mfasta' field in result")
        output_path.write_text(mfasta, encoding="utf-8")
        return output_path


def predict_sequences(
    pdb_path: str | Path,
    api_key: str | None = None,
    base_url: str | None = None,
    use_self_hosted: bool = False,
    **kwargs,
) -> dict:
    """
    Convenience function to predict sequences from a PDB file.

    Args:
        pdb_path: Path to PDB file
        api_key: NVIDIA API key (or set NVIDIA_API_KEY env var)
        base_url: Custom base URL
        use_self_hosted: Use self-hosted endpoint
        **kwargs: Additional arguments passed to predict_from_file

    Returns:
        Prediction result dictionary
    """
    api_key = get_api_key(api_key)
    client = ProteinMPNNClient(
        api_key=api_key,
        base_url=base_url,
        use_self_hosted=use_self_hosted,
    )
    return client.predict_from_file(pdb_path, **kwargs)
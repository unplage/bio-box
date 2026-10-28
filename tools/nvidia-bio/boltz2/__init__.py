"""Boltz-2 client for NVIDIA NIM API."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

# Add parent directory to path for common module
sys.path.insert(0, str(Path(__file__).parent.parent))

from common.client import NVIDIAClient, NVIDIAAPIError, get_api_key


class Boltz2Client(NVIDIAClient):
    """Client for Boltz-2 NIM API."""

    DEFAULT_HOSTED_URL = "https://health.api.nvidia.com/v1"
    DEFAULT_SELF_HOSTED_URL = "http://localhost:8000"
    PREDICT_ENDPOINT = "/biology/mit/boltz2/predict"
    HEALTH_ENDPOINT = "/v1/health/ready"

    def __init__(
        self,
        api_key: str | None = None,
        base_url: str | None = None,
        use_self_hosted: bool = False,
        timeout: float = 600.0,
        max_retries: int = 3,
    ):
        """
        Initialize Boltz-2 client.

        Args:
            api_key: NVIDIA API key (or set NVIDIA_API_KEY env var)
            base_url: Custom base URL (overrides use_self_hosted)
            use_self_hosted: Use self-hosted endpoint (http://localhost:8000)
            timeout: Request timeout in seconds (longer for structure prediction)
            max_retries: Maximum retry attempts
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
        polymers: list[dict],
        recycling_steps: int = 3,
        sampling_steps: int = 50,
        diffusion_samples: int = 1,
        step_scale: float = 1.638,
        output_format: str = "mmcif",
        msa: dict | None = None,
        ligands: list[dict] | None = None,
        templates: list[dict] | None = None,
        constraints: list[dict] | None = None,
        properties: list[str] | None = None,
    ) -> dict:
        """
        Predict protein structure(s).

        Args:
            polymers: List of polymer dictionaries, each with:
                - id: Chain identifier (e.g., "A", "B")
                - molecule_type: "protein", "rna", "dna"
                - sequence: Amino acid / nucleotide sequence
            recycling_steps: Number of recycling steps (default 3)
            sampling_steps: Number of diffusion sampling steps (default 50)
            diffusion_samples: Number of diffusion samples (default 1)
            step_scale: Step scale factor (default 1.638)
            output_format: Output format - "mmcif" or "pdb" (default "mmcif")
            msa: Optional MSA configuration dict
            ligands: Optional list of ligand dictionaries
            templates: Optional list of template dictionaries
            constraints: Optional constraint dictionaries
            properties: List of properties to compute (e.g., ["affinity", "pae", "pde"])

        Returns:
            Prediction result dictionary with structures, confidence scores, etc.
        """
        payload = {
            "polymers": polymers,
            "recycling_steps": recycling_steps,
            "sampling_steps": sampling_steps,
            "diffusion_samples": diffusion_samples,
            "step_scale": step_scale,
            "output_format": output_format,
        }

        if msa is not None:
            payload["msa"] = msa
        if ligands is not None:
            payload["ligands"] = ligands
        if templates is not None:
            payload["templates"] = templates
        if constraints is not None:
            payload["constraints"] = constraints
        if properties is not None:
            payload["properties"] = properties

        response = self.post(self.PREDICT_ENDPOINT, json_data=payload)
        return response.json()

    def predict_protein(
        self,
        sequence: str,
        chain_id: str = "A",
        recycling_steps: int = 3,
        sampling_steps: int = 50,
        diffusion_samples: int = 1,
        step_scale: float = 1.638,
        output_format: str = "mmcif",
        **kwargs,
    ) -> dict:
        """Convenience method for single protein structure prediction."""
        polymers = [{
            "id": chain_id,
            "molecule_type": "protein",
            "sequence": sequence,
        }]
        return self.predict(polymers=polymers, **kwargs)

    def predict_complex(
        self,
        chains: list[dict],
        **kwargs,
    ) -> dict:
        """Predict structure for a multi-chain complex."""
        return self.predict(polymers=chains, **kwargs)

    def predict_with_ligand(
        self,
        protein_sequence: str,
        ligand_ccd: str,
        ligand_smiles: str | None = None,
        protein_chain_id: str = "A",
        ligand_chain_id: str = "B",
        **kwargs,
    ) -> dict:
        """Predict protein-ligand complex structure."""
        polymers = [{"id": protein_chain_id, "molecule_type": "protein", "sequence": protein_sequence}]
        ligands = [{"id": ligand_chain_id, "ccd": ligand_ccd}]
        if ligand_smiles:
            ligands[0]["smiles"] = ligand_smiles
        return self.predict(polymers=polymers, ligands=ligands, **kwargs)

    

    def save_structure(self, result: dict, output_path: str | Path, index: int = 0) -> Path:
        """Save predicted structure to file."""
        output_path = Path(output_path)
        structures = result.get("structures", [])
        if not structures or index >= len(structures):
            raise ValueError(f"No structure at index {index}")
        structure = structures[index]
        content = structure.get("structure", "")
        if not content:
            raise ValueError(f"No structure content at index {index}")
        output_path.write_text(content, encoding="utf-8")
        return output_path

    def save_all_structures(self, result: dict, output_dir: str | Path) -> list[Path]:
        """Save all predicted structures to directory."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        structures = result.get("structures", [])
        paths = []
        for i, structure in enumerate(structures):
            content = structure.get("structure", "")
            if content:
                fmt = structure.get("format", "mmcif")
                ext = ".cif" if fmt == "mmcif" else ".pdb"
                path = output_dir / f"structure_{i}{ext}"
                path.write_text(content, encoding="utf-8")
                paths.append(path)
        return paths


def predict_structure(
    sequence: str | None = None,
    chains: list[dict] | None = None,
    api_key: str | None = None,
    base_url: str | None = None,
    use_self_hosted: bool = False,
    **kwargs,
) -> dict:
    """
    Predict structure for a protein or complex.

    Args:
        sequence: Protein amino acid sequence
        chains: List of chain dicts for multi-chain/ligand (alternative to sequence)
        api_key: NVIDIA API key (or set NVIDIA_API_KEY env var)
        base_url: Custom base URL
        use_self_hosted: Use self-hosted endpoint
        **kwargs: Additional arguments passed to predict/predict_protein

    Returns:
        Prediction result dictionary
    """
    if sequence is None and chains is None:
        raise ValueError("Either sequence or chains must be provided")

    api_key = get_api_key(api_key)
    client = Boltz2Client(
        api_key=api_key,
        base_url=base_url,
        use_self_hosted=use_self_hosted,
    )

    if chains is not None:
        return client.predict(polymers=chains, **kwargs)
    else:
        assert sequence is not None
        return client.predict_protein(sequence, **kwargs)
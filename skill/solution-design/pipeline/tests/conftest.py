"""Shared fixtures: load STEAP1 8UCD template, harmonize, reference-number.

Run with: python -m pytest tests
"""
import sys
from pathlib import Path

import yaml
import pytest

PIPELINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE))

from tmkit.pdbio import parse_pdb, harmonize_chains, reference_numbering  # noqa: E402
from tmkit.mask import build_mask  # noqa: E402


@pytest.fixture(scope='session')
def cfg():
    return yaml.safe_load((PIPELINE / 'config.steap1.yaml').read_text())


@pytest.fixture(scope='session')
def target(cfg):
    t = parse_pdb(PIPELINE / 'pdb' / '8UCD.pdb',
                  chain_ids=cfg['chains'])
    harmonize_chains(t, cfg['chains'])
    for c in cfg['chains']:
        reference_numbering(t, c, cfg['reference_sequences'][c])
    return t


@pytest.fixture(scope='session')
def ref_seq(cfg):
    return cfg['reference_sequences']['A']


@pytest.fixture(scope='session')
def mask(target, cfg):
    return build_mask(target, cfg['chains'], cfg)


TM = [[71, 91], [119, 139], [164, 184], [218, 238], [258, 278], [291, 311]]
ECD = [[92, 118], [185, 217], [279, 290]]
"""Membrane geometry for lipid-facing residue selection.

Two ways to define the membrane slab:

* explicit normal / center / thickness::

    membrane:
      normal: [0, 0, 1]     # or "auto" (PCA of TM C-alphas)
      thickness: 30.0        # hydrophobic core (A)
      headgroup: 15.0        # extra interface/headgroup zone on each side
      center: auto           # or a number (projection along the normal)
      exclude_aromatic_belt: true
      restrict_to_core: false

* an external OPM/PPM transformed structure (membrane plane is z=0, normal +z)::

    membrane_opm_pdb: /path/to/opm_transformed.pdb

When neither is given, no membrane restriction is applied and the mask keeps
its SASA-only behaviour.
"""
import numpy as np

from .pdbio import ca_coords

AROMATIC = {'TRP', 'PHE', 'TYR', 'W', 'F', 'Y'}


def _unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    return v / n if n > 0 else np.array([0.0, 0.0, 1.0])


def _tm_coords(target, chains, tm):
    pts = []
    for cid in chains:
        coords, refs = ca_coords(target, cid)
        for coord, refn in zip(coords, refs):
            if not tm or any(a <= refn <= b for a, b in tm):
                pts.append(coord)
    return np.asarray(pts, dtype=float) if pts else np.zeros((0, 3))


def _auto_normal(pts):
    """Principal axis of the TM bundle = membrane normal."""
    if len(pts) < 3:
        return np.array([0.0, 0.0, 1.0])
    centered = pts - pts.mean(axis=0)
    try:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
    except np.linalg.LinAlgError:
        return np.array([0.0, 0.0, 1.0])
    return _unit(vh[0])


def membrane_geometry(target, chains, config):
    """Return a geometry dict, or None when membrane handling is disabled."""
    opm = config.get('membrane_opm_pdb') or config.get('opm_pdb')
    mem = config.get('membrane')
    if opm:
        mem = dict(mem or {})
        mem.setdefault('normal', [0, 0, 1])
        mem.setdefault('center', 0.0)
    if mem is False or mem is None and not opm and not config.get('membrane_auto'):
        return None
    if mem is None:
        mem = {}
    if not isinstance(mem, dict):
        mem = {}
    tm = [(int(a), int(b)) for a, b in (config.get('tm_segments') or [])]
    pts = _tm_coords(target, chains, tm)

    normal_spec = mem.get('normal', 'auto' if config.get('membrane_auto') else None)
    if normal_spec is None or (isinstance(normal_spec, str)
                               and normal_spec.lower() == 'auto'):
        normal = _auto_normal(pts)
    else:
        normal = _unit(normal_spec)

    center_spec = mem.get('center', 'auto')
    if isinstance(center_spec, str) and center_spec.lower() == 'auto':
        center = float(np.dot(pts, normal).mean()) if len(pts) else 0.0
    else:
        center = float(center_spec)

    core = float(mem.get('thickness', 30.0)) / 2.0
    head = float(mem.get('headgroup', 15.0))
    return {
        'normal': normal,
        'center': center,
        'core': core,
        'head': head,
        'exclude_aromatic_belt': bool(mem.get('exclude_aromatic_belt', True)),
        'restrict_to_core': bool(mem.get('restrict_to_core', False)),
    }


def signed_distance(coord, geom):
    return float(np.dot(np.asarray(coord, dtype=float), geom['normal'])
                 - geom['center'])


def membrane_exclusion(coord, geom, resname):
    """Return an exclusion reason, or None if the residue may be designed."""
    d = abs(signed_distance(coord, geom))
    if d > geom['core'] + geom['head']:
        return 'outside_membrane'
    if geom['restrict_to_core'] and d > geom['core']:
        return 'headgroup'
    if (geom['exclude_aromatic_belt'] and d > geom['core']
            and resname.upper() in AROMATIC):
        return 'aromatic_belt'
    return None

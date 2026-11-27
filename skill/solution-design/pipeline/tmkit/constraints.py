"""Design constraints shared by the mask and evaluation layers.

Constraints are expressed in reference numbering (1-based) and can be given
globally (homomers) or per chain (heteromers).  The builder turns a target
config into:

* fixed regions (ECD by default, always enforced when ``fix_ecd`` is true),
* fixed residues (user seeds + protection pads + pockets + disulfides + PTMs
  + epitopes + pore residues),
* interface handling (all fixed, core fixed / edge free, or none),
* membrane/pocket/pore shells used by the report.

Nothing here knows about STEAP1; the same rules work for a GPCR, a single-pass
receptor, a beta-barrel or a soluble control.
"""
import re

import numpy as np

from .pdbio import AMINO3, ca_coords, hetero_residues


def symmetry_groups(config, chain_ids):
    """Return explicit chain-equivalence groups for ties / symmetry metrics.

    ``symmetry.groups`` wins.  Without it, identical reference sequences are
    treated as one homomer group; heteromers get no implicit groups.  This
    replaces the old "designable sets happen to match" heuristic, which tied
    unrelated positions together.
    """
    sym = config.get('symmetry') or {}
    groups = sym.get('groups') if isinstance(sym, dict) else None
    if groups:
        out = []
        for group in groups:
            if isinstance(group, str):
                group = [group]
            keep = [c for c in group if c in chain_ids]
            if len(keep) > 1:
                out.append(keep)
        return out
    refs = config.get('reference_sequences') or {}
    seqs = [refs.get(c) for c in chain_ids]
    if len(chain_ids) > 1 and all(s and s == seqs[0] for s in seqs):
        return [list(chain_ids)]
    return []


def _per_chain(value, chain_ids):
    """Normalize a config value to {chain: set(refnums)}.

    Accepts a list of reference positions (applied to all chains), a mapping
    ``chain -> list``, or ``None``.  PDB chain ids are single letters but the
    helper does not assume that.
    """
    out = {c: set() for c in chain_ids}
    if value is None:
        return out
    if isinstance(value, dict):
        for c in chain_ids:
            vals = value.get(c, value.get('*'))
            if vals:
                out[c] = set(int(x) for x in vals)
        return out
    vals = {int(x) for x in value}
    for c in chain_ids:
        out[c] = set(vals)
    return out


def _ref_of(target, chain_id, pdb_resnum):
    return target.reference_map.get(chain_id, {}).get(pdb_resnum, pdb_resnum)


def detect_disulfides(target, chain_ids, cutoff=2.5):
    """Return reference-numbered Cys pairs whose SG atoms are bonded.

    Only residues with an SG atom are considered, so free cysteines are not
    silently fixed.  The result is a list of ``(chain, refn, chain, refn)``
    tuples; each Cys is fixed by the mask.
    """
    sg = []
    for cid in chain_ids:
        for residue in target.chains[cid].residues:
            if residue.resname != 'CYS':
                continue
            try:
                atom = residue['SG']
            except KeyError:
                continue
            sg.append((cid, _ref_of(target, cid, residue.id[1]),
                       atom.get_coord()))
    pairs = []
    cutoff2 = float(cutoff) ** 2
    for i in range(len(sg)):
        for j in range(i + 1, len(sg)):
            if ((sg[i][2] - sg[j][2]) ** 2).sum() <= cutoff2:
                pairs.append((sg[i][0], sg[i][1], sg[j][0], sg[j][1]))
    return pairs


def detect_ptm_sequons(seq, motifs=None):
    """Return 1-based positions in N-X-S/T sequons (X != P).

    ``motifs`` may contain additional regular expressions; every match's span
    is returned.  This is a conservative "potential PTM" annotation, not a
    predictor.
    """
    positions = set()
    for i in range(len(seq) - 2):
        if seq[i] == 'N' and seq[i + 1] != 'P' and seq[i + 2] in 'ST':
            positions.add(i + 1)
    for pattern in motifs or []:
        for match in re.finditer(pattern, seq):
            positions.update(range(match.start() + 1, match.end() + 1))
    return sorted(positions)


def expand_ca_shell(target, chain_ids, seeds, dist):
    """Grow a reference-numbered seed set with CA neighbours per chain."""
    seeds = {int(x) for x in seeds or []}
    out = {c: set() for c in chain_ids}
    if not seeds:
        return out
    cutoff2 = float(dist) ** 2
    for cid in chain_ids:
        coords, refs = ca_coords(target, cid)
        seed_coords = [coord for ref, coord in zip(refs, coords)
                       if ref in seeds]
        if not seed_coords:
            continue
        seed_arr = np.asarray(seed_coords)
        for ref, coord in zip(refs, coords):
            if ref in seeds or ((seed_arr - coord) ** 2).sum(1).min() <= cutoff2:
                out[cid].add(ref)
    return out


def _parse_pockets(config):
    pockets = []
    if config.get('pockets'):
        for item in config['pockets']:
            if not isinstance(item, dict):
                raise ValueError('pockets entries must be mappings with '
                                 'name/residues/shell/fix_shell')
            pockets.append({
                'name': item.get('name', 'pocket'),
                'residues': [int(x) for x in item.get('residues', [])],
                'shell': float(item.get('shell', config.get(
                    'pocket_shell_dist', 8.0))),
                'fix_shell': bool(item.get('fix_shell', True)),
            })
    elif config.get('pocket_residues'):
        pockets.append({
            'name': 'pocket',
            'residues': [int(x) for x in config['pocket_residues']],
            'shell': float(config.get('pocket_shell_dist', 8.0)),
            'fix_shell': bool(config.get('fix_pocket_shell', True)),
        })
    return pockets


def _parse_disulfides(config, target, chain_ids):
    spec = config.get('disulfides', 'auto')
    if spec is False or spec is None:
        return []
    if spec is True or (isinstance(spec, str) and spec.lower() == 'auto'):
        return detect_disulfides(
            target, chain_ids,
            cutoff=float(config.get('disulfide_sg_dist', 2.5)))
    pairs = []
    for item in spec:
        if not isinstance(item, (list, tuple)) or len(item) != 4:
            raise ValueError('disulfides entries must be [chainA, refA, '
                             'chainB, refB]')
        pairs.append((str(item[0]), int(item[1]), str(item[2]), int(item[3])))
    return pairs


def _parse_ptms(config, target, chain_ids):
    spec = config.get('ptms', 'auto')
    if spec is False or spec is None:
        return _per_chain(None, chain_ids)
    if spec is True or (isinstance(spec, str) and spec.lower() == 'auto'):
        motifs = config.get('ptm_motifs') or []
        out = {}
        for cid in chain_ids:
            seq = ''.join(AMINO3.get(r.resname, 'X')
                          for r in target.chains[cid].residues)
            refs = [_ref_of(target, cid, r.id[1])
                    for r in target.chains[cid].residues]
            hits = {refs[i - 1] for i in detect_ptm_sequons(seq, motifs)
                    if 0 <= i - 1 < len(refs)}
            out[cid] = hits
        return out
    return _per_chain(spec, chain_ids)


def _parse_cofactors(config, target, chain_ids):
    """Fixed residues around ligands/cofactors.

    A spec is either a residue-name string (``HEM``) or a mapping with
    ``resnames``, optional ``chain``, explicit ``residues``, ``shell`` and
    ``fix_shell``.  Contacts are reference-numbered and protein-chain scoped.
    """
    out = []
    for i, spec in enumerate(config.get('cofactors') or []):
        if isinstance(spec, str):
            spec = {'resnames': [spec]}
        if not isinstance(spec, dict):
            raise ValueError(
                f'cofactors[{i}] must be a residue-name string or a mapping')
        resnames = spec.get('resnames')
        if resnames is None and spec.get('resname'):
            resnames = [spec['resname']]
        if isinstance(resnames, str):
            resnames = [resnames]
        resnames = [str(x).upper() for x in (resnames or [])]
        shell = float(spec.get('shell', config.get('cofactor_shell', 4.5)))
        fix_shell = bool(spec.get('fix_shell', True))
        ligand_chains = [spec['chain']] if spec.get('chain') else None
        ligands = hetero_residues(target.structure, resnames,
                                  ligand_chains) if resnames else []
        contacts = set()
        for cid, nums in _per_chain(spec.get('residues'), chain_ids).items():
            contacts.update((cid, int(n)) for n in nums)
        if ligands and fix_shell:
            coords = np.vstack([lig['coords'] for lig in ligands])
            for cid in chain_ids:
                for residue in target.chains[cid].residues:
                    try:
                        ca = residue['CA'].get_coord()
                    except KeyError:
                        continue
                    if np.min(np.linalg.norm(coords - ca, axis=1)) <= shell:
                        contacts.add((cid, _ref_of(
                            target, cid, residue.id[1])))
        out.append({
            'name': spec.get('name', f'cofactor_{i}'),
            'resnames': resnames,
            'contacts': sorted(contacts),
            'ligands': ligands,
        })
    return out


def build_constraints(target, chain_ids, config, topology):
    """Build the normalized constraint bundle consumed by ``build_mask``."""
    fixed_regions = config.get('fixed_regions')
    if fixed_regions is None:
        fixed_regions = ['ECD'] if config.get('fix_ecd', True) else []
    fixed_regions = {str(x).upper() for x in fixed_regions}

    protect_pad = int(config.get('protect_pad', 3))
    seed_fixed = _per_chain(config.get('fixed_residues'), chain_ids)
    fixed_reasons = {c: {} for c in chain_ids}

    def fix(cid, refn, reason):
        if not (1 <= int(refn)):
            return
        fixed_reasons[cid][int(refn)] = reason

    for cid in chain_ids:
        for refn in seed_fixed[cid]:
            for d in range(-protect_pad, protect_pad + 1):
                fix(cid, refn + d, 'protected')

    pockets = _parse_pockets(config)
    pocket_shell = {c: set() for c in chain_ids}
    for pocket in pockets:
        for cid in chain_ids:
            pocket_shell[cid].update(pocket['residues'])
        if pocket['fix_shell']:
            shell = expand_ca_shell(target, chain_ids, pocket['residues'],
                                    pocket['shell'])
            for cid in chain_ids:
                pocket_shell[cid].update(shell[cid])
                for refn in shell[cid]:
                    if refn not in fixed_reasons[cid]:
                        fix(cid, refn, 'pocket')
        else:
            for cid in chain_ids:
                for refn in pocket['residues']:
                    if refn not in fixed_reasons[cid]:
                        fix(cid, refn, 'pocket_seed')

    disulfides = _parse_disulfides(config, target, chain_ids)
    for ca, ra, cb, rb in disulfides:
        if ca in fixed_reasons:
            fix(ca, ra, 'disulfide')
        if cb in fixed_reasons:
            fix(cb, rb, 'disulfide')

    ptms = _parse_ptms(config, target, chain_ids)
    for cid in chain_ids:
        for refn in ptms[cid]:
            if refn not in fixed_reasons[cid]:
                fix(cid, refn, 'ptm')

    cofactors = _parse_cofactors(config, target, chain_ids)
    for cid, refn in (contact for cf in cofactors
                      for contact in cf['contacts']):
        fix(cid, refn, 'cofactor')

    epitopes = _per_chain(config.get('epitopes'), chain_ids)
    if config.get('fix_epitopes', True):
        for cid in chain_ids:
            for refn in epitopes[cid]:
                if refn not in fixed_reasons[cid]:
                    fix(cid, refn, 'epitope')

    pore = _per_chain(config.get('pore_residues'), chain_ids)
    pore_shell = {c: set() for c in chain_ids}
    if config.get('pore_shell'):
        grown = expand_ca_shell(target, chain_ids,
                                {r for s in pore.values() for r in s},
                                float(config['pore_shell']))
        for cid in chain_ids:
            pore_shell[cid].update(grown[cid])
    else:
        for cid in chain_ids:
            pore_shell[cid].update(pore[cid])
    if config.get('fix_pore', True):
        for cid in chain_ids:
            for refn in pore_shell[cid]:
                if refn not in fixed_reasons[cid]:
                    fix(cid, refn, 'pore')

    interface_mode = config.get('interface_mode')
    if interface_mode is None:
        interface_mode = 'fixed' if config.get('fix_interface_contacts', True) \
            else 'none'
    interface_mode = str(interface_mode).lower()
    if interface_mode not in ('fixed', 'core', 'none'):
        raise ValueError('interface_mode must be fixed, core or none')
    interface_counts = {}
    interface_core = set()
    if interface_mode != 'none' and len(chain_ids) > 1:
        from . import stability
        interface_counts = stability.interface_contact_counts(
            target, chain_ids, dist=float(config.get('interface_contact_dist',
                                                     6.0)))
        core_threshold = int(config.get('interface_core_contacts', 3))
        if interface_mode == 'core':
            for key, count in interface_counts.items():
                if count >= core_threshold:
                    interface_core.add(key)
        else:
            interface_core = set(interface_counts)

    return {
        'fixed_regions': fixed_regions,
        'fix_ecd': bool(config.get('fix_ecd', True)),
        'fixed_reasons': fixed_reasons,
        'pockets': pockets,
        'pocket_shell': pocket_shell,
        'disulfides': disulfides,
        'ptms': ptms,
        'epitopes': epitopes,
        'pore_residues': pore,
        'pore_shell': pore_shell,
        'cofactors': cofactors,
        'interface_mode': interface_mode,
        'interface_counts': interface_counts,
        'interface_core': interface_core,
        'interface_contact_dist': float(config.get('interface_contact_dist',
                                                   6.0)),
    }


def fixed_positions(constraints, chain_id):
    return set(constraints['fixed_reasons'].get(chain_id, {}))

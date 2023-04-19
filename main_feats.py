import sys
import os
import re
import ase
import moltop
import numpy as np
import pandas as pd
import networkx as nx
import libarvo as la
import morfeus as mf

from ast import literal_eval

# Local Packages
from moltop import adjacency as ca
from moltop import descriptors as ctd

from utils.get_3d import get_3d


def mol_from_nbo(log):

    try:
        block = re.search(
            "\s{23}Input\sorientation:.*?Z\n\s-{69}\n(.*?)\n\s-{69}", log, re.DOTALL
        ).group(1)
        mol_tab = np.asarray(
            [
                np.asarray([float(val) for val in line.split()])
                for line in block.split("\n")
            ]
        )
        mol = ase.Atoms(numbers=mol_tab[:, 1], positions=mol_tab[:, [3, 4, 5]])
    except AttributeError:

        block = re.search(
            "\s Symbolic Z-matrix:\n\sCharge.*?=\s[0-9]+\n(.*?)\n\s\n", log, re.DOTALL
        ).group(1)
        mol_tab = np.asarray(
            [np.asarray([val for val in line.split()]) for line in block.split("\n")]
        )
        mol = ase.Atoms(symbols=mol_tab[:, 0], positions=mol_tab[:, -3:].astype(float))
    return mol


def normalize(vec):
    """Normalize vector
    :param vec: A N-Dim vector
    :return norm_vec: The normalized vector
    """
    norm = np.linalg.norm(vec)
    if norm == 0:
        norm = np.finfo(vec.dtype).eps
    return vec / norm


def rotation_matrix_from_vectors(vec1, vec2):
    """Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (
        vec2 / np.linalg.norm(vec2)
    ).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s**2))

    return rotation_matrix


def align_w_origin(mol, inds, m_c, symm=False, orient=True):

    inds = inds.copy() - 1
    pos = mol.get_positions()

    pos -= m_c

    assert len(inds) == 2, f"'inds' length must be 2, currently {len(inds)}."

    # Rotate in place along z-axis

    dot_vec = normalize(pos[inds[0]] + pos[inds[1]])

    if all(dot_vec == np.array([0.0, 0.0, 1.0])):
        pos *= np.array([1.0, 1.0, -1.0])

    elif not all(dot_vec == np.array([0.0, 0.0, -1.0])):
        r_z = rotation_matrix_from_vectors(
            dot_vec,
            np.array([0.0, 0.0, -1.0]),
        )
        pos = np.asarray([r_z @ x for x in pos])

    # If not already, rotate in place along y-axis
    cross_vec = normalize(np.cross(pos[inds[0]], pos[inds[1]]))

    if all(cross_vec == np.array([0.0, -1.0, 0.0])):
        pos *= np.array([1.0, -1.0, 1.0])

    elif not all(cross_vec == np.array([0.0, 1.0, 0.0])):
        r_z = rotation_matrix_from_vectors(
            cross_vec,
            np.array([0.0, 1.0, 0.0]),
        )
        pos = np.asarray([r_z @ x for x in pos])

    mol.set_positions(pos)

    if orient:

        com = mol.get_center_of_mass()
        print(str(com), file=sys.stderr)
        if com[0] < 0:
            pos *= np.array([-1.0, -1.0, 1.0])
            mol.set_positions(pos)

    if symm:

        com = mol.get_center_of_mass()

        inflect = (com / np.abs(com)) * np.array([-1.0, 1.0, 1.0])
        inflect[2] = 1.0

        pos *= inflect
        mol.set_positions(pos)

    corr_inds = inds[pos[inds][:, 0].argsort()]

    return mol, corr_inds + 1


def get_nbo(log, inds):
    block = re.search("\sSummary.*?-{71}\n(.*?)\n\s[=]+", log, re.DOTALL).group(1)
    cnbos = np.asarray([float(x.split()[2]) for x in block.split("\n")])[inds - 1]
    block = re.search(
        "\sMolecular unit.*?\(.*?\)\s*\n(.*?)\n\s{7}-{31}", log, re.DOTALL
    ).group(1)

    enbos = []
    onbos = []
    for ind in inds:
        b_o = 0
        b_e = None
        for line in block.split("\n"):
            if line[7:11] == " LP " and int(line[19:23]) == ind:
                o = float(line[40:48])
                e = float(line[52:60])
                if o > b_o:
                    b_o = o
                    b_e = e
        assert b_e is not None
        onbos += (b_o,)
        enbos += (b_e,)

    return np.concatenate((enbos, onbos, cnbos))


def get_geo(mol, inds, m_c):

    inds = inds.copy() - 1
    pos = mol.get_positions()
    assert len(inds) == 2, f"Only 2 indices allowed. Got {len(inds)}."

    pa = m_c
    pb = pos[inds[0]]
    pc = pos[inds[1]]

    ab = pb - pa
    ac = pc - pa

    cos_ang = np.dot(ab, ac) / (np.linalg.norm(ab) * np.linalg.norm(ac))
    ang = np.arccos(cos_ang)

    return np.array([abs(ang)])


def get_mo(log):
    # Adapted from https://github.com/doyle-lab-ucla/auto-qchem
    float_or_int_regex = "[-+]?[0-9]*\.[0-9]+|[0-9]+"

    block = re.search(
        "Population.*?SCF [Dd]ensity.*?(\sAlph.*?)\n\s*Molecular", log, re.DOTALL
    ).group(1)
    energies = [
        re.findall(f"({float_or_int_regex})", s_part)
        for s_part in block.split("Alpha virt.", 1)
    ]
    occupied_energies, unoccupied_energies = [map(float, e) for e in energies]
    homo, lumo = max(occupied_energies), min(unoccupied_energies)

    block = re.search("Dipole\smoment.*?(Tot=.*?)\n", log, re.DOTALL).group(1)
    dipole = re.findall(f"({float_or_int_regex})", block)[0]

    return np.array([homo, lumo, float(dipole)])


def get_top(mol, inds, file):

    inds = inds.copy() - 1

    tops = pd.DataFrame()
    try:
        G = ca.g_from_mol(mol, mult=1.05, radii=ase.data.covalent_radii)
    except AssertionError:
        try:
            G = ca.g_from_mol(mol, mult=1.1, radii=ase.data.covalent_radii)
        except AssertionError:
            G = ca.g_from_mol(mol, mult=1.2, radii=ase.data.covalent_radii)

    if np.isnan(inds[0]) or np.isnan(inds[1]):
        t1 = 0
        t2 = 0
        t3 = 0
        t4 = 0
    else:

        t1 = ctd.local_top_state(G, inds[0])
        t2 = ctd.local_top_state(G, inds[1])
        t3 = ctd.el_top_state(G, inds[0])
        t4 = ctd.el_top_state(G, inds[1])

    bos = moltop.properties.get_bos(file, file_type="gaussian")

    nx.set_edge_attributes(G, 0, "bond_order")
    for edge in G.edges():
        try:
            G.edges[edge]["bond_order"] = bos[edge]
        except KeyError:
            G.edges[edge]["bond_order"] = 0

    g_ke = moltop.utils.no_terminal(G)
    g_an = moltop.utils.no_hydrogen(G)

    na = G.number_of_nodes()
    nb = G.number_of_edges()

    tops = pd.concat(
        [
            tops,
            pd.DataFrame(
                {
                    "na": [na],
                    "na_an": [g_an.number_of_nodes()],
                    "nb": [nb],
                    "bpa": [nb / na],
                    "estrada": [nx.estrada_index(g_an)],
                    "wiener": [nx.wiener_index(g_an)],
                    "global_eff": [nx.global_efficiency(g_an)],
                    "balaban": [ctd.balaban_j(g_an)],
                    "hosoya": [ctd.hosoya_z(g_an)],
                    "zagreb1": [ctd.zagreb(g_an, 1)],
                    "zagreb2": [ctd.zagreb(g_an, 2)],
                    "global_simple": [ctd.global_simple(g_an, bo_label="bond_order")],
                    "crest_flex": [ctd.crest_flex(g_ke, bo_label="bond_order")],
                    "kier_a": [ctd.kier_alpha(g_an, mode="a")],
                    "kier_b": [ctd.kier_alpha(g_an, mode="b")],
                    "kier_al": [ctd.kier_alpha(G, mode="legacy")],
                    "0k": [ctd.kier_mkappa(g_an, 0, alpha=False)],
                    "1k": [ctd.kier_mkappa(g_an, 1, alpha=False)],
                    "2k": [ctd.kier_mkappa(g_an, 2, alpha=False)],
                    "3k": [ctd.kier_mkappa(g_an, 3, alpha=False)],
                    "1ka": [ctd.kier_mkappa(g_an, 1, alpha=True, mode="a")],
                    "2ka": [ctd.kier_mkappa(g_an, 2, alpha=True, mode="a")],
                    "3ka": [ctd.kier_mkappa(g_an, 3, alpha=True, mode="a")],
                    "1kb": [ctd.kier_mkappa(g_an, 1, alpha=True, mode="b")],
                    "2kb": [ctd.kier_mkappa(g_an, 2, alpha=True, mode="b")],
                    "3kb": [ctd.kier_mkappa(g_an, 3, alpha=True, mode="b")],
                    "1kal": [ctd.kier_mkappa(G, 1, alpha=True, mode="legacy")],
                    "2kal": [ctd.kier_mkappa(G, 2, alpha=True, mode="legacy")],
                    "3kal": [ctd.kier_mkappa(G, 3, alpha=True, mode="legacy")],
                    "k_phi": [ctd.kier_phi(g_an, alpha=False)],
                    "k_phia": [ctd.kier_phi(g_an, alpha=True, mode="a")],
                    "k_xia": [ctd.kier_xi(g_an, mode="a")],
                    "k_phib": [ctd.kier_phi(g_an, alpha=True, mode="b")],
                    "k_xib": [ctd.kier_xi(g_an, mode="b")],
                    "k_phial": [ctd.kier_phi(G, alpha=True, mode="legacy")],
                    "k_xial": [ctd.kier_xi(G, mode="legacy")],
                    "redu": [ctd.redundancy(g_an)],
                    "0chi": [ctd.hall_kier_mchi(G, 0, valence=False)],
                    "1chi": [ctd.hall_kier_mchi(G, 1, valence=False)],
                    "2chi": [ctd.hall_kier_mchi(G, 2, valence=False)],
                    "3chi": [ctd.hall_kier_mchi(G, 3, valence=False)],
                    "4chi": [ctd.hall_kier_mchi(G, 4, valence=False)],
                    "5chi": [ctd.hall_kier_mchi(G, 5, valence=False)],
                    "0chiv": [ctd.hall_kier_mchi(G, 0, valence=True)],
                    "1chiv": [ctd.hall_kier_mchi(G, 1, valence=True)],
                    "2chiv": [ctd.hall_kier_mchi(G, 2, valence=True)],
                    "3chiv": [ctd.hall_kier_mchi(G, 3, valence=True)],
                    "4chiv": [ctd.hall_kier_mchi(G, 4, valence=True)],
                    "5chiv": [ctd.hall_kier_mchi(G, 5, valence=True)],
                    "T": [ctd.global_top_state(G)],
                    "Tnbo1": [t1],
                    "Tnbo2": [t2],
                    "Tenbo1": [t3],
                    "Tenbo2": [t4],
                }
            ),
        ]
    )

    tops_arr = tops.reset_index(drop=True).to_numpy()
    assert len(tops_arr) == 1, f"Array of length {len(tops_arr)}. Should be 1."

    return tops_arr[0]


def get_steric(mol, bondi_scale=1.17, b_sph_r=3.5, density=1e-4):
    pos = mol.get_positions()
    els = mol.get_chemical_symbols()

    h_mask = [False if an == "H" else True for an in els]
    pos_nh = pos[h_mask]

    # Total
    radii = np.asarray(mf.utils.get_radii(els, radii_type="bondi", scale=bondi_scale))
    t_vol, t_sur = la.molecular_vs(pos, radii)
    t_ova = t_sur / (4 * np.pi * ((3 * t_vol) / (4 * np.pi)) ** (2 / 3))

    # Total nH
    radii_nh = np.asarray(
        mf.utils.get_radii(els, radii_type="bondi", scale=bondi_scale)
    )[h_mask]
    t_vol_nh, t_sur_nh = la.molecular_vs(pos_nh, radii_nh)
    t_ova_nh = t_sur_nh / (4 * np.pi * ((3 * t_vol_nh) / (4 * np.pi)) ** (2 / 3))

    # Buried
    b_radii = np.concatenate((np.array([b_sph_r]), radii))
    b_pos = np.concatenate((np.array([[0.0, 0.0, 0.0]]), pos))
    b_t_vol, _ = la.molecular_vs(b_pos, b_radii)
    b_vol = ((4 / 3) * np.pi * b_sph_r**3) - (b_t_vol - t_vol)

    # Buried nH
    b_radii_nh = np.concatenate((np.array([b_sph_r]), radii_nh))
    b_pos_nh = np.concatenate((np.array([[0.0, 0.0, 0.0]]), pos_nh))
    b_t_vol_nh, _ = la.molecular_vs(b_pos_nh, b_radii_nh)
    b_vol_nh = ((4 / 3) * np.pi * b_sph_r**3) - (b_t_vol_nh - t_vol_nh)

    tvs = np.array([t_vol, t_sur, t_ova, t_vol_nh, t_sur_nh, t_ova_nh, b_vol, b_vol_nh])

    max_radius = np.max(
        np.linalg.norm(pos, axis=1)
        + mf.utils.get_radii(els, radii_type="bondi", scale=bondi_scale)
    )

    b_els = np.concatenate((np.array(["Cu"]), els))

    # Buried Octants
    b_bv = mf.BuriedVolume(
        elements=b_els,
        coordinates=b_pos,
        metal_index=1,
        radius=b_sph_r,
        radii_type="bondi",
        radii_scale=bondi_scale,
        density=density,
        include_hs=False,
    )
    b_octants = b_bv.octant_analysis().octants
    b_octants_a = np.array([float(b_octants["buried_volume"][i]) for i in range(8)])

    # Octants
    bv = mf.BuriedVolume(
        elements=b_els,
        coordinates=b_pos,
        metal_index=1,
        radius=max_radius + 0.1,
        radii_type="bondi",
        radii_scale=bondi_scale,
        density=density,
        include_hs=True,
    )
    octants = bv.octant_analysis().octants
    octants_a = np.array([float(octants["buried_volume"][i]) for i in range(8)])

    # Octants nH
    bv_nh = mf.BuriedVolume(
        elements=b_els,
        coordinates=b_pos,
        metal_index=1,
        radius=max_radius + 0.1,
        radii_type="bondi",
        radii_scale=bondi_scale,
        density=density,
        include_hs=False,
    )
    octants_nh = bv_nh.octant_analysis().octants
    octants_nh_a = np.array([float(octants_nh["buried_volume"][i]) for i in range(8)])

    return np.concatenate((tvs, b_octants_a, octants_a, octants_nh_a))


file = sys.argv[1]

if sys.argv[2] == "1":
    ind_lut = pd.read_csv("ligs/csd_ligs.csv")
    name = os.path.basename(file)[:-13] + ".xyz"
    # name = os.path.basename(file)[:-15] + ".xyz"

if sys.argv[2] == "2":
    ind_lut = pd.read_csv("ligs/lit_ligs.tsv", sep="\t")
    name = os.path.basename(file)

print(f"Featurize: '{name}'", file=sys.stderr)

with open(file, "r") as f:
    log = f.read()
    f.close()

# Start
# Get mol
mol = mol_from_nbo(log)
# Get indices

try:
    q = ind_lut.query("Filename == @name")
    p_inds = q[["Biting_Atom_1", "Biting_Atom_2", "Metal_coord"]].to_numpy()

    assert len(p_inds) == 1
    c_smi = list(q["Smiles_Regen"])[0]
    inds = p_inds[0][[0, 1]].copy().astype(int) + 1
    m_c = np.asarray(literal_eval(p_inds[0][2])).astype(float)

except KeyError:
    q = ind_lut.query("Filename == @name")
    p_inds = q[["Biting_Atom_1", "Biting_Atom_2"]].to_numpy()

    assert len(p_inds) == 1
    c_smi = list(q["Smiles_Regen"])[0]
    inds = p_inds[0].copy() + 1
    m_c = np.array([0.0, 0.0, 0.0])

# Align Mol
mol, corr_inds = align_w_origin(mol, inds, m_c)

# DEBUG/Write xyz
# ase.io.write(f"{name[:-3]}xyz", mol) # Generate xyz
# get_3d()
# exit()

### FEATS ###

# Get NBOs
nbos = get_nbo(log, corr_inds)
assert len(nbos) == 6, f"'get_nbo': {len(nbos)} values found. Should be 6."
# print(nbos)

# Get geom
geos = get_geo(mol, corr_inds, m_c)
assert len(geos) == 1, f"'get_geo': {len(geos)} values found. Should be 1."
# print(geos)

# Get MO
mos = get_mo(log)
assert len(mos) == 3, f"'get_mo': {len(mos)} values found. Should be 3."
# print(mos)

# Get TOP
tops = get_top(mol, corr_inds, file)
assert len(tops) == 54, f"'get_top': {len(tops)} values found. Should be 54."
# print(tops)

# Get Steric
sterics = get_steric(mol)
assert len(sterics) == 32, f"'get_steric': {len(sterics)} values found. Should be 32."
# print(sterics)

full_arr = np.concatenate((nbos, geos, mos, tops, sterics))

# name,c_smiles,nbo1,nbo2,onbo1,onbo2,cnbo1,cnbo2,angle,homo,lumo,dipole,na,na_an,nb,bpa,estrada,wiener,global_eff,balaban,hosoya,zagreb1,zagreb2,global_simple,crest_flex,kier_a,kier_b,kier_al,0k,1k,2k,3k,1ka,2ka,3ka,1kb,2kb,3kb,1kal,2kal,3kal,k_phi,k_phia,k_xia,k_phib,k_xib,k_phial,k_xial,redu,0chi,1chi,2chi,3chi,4chi,5chi,0chiv,1chiv,2chiv,3chiv,4chiv,5chiv,T,Tnbo1,Tnbo2,Tenbo1,Tenbo2,t_vol,t_sur,t_ova,t_vol_nh,t_sur_nh,t_ova_nh,b_vol,b_vol_nh,bonh1,bonh2,bonh3,bonh4,bonh5,bonh6,bonh7,bonh8,o1,o2,o3,o4,o5,o6,o7,o8,onh1,onh2,onh3,onh4,onh5,onh6,onh7,onh8
print(name + "," + c_smi + "," + ",".join(full_arr.astype(str)))

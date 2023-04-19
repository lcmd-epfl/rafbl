import ase
import sys

import networkx as nx
import numpy as np
import scipy as sp

from ase import neighborlist
from ase.io.extxyz import read_xyz
from ase.io.gaussian import read_gaussian_out
from sklearn.neighbors import KDTree as skKDTree

try:
    import pymatgen
    from pymatgen.analysis.chemenv.connectivity.connectivity_finder import (
        ConnectivityFinder,
    )
    from pymatgen.analysis.chemenv.connectivity.structure_connectivity import (
        StructureConnectivity,
    )
    from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import (
        MultiWeightsChemenvStrategy,
        SimplestChemenvStrategy,
    )
    from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import (
        LocalGeometryFinder,
    )
    from pymatgen.analysis.chemenv.coordination_environments.structure_environments import (
        LightStructureEnvironments,
    )
    from pymatgen.core import Lattice, Molecule, Structure
    from pymatgen.ext.matproj import MPRester
except ImportError:
    pymatgen = None

try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    rdkit = None


def get_spread(z, coords, radii=ase.data.covalent_radii, pad=5):
    libx, liby, libz = coords.argmin(0)
    uibx, uiby, uibz = coords.argmax(0)
    sx = (
        np.abs(-coords[libx][0] + radii[z[libx]] + coords[uibx][0] - radii[z[uibx]])
        + pad
    )
    sy = (
        np.abs(-coords[liby][1] + radii[z[liby]] + coords[uiby][1] - radii[z[uiby]])
        + pad
    )
    sz = (
        np.abs(-coords[libz][2] + radii[z[libz]] + coords[uibz][2] - radii[z[uibz]])
        + pad
    )
    return sx, sy, sz


def check_symmetric(am, tol=1e-8):
    return sp.linalg.norm(am - am.T, sp.Inf) < tol


def check_connected(am, tol=1e-8):
    sums = am.sum(axis=1)
    lap = np.diag(sums) - am
    eigvals, eigvects = np.linalg.eig(lap)
    return len(np.where(abs(eigvals) < tol)[0]) < 2


def check_am(am):
    if check_connected(am) and check_symmetric(am):
        return True
    else:
        print(f"Connected: {check_connected(am)}; Symmetric: {check_symmetric(am)}", file=sys.stderr)
        return False


def get_cutoffs(z, radii=ase.data.covalent_radii, mult=1.05):
    return [radii[zi] * mult for zi in z]


def ase_mol_to_pymatgen_structure(mol, radii=ase.data.vdw_radii, pad=5):
    a, b, c = get_spread(mol.get_atomic_numbers(), mol.get_positions(), radii, pad)
    lattice = Lattice.from_parameters(a=a, b=b, c=c, alpha=90, beta=90, gamma=90)
    structure = Structure(
        lattice,
        mol.get_atomic_numbers(),
        mol.get_positions(),
        coords_are_cartesian=True,
    )
    spread = max([a, b, c])
    return structure, spread


def run_chemenv(structure, spread):
    lgf = LocalGeometryFinder()
    lgf.setup_structure(structure)
    se = lgf.compute_structure_environments(
        maximum_distance_factor=1.4,
        only_cations=False,
        voronoi_distance_cutoff=spread,
    )
    strategy = SimplestChemenvStrategy(distance_cutoff=1.4, angle_cutoff=0.3)
    lse = LightStructureEnvironments.from_structure_environments(
        strategy=strategy, structure_environments=se
    )
    cf = ConnectivityFinder()
    sc = cf.get_structure_connectivity(lse)
    return sc


def radius_adjacency(
    z, coords, radii=ase.data.covalent_radii, mult=1.05, pruning_mult=3
):
    tree = skKDTree(coords)
    am = np.zeros((len(z), len(z)), dtype=int)
    z_radii = np.array(
        [radii[zi] if zi > 1 else radii[zi] * mult for zi in z], dtype=float
    )
    idxs = tree.query_radius(coords, radii[z] * mult * pruning_mult)
    for i, zi in enumerate(z):
        for j in idxs[i]:
            if j > i:
                if np.linalg.norm(coords[i] - coords[j]) <= mult * (
                    z_radii[i] + z_radii[j]
                ):
                    am[i, j] = am[j, i] = 1
    return am


def am_from_mol(
    mol,
    radii=ase.data.covalent_radii,
    atomids=None,
    mode="default",
    mult=1.05,
):
    am = np.zeros(
        (len(mol.get_atomic_numbers()), len(mol.get_atomic_numbers())), dtype=int
    )
    if mode == "default":
        am = radius_adjacency(
            mol.get_atomic_numbers(), mol.get_positions(), radii, mult
        )
    if mode == "ase":
        cutoff = get_cutoffs(mol.get_atomic_numbers(), radii, mult)
        nl = neighborlist.NeighborList(cutoff, self_interaction=False, bothways=True)
        nl.update(mol)
        am = nl.get_connectivity_matrix(sparse=False)
    if mode == "chemenv":
        if pymatgen is None:
            raise ImportError(
                "pymatgen must be installed for chemenv connectivity to be used."
            )
        structure, spread = ase_mol_to_pymatgen_structure(mol)
        sc = run_chemenv(structure, spread)
        for node, adj_dict in sc._graph.adjacency():
            am[node, [*adj_dict]] = 1
    if atomids is not None:
        z = mol.get_atomic_numbers()[atomids]
        coords = mol.get_positions()[atomids]
        am = am[atomids, atomids]
    else:
        z = mol.get_atomic_numbers()
        coords = mol.get_positions()
    return am


def am_from_file(
    filename,
    radii=ase.data.covalent_radii,
    atomids=None,
    mode="default",
    mult=1.05,
):
    if filename[-3:] == "xyz":
        mol = next(read_xyz(open(filename, "r")))
    elif filename[-3:] == "log":
        mol = read_gaussian_out(open(filename, "r"))
    else:
        raise Exception(f"No filename termination understood for {filename[-3:]}")
    am = am_from_mol(mol, radii, atomids, mode, mult)
    return am


def am_to_g(mol, atomids, am):
    assert check_am(am)
    if atomids is not None:
        z = mol.get_atomic_numbers()[atomids]
        coords = mol.get_positions()[atomids]
        am = am[atomids, atomids]
    else:
        z = mol.get_atomic_numbers()
        coords = mol.get_positions()
    # G = nx.from_numpy_matrix(am, create_using=nx.MultiGraph)
    G = nx.from_numpy_matrix(am, create_using=nx.Graph)
    an_dict = {i: z[i] for i in range(len(z))}
    coord_dict = {i: coords[i] for i in range(len(z))}
    nx.set_node_attributes(G, an_dict, "atomic_number")
    nx.set_node_attributes(G, coord_dict, "coordinates")
    ds = np.zeros((len(G.edges())))
    cs = np.zeros_like(ds)
    for i, edge in enumerate(G.edges()):
        ds[i] = np.linalg.norm(coords[edge[0]] - coords[edge[1]])
        cs[i] = z[edge[0]] * z[edge[1]]
    b_dict = nx.edge_betweenness_centrality(G, normalized=False)
    d_dict = {edge: d for edge, d in zip(b_dict.keys(), ds)}
    c_dict = {edge: c for edge, c in zip(b_dict.keys(), cs)}
    nx.set_edge_attributes(G, b_dict, "betweenness")
    nx.set_edge_attributes(G, d_dict, "distance")
    nx.set_edge_attributes(G, c_dict, "coulomb_energy")
    G.graph["graph_type"] = "plerogram"
    return G


def g_from_mol(
    mol,
    radii=ase.data.covalent_radii,
    atomids=None,
    mode="default",
    mult=1.05,
):
    am = am_from_mol(mol, radii, atomids, mode, mult)
    g = am_to_g(mol, atomids, am)
    return g


def g_from_file(
    filename,
    radii=ase.data.covalent_radii,
    atomids=None,
    mode="default",
    mult=1.05,
):
    if filename[-3:] == "xyz":
        mol = next(read_xyz(open(filename, "r")))
    elif filename[-3:] == "log":
        mol = read_gaussian_out(open(filename, "r"))
    else:
        raise Exception(f"No filename termination understood for {filename[-3:]}")
    g = g_from_mol(mol, radii, atomids, mode, mult)
    return g


def g_from_smiles(
    smiles,
    noh=False,
):
    G = nx.Graph()
    if rdkit is None:
        raise ImportError("rdkit must be installed for SMILES handling to be used.")
    m = Chem.MolFromSmiles(smiles)
    if not noh:
        m = Chem.AddHs(m)
        params = AllChem.srETKDGv3()
        params.useSmallRingTorsions = True
        AllChem.EmbedMolecule(m, params=params)
    pt = Chem.rdchem.GetPeriodicTable()
    z = []    
    for atom in m.GetAtoms():
        an = atom.GetAtomicNum()
        G.add_node(
            atom.GetIdx(),
            atomic_number=an,
            dg=atom.GetDegree(),
            charge=atom.GetFormalCharge(),
            hyb=atom.GetHybridization(),
            val=pt.GetNOuterElecs(an),
            is_aromatic=atom.GetIsAromatic(),
            Hs=atom.GetTotalNumHs(),
            rad=atom.GetNumRadicalElectrons,
            ring=atom.IsInRing(),
            chirality=atom.GetChiralTag(),
        )
        z.append(an)
    for bond in m.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_order=bond.GetBondType()
        )
    coords = m.GetConformer().GetPositions()
    coord_dict = {i: coords[i] for i, _ in enumerate(m.GetAtoms())}
    nx.set_node_attributes(G, coord_dict, "coordinates")
    ds = np.zeros((len(G.edges())))
    cs = np.zeros_like(ds)
    for i, edge in enumerate(G.edges()):
        ds[i] = np.linalg.norm(coords[edge[0]] - coords[edge[1]])
        cs[i] = z[edge[0]] * z[edge[1]]
    b_dict = nx.edge_betweenness_centrality(G, normalized=False)
    d_dict = {edge: d for edge, d in zip(b_dict.keys(), ds)}
    c_dict = {edge: c for edge, c in zip(b_dict.keys(), cs)}
    nx.set_edge_attributes(G, b_dict, "betweenness")
    nx.set_edge_attributes(G, d_dict, "distance")
    nx.set_edge_attributes(G, c_dict, "coulomb_energy")
    G.graph["graph_type"] = "plerogram"
    return G

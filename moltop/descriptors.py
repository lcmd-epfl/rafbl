import ase
import warnings
import networkx as nx
import numpy as np
from scipy.stats import gmean
from itertools import combinations, compress
from collections import Counter
from .utils import *
from .data import *


def path_finder(g, u, n, excludeSet=None):  # deprecated
    if excludeSet is None:
        excludeSet = set([u])
    else:
        excludeSet.add(u)
    if n == 0:
        return [[u]]
    paths = [
        [u] + path
        for neighbor in g.neighbors(u)
        for path in path_finder(g, neighbor, n - 1, excludeSet)
    ]
    excludeSet.remove(u)
    return paths


def weight_paths(g, paths, radii=ase.data.covalent_radii):  # deprecated
    t = 0
    for path in paths:
        if len(path) == 2:
            wt = 0
            for na in path:
                wt += radii[g.nodes[na]["atomic_number"]] / radii[6]
            t += wt / 2
        if len(path) == 3:
            wt = 0
            for na in path:
                wt += radii[g.nodes[na]["atomic_number"]] / radii[6]
            t += wt / 3
        if len(path) == 4:
            wt = 0
            for na in path:
                wt += radii[g.nodes[na]["atomic_number"]] / radii[6]
            t += wt / 4
    return t


def paths_finder(g, n):  # deprecated
    paths = []
    for na in g.nodes():
        paths.extend(path_finder(g, na, n))
    sorted_paths = [
        tuple(sorted(path[:: len(path) - 1]) + [path[1]])
        for path in paths
        if len(set(path)) == len(path)
    ]
    unique_paths = set(sorted_paths)
    return unique_paths


def _paths_finder_rev(g, n):
    def recu_path(g, na, n, sub_paths=[], path=None):
        if path is None:
            path = [na]
        for neighbor in g.neighbors(path[-1]):
            if neighbor not in path:
                if n - 1 == 0:
                    if neighbor > na:
                        sub_paths += [path + [neighbor]]
                else:
                    recu_path(g, na, n - 1, sub_paths, path + [neighbor])
        return sub_paths

    paths = []
    if n < 1:
        for na in g.nodes():
            paths.append([na])
    else:
        for na in g.nodes():
            paths.extend(recu_path(g, na, n, []))
    return paths


def _cluster_finder(g, n):
    # TODO
    warnings.warn("Cluster implementation is wrong as is now.")
    clusters = []
    if n > 2:
        for na in g.nodes():
            neighs = list(g.neighbors(na))
            cluster = [[na] + list(x) for x in combinations(neighs, n)]
            clusters.extend(cluster)
    return clusters


def _ring_finder(g, n):
    return [x for x in nx.cycle_basis(g) if len(x) == n]


def hall_kier_a(g, radii=ase.data.covalent_radii, mode="a"):  # deprecated
    A = 0
    if mode == "a":
        for na in g.nodes():
            A += radii[g.nodes[na]["atomic_number"]] / radii[6]
    elif mode == "b":
        for na in g.nodes():
            tot_dist = []
            for edge in g.edges(na):
                tot_dist.append(g.edges[edge]["distance"])
            A += np.mean(tot_dist) / 1.535  # Taken from Wiki, Sp3 C-C bond
    else:
        raise NotImplementedError
    return A


def kna(g, n, radii=ase.data.covalent_radii, mode="a"):  # deprecated
    A = hall_kier_a(g, radii, mode=mode)
    num = (A - (n - 1)) * ((A - n) ** 2)
    if n == 1:
        p = weight_paths(g, g.edges(), radii)
    elif n == 3:
        assert A > 2, f"Needs at least 3 atoms, got '{A}'."
        if A % 2 == 0:
            num = (A - 3) * ((A - 2) ** 2)
        else:
            num = (A - 1) * ((A - 3) ** 2)
        unique_paths = _paths_finder_rev(g, 3)
        p = weight_paths(g, unique_paths, radii)
        return num / ((p) ** 2)
    else:
        unique_paths = paths_finder(g, n)
        p = weight_paths(g, unique_paths, radii)
    return num / ((p) ** 2)


def kn(g, n):  # deprecated
    A = g.number_of_nodes()
    B = g.number_of_edges()
    num = (A - (n - 1)) * ((A - n) ** 2)
    if n == 1:
        p = B
    elif n == 3:
        assert A > 2, f"Needs at least 3 atoms, got '{A}'."
        if A % 2 == 0:
            num = (A - 3) * ((A - 2) ** 2)
        else:
            num = (A - 1) * ((A - 3) ** 2)
        unique_paths = _paths_finder_rev(g, 3)
        p = len(unique_paths)
        return num / ((p) ** 2)
    else:
        unique_paths = paths_finder(g, n)
        p = len(unique_paths)
    return num / ((p) ** 2)


def kier(g):  # deprecated
    A = g.number_of_nodes()
    return kn(g, 1) * kn(g, 2) / A


def hall_kier(g):  # deprecated
    A = g.number_of_nodes()
    return kna(g, 1) * kna(g, 2) / A


def kier_alpha(g, radii=ase.data.covalent_radii, mode="a"):
    k_alpha = 0

    def _plero_warn(g):
        if g.graph["graph_type"] == "plerogram":
            warnings.warn(
                "The graph type is 'plerogram'. This may lead to undesired results for Kier descriptors."
            )

    if mode == "a":
        _plero_warn(g)
        for na in g.nodes():
            k_alpha += (radii[g.nodes[na]["atomic_number"]] / radii[6]) - 1
    elif mode == "b":
        _plero_warn(g)
        for edge in g.edges():
            k_alpha += (
                g.edges[edge]["distance"] / 1.535
            ) - 1  # Taken from Wiki, Sp3 C-C bond
    elif mode == "legacy":
        warnings.warn("Legacy mode. Use only for uncharged/non-radical molecules.")
        assert (
            g.graph["graph_type"] == "plerogram"
        ), f"Plerogram required to compute hybridization."
        for na in g.nodes():
            an = g.nodes[na]["atomic_number"]
            if an != 1:
                try:
                    k_alpha += (
                        kier_radii[(an, len(list(g.neighbors(na))))]
                        / kier_radii[(6, 4)]
                        - 1
                    )
                except KeyError:
                    warnings.warn(
                        f"Atomic number '{an}' not tabulated. Using sp3 carbon."
                    )
    else:
        raise NotImplementedError(f"No mode '{mode}'.")

    return k_alpha


def molecular_shannon_i(g):
    atom_types = []
    a = g.number_of_nodes()
    for na in g.nodes():
        atom_types.append(
            tuple([g.nodes[na]["atomic_number"], len(list(g.neighbors(na)))])
        )
    freqs = dict([(k, v / a) for k, v in Counter(atom_types).items()])
    i = 0
    for _, rho_i in freqs.items():
        i -= rho_i * np.log10(rho_i)
    return i


def kier_mkappa(g, m, alpha=False, mode="a"):
    if alpha:
        alf = kier_alpha(g, mode=mode)
    else:
        alf = 0
    if mode == "legacy":
        g = no_hydrogen(g)
    a = g.number_of_nodes()
    if m == 0:
        return molecular_shannon_i(g) * a
    elif m == 1:
        num = (a + alf - 0) * (a + alf - 1) ** 2
        p = g.number_of_edges()
    elif m == 2:
        num = (a + alf - 1) * (a + alf - 2) ** 2
        unique_paths = _paths_finder_rev(g, m)
        p = len(unique_paths)
    elif m == 3:
        assert a > 2, f"Needs at least 3 atoms, got '{a}'."
        # TODO
        warnings.warn("3K may not works with cyclopropanes.")
        if a % 2 == 0:
            num = (a + alf - 3) * ((a + alf - 2) ** 2)
        else:
            num = (a + alf - 1) * ((a + alf - 3) ** 2)
        unique_paths = _paths_finder_rev(g, m)
        p = len(unique_paths)
    else:
        raise NotImplementedError("Invalid 'm', '{m}'.")
    return num / ((p + alf) ** 2)


def kier_phi(g, alpha=False, mode="a"):
    if mode == "legacy":
        num_a = no_hydrogen(g).number_of_nodes()
    else:
        num_a = g.number_of_nodes()
    return (
        kier_mkappa(g, 1, alpha=alpha, mode=mode)
        * kier_mkappa(g, 2, alpha=alpha, mode=mode)
        / num_a
    )


def kier_xi(g, mode="a"):
    return (
        2 * kier_mkappa(g, 1, alpha=True, mode=mode)
        - kier_mkappa(g, 0)
        - kier_mkappa(g, 3, alpha=True, mode=mode)
    )


def redundancy(g):
    a = g.number_of_nodes()
    return 1 - (molecular_shannon_i(g) / np.log10(a))


def balaban_j(g):
    n = g.number_of_nodes()
    m = g.number_of_edges()
    gamma = m - n + nx.number_connected_components(g)
    prefactor = m / (gamma + 1)
    dict_length = dict(nx.all_pairs_shortest_path_length(g))
    summation = 0
    for (i, j) in g.edges():
        di = sum([dict_length[i][k] for k in g.nodes()])
        dj = sum([dict_length[j][k] for k in g.nodes()])
        summation += (di * dj) ** -0.5
    return prefactor * summation


def hosoya_z(g):
    return len(nx.algorithms.approximation.maximum_independent_set(g))


def zagreb(g, index=1):
    z = 0
    if index == 1:
        for na in g.nodes():
            z += len(list(g.neighbors(na))) ** 2
    elif index == 2:
        for na in g.nodes():
            nei_na = len(list(g.neighbors(na)))
            for nb in g.neighbors(na):
                nei_nb = len(list(g.neighbors(nb)))
                if na < nb:
                    z += nei_na * nei_nb
    else:
        raise ValueError("Only 1 and 2")
    return z


def local_simple(g, source, target, bo_label="bond_order", normalized=False):
    if not nx.get_edge_attributes(g, bo_label) or bo_label == "":
        nx.set_edge_attributes(g, 0.0, bo_label)
        warnings.warn("'local_simple' No bond orders founds.")
    paths = nx.all_shortest_paths(g, source, target)
    local_simples = []
    for path in paths:
        spn = len(path) - 1
        assert spn != 0
        nrb = 0
        bx = 0
        by = 0
        cycles = nx.cycle_basis(g)
        for i, n in enumerate(path):
            nn = len(list(g.neighbors(n)))
            if nn == 3:
                bx += 1
            elif nn == 4:
                by += 1
            if i != 0:
                iring = any([np.isin([path[i - 1], path[i]], x).all() for x in cycles])
                mbond = g.edges[path[i - 1], path[i]][bo_label] > float(4 / 3)
                if iring:
                    nrb += 1
                elif mbond:
                    nrb += 1
        if normalized:
            local_simples.append(((nrb + 0.75 * bx + 0.5 * by) / 2) / spn)
        else:
            local_simples.append(spn - ((nrb + 0.75 * bx + 0.5 * by) / 2))
    return np.mean(local_simples)


def global_simple(g, bo_label="bond_order", normalized=False):
    global_simples = []
    for c in combinations(g.nodes(), 2):
        global_simples.append(
            local_simple(g, c[0], c[1], bo_label=bo_label, normalized=normalized)
        )
    return np.mean(global_simples)


def radius_simple(g, n, radius, bo_label="bond_order", normalized=False):
    sg = nx.ego_graph(g, n, radius)
    return global_simple(sg, bo_label=bo_label, normalized=normalized)


def crest_flex(g, bo_label="bond_order"):
    if g.graph["graph_type"] != "kenogram":
        warnings.warn("Per definition 'crest_flex' works on kenograms.")
    if not nx.get_edge_attributes(g, bo_label) or bo_label == "":
        nx.set_edge_attributes(g, 0.0, bo_label)
        warnings.warn("'crest_flex' No bond orders founds.")
    cycles = nx.cycle_basis(g)
    edges = list(g.edges())
    m = len(edges)
    av2 = 0.0
    for edge in edges:
        assert len(edge) == 2
        cns = np.array([len(list(g.neighbors(n))) for n in edge])
        hybf = 1.0
        hybf *= 0.5 if g.nodes[edge[0]]["atomic_number"] == 6 and cns[0] < 4 else 1.0
        hybf *= 0.5 if g.nodes[edge[1]]["atomic_number"] == 6 and cns[1] < 4 else 1.0
        doublef = 1.0 - np.exp(-4.0 * (g.edges[edge][bo_label] - 2.0) ** 6)
        branch = 2.0 / np.sqrt(np.prod(cns))
        iring = [
            np.isin([edge[0], edge[1]], x).any() for x in cycles
        ]  # Replace by .all() for better implementation
        try:
            k = min([len(x) for x in list(compress(cycles, iring))])
        except ValueError:
            k = 0
        ringf = 0.5 * (1.0 - np.exp(-0.06 * k)) if k > 0 else 1.0
        val = branch * ringf * doublef * hybf
        av2 += val**2
        av2 = np.sqrt(av2 / m) if m > 0 else av2
    return av2


def _guess_valence(an):
    val = an
    for period in [2, 8, 8, 18, 18, 32, 32]:
        if val - period <= 0:
            break
        else:
            val -= period
    return val


def hall_kier_mchi(g, m, valence=False, subgraphs="p", valence_label="val"):
    # ring type is abbreviated 'r' and not 'ch' for simplicity
    assert (
        g.graph["graph_type"] == "plerogram"
    ), f"Plerogram required to compute valence."
    paths = []
    if valence:
        if not nx.get_node_attributes(g, valence_label):
            nx.set_node_attributes(
                g,
                {na: _guess_valence(g.nodes[na]["atomic_number"]) for na in g.nodes()},
                valence_label,
            )
            warnings.warn("Guessing valences")
    g_noh = no_hydrogen(g)
    if m == -1:
        paths += [list(g_noh.nodes())]
    else:
        if "p" in subgraphs:
            paths += _paths_finder_rev(g_noh, m)
        if "c" in subgraphs:
            paths += _cluster_finder(g_noh, m)
        if "r" in subgraphs:
            paths += _ring_finder(g_noh, m)
    cs = []
    for path in paths:
        if valence:
            deltas = []
            for na in path:
                an = g_noh.nodes[na]["atomic_number"]
                val = g_noh.nodes[na][valence_label]
                n_h = np.sum(
                    [
                        g.nodes[neigh]["val"]
                        for neigh in g.neighbors(na)
                        if g.nodes[neigh]["atomic_number"] == 1
                    ]
                )
                if an > 10:
                    delta = (val - n_h) / (an - val - 1)
                else:
                    delta = val - n_h
                deltas.append(delta)
            cs.append(np.prod(deltas) ** -0.5)
        else:
            cs.append(np.prod([len(list(g_noh.neighbors(na))) for na in path]) ** -0.5)
    return np.sum(cs)


def _single_top_state(g, source, target):
    if source == target:
        return g.nodes[source]["delta_v"]
    else:
        paths = nx.all_simple_paths(g, source, target)
        tij = 0
        for path in paths:
            tij += gmean([g.nodes[na]["delta_v"] for na in path]) / len(path)
        return tij


def local_top_state(g, ai, valence_label="val"):
    assert (
        g.graph["graph_type"] == "plerogram"
    ), f"Plerogram required to compute valence."
    if not nx.get_node_attributes(g, valence_label):
        nx.set_node_attributes(
            g,
            {na: _guess_valence(g.nodes[na]["atomic_number"]) for na in g.nodes()},
            valence_label,
        )
        warnings.warn("Guessing valences")
    deltas = {}
    for na in g.nodes():
        an = g.nodes[na]["atomic_number"]
        val = g.nodes[na][valence_label]
        n_h = np.sum(
            [
                g.nodes[neigh]["val"]
                for neigh in g.neighbors(na)
                if g.nodes[neigh]["atomic_number"] == 1
            ]
        )
        if an > 10:
            delta = (val - n_h) / (an - val - 1)
        else:
            delta = val - n_h
        deltas[na] = delta
    nx.set_node_attributes(g, deltas, "delta_v")
    g_noh = no_hydrogen(g)
    tis = []
    for c in g_noh.nodes():
        tis.append(_single_top_state(g_noh, ai, c))
    return np.sum(tis)


def global_top_state(g, valence_label="val"):
    assert (
        g.graph["graph_type"] == "plerogram"
    ), f"Plerogram required to compute valence."
    if not nx.get_node_attributes(g, valence_label):
        nx.set_node_attributes(
            g,
            {na: _guess_valence(g.nodes[na]["atomic_number"]) for na in g.nodes()},
            valence_label,
        )
        warnings.warn("Guessing valences")
    deltas = {}
    for na in g.nodes():
        an = g.nodes[na]["atomic_number"]
        val = g.nodes[na][valence_label]
        n_h = np.sum(
            [
                g.nodes[neigh]["val"]
                for neigh in g.neighbors(na)
                if g.nodes[neigh]["atomic_number"] == 1
            ]
        )
        if an > 10:
            delta = (val - n_h) / (an - val - 1)
        else:
            delta = val - n_h
        deltas[na] = delta
    nx.set_node_attributes(g, deltas, "delta_v")
    g_noh = no_hydrogen(g)
    tis = []
    for c in g_noh.nodes():
        tis.append(_single_top_state(g_noh, c, c))
    for c in combinations(g_noh.nodes(), 2):
        tis.append(_single_top_state(g_noh, c[0], c[1]))
    return np.sum(tis)


def _get_n_quantum(an):
    if an < 1:
        raise ValueError("Atomic number cannot be negative or zero.")
    elif an <= 2:
        return 1
    elif an <= 10:
        return 2
    elif an <= 18:
        return 3
    elif an <= 36:
        return 4
    elif an <= 54:
        return 5
    elif an <= 86:
        return 6
    elif an <= 118:
        return 7
    else:
        return 8  # maximum principal quantum number set to this value


def el_top_state(g, ai, valence_label="val"):
    assert (
        g.graph["graph_type"] == "plerogram"
    ), f"Plerogram required to compute valence."
    if not nx.get_node_attributes(g, valence_label):
        nx.set_node_attributes(
            g,
            {na: _guess_valence(g.nodes[na]["atomic_number"]) for na in g.nodes()},
            valence_label,
        )
        warnings.warn("Guessing valences")
    deltas = {}
    for na in g.nodes():
        an = g.nodes[na]["atomic_number"]
        val = g.nodes[na][valence_label]
        n_h = np.sum(
            [
                g.nodes[neigh]["val"]
                for neigh in g.neighbors(na)
                if g.nodes[neigh]["atomic_number"] == 1
            ]
        )
        if an > 10:
            delta = val - n_h
            # delta = (val - n_h) / (an - val - 1)
        else:
            delta = val - n_h
        deltas[na] = delta
    nx.set_node_attributes(g, deltas, "delta_v")
    g_noh = no_hydrogen(g)
    nx.set_node_attributes(
        g_noh, {na: len(list(g_noh.neighbors(na))) for na in g_noh.nodes()}, "delta"
    )
    el_states = {}
    for na in g_noh.nodes():
        d = g_noh.nodes[na]["delta"]
        d_v = g_noh.nodes[na]["delta_v"]
        nq = _get_n_quantum(g_noh.nodes[na]["atomic_number"])
        el_states[na] = (((2 / nq) ** 2 * d_v) + 1) / d
    nx.set_node_attributes(g_noh, el_states, "el_state")
    si = 0
    for na in g_noh.nodes():
        if na == ai:
            si += g_noh.nodes[na]["el_state"]
        else:
            spl = nx.shortest_path_length(g_noh, ai, na) + 1
            si += (g_noh.nodes[ai]["el_state"] - g_noh.nodes[na]["el_state"]) / spl**2
    return si

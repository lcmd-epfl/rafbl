def no_hydrogen(g, cp=True):

    non_hydrogens = []
    for node in g.nodes():
        if g.nodes[node]["atomic_number"] != 1:
            non_hydrogens.append(node)

    g_sub = g.subgraph(non_hydrogens).copy() if cp else g.subgraph(non_hydrogens)
    g_sub.graph["graph_type"] = "anydrogram"

    return g_sub


def no_terminal(g, cp=True):

    non_terminals = []
    for node in g.nodes():
        if len(list(g.neighbors(node))) > 1:
            non_terminals.append(node)

    g_sub = g.subgraph(non_terminals).copy() if cp else g.subgraph(non_terminals)
    g_sub.graph["graph_type"] = "kenogram"

    return g_sub

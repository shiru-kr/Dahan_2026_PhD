# this script computes the following multiplex graph metrics:
# Multilayer communities + modularity Q - nodal per layer: (N,), global: float.
# Flexibility - nodal: (N,L), global: float.
# Persistence - global: float.
# Alligance - nodal: (N,), global: float.

import numpy as np
import networkx as nx
import community as community_louvain


##### Preprocessing matrices: #####

def threshold_binarize(layers, threshold=0.10, binarize=True, eps=1e-12):
    """
    Proportional threshold (by density) and optionally binarize layers.

    Keeps the top `threshold` proportion of edges (by weight) in each layer. Optionally binarizes the retained edges.

    Parameters
    ----------
    layers : list of (N, N) np.ndarray
        Weighted adjacency matrices (symmetric). Assumed non-negative.
    threshold : float
        Proportion of strongest edges to keep, in (0, 1]. Example: 0.15 keeps top 15%.
    binarize : bool
        If True, return binary matrices (0/1). If False, keep original weights for retained edges.
    eps : float
        Values <= eps are treated as zero (absent) when building the threshold distribution.

    Returns
    -------
    out_layers : list of (N, N) np.ndarray
        Thresholded (and optionally binarized) adjacency matrices.
    """
    if not (0 < threshold <= 1):
        raise ValueError("threshold must be in (0, 1]. For example, 0.15 keeps top 15% edges.")

    out_layers = []

    for A in layers:
        A = np.asarray(A)

        # Upper triangle weights (exclude diagonal), use only positive edges
        ut = np.triu(A, k=1)
        w = ut[ut > eps]

        if w.size == 0:
            # No edges: return all zeros
            out_layers.append(np.zeros_like(A, dtype=float))
            continue

        # Weight cutoff for top `threshold` proportion
        cutoff = np.percentile(w, 100 - 100 * threshold)

        # Mask of retained edges
        mask = A >= cutoff

        if binarize:
            B = mask.astype(float)
        else:
            B = np.zeros_like(A, dtype=float)
            B[mask] = A[mask]

        # Enforce symmetry and zero diagonal
        B = np.maximum(B, B.T)
        np.fill_diagonal(B, 0.0)

        out_layers.append(B)

    return out_layers


##### Community detection and modulartrity functions: #####
# These functions assumes each layer has already been preprocessed:
#    - Binary: elements are 0/1
#    - Weighted: elements are between 0 to 1

def supra_g(layers, omega=0.8, eps=1e-12):
    """
    Build supra-adjacency matrix from layers with inter-layer coupling omega.

    Parameters
    ----------
    layers : list of (N, N) np.ndarray
        Symmetric, non-negative adjacency matrices (one per layer).
    omega : float
        Inter-layer coupling weight between identical nodes across layers.

    Returns
    -------
    G : nx.Graph
        Supra-adjacency graph.
    """
    L = len(layers)
    N = layers[0].shape[0]
    G = nx.Graph()

    # Intralayer edges
    for ell in range(L):
        A = layers[ell]
        for i in range(N):
            for j in range(i+1, N):
                w = A[i, j]
                if w > eps:  # Use eps from global params
                    u = i + ell * N
                    v = j + ell * N
                    G.add_edge(u, v, weight=float(w))

    # Interlayer coupling
    for i in range(N):
        for ell1 in range(L):
            for ell2 in range(ell1+1, L):
                u = i + ell1 * N
                v = i + ell2 * N
                G.add_edge(u, v, weight=float(omega))

    return G


def multilayer_communities(
    G,
    layers,
    resolution=1.5,
    consensus=True,
    n_runs=500,
    consensus_threshold=0.5,
    random_state=0,
    weight="weight",
):
    """
    Detects community structure in a multilayer supra-graph using Louvain modularity optimisation.

    Parameters
    ----------
    G : nx.Graph
        Supra-adjacency graph containing all intra- and inter-layer connections.
        Nodes must be indexed consecutively from 0 to N*L - 1.
    layers : list of (N, N) np.ndarray
        List of adjacency matrices defining the individual network layers.
        Used only to infer the number of nodes (N) and layers (L).
    resolution (gamma) : float, optional
        Louvain resolution parameter controlling community granularity.
        Higher values produce smaller communities.
    consensus : bool, optional
        If False, performs a single Louvain partition.
        If True, performs repeated Louvain runs and derives a consensus partition
        from the supra-node agreement matrix.
    n_runs : int, optional
        Number of Louvain repetitions used when consensus=True.
    consensus_threshold : float, optional
        Minimum co-assignment frequency required to retain an edge in the
        agreement matrix during consensus clustering.
    random_state : int, optional
        Seed used for reproducibility in single-run mode and in the final
        consensus reclustering step.
    weight : str, optional
        Edge attribute used as graph weight (see community_louvain docs).

    Returns
    -------
    comm_cons : np.ndarray, shape (N, L)
        Community assignment matrix where rows correspond to nodes and columns
        correspond to layers.
    Q_cons : float
        Louvain modularity of the returned partition evaluated on the original
        supra-graph.

    Notes
    -----
    In consensus mode, Louvain is repeated across seeds 0..n_runs-1 (reproducible). A dense
    agreement matrix is constructed from the frequency with which supra-nodes
    are assigned to the same community across runs. After thresholding weak
    co-assignment frequencies, Louvain is applied once more to this agreement
    graph to obtain the final consensus partition.
    """
    L = len(layers)
    N = layers[0].shape[0]
    NL = N * L
    assert G.number_of_nodes() == NL, "Supra-graph size does not match N*L."
    if set(G.nodes()) != set(range(NL)):
        raise ValueError("Supra-graph nodes must be exactly 0..NL-1.")

    def one_run(seed):
        part = community_louvain.best_partition(
            G, weight=weight, resolution=resolution, random_state=seed
        )
        # ensure every supra-node has a label
        if set(part.keys()) != set(range(NL)):
            raise ValueError("Partition missing nodes (check supra-graph construction).")
        labels_1d = np.array([part[i] for i in range(NL)], dtype=int)
        comm = labels_1d.reshape(L, N).T
        Q = community_louvain.modularity(part, G, weight=weight)
        return comm, Q

    # single run
    if not consensus:
        return one_run(random_state)

    # build dense agreement matrix
    agreement = np.zeros((NL, NL), dtype=float)

    for seed in range(n_runs):
        comm, _ = one_run(seed)
        labels = comm.T.reshape(NL)

        for c in np.unique(labels):
            nodes_c = np.where(labels == c)[0]
            if nodes_c.size > 1:
                agreement[np.ix_(nodes_c, nodes_c)] += 1.0

    # remove diagonal and normalise
    np.fill_diagonal(agreement, 0.0)
    agreement /= float(n_runs)

    # threshold weak co-assignments
    agreement[agreement < consensus_threshold] = 0.0

    # recluster agreement graph to get consensus partition
    G_cons = nx.from_numpy_array(agreement)
    part_cons = community_louvain.best_partition(
        G_cons, weight=weight, resolution=resolution, random_state=random_state
    )

    if set(part_cons.keys()) != set(range(NL)):
        raise ValueError("Consensus partition missing nodes (unexpected).")

    labels_cons_1d = np.array([part_cons[i] for i in range(NL)], dtype=int)
    comm_cons = labels_cons_1d.reshape(L, N).T

    # modularity of the consensus partition on ORIGINAL supra-graph
    part_dict = {i: int(labels_cons_1d[i]) for i in range(NL)}
    Q_cons = community_louvain.modularity(part_dict, G, weight=weight)

    return comm_cons, Q_cons



def flexibility(comm):
    '''
    The flexibility of each node is calculated as the number of times that it changes community assignment, 
    normalized by the total possible number of changes. 
    In categorical multilayer networks, community assignment changes are possible between any pairs of layers.
    
    Parameters
    ----------
    comm : (N, L) array of community labels

    Returns
    -------
    F : numpy array of shape (N,)
        The flexibility score of each node range [0, 1]. 0 - node does not change community across layers. 1 - node changes community across all layers.
    GF : float
        Global flexibility measure for the network
    '''
    
    N, L = comm.shape
    F = np.zeros(N, dtype=float)
    total_pairs = L * (L - 1) / 2  # possible layer pairs

    if total_pairs == 0:
        return F  # all zeros if only 1 layer

    for i in range(N):
        changes = 0
        for layer1 in range(L):
            for layer2 in range(layer1 + 1, L):
                if comm[i, layer1] != comm[i, layer2]:
                    changes += 1
        F[i] = changes / total_pairs
    GF = F.mean()

    return F, GF


def persistence(comm):
    """
    Global persistence (categorical):

    The persistence (global) of a network is the normalized count of non-changes
    across nodes and unordered layer pairs. It ranges from 0 (all changes) to 1
    (no changes).

    Parameters
    ----------
    comm : (N, L) array of community labels

    Returns
    -------
    GP : float
        Global persistence measure in [0, 1].
    """
    N, L = comm.shape
    total_pairs = L * (L - 1) / 2  # possible layer pairs

    if total_pairs == 0:
        return 0.0  # only 1 layer

    same_count = 0  # total number of non-changes across all nodes and pairs

    for i in range(N):
        for layer1 in range(L):
            for layer2 in range(layer1 + 1, L):
                if comm[i, layer1] == comm[i, layer2]:
                    same_count += 1

    GP = same_count / (N * total_pairs)
    return GP


def multiplex_participation(G, comm):
    """
    Multilayer (supra-graph) participation coefficient.

    Parameters
    ----------
    G : nx.Graph
        Supra-graph with weights.
    comm : np.ndarray (N, L)
        Community labels from multilayer detection.

    Returns
    -------
    P_node : (N,)
        Participation coefficient per physical node.
    """

    N, L = comm.shape
    NL = N * L

    # flatten community labels to supra ordering
    labels = comm.T.reshape(NL)

    # total degree per supra node
    k_total = dict(G.degree(weight="weight"))

    # participation per supra node
    P_supra = np.zeros(NL)

    for u in range(NL):
        k_i = k_total[u]
        if k_i == 0:
            continue

        # degree to each community
        k_im = {}
        for v, data in G[u].items():
            w = data.get("weight", 1.0)
            m = labels[v]
            k_im[m] = k_im.get(m, 0.0) + w

        P_supra[u] = 1 - sum((k_im[m] / k_i) ** 2 for m in k_im)

    # reshape to (N, L)
    P_layer = P_supra.reshape(L, N).T

    # aggregate across layers
    P = P_layer.mean(axis=1)
  
    return P



def foo(x):
    return x * 2           

if __name__ == "__main__":
    # quick test, not used when imported
    print(foo(3))


import pandas as pd
import networkx as nx
import pickle
import sys


def read_file(path, delimiter):
    try:
        # Read the TSV file
        df = pd.read_csv(path, sep=delimiter)
        return df
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None


def export_graph_to_pickle(graph, filename):
    with open(filename, "wb") as f:
        pickle.dump(graph, f)


def import_graph_from_pickle(filename):
    with open(filename, "rb") as f:
        return pickle.load(f)


def create_interactome(data, columns):
    G = nx.DiGraph()
    for index, row in data.iterrows():
        node1 = row[columns[0]]
        node2 = row[columns[1]]
        experiment = row[columns[2]]
        if experiment != 0:
            if not G.has_node(node1):
                G.add_node(node1)
            if not G.has_node(node2):
                G.add_node(node2)

            if not G.has_edge(node1, node2) or not G.has_edge(node2, node1):
                G.add_edge(node1, node2, weight=experiment)
                G.add_edge(node2, node1, weight=experiment)
    return G


def bfs_k_hop(G: nx.DiGraph, node_set, k):
    k_hop_neighbors_set = set()

    # Iterate over each node in the node set
    for node in node_set:
        # Check if the node exists in the graph
        if node in G:
            # Perform BFS with depth limit k
            bfs_edges = nx.bfs_edges(G, node, depth_limit=k)

            # Extract nodes from the BFS traversal
            bfs_nodes = set([node for edge in bfs_edges for node in edge])
            bfs_nodes.add(node)  # Include the starting node itself

            # Add the BFS nodes to the global set
            k_hop_neighbors_set.update(bfs_nodes)

    return k_hop_neighbors_set


def main():
    pathway_path = "string_interactions_short.tsv"
    stringdb_path = "9606.protein.links.detailed.v12.0.txt"

    # read in files
    df_pathway = read_file(pathway_path, "\t")
    df_stringdb = read_file(stringdb_path, " ")

    columns = ["protein1", "protein2", "experimental"]
    G: nx.DiGraph = import_graph_from_pickle("interactome-experimental.pickle")
    # G: nx.DiGraph = create_interactome(df_stringdb, columns)
    # export_graph_to_pickle(G, "interactome-experimental.pickle")
    print("Interactome")
    print("nodes:" + str(len(G.nodes())))
    print("edges:" + str(len(G.edges()) / 2))
    print()

    columns = [
        "node1_string_id",
        "node2_string_id",
        "experimentally_determined_interaction",
    ]
    # D: nx.DiGraph = import_graph_from_pickle("pathway.pickle")
    D: nx.DiGraph = create_interactome(df_pathway, columns)
    # export_graph_to_pickle(D, "pathway.pickle")
    print("Pathway")
    print("nodes:" + str(len(D.nodes())))
    print("edges:" + str(len(D.edges()) / 2))
    print()
    # check how many edges in pathway is missing in network
    count = 0
    remove_edges = []
    for edge in D.edges():
        if not G.has_edge(edge[0], edge[1]) or not G.has_edge(edge[1], edge[0]):
            # remove the edge?
            remove_edges.append([edge[0], edge[1]])
            count += 1
    print("edges in pathway missing in interactome: ", count)

    for edge in remove_edges:
        if D.has_edge(edge[0], edge[1]):
            D.remove_edge(edge[0], edge[1])
        if D.has_edge(edge[1], edge[0]):
            D.remove_edge(edge[1], edge[0])

    # trying to subset the network
    print()
    print("Subnetwork generation")
    print()

    k = 2

    # Get the k-hop neighbors from the small nodes
    neighbors = bfs_k_hop(G, set(D.nodes()), k)

    # Create a subgraph containing only the k-hop neighbors and their interactions
    k_hop_subgraph = G.subgraph(neighbors)

    print("nodes:" + str(len(k_hop_subgraph.nodes())))
    print("edges:" + str(len(k_hop_subgraph.edges())/2))



    interaction = []
    for index, row in df_stringdb.iterrows():
        node1 = row["protein1"]
        node2 = row["protein2"]
        combined = '_'.join(sorted([node1, node2]))
        interaction.append(combined)
    df_stringdb["interaction"] = interaction

if __name__ == "__main__":
    main()

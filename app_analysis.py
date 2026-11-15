import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms.community import greedy_modularity_communities, modularity

# ==========================================
# STEP 1: Load and Preprocess Real-World APP PPI Data
# ==========================================
script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, 'app_interactions.tsv')

df = pd.read_csv(file_path, sep='\t')

# Filter for Highest Confidence (0.9) and remove the ego-node (APP)
# to eliminate artificial star-graph sampling bias
df_filtered = df[(df['combined_score'] >= 0.9) & 
                 (df['#node1'] != 'APP') & 
                 (df['node2'] != 'APP')]

G = nx.from_pandas_edgelist(df_filtered, source='#node1', target='node2')

N = G.number_of_nodes()
L = G.number_of_edges()
avg_clustering = nx.average_clustering(G)
er_clustering = (2 * L) / (N * (N - 1)) if N > 1 else 0

print("--- APP Network Basic Statistics ---")
print(f"Total initial edges (raw): {len(df)}")
print(f"Number of nodes (N): {N}")
print(f"Number of edges (L): {L}")
print(f"Average Degree (<k>): {(2 * L / N):.2f}")
print(f"Average Clustering Coefficient (<C>): {avg_clustering:.3f}")
print(f"Equivalent Erdos-Renyi Clustering (C_rand): {er_clustering:.4f}")
print(f"Small-World Clustering Ratio (<C> / C_rand): {(avg_clustering / er_clustering):.1f}x\n")

# ==========================================
# STEP 2: Degree Distribution & Power-Law Fit (CCDF)
# ==========================================
degrees = np.array([d for n, d in G.degree() if d > 0])
unique_k, counts = np.unique(degrees, return_counts=True)
pk = counts / len(degrees)
ccdf = [np.sum(degrees >= k) / len(degrees) for k in unique_k]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Left plot: Raw P(k) showing fat-tail noise
ax1.loglog(unique_k, pk, 'o', color='royalblue', alpha=0.7)
ax1.set_title('Raw Degree Distribution P(k) (APP Network)')
ax1.set_xlabel('Degree (k)')
ax1.set_ylabel('P(k)')
ax1.grid(True, which="both", ls="--", alpha=0.5)

# Right plot: CCDF with Power-Law Fit
ax2.loglog(unique_k, ccdf, 'o-', color='crimson', markersize=5, alpha=0.8, label='Empirical CCDF')

# Fit power-law to the tail (k >= 10)
mask = unique_k >= 10
if np.sum(mask) > 2:
    log_k = np.log10(unique_k[mask])
    log_ccdf = np.log10(np.array(ccdf)[mask])
    slope, intercept = np.polyfit(log_k, log_ccdf, 1)
    gamma = -slope + 1
    print(f"Estimated Power-Law Exponent (gamma): {gamma:.2f}\n")
    ax2.loglog(unique_k[mask], 10**(intercept + slope * log_k), 'k--', lw=2, 
               label=f'Power-law tail fit ($\\gamma \\approx {gamma:.2f}$)')

ax2.set_title('Cumulative Degree Distribution (CCDF)')
ax2.set_xlabel('Degree (k)')
ax2.set_ylabel('P(K >= k)')
ax2.legend()
ax2.grid(True, which="both", ls="--", alpha=0.5)

plt.tight_layout()
plt.savefig(os.path.join(script_dir, 'degree_distribution.png'), dpi=300)
plt.show()

# ==========================================
# STEP 3: Hub Identification (Centrality Analysis)
# ==========================================
deg_centrality = nx.degree_centrality(G)
bet_centrality = nx.betweenness_centrality(G)

top_degree = sorted(deg_centrality.items(), key=lambda x: x[1], reverse=True)[:10]
top_betweenness = sorted(bet_centrality.items(), key=lambda x: x[1], reverse=True)[:10]

print("--- Top 10 Hubs by Degree Centrality (Local Connectors) ---")
for rank, (node, val) in enumerate(top_degree, 1):
    print(f"{rank}. {node}: {val:.4f} (Degree: {G.degree(node)})")

print("\n--- Top 10 Hubs by Betweenness Centrality (Network Bridges) ---")
for rank, (node, val) in enumerate(top_betweenness, 1):
    print(f"{rank}. {node}: {val:.4f} (Degree: {G.degree(node)})")

# ==========================================
# STEP 4: Community Detection & Modularity
# ==========================================
largest_cc = max(nx.connected_components(G), key=len)
G_lcc = G.subgraph(largest_cc).copy()

communities = greedy_modularity_communities(G_lcc)
mod_score = modularity(G_lcc, communities)

print(f"\n--- Community Detection (Greedy Modularity on LCC) ---")
print(f"LCC Nodes: {G_lcc.number_of_nodes()}, LCC Edges: {G_lcc.number_of_edges()}")
print(f"Number of detected communities: {len(communities)}")
print(f"Modularity Score (Q): {mod_score:.3f}")
for i, comm in enumerate(communities[:5], 1):
    sample_proteins = list(comm)[:6]
    print(f"Community {i} ({len(comm)} nodes) - Sample proteins: {sample_proteins}")

# ==========================================
# STEP 5: Network Visualization (Communities & Hubs)
# ==========================================
plt.figure(figsize=(12, 10))
pos = nx.spring_layout(G_lcc, seed=42, k=0.18)

palette = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']
node_colors = []
for node in G_lcc.nodes():
    color_assigned = '#cccccc'
    for i, comm in enumerate(communities[:5]):
        if node in comm:
            color_assigned = palette[i]
            break
    node_colors.append(color_assigned)

node_sizes = [G_lcc.degree(n) * 3.5 for n in G_lcc.nodes()]

nx.draw_networkx_edges(G_lcc, pos, alpha=0.12, edge_color='gray')
nx.draw_networkx_nodes(G_lcc, pos, node_color=node_colors, node_size=node_sizes, alpha=0.85)

# Label top 8 degree hubs
top_hub_names = [node for node, _ in top_degree[:8] if node in G_lcc]
labels = {node: node for node in top_hub_names}
nx.draw_networkx_labels(G_lcc, pos, labels, font_size=10, font_weight='bold', 
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

plt.title(f'Human APP (Alzheimer\'s) Interaction Network Communities (Modularity Q = {mod_score:.2f})', fontsize=14)
plt.axis('off')
plt.tight_layout()
plt.savefig(os.path.join(script_dir, 'network_communities.png'), dpi=300)
plt.show()
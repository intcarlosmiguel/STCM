import networkx as nx
import graph_tool.all as gt
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns
import pandas as pd
from pybib.calc import *
import os
import matplotlib.font_manager as fm
from matplotlib.gridspec import GridSpec
from pybib.template import *
from pybib.models import *
import graph_tool.all as gt

def transform_faixa(a,coluna,quartis):
    faixas = []
    for i in a[coluna]:
        for j in range(len(quartis)):
            if(i<quartis[j]):
                break
        faixas.append(j-1)
    faixas = np.array(faixas)
    a[coluna+"Faixas"] = faixas
    return a


#!pip freeze > requirements.txt
def init_polymod():
    polymod = pd.read_csv("./POLYMOD/2008_Mossong_POLYMOD_contact_common.csv")
    polymod_ = pd.read_csv("./POLYMOD/2008_Mossong_POLYMOD_hh_common.csv")
    polymod_ids = pd.read_csv("./POLYMOD/2008_Mossong_POLYMOD_participant_common.csv")

    contatos = polymod[['part_id','cont_id','cnt_age_exact','cnt_age_est_min','cnt_age_est_max']]
    contatos = contatos.merge(polymod_ids[['part_id','part_age']], on='part_id', how='inner')
    df = polymod[['phys_contact','duration_multi']]
    df = df.dropna()
    merged_df = contatos.dropna(subset=['cnt_age_exact', 'cnt_age_est_min'], how='all')[['part_id','part_age','cnt_age_exact']]
    quartis = [0,20,30,50,70,10000]
    part = np.unique(merged_df['part_id'])
    merged_df = transform_faixa(merged_df,"part_age",quartis)
    merged_df = transform_faixa(merged_df,"cnt_age_exact",quartis)
    polymod_ids = transform_faixa(polymod_ids,"part_age",quartis)
    pivot_df = merged_df.pivot_table(index='part_id', columns='cnt_age_exactFaixas', aggfunc='size', fill_value=0)
    faixas = polymod_ids[polymod_ids['part_id'].isin(part)]['part_ageFaixas'].values
    contatos = pivot_df.values
    polymod = polymod.dropna(subset=['duration_multi'])
    connections = polymod[["part_id","cont_id","cnt_age_exact","duration_multi"]]
    connections = pd.merge(polymod_ids[["part_id","part_age"]], connections, on='part_id', how='inner')
    quartis = [0,20,30,50,70,10000]
    connections = transform_faixa(connections,"part_age",quartis)
    connections = transform_faixa(connections,"cnt_age_exact",quartis)
    connections = connections[["part_id","part_ageFaixas","cnt_age_exactFaixas","duration_multi"]]
    degree = np.sum(contatos,axis = 1)
    degree = degree[degree > 0]
    POLYMOD_DEGREE = np.mean(degree)
    connections["id"] = np.arange(len(connections))+8002
    return connections, POLYMOD_DEGREE
# Create an empty directed graph
G_polymod = nx.Graph()
connections, POLYMOD_DEGREE = init_polymod()
participant_categories = dict()
# Add nodes with their categorie
unique_participants = np.unique(connections['part_id']).tolist()
for i in zip(connections['part_id'],connections['part_ageFaixas']):
    participant_categories[i[0]] = i[1]
unique_participants += np.unique(connections['id']).tolist()
for i in zip(connections['id'],connections['cnt_age_exactFaixas']):
    participant_categories[i[0]] = i[1]

# Add nodes with their category attributes
for participant in unique_participants:
    G_polymod.add_node(participant, category=participant_categories[participant])

# Add edges from connections DataFrame
for _, row in connections.iterrows():
    G_polymod.add_edge(row['part_id'], row['id'], )

# Print basic network statistics
print(f"Number of nodes: {G_polymod.number_of_nodes()}")
print(f"Number of edges: {G_polymod.number_of_edges()}")




# Create a new graph-tool Graph
g = gt.Graph(directed=False)

# Add vertex properties
vprop_category = g.new_vertex_property("int")

# Add vertices first and create a mapping
vertex_map = {i: g.add_vertex() for i in G_polymod.nodes()}
for v in G_polymod.nodes():
    new_v = int(vertex_map[v])
    vertex_map[int(v)] = new_v
    vprop_category[new_v] = participant_categories[int(v)]

# Add edges
for e in G_polymod.edges():
    g.add_edge(vertex_map[e[0]], vertex_map[e[1]])

# Add the vertex property to the graph
g.vertex_properties["category"] = vprop_category
state = gt.BlockState(g, b=g.vp["category"], deg_corr=True)
num_category = 5
N = int(30000)
K = np.zeros(N).astype(int)
freq = np.loadtxt('./input/POLYMOD/faixas.txt')
freq = np.cumsum(freq)
p = np.random.rand(N)
stratified = np.searchsorted(freq, p)

histogram = np.linspace(0, 1, 101)
H = np.zeros((len(histogram),num_category,num_category)).T
results = []  

for i in range(len(stratified)):
    distribution = np.loadtxt(f'./input/POLYMOD/distribution_{stratified[i]}_out.txt')
    K[i] = empiric_distribution(distribution, len(distribution)-1)
print(len(K),len(freq),len(p))
for i in tqdm(range(100)):
    gt.seed_rng(42*N+i)
    Grafo = gt.generate_sbm(stratified, gt.adjacency(state.get_bg(), state.get_ers()).T, K, directed=False)

    gt.remove_parallel_edges(Grafo)
    gt.remove_self_loops(Grafo)
    # Remove vértices com grau 0
    deg = Grafo.get_total_degrees(Grafo.get_vertices())
    v_filter = Grafo.new_vertex_property('bool')
    v_filter.a = deg > 0
    stratified = stratified[deg > 0]
    Grafo.set_vertex_filter(v_filter)
    Grafo.purge_vertices()
    Grafo = convert_graph(Grafo)
    contatos_sucessores = []
    contatos_predecessores = []
    for i in Grafo.nodes():
        # Initialize histograms for this node i
        hist_sucessores = np.zeros(num_category, dtype=int)

        hist_predecessores = np.zeros(num_category, dtype=int)

        if nx.is_directed(Grafo):
            for j in Grafo.successors(i): # Get nodes pointed TO
                j_cat = int(stratified[j])
                
                hist_sucessores[j_cat] += 1
            contatos_sucessores.append(list(hist_sucessores))
            for j in Grafo.predecessors(i):
                j_cat = int(stratified[j])
                
                hist_predecessores[j_cat] += 1
            contatos_predecessores.append(list(hist_predecessores))
        else:
            for j in Grafo.neighbors(i): # Get nodes pointed TO
                j_cat = int(stratified[j])
                
                hist_sucessores[j_cat] += 1
            contatos_sucessores.append(list(hist_sucessores))
    
    contatos_sucessores = np.array(contatos_sucessores)
    degree_in = np.sum(contatos_sucessores,axis = 1)
    if nx.is_directed(Grafo):
        contatos_predecessores = np.array(contatos_predecessores)
        degree_out = np.sum(contatos_predecessores,axis = 1)

    H = calculate_histogram_undirected(contatos_sucessores, stratified, degree_in, num_category, H)
    #get_parameters(Grafo, results)
for i in range(num_category):
    np.savetxt(f'./output/DCSBM/Matrix_DCSBM_POLYMOD_{N}_{i}.txt',H[i], fmt = "%f")

import networkx as nx
import graph_tool.all as gt
import numpy as np
from tqdm import tqdm
from pybib.models import *
import sys

def create_network_DCSBM(g1, N,seed):
    """
    Gera um grafo sintético usando o Modelo de Blocos Estocástico (SBM).

    Utiliza um grafo de entrada como modelo para sua estrutura de comunidades.
    O grafo gerado terá N * |V_template| vértices.

    Args:
        g1 (gt.Graph): O grafo modelo, que deve ter uma propriedade
                               de vértice interna chamada "value" para as categorias.
        N (int): O número de cópias do modelo a serem usadas para criar
                 a base estatística do novo grafo.

    Returns:
        tuple: Uma tupla contendo:
               - u (gt.Graph): O novo grafo sintético gerado.
               - value_u (gt.VertexPropertyMap): Um mapa de propriedades contendo
                                                 a categoria (bloco) de cada nó em 'u'.
    """
    # --- 1. Preparação do Modelo ---
    # Garante que estamos trabalhando com o maior componente conectado do grafo modelo.
    gt.seed_rng(42*N+seed)
    comp, hist = gt.label_components(g1, directed=False)
    if not hist.size:
        raise ValueError("O grafo modelo está vazio.")
        
    largest_comp_idx = np.argmax(hist)
    v_filter = g1.new_vertex_property('bool')
    v_filter.a = (comp.a == largest_comp_idx)
    g1 = gt.Graph(gt.GraphView(g1, vfilt=v_filter), prune=True)

    # --- 2. Construção do Super-Grafo Estatístico ---
    g_final = gt.Graph(directed=g1.is_directed())
    
    # Loop N vezes para criar um grafo com N cópias, corrigindo o bug de N+1.
    for _ in range(N):
        # CORREÇÃO CRÍTICA: Adicionado internal_props=True para copiar "value".
        g_final = gt.graph_union(g_final, g1, include=True, internal_props=True)

    # --- 3. Extração dos Parâmetros do Modelo ---
    # Acessa a propriedade de categoria que foi corretamente propagada.
    value = g_final.vp["value"]
    
    # Cria o estado do bloco para medir as estatísticas da comunidade.
    state = gt.BlockState(g_final, b=value, deg_corr=True)
    
    # --- 4. Geração do Grafo Sintético 'u' ---
    # Gera um novo grafo 'u' com o mesmo número de nós e a mesma estrutura de blocos.
    u = gt.generate_sbm(state.b.a, gt.adjacency(state.get_bg(), state.get_ers()).T,
                        g_final.degree_property_map("out").a,
                        g_final.degree_property_map("in").a, directed=True)
    
    # --- 5. Atribuição das Categorias ao Novo Grafo ---
    value_u = u.new_vertex_property("int")
    
    # FORMA VETORIZADA (MUITO RÁPIDA): Copia as categorias diretamente.
    value_u.a = state.b.a
    
    # Limpeza do grafo gerado.
    gt.remove_parallel_edges(u)
    gt.remove_self_loops(u)
    # Remove vértices com grau 0
    deg = u.get_total_degrees(u.get_vertices())
    v_filter = u.new_vertex_property('bool')
    v_filter.a = deg > 0
    u.set_vertex_filter(v_filter)
    u.purge_vertices()


    # A função agora retorna tanto o grafo quanto suas categorias.
    return u, np.array(list(value_u))


g_real = gt.collection.data["polblogs"]
num_category = 2

# Get size from command line arguments
if len(sys.argv) != 2:
    print("Usage: python DCSBM.py <size>")
    sys.exit(1)
size = int(sys.argv[1])

# --- ETAPA 1: Encontrar os parâmetros da rede real ---
histogram = np.linspace(0, 1, 101)
H = np.zeros((len(histogram),num_category,num_category)).T
H_in = np.zeros((len(histogram),num_category,num_category)).T
results = []    
for i in tqdm(range(100)):

    g_sintetico, categories = create_network_DCSBM(g_real, size,i)
    G = convert_graph(g_sintetico)
    
    # Remove self-loops
    contatos_sucessores = []
    G.remove_edges_from(nx.selfloop_edges(G))

    contatos_predecessores = []

    # Iterate through each node in the graph ONCE
    for i in G.nodes():
        # Initialize histograms for this node i
        hist_sucessores = np.zeros(2, dtype=int)
        hist_predecessores = np.zeros(2, dtype=int)
        for j in G.successors(i): # Get nodes pointed TO
            j_cat = int(categories[j])
            
            hist_sucessores[j_cat] += 1
        for j in G.predecessors(i): 
            j_cat = int(categories[j])
            
            hist_predecessores[j_cat] += 1
        contatos_sucessores.append(list(hist_sucessores))
        contatos_predecessores.append(list(hist_predecessores))
    
    contatos_sucessores = np.array(contatos_sucessores)
    contatos_predecessores = np.array(contatos_predecessores)
    degree_in = np.sum(contatos_sucessores,axis = 1)
    degree_out = np.sum(contatos_predecessores,axis = 1)

    H,H_in = calculate_histogram(contatos_predecessores,contatos_sucessores,categories,degree_in,degree_out,num_category,H,H_in)
    get_parameters(G, results)

for i in range(2):
    np.savetxt(f'./output/DCSBM/Matrix_DCSBM_{size}_{i}.txt',H[i], fmt = "%f")
    np.savetxt(f'./output/DCSBM/Matrix_DCSBM_in_{size}_{i}.txt',H_in[i], fmt = "%f")
results = np.array(results)
#results = np.mean(results,axis = 0)
np.savetxt(f'./output/DCSBM/DCSBM_{size}.txt',results, fmt = "%f")
print(f'Tamanho: {size} concluído') 

from collections import defaultdict
from scipy.stats import binom
import numpy as np
import networkx as nx
from tqdm import tqdm
#import graph_tool.all as gt
import collections
#from calc import empiric_distribution



def CM(is_directed,N,total):

    freq = []
    k_in = []
    k_out = []

    num_category = 0
    if(not is_directed):
        freq = np.loadtxt('./input/POLYMOD/faixas.txt')
        freq = np.cumsum(freq)
        num_category = len(freq)
    else :
        freq = np.loadtxt('./input/blogs2/faixas.txt')
        k_in = np.loadtxt('./input/blogs2/in.txt')
        k_out = np.loadtxt('./input/blogs2/out.txt')
        num_category = 2
    values = []
    histogram = np.linspace(0, 1, 101)
    H = np.zeros((len(histogram),num_category,num_category)).T
    H_in = np.zeros((len(histogram),num_category,num_category)).T
    if(not is_directed):
        del H_in
    for seed in range(total):
        np.random.seed(42+seed)
        p = np.random.rand(N)
        if(not is_directed):
            contatos = []
            G = generate_graph(is_directed,N,seed,freq,p)
            for i in G.nodes():
                hist = list(np.zeros(5).astype(int))
                lista = [int(categories[j].astype(int)) for j in list(G.neighbors(i))]
                for j in lista:
                    hist[j] += 1
                contatos.append(list(hist))
            contatos = np.array(contatos)
            degree = np.sum(contatos,axis = 1)
            contatos = contatos[degree > 0]
            idade_CM = idade_CM[degree > 0]
            for category in range(5):
                contagem = contatos[idade_CM == category]
                B = contagem/(np.sum(contagem,axis = 1)[:, np.newaxis])
                for i in range(5):
                    for value in B.T[i]:
                        for j in range(len(histogram)-1):
                            if(value <= histogram[j+1]):
                                H[category][i][j] += 1
                                break
            
        else:
            categories = freq.tolist()*N
            G = generate_graph(is_directed,N,seed,freq,categories,k_in,k_out)

            contatos_sucessores = []
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

            for cat in range(num_category):
                contagem = contatos_sucessores[categories == cat]
                contagem = contagem[degree_in[categories == cat] > 0]
                B = contagem/(np.sum(contagem,axis = 1)[:, np.newaxis])
                for i in range(num_category):
                    for value in B.T[i]:
                        for j in range(len(histogram)-1):
                            if(value <= histogram[j+1]):
                                H[cat][i][j] += 1
                                break
                contagem = contatos_predecessores[categories == cat]
                contagem = contagem[degree_out[categories == cat] > 0]
                B = contagem/(np.sum(contagem,axis = 1)[:, np.newaxis])
                for i in range(num_category):
                    for value in B.T[i]:
                        for j in range(len(histogram)-1):
                            if(value <= histogram[j+1]):
                                H_in[cat][i][j] += 1
                                break

        G = nx.Graph(G)

        # Calculate average degree
        average_degree = np.mean([d for n, d in G.degree()])

        # Calculate total clustering
        total_clustering = nx.transitivity(G)

        # Calculate average clustering
        average_clustering = nx.average_clustering(G)
        # Calculate diameter
        shortest_paths = dict(nx.shortest_path_length(G))
        paths = []
        diameter = 0
        maximo = 0
        for node in shortest_paths:
            path = np.array(list(shortest_paths[node].values()))
            path = path[path != 0]
            paths += path.tolist()
            if(len(path) == 0):
                continue
            maximo = np.max(path)
            if(maximo > diameter):
                diameter = maximo
        average_path_length = np.mean(paths)
        del paths 
        del path 
        del shortest_paths
        # Calculate average shortest path length

        # Calculate degree-clustering correlation
        degrees = np.array([d for n, d in G.degree()])
        clustering_coeffs = np.array([nx.clustering(G, n) for n in G.nodes()])
        correlation = np.corrcoef(degrees, clustering_coeffs)[0, 1]

        # Save results in values
        values.append([average_degree, total_clustering, average_clustering, diameter, average_path_length, correlation])
        del G 
        print(f"Rede: {seed+1}")
    if(not is_directed):
        np.savetxt(f'./output/POLYMOD/CM_{N}.txt',values, fmt = "%f")
        for i in range(5):
            np.savetxt(f'./output/POLYMOD/Matrix_CM_{N}_{i}.txt',H[i], fmt = "%f")
    else:
        np.savetxt(f'./output/blogs/CM_{N}.txt',values, fmt = "%f")
        for i in range(2):
            np.savetxt(f'./output/blogs/Matrix_CM_{N}_{i}.txt',H[i], fmt = "%f")
            np.savetxt(f'./output/blogs/Matrix_CM_in_{N}_{i}.txt',H_in[i], fmt = "%f")

def find_network_parameters(graph: gt.Graph):
    """
    Analisa uma rede direcionada com uma partição de blocos conhecida e extrai
    os parâmetros para um modelo DC-SBM.

    Args:
        graph: O grafo do graph-tool. Espera-se que ele tenha uma propriedade
               de vértice chamada 'value' contendo a partição de blocos.

    Returns:
        Um dicionário contendo os parâmetros do modelo:
        - 'propensity_matrix' (e): Matriz de afinidade entre blocos.
        - 'block_proportions': A proporção de nós em cada bloco.
        - 'in_degrees_by_block': Um dicionário mapeando ID do bloco para uma lista de graus de entrada.
        - 'out_degrees_by_block': Um dicionário mapeando ID do bloco para uma lista de graus de saída.
    """
    print("Iniciando a extração de parâmetros do grafo original...")

    # A partição é dada pela propriedade de vértice 'value' no dataset polblogs
    b_original = graph.vp.value

    # 1. Estimar a matriz de propensão 'e'
    # Criamos um BlockState com a partição conhecida e deg_corr=True
    state = gt.BlockState(graph, b=b_original, deg_corr=True)
    e_matrix = state.get_matrix().toarray()
    
    # 2. Calcular a proporção de cada bloco
    block_ids, counts = np.unique(b_original.a, return_counts=True)
    block_proportions = counts / graph.num_vertices()

    # 3. Coletar as distribuições de grau (entrada e saída) para cada bloco
    in_degs_by_block = collections.defaultdict(list)
    out_degs_by_block = collections.defaultdict(list)

    for v in graph.vertices():
        block_id = b_original[v]
        in_degs_by_block[block_id].append(v.in_degree())
        out_degs_by_block[block_id].append(v.out_degree())

    # Empacotar tudo em um dicionário para retorno
    parameters = {
        "propensity_matrix": e_matrix,
        "block_proportions": block_proportions,
        "in_degrees_by_block": dict(in_degs_by_block),
        "out_degrees_by_block": dict(out_degs_by_block),
        "block_ids": block_ids
    }
    
    print("Extração de parâmetros concluída.")
    return parameters

def convert_graph(g):
    G_nx = nx.DiGraph()

    # 4.2: Adicionando os nós e suas propriedades
    # Iteramos sobre os vértices do grafo graph-tool.
    # O 'graph-tool' usa 'Vertex' objects, enquanto o 'networkx' geralmente usa inteiros
    # ou strings como identificadores de nó. Vamos usar o índice do vértice (um inteiro).
    for v in g.vertices():
        node_id = int(v)
        # Criamos um dicionário para armazenar as propriedades do nó
        props = {}
        for key, prop_map in g.vertex_properties.items():
            props[key] = prop_map[v]
            
        G_nx.add_node(node_id, **props)

    # 4.3: Adicionando as arestas
    # Iteramos sobre as arestas do grafo graph-tool.
    # É importante converter os descritores de vértice (v1, v2) para inteiros.
    for e in g.edges():
        source_id = int(e.source())
        target_id = int(e.target())
        # No momento, não temos propriedades de aresta, mas se tivéssemos,
        # o processo seria similar ao dos nós.
        G_nx.add_edge(source_id, target_id)
    return G_nx

# ==============================================================================
# FUNÇÃO CORRIGIDA
# ==============================================================================
def generate_synthetic_network(params: dict, N: int, seed: int):
    """
    Gera uma nova rede direcionada de tamanho N usando os parâmetros de um modelo DC-SBM.

    Args:
        params: Dicionário de parâmetros retornado pela função find_network_parameters.
        N: O número de nós desejado para a nova rede.
        seed: Semente para reprodutibilidade.

    Returns:
        g_new (gt.Graph): Um novo grafo sintético do graph-tool.
        b_new_vals (np.ndarray): Um array NumPy com a atribuição de bloco para cada nó.
    """

    # 1. Criar a nova partição de blocos (como um array NumPy)
    np.random.seed(42 + seed)
    b_new_vals = np.random.choice(
        params["block_ids"],
        size=N,
        p=params["block_proportions"]
    )
    
    # 2. Gerar as novas sequências de grau (entrada e saída)
    theta_in_new = np.zeros(N)
    theta_out_new = np.zeros(N)

    for i in range(N):
        assigned_block = b_new_vals[i]
        
        in_deg_pool = params["in_degrees_by_block"][assigned_block]
        theta_in_new[i] = np.random.choice(in_deg_pool)

        out_deg_pool = params["out_degrees_by_block"][assigned_block]
        theta_out_new[i] = np.random.choice(out_deg_pool)

    # 3. Gerar o grafo
    g_new = gt.generate_sbm(
        b_new_vals,
        params["propensity_matrix"],
        out_degs=theta_out_new,
        in_degs=theta_in_new
    )
    
    return g_new, b_new_vals



def calculate_histogram(contatos_predecessores,contatos_sucessores,categories,degree_in,degree_out,num_category,H,H_in):
    histogram = np.linspace(0, 1, 101)
    for cat in range(num_category):
        contagem = contatos_sucessores[categories == cat]
        contagem = contagem[degree_in[categories == cat] > 0]
        B = contagem/(np.sum(contagem,axis = 1)[:, np.newaxis])
        for i in range(num_category):
            for value in B.T[i]:
                for j in range(len(histogram)-1):
                    if(value <= histogram[j+1]):
                        H[cat][i][j] += 1
                        break
        contagem = contatos_predecessores[categories == cat]
        contagem = contagem[degree_out[categories == cat] > 0]
        B = contagem/(np.sum(contagem,axis = 1)[:, np.newaxis])
        for i in range(num_category):
            for value in B.T[i]:
                for j in range(len(histogram)-1):
                    if(value <= histogram[j+1]):
                        H_in[cat][i][j] += 1
                        break
    return H,H_in

def fit_bccm_directed(edge_list: list[tuple], partition: dict[object, int]):
    """
    Ajusta um Block-Constrained Configuration Model (BCCM) a uma rede DIRECIONADA
    observada e a uma partição de nós pré-definida.

    Args:
        edge_list (list[tuple]): Uma lista de tuplas (origem, destino) representando
                                 as arestas direcionadas da rede.
        partition (dict): Um dicionário mapeando cada nó ao seu índice de bloco.

    Returns:
        tuple: Uma tupla contendo:
               - B_omega (np.ndarray): A matriz B x B de propensões estimadas (não simétrica).
               - metrics (dict): Um dicionário com 'log_likelihood', 'aic', 'bic'.
    """
    print("--- Iniciando o ajuste do modelo BCCM para rede DIRECIONADA ---")

    # --- Passo 0: Pré-processamento ---
    nodes = set(partition.keys())
    num_nodes = len(nodes)
    num_edges = len(edge_list)
    
    block_labels = sorted(list(set(partition.values())))
    block_map = {label: i for i, label in enumerate(block_labels)}
    num_blocks = len(block_labels)
    node_to_block_idx = {node: block_map[partition[node]] for node in nodes}

    print(f"Rede com {num_nodes} nós, {num_edges} arestas e {num_blocks} blocos.")

    # --- Passo 1.1: Calcular Graus de Entrada e Saída (k_in, k_out) ---
    out_degrees = defaultdict(int)
    in_degrees = defaultdict(int)
    for u, v in edge_list:
        out_degrees[u] += 1
        in_degrees[v] += 1

    # --- Passo 1.2: Calcular a Matriz de Arestas Observadas (A_block) ---
    # Agora A_block[α, β] conta as arestas que vão do bloco α para o bloco β.
    A_block = np.zeros((num_blocks, num_blocks), dtype=np.float64)
    for u, v in edge_list:
        b_u = node_to_block_idx[u]
        b_v = node_to_block_idx[v]
        A_block[b_u, b_v] += 1
    
    print("\nMatriz de Arestas Observadas (A_block):\n", A_block)

    # --- Passo 2: Calcular a Matriz de Arestas Potenciais (Ξ_block) ---
    # Ξ_{αβ} = (soma dos k_out no bloco α) * (soma dos k_in no bloco β)
    
    K_out_block = np.zeros(num_blocks, dtype=np.float64)
    for node, deg in out_degrees.items():
        block_idx = node_to_block_idx.get(node)
        if block_idx is not None:
            K_out_block[block_idx] += deg

    K_in_block = np.zeros(num_blocks, dtype=np.float64)
    for node, deg in in_degrees.items():
        block_idx = node_to_block_idx.get(node)
        if block_idx is not None:
            K_in_block[block_idx] += deg
            
    # O produto externo calcula todas as combinações Ξ_{αβ} de uma vez.
    Xi_block = np.outer(K_out_block, K_in_block)

    # Correção para auto-loops (nós no mesmo bloco)
    # A fórmula exata é complexa. A aproximação acima é padrão na literatura.
    # Para grafos simples, subtraímos as conexões de um nó para si mesmo.
    for i in range(num_blocks):
        sum_kout_kin_prod = 0
        nodes_in_block_i = [n for n, b_idx in node_to_block_idx.items() if b_idx == i]
        for node in nodes_in_block_i:
            sum_kout_kin_prod += out_degrees[node] * in_degrees[node]
        Xi_block[i, i] -= sum_kout_kin_prod


    print("\nMatriz de Arestas Potenciais (Xi_block):\n", Xi_block)

    # --- Passo 3: Estimar a Matriz de Propensão (B_omega) ---
    # A fórmula permanece a mesma, mas agora aplicada a matrizes não simétricas.
    with np.errstate(divide='ignore', invalid='ignore'):
        p_block = np.divide(A_block, Xi_block)
        p_block[np.isnan(p_block)] = 0

    with np.errstate(divide='ignore'):
        B_omega = -np.log(1 - p_block)

    B_omega[np.isinf(B_omega)] = 1e9
    B_omega[np.isnan(B_omega)] = 0
    
    print("\nMatriz de Propensão Estimada (B_omega):\n", B_omega)

    # --- Passo 4: Calcular Qualidade do Ajuste ---
    # A lógica da verossimilhança binomial ainda se aplica.
    log_likelihood = 0
    for i in range(num_blocks):
        for j in range(num_blocks):
            num_successes = A_block[i, j]
            num_trials = Xi_block[i, j]
            if num_trials > 0:
                prob_success = num_successes / num_trials
                log_likelihood += binom.logpmf(k=num_successes, n=num_trials, p=prob_success)

    # O número de parâmetros agora é B*B, pois a matriz não é simétrica.
    num_params = num_blocks * num_blocks

    aic = 2 * num_params - 2 * log_likelihood
    bic = np.log(num_edges) * num_params - 2 * log_likelihood

    metrics = {'log_likelihood': log_likelihood, 'aic': aic, 'bic': bic}
    
    print("\n--- Métricas de Qualidade do Ajuste ---")
    for key, value in metrics.items():
        print(f"{key.upper()}: {value:.4f}")

    return B_omega, metrics


def generate_bccm_network_directed_optimized(
    num_nodes: int,
    partition: dict[int, int],
    out_degrees: dict[int, int],
    in_degrees: dict[int, int],
    B_omega: np.ndarray,
    allow_self_loops: bool = False,
    allow_multi_edges: bool = False,
) -> list[tuple]:
    """
    Gera uma rede sintética DIRECIONADA a partir de um BCCM especificado.
    Versão otimizada usando NumPy para operações vetorizadas.

    Args:
        num_nodes (int): O número total de nós.
        partition (dict): Mapeamento nó -> bloco.
        out_degrees (dict): Mapeamento nó -> grau de saída desejado.
        in_degrees (dict): Mapeamento nó -> grau de entrada desejado.
        B_omega (np.ndarray): Matriz B x B de propensões (origem -> destino).
        allow_self_loops (bool): Permite a formação de arestas (u, u).
        allow_multi_edges (bool): Permite múltiplas arestas entre os mesmos nós.

    Returns:
        list[tuple]: Lista de arestas direcionadas (origem, destino).
    """
    print("--- Iniciando a geração de rede (versão otimizada com NumPy) ---")

    # --- Passo 0: Validação e Inicialização ---
    if sum(out_degrees.values()) != sum(in_degrees.values()):
        raise ValueError("A soma dos graus de saída deve ser igual à soma dos graus de entrada.")
    
    num_edges = sum(out_degrees.values())

    # Criar arrays NumPy de stubs
    out_stubs = np.array([node for node, deg in out_degrees.items() for _ in range(deg)], dtype=np.int32)
    in_stubs = np.array([node for node, deg in in_degrees.items() for _ in range(deg)], dtype=np.int32)
    
    np.random.shuffle(out_stubs)
    np.random.shuffle(in_stubs)
    
    # --- Passo 1: Pré-calcular vetores de mapeamento para eficiência ---
    # Mapeia cada nó ao seu bloco
    node_to_block = np.zeros(num_nodes, dtype=np.int8)
    for node, block in partition.items():
        node_to_block[node] = block
        
    # Mapeia cada stub na lista `in_stubs` ao seu respectivo bloco
    in_stub_blocks = node_to_block[in_stubs]

    print(f"Gerando uma rede com {num_nodes} nós e {num_edges} arestas direcionadas.")
    
    edge_list = []
    # Usar um set para checar duplicatas é muito mais rápido que 'in list'
    generated_edges = set() if not allow_multi_edges else None

    # --- Passo 2 & 3: Amostragem Sequencial Otimizada ---
    num_stubs_to_pair = len(out_stubs)
    for i in range(num_stubs_to_pair):
        # 1. Escolher o próximo stub de saída.
        # Processamos em ordem após o shuffle, o que é equivalente a pop(0) mas mais rápido.
        stub_u = out_stubs[i]
        
        # 2. Calcular os pesos para TODOS os stubs de entrada restantes de forma vetorizada.
        stub_u_block = node_to_block[stub_u]
        
        # Pega a linha da matriz de propensão correspondente ao bloco do stub de saída.
        # Ex: Se stub_u está no bloco 2, pegamos B_omega[2, :]
        propensity_row = B_omega[stub_u_block, :]
        
        # Usa o vetor `in_stub_blocks` para atribuir a cada stub de entrada a propensão correta.
        # Esta é a operação chave: `propensity_row[in_stub_blocks]` cria um vetor de pesos
        # com o mesmo tamanho de `in_stubs`.
        weights = propensity_row[in_stub_blocks]

        # 3. Filtrar candidatos inválidos (self-loops e multi-arestas)
        # Criar uma máscara booleana é muito mais rápido que construir listas.
        valid_mask = np.ones_like(in_stubs, dtype=bool)

        if not allow_self_loops:
            valid_mask[in_stubs == stub_u] = False
            
        # A checagem de multi-arestas ainda é a parte mais lenta.
        # Uma otimização é só fazer essa checagem se o número de candidatos for pequeno
        # ou se a rede for esparsa. Para simplificar, faremos sempre.
        if not allow_multi_edges:
            # Esta parte não é facilmente vetorizável, mas podemos torná-la mais eficiente
            # limitando a busca.
            invalid_indices = [idx for idx, stub_v in enumerate(in_stubs) if (stub_u, stub_v) in generated_edges]
            if invalid_indices:
                valid_mask[invalid_indices] = False
        
        # Aplicar a máscara aos pesos. Peso 0 para candidatos inválidos.
        weights[~valid_mask] = 0

        # Se a soma dos pesos for zero, não há candidatos válidos.
        # Isso é raro, mas pode acontecer.
        sum_weights = np.sum(weights)
        if sum_weights <= 0:
            # print(f"Aviso: Stub de saída {stub_u} não encontrou parceiro válido.")
            continue

        # Normalizar os pesos para que somem 1 (necessário para np.random.choice)
        probabilities = weights / sum_weights
        
        # 4. Sortear UM stub de entrada com base nas probabilidades calculadas.
        chosen_idx = np.random.choice(len(in_stubs), p=probabilities)
        stub_v = in_stubs[chosen_idx]
        
        # 5. Formar a aresta e registrar
        new_edge = (stub_u, stub_v)
        edge_list.append(new_edge)
        if not allow_multi_edges:
            generated_edges.add(new_edge)
        
        # 6. Remover o stub de entrada escolhido da consideração futura.
        # Em vez de `pop` (lento), trocamos o stub escolhido com o último
        # e depois encurtamos a "visão" dos arrays em 1.
        last_idx = len(in_stubs) - 1
        # Trazer o último stub para a posição do escolhido
        in_stubs[chosen_idx] = in_stubs[last_idx]
        in_stub_blocks[chosen_idx] = in_stub_blocks[last_idx]
        # Encurtar os arrays
        in_stubs = in_stubs[:last_idx]
        in_stub_blocks = in_stub_blocks[:last_idx]

        if (i + 1) % 1000 == 0:
            print(f"  ... {i + 1} de {num_edges} arestas geradas.")

    print(f"\n--- Geração Concluída. Total de {len(edge_list)} arestas criadas. ---")
    if len(in_stubs) > 0:
         print(f"Aviso: {len(in_stubs)} stubs não puderam ser pareados.")

    return edge_list

def get_parameters(G, results):
    average_degree = np.mean([d for n, d in G.degree()])

    # Calculate total clustering
    total_clustering = nx.transitivity(G)

    # Calculate average clustering
    average_clustering = nx.average_clustering(G)
    # Calculate diameter
    shortest_paths = dict(nx.shortest_path_length(G))
    paths = []
    diameter = 0
    maximo = 0
    for node in shortest_paths:
        path = np.array(list(shortest_paths[node].values()))
        path = path[path != 0]
        paths += path.tolist()
        if(len(path) == 0):
            continue
        maximo = np.max(path)
        if(maximo > diameter):
            diameter = maximo
    average_path_length = np.mean(paths)
    del paths 
    del path 
    del shortest_paths
    # Calculate average shortest path length

    # Calculate degree-clustering correlation
    degrees = np.array([d for n, d in G.degree()])
    clustering_coeffs = np.array([nx.clustering(G, n) for n in G.nodes()])
    if(np.std(degrees) == 0 or np.std(clustering_coeffs) == 0):
        correlation = 0.0
    else:
        correlation = np.corrcoef(degrees, clustering_coeffs)[0, 1]

    # Save results in values
    results.append([average_degree, total_clustering, average_clustering, diameter, average_path_length, correlation])

def generate_graph(is_directed,N,seed,freq,categories,k_in = None,k_out = None):
    if(not is_directed):
        categories = np.searchsorted(freq, categories)
        K = np.zeros(N).astype(int)
        for i in range(len(categories)):
            distribution = np.loadtxt(f'./input/POLYMOD/distribution_{categories[i]}_out.txt')
            K[i] = empiric_distribution(distribution, len(distribution)-1)
        if(np.sum(K)%2 != 0):
            K[0] += 1
        G = nx.configuration_model(K,seed =42+seed)
        return G
    else:
        k_out = k_out.tolist()*N
        k_in = k_in.tolist()*N 
        categories = freq.tolist()*N

        k_out = np.array(k_out).astype(int)
        categories = np.array(categories).astype(int)
        k_out = np.sum(k_out,axis = 1)
        k_in = np.array(k_in).astype(int)
        k_in = np.sum(k_in,axis = 1)
        G = nx.directed_configuration_model(k_out,k_in,seed=42+seed)
        return G
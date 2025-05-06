import numpy as np
import networkx as nx
from calc import *
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
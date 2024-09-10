import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def load_gml(file_path):
    labels = {}
    edges = []
    is_edge = False
    edges = []
    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        for line in lines:
            if(not is_edge):
                if('id' in line):
                    try:
                        id = int(line.split('id ')[1])
                        df[id] = ''
                    except:
                        continue
                if('value' in line):
                    try:
                        df[id] = int(line.split('value ')[1])
                    except:
                        continue
            else:
                if('source' in line):
                    try:
                        id = int(line.split('source ')[1])
                        aux = [0]*2
                        aux[0]= id
                    except:
                        continue
                if('target' in line):
                    try:
                        id = int(line.split('target ')[1])
                        aux[1]= id
                        edges.append(aux)
                    except:
                        continue
            if('edge ' in line):
                is_edge = True
    G = nx.DiGraph()
    G.add_edges_from(edges)
    for no, peso in df.items():
        G.add_node(no, weight=peso)
    return G
def apagar_nos_com_grau_zero(grafo):
    """
    Remove todos os nós com grau zero de um grafo.
    
    Parâmetros:
    - grafo: um grafo do networkx.
    
    Retorna:
    - O grafo sem os nós com grau zero.
    """
    # Identifica os nós com grau zero
    nos_grau_zero = [no for no, grau in grafo.degree() if grau == 0]
    
    # Remove os nós com grau zero
    grafo.remove_nodes_from(nos_grau_zero)
    
    return grafo

def calcular_graus_por_categoria(grafo, categoria):
    """
    Calcula o grau de entrada (in-degree) e o grau de saída (out-degree) para nós de uma categoria específica.
    
    Parâmetros:
    - grafo: um grafo direcionado (DiGraph) do networkx.
    - categoria: o valor do peso da categoria (0 ou 1).
    
    Retorna:
    - Um dicionário com o nó como chave e uma tupla (in-degree, out-degree) como valor.
    """
    graus_in = []
    graus_out = []
    
    # Itera sobre todos os nós no grafo
    for no in grafo.nodes(data=True):
        id_no = no[0]
        peso_no = no[1].get('weight', None)  # Pega o peso do nó
        
        # Verifica se o nó pertence à categoria especificada
        if peso_no == categoria:
            # Calcula o in-degree e o out-degree do nó
            in_degree = grafo.in_degree(id_no)
            out_degree = grafo.out_degree(id_no)
            
            # Armazena os graus de entrada e saída no dicionário
            graus_in.append(in_degree)
            graus_out.append(out_degree)
    
    return np.array(graus_in),np.array(graus_out)
def histogram(x):
    arr = np.arange(1,np.max(x))
    hist = np.zeros(len(arr))
    for i in x:
        for j in range(len(arr)):
            if(i < arr[j]):
                hist[j] += 1
                break
    arr = arr[hist!=0]
    hist = hist[hist!=0]
    return arr,hist
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pandas as pd
from sklearn.linear_model import LinearRegression
import os

def process_files_to_dataframe(directory):
    # Prepare an empty list to store the data
    data = []
    
    # Iterate through files in the directory
    for file_name in os.listdir(directory):
            # Extract components from the file name
        parts = file_name.replace("resultados_", "").replace(".txt", "").split("_")
        base, network_size = parts[0], parts[1]
        probability = parts[2] if len(parts) > 2 else 'N/A'
        # Append the extracted data
        data.append({
            "file_name": file_name,
            "network_size": int(network_size),
            "probability": float(probability) if probability != 'N/A' else 'N/A'    
        })
    
    # Convert the data into a DataFrame
    df = pd.DataFrame(data)
    return df


def load_gml(file_path):
    df = {}
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
    arr = np.arange(np.max(x))
    hist = np.zeros(len(arr))
    for i in x:
        for j in range(len(arr)):
            if(i <= arr[j]):
                hist[j] += 1
                break
    arr = arr[hist!=0]
    hist = hist[hist!=0]
    return arr,hist
def contar_ligacoes_in_por_categoria(grafo):
    """
    Cria um vetor onde cada coluna representa o número de ligações de entrada 
    que o nó tem com cada categoria de nós (peso 0 ou peso 1).
    
    Parâmetros:
    - grafo: Um grafo direcionado (DiGraph) do networkx onde cada nó tem peso 0 ou 1.
    
    Retorna:
    - Um dicionário onde a chave é o nó e o valor é uma tupla (ligacoes_com_peso_0, ligacoes_com_peso_1).
    """
    resultado = {}
    
    # Itera sobre todos os nós no grafo
    for no in grafo.nodes:
        ligacoes_com_peso_0 = 0
        ligacoes_com_peso_1 = 0
        
        # Pega os nós que têm arestas direcionadas para o nó atual (in-degree)
        nos_de_entrada = grafo.predecessors(no)
        
        # Conta quantas ligações de entrada vêm de nós com peso 0 ou 1
        for predecessor in nos_de_entrada:
            peso_predecessor = grafo.nodes[predecessor].get('weight', None)
            if peso_predecessor == 0:
                ligacoes_com_peso_0 += 1
            elif peso_predecessor == 1:
                ligacoes_com_peso_1 += 1
        
        # Armazena o resultado no dicionário
        resultado[no] = (ligacoes_com_peso_0, ligacoes_com_peso_1)
    
    return np.array(list(resultado.values()))
def generate_distribution_byfaixas(contagem,faixas):

    graus = np.sum(contagem,axis = 1)
    contagem = contagem[graus>0]
    faixas = faixas[graus>0]
    graus = graus[graus>0]
    B = contagem/(np.sum(contagem,axis = 1)[:, np.newaxis])
    x  = np.copy(B)
    B = [np.mean(B[faixas == i,:],axis = 0) for i in np.unique(faixas)]
    return B,x

def empiric_distribution(distribution, limit):
    r = np.random.rand()  # Gera número aleatório uniforme em [0, 1)
    n = 0
    while r > distribution[n]:
        n += 1
        if n > limit:
            return 0
    return n

def contar_ligacoes_out_por_categoria(grafo):
    """
    Cria um vetor onde cada coluna representa o número de ligações de saída 
    que o nó tem com cada categoria de nós (peso 0 ou peso 1).
    
    Parâmetros:
    - grafo: Um grafo direcionado (DiGraph) do networkx onde cada nó tem peso 0 ou 1.
    
    Retorna:
    - Um dicionário onde a chave é o nó e o valor é uma tupla (ligacoes_com_peso_0, ligacoes_com_peso_1).
    """
    resultado = {}
    
    # Itera sobre todos os nós no grafo
    for no in grafo.nodes:
        ligacoes_com_peso_0 = 0
        ligacoes_com_peso_1 = 0
        
        # Pega os nós que têm arestas direcionadas a partir do nó atual (out-degree)
        nos_de_saida = grafo.successors(no)
        
        # Conta quantas ligações de saída vão para nós com peso 0 ou 1
        for sucessor in nos_de_saida:
            peso_sucessor = grafo.nodes[sucessor].get('weight', None)
            if peso_sucessor == 0:
                ligacoes_com_peso_0 += 1
            elif peso_sucessor == 1:
                ligacoes_com_peso_1 += 1
        
        # Armazena o resultado no dicionário
        resultado[no] = (ligacoes_com_peso_0, ligacoes_com_peso_1)
    
    return np.array(list(resultado.values()))
def cumulative_distribution(x):
    x = np.sort(x)
    y = np.arange(1,len(x)+1)/len(x)
    return x,y

def get_faixa(idade):
    faixas = np.array([ 20,30,50,70,1e8])
    indices = np.searchsorted(faixas,idade,side='right')
    indices= np.where(indices < len(faixas),indices,None)
    return indices
def LM(x,y,axs,color = 'reds'):
    X = x.reshape(-1,1)
    regressor = LinearRegression()
    regressor.fit(X, y)
    coeficientes = regressor.coef_[0]
    intercepto = regressor.intercept_

    # Calcular o coeficiente de determinação (R²)
    r2 = regressor.score(X, y)
    if(axs != -1):
        axs.plot(x,x*coeficientes + intercepto,c = color,label = 'R² :{:.3f}, {:.3f}x + {:.3f}'.format(r2,coeficientes,intercepto))
    return coeficientes,intercepto,r2
# MCE — Modelos de Contato Estruturado / Structured Contact Models

> **Português** | [English version below](#english-version)

---

## Visão Geral

Este repositório implementa e compara diferentes **modelos de geração de redes de contato estruturado** com aplicação em simulações epidemiológicas. O projeto utiliza dados empíricos de contato humano (POLYMOD) e redes sociais reais (Political Blogs), e os emprega para calibrar três famílias distintas de modelos generativos:

| Modelo | Descrição |
|---|---|
| **MCE / STCM** (C) | *Stochastic Configuration Model* — geração de redes com estrutura de blocos por categorias de idade ou político usando a biblioteca igraph |
| **DC-SBM** (Python) | Degree-Corrected Stochastic Block Model — geração de redes sintéticas baseadas em blocos com correção de grau, usando `graph-tool` |
| **gHypEG** (R) | Generalized Hypergeometric Ensemble of Graphs — ajuste de um modelo probabilístico ao grafo observado e geração de amostras sintéticas |

---

## Estrutura do Repositório

```
MCE/
├── main.c                    # Ponto de entrada do programa C (MCE/LCM)
├── compile.sh                # Script de compilação auxiliar
├── roda.sh                   # Script de execução do programa C
├── DCSBM.py                  # Script principal do modelo DC-SBM (rede polblogs)
├── DCSBM_POLYMOD.py          # Script DC-SBM adaptado para dados POLYMOD
├── ghype.r                   # Script gHypEG para redes tipo polblogs
├── ghype_polymod.r            # Script gHypEG para dados POLYMOD
├── transform.py              # Utilitário de conversão de PDFs para imagens PNG
├── requirements.txt          # Dependências Python
│
├── bib/                      # Biblioteca C (cabeçalhos)
│   ├── STCM.h                # Implementação principal do Stochastic Configuration Model
│   ├── calc.h                # Funções matemáticas auxiliares em C (todas funções `static`)
│   ├── define.h              # Constantes epidemiológicas e parâmetros globais
│   └── mtwister.h            # Gerador de números aleatórios Mersenne Twister (64-bit)
│
├── pybib/                    # Biblioteca Python
│   ├── models.py             # Funções de geração e ajuste de modelos de rede
│   ├── calc.py               # Funções de análise e cálculo de propriedades de redes
│   ├── cm.py                 # Modelo de Configuração (Configuration Model)
│   └── template.py           # Templates e utilitários de plotagem
│
├── POLYMOD/                  # Dados brutos do estudo POLYMOD (2008, Mossong et al.)
│   ├── 2008_Mossong_POLYMOD_contact_common.csv
│   ├── 2008_Mossong_POLYMOD_participant_common.csv
│   ├── 2008_Mossong_POLYMOD_hh_common.csv
│   └── 2008_Mossong_POLYMOD_dictionary.xls
│
├── input/                    # Dados de entrada pré-processados, por modelo
│   ├── POLYMOD/              # Distribuições de grau e probabilidades (POLYMOD)
│   ├── polblogs/             # Distribuições de grau e probabilidades (polblogs)
│   ├── gHypEG/               # Arquivos de nós e arestas para o modelo gHypEG
│   ├── ComesF/               # Dados da rede ComesF
│   └── test/                 # Dados para testes unitários do modelo C
│
├── output/                   # Resultados gerados pelos modelos
│   ├── DCSBM/                # Matrizes e métricas do modelo DC-SBM
│   ├── POLYMOD/              # Matrizes e métricas do modelo MCE (POLYMOD)
│   ├── blogs/                # Matrizes e métricas do modelo MCE (polblogs)
│   └── gHypEG/               # Matrizes e métricas do modelo gHypEG
│
├── analysis.ipynb            # Notebook Jupyter para análise e visualização dos resultados
├── img/                      # Figuras e imagens geradas
└── font/                     # Fontes utilizadas nas visualizações
```

---

## Dependências

### Programa C (MCE/LCM)

- **GCC** (com suporte a OpenMP: `-fopenmp`)
- **igraph** (≥ 0.10): biblioteca de análise de redes ([igraph.org](https://igraph.org/c/))
- **libm** (padrão C)
- **libstdc++** (padrão)

### Scripts Python

Instale as dependências via pip:

```bash
pip install -r requirements.txt
```

Dependências principais:
- `networkx` — manipulação de grafos
- `graph-tool` — geração eficiente de redes SBM (**instalação separada**, ver [graph-tool.skewed.de](https://graph-tool.skewed.de))
- `numpy`, `scipy`, `pandas`
- `matplotlib`, `seaborn`
- `tqdm`

> **Nota:** O `graph-tool` **não** está disponível via pip. Instale-o pelo gerenciador de pacotes do sistema ou via conda.

### Scripts R

```r
install.packages(c("ghypernet", "igraph"))
```

---

## Compilação e Execução

### 1. Programa C — MCE (Stochastic Configuration Model — STCM)

#### Compilação

```bash
gcc -I/caminho/para/igraph/build/include \
    -I/caminho/para/igraph/include \
    main.c -o main \
    -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
```

Substitua `/caminho/para/igraph/` pelo caminho real de instalação da biblioteca igraph em seu sistema. O arquivo `roda.sh` contém um exemplo de comando de compilação configurado para o ambiente de desenvolvimento original.

#### Execução

```bash
./main <seed> <N> <is_undirect> <redes> <prob>
```

| Parâmetro | Tipo | Descrição |
|---|---|---|
| `seed` | inteiro | Semente para o gerador de números aleatórios |
| `N` | inteiro | Número de nós da rede |
| `is_undirect` | 0 ou 1 | `1` para grafo não-direcionado (POLYMOD), `0` para direcionado (polblogs) |
| `redes` | inteiro | Número de redes a serem geradas |
| `prob` | inteiro (%) | Probabilidade de conexão triangular (0 a 100) |

**Exemplo:**

```bash
./main 6437 30000 1 100 0
```

Gera 100 redes não-direcionadas com 30.000 nós, semente 6437 e probabilidade 0%. Internamente invoca `generate_stochastic_configuration_model()`.

#### Arquivos de Saída (Programa C)

Os resultados são salvos de forma incremental (append) durante a execução:

- **Modo não-direcionado (POLYMOD):**
  - `output/POLYMOD/matrix/out_matrix_<N>_<prob>.txt` — matriz de conexões de saída por faixa etária
  
- **Modo direcionado (polblogs):**
  - `output/blogs/matrix/out_matrix_<N>_<prob>.txt` — conexões de saída
  - `output/blogs/matrix/in_matrix_<N>_<prob>.txt` — conexões de entrada

---

### 2. Scripts Python — DC-SBM

#### `DCSBM.py` — Modelo para a rede Political Blogs

Gera 100 redes sintéticas usando o DC-SBM calibrado na rede `polblogs` (embutida no `graph-tool`).

```bash
python DCSBM.py <tamanho>
```

- `<tamanho>`: fator multiplicador de escala da rede (número de cópias do grafo modelo).

**Arquivos de Saída:**

Salvos em `output/DCSBM/`:
- `Matrix_DCSBM_<tamanho>_<categoria>.txt` — histograma de proporções de conexão de saída por categoria
- `Matrix_DCSBM_in_<tamanho>_<categoria>.txt` — histograma de proporções de conexão de entrada por categoria
- `DCSBM_<tamanho>.txt` — métricas das redes geradas (grau médio, clustering, diâmetro, etc.)

#### `DCSBM_POLYMOD.py` — Modelo para dados POLYMOD

Gera 100 redes sintéticas não-direcionadas com 30.000 nós, calibradas nos dados de contato do estudo POLYMOD, estratificadas em 5 faixas etárias (0–20, 20–30, 30–50, 50–70, 70+).

Requer os arquivos de entrada em `input/POLYMOD/`:
- `faixas.txt` — distribuição de frequência das faixas etárias
- `distribution_<i>_out.txt` — distribuição de grau para a faixa `i`

```bash
python DCSBM_POLYMOD.py
```

**Arquivos de Saída:**

Salvos em `output/DCSBM/`:
- `Matrix_DCSBM_POLYMOD_30000_<categoria>.txt` — histograma de proporções de conexão por categoria etária

---

### 3. Scripts R — gHypEG

#### `ghype.r` — Modelo para redes tipo polblogs

Ajusta o modelo gHypEG à rede de Blogs Políticos e gera 100 redes sintéticas.

```bash
Rscript ghype.r <tamanho>
```

Requer arquivos de entrada em `input/gHypEG/`:
- `nodes_<tamanho>.txt` — lista de nós com categorias
- `edges_<tamanho>.txt` — lista de arestas

**Arquivos de Saída:**

Salvos em `output/gHypEG/`:
- `metrics/gHypEG_<tamanho>.txt` — métricas das redes geradas
- `matrix/Matrix_gHypEG_out_<categoria>_<tamanho>.txt` — histogramas de conexões de saída
- `matrix/Matrix_gHypEG_in_<categoria>_<tamanho>.txt` — histogramas de conexões de entrada

#### `ghype_polymod.r` — Modelo para dados POLYMOD

Ajusta o modelo gHypEG aos dados POLYMOD (não-direcionado).

```bash
Rscript ghype_polymod.r
```

Requer arquivos de entrada em `input/gHypEG/`:
- `polymod_nodes.txt`
- `polymod_edges.txt`

**Arquivos de Saída:**

Salvos em `output/gHypEG/`:
- `metrics/gHypEG_polymod.txt`
- `matrix/Matrix_gHypEG_out_<categoria>_polymod.txt`
- `matrix/Matrix_gHypEG_in_<categoria>_polymod.txt`

---

### 4. Utilitário — `transform.py`

Converte todos os arquivos PDF de uma pasta para imagens PNG redimensionadas.

```bash
python transform.py
```

- **Entrada:** `img/pdfs/` — pasta com arquivos PDF
- **Saída:** `img/pngs/` — imagens PNG com largura padronizada de 600px

---

### 5. Análise — `analysis.ipynb`

Notebook Jupyter para carregamento, análise estatística e visualização dos resultados gerados pelos três modelos. Execute com:

```bash
jupyter notebook analysis.ipynb
```

---

## Constantes Epidemiológicas (`bib/define.h`)

O modelo C utiliza os seguintes parâmetros de compartimento SEIR/SEIHR:

| Parâmetro | Valor | Descrição |
|---|---|---|
| `beta1` | 0.5 | Taxa de transmissão (fase 1) |
| `beta2` | 0.41 | Taxa de transmissão (fase 2) |
| `sigma` | 1/5.1 | Taxa de progressão da exposição |
| `gamma_A` | 1/7 | Taxa de recuperação (assintomático) |
| `gamma_I` | 1/7 | Taxa de recuperação (sintomático) |
| `gamma_H` | 1/12.1 | Taxa de recuperação (hospitalizado) |
| `delta` | 1/13.7 | Taxa de mortalidade |
| `recupera` | 1/40 | Taxa de perda de imunidade |
| `dias` | 465 | Duração da simulação (dias) |

---

## Dados

Os dados brutos do estudo **POLYMOD** (Mossong et al., 2008) estão disponíveis na pasta `POLYMOD/` e correspondem a um levantamento de contatos humanos em oito países europeus. O conjunto de dados da rede de **Political Blogs** (Adamic & Glance, 2005) está disponível através do `graph-tool` e do arquivo `polblogs.zip`.

---

---

# English Version

## Overview

This repository implements and compares different **structured contact network generation models** with applications in epidemiological simulations. The project uses empirical human contact data (POLYMOD) and real-world social networks (Political Blogs) to calibrate three distinct families of generative models:

| Model | Description |
|---|---|
| **MCE / STCM** (C) | Stochastic Configuration Model — generates block-structured networks by age or political categories using the igraph library |
| **DC-SBM** (Python) | Degree-Corrected Stochastic Block Model — generates synthetic block-structured networks with degree correction, using `graph-tool` |
| **gHypEG** (R) | Generalized Hypergeometric Ensemble of Graphs — fits a probabilistic model to an observed graph and generates synthetic samples |

---

## Repository Structure

```
MCE/
├── main.c                    # Entry point for the C program (MCE/LCM)
├── compile.sh                # Auxiliary compilation script
├── roda.sh                   # Execution script for the C program
├── DCSBM.py                  # Main script for the DC-SBM model (polblogs network)
├── DCSBM_POLYMOD.py          # DC-SBM script adapted for POLYMOD data
├── ghype.r                   # gHypEG script for polblogs-type networks
├── ghype_polymod.r           # gHypEG script for POLYMOD data
├── transform.py              # Utility to convert PDFs to PNG images
├── requirements.txt          # Python dependencies
│
├── bib/                      # C library (header files)
│   ├── STCM.h                # Core implementation of the Stochastic Configuration Model
│   ├── calc.h                # Auxiliary mathematical functions in C (all functions declared `static`)
│   ├── define.h              # Epidemiological constants and global parameters
│   └── mtwister.h            # 64-bit Mersenne Twister random number generator
│
├── pybib/                    # Python library
│   ├── models.py             # Network model generation and fitting functions
│   ├── calc.py               # Network property analysis and calculation functions
│   ├── cm.py                 # Configuration Model implementation
│   └── template.py           # Plotting templates and utilities
│
├── POLYMOD/                  # Raw data from the POLYMOD study (2008, Mossong et al.)
│   ├── 2008_Mossong_POLYMOD_contact_common.csv
│   ├── 2008_Mossong_POLYMOD_participant_common.csv
│   ├── 2008_Mossong_POLYMOD_hh_common.csv
│   └── 2008_Mossong_POLYMOD_dictionary.xls
│
├── input/                    # Pre-processed input data, organized by model
│   ├── POLYMOD/              # Degree distributions and probabilities (POLYMOD)
│   ├── polblogs/             # Degree distributions and probabilities (polblogs)
│   ├── gHypEG/               # Node and edge files for the gHypEG model
│   ├── ComesF/               # ComesF network data
│   └── test/                 # Test data for C model unit tests
│
├── output/                   # Results generated by the models
│   ├── DCSBM/                # DC-SBM model matrices and metrics
│   ├── POLYMOD/              # MCE model matrices and metrics (POLYMOD)
│   ├── blogs/                # MCE model matrices and metrics (polblogs)
│   └── gHypEG/               # gHypEG model matrices and metrics
│
├── analysis.ipynb            # Jupyter notebook for analysis and visualization
├── img/                      # Generated figures and images
└── font/                     # Fonts used in visualizations
```

---

## Dependencies

### C Program (MCE/LCM)

- **GCC** (with OpenMP support: `-fopenmp`)
- **igraph** (≥ 0.10): network analysis library ([igraph.org](https://igraph.org/c/))
- **libm** (standard C)
- **libstdc++** (standard)

### Python Scripts

Install dependencies via pip:

```bash
pip install -r requirements.txt
```

Key dependencies:
- `networkx` — graph manipulation
- `graph-tool` — efficient SBM network generation (**separate installation**, see [graph-tool.skewed.de](https://graph-tool.skewed.de))
- `numpy`, `scipy`, `pandas`
- `matplotlib`, `seaborn`
- `tqdm`

> **Note:** `graph-tool` is **not** available via pip. Install it through your system package manager or conda.

### R Scripts

```r
install.packages(c("ghypernet", "igraph"))
```

---

## Compilation and Execution

### 1. C Program — MCE (Stochastic Configuration Model — STCM)

#### Compilation

```bash
gcc -I/path/to/igraph/build/include \
    -I/path/to/igraph/include \
    main.c -o main \
    -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
```

Replace `/path/to/igraph/` with the actual installation path of the igraph library on your system. The `roda.sh` file contains a compilation command example configured for the original development environment.

#### Execution

```bash
./main <seed> <N> <is_undirect> <redes> <prob>
```

| Parameter | Type | Description |
|---|---|---|
| `seed` | integer | Random number generator seed |
| `N` | integer | Number of network nodes |
| `is_undirect` | 0 or 1 | `1` for undirected graph (POLYMOD), `0` for directed (polblogs) |
| `redes` | integer | Number of networks to generate |
| `prob` | integer (%) | Triangular connection probability (0 to 100) |

**Example:**

```bash
./main 6437 30000 1 100 0
```

Generates 100 undirected networks with 30,000 nodes, seed 6437 and probability 0%. Internally calls `generate_stochastic_configuration_model()`.

#### Output Files (C Program)

Results are saved incrementally (append mode) during execution:

- **Undirected mode (POLYMOD):**
  - `output/POLYMOD/matrix/out_matrix_<N>_<prob>.txt` — outgoing connection matrix by age group

- **Directed mode (polblogs):**
  - `output/blogs/matrix/out_matrix_<N>_<prob>.txt` — outgoing connections
  - `output/blogs/matrix/in_matrix_<N>_<prob>.txt` — incoming connections

---

### 2. Python Scripts — DC-SBM

#### `DCSBM.py` — Model for the Political Blogs Network

Generates 100 synthetic networks using the DC-SBM calibrated on the `polblogs` network (embedded in `graph-tool`).

```bash
python DCSBM.py <size>
```

- `<size>`: scale multiplier factor (number of copies of the template graph).

**Output Files:**

Saved to `output/DCSBM/`:
- `Matrix_DCSBM_<size>_<category>.txt` — histogram of outgoing connection proportions per category
- `Matrix_DCSBM_in_<size>_<category>.txt` — histogram of incoming connection proportions per category
- `DCSBM_<size>.txt` — network metrics (average degree, clustering, diameter, etc.)

#### `DCSBM_POLYMOD.py` — Model for POLYMOD Data

Generates 100 synthetic undirected networks with 30,000 nodes, calibrated on POLYMOD contact data, stratified into 5 age groups (0–20, 20–30, 30–50, 50–70, 70+).

Requires input files in `input/POLYMOD/`:
- `faixas.txt` — age group frequency distribution
- `distribution_<i>_out.txt` — degree distribution for age group `i`

```bash
python DCSBM_POLYMOD.py
```

**Output Files:**

Saved to `output/DCSBM/`:
- `Matrix_DCSBM_POLYMOD_30000_<category>.txt` — histogram of connection proportions per age category

---

### 3. R Scripts — gHypEG

#### `ghype.r` — Model for polblogs-type Networks

Fits the gHypEG model to the Political Blogs network and generates 100 synthetic networks.

```bash
Rscript ghype.r <size>
```

Requires input files in `input/gHypEG/`:
- `nodes_<size>.txt` — node list with categories
- `edges_<size>.txt` — edge list

**Output Files:**

Saved to `output/gHypEG/`:
- `metrics/gHypEG_<size>.txt` — metrics of generated networks
- `matrix/Matrix_gHypEG_out_<category>_<size>.txt` — outgoing connection histograms
- `matrix/Matrix_gHypEG_in_<category>_<size>.txt` — incoming connection histograms

#### `ghype_polymod.r` — Model for POLYMOD Data

Fits the gHypEG model to POLYMOD data (undirected).

```bash
Rscript ghype_polymod.r
```

Requires input files in `input/gHypEG/`:
- `polymod_nodes.txt`
- `polymod_edges.txt`

**Output Files:**

Saved to `output/gHypEG/`:
- `metrics/gHypEG_polymod.txt`
- `matrix/Matrix_gHypEG_out_<category>_polymod.txt`
- `matrix/Matrix_gHypEG_in_<category>_polymod.txt`

---

### 4. Utility — `transform.py`

Converts all PDF files in a folder to resized PNG images.

```bash
python transform.py
```

- **Input:** `img/pdfs/` — folder with PDF files
- **Output:** `img/pngs/` — PNG images with standardized 600px width

---

### 5. Analysis — `analysis.ipynb`

Jupyter notebook for loading, statistical analysis, and visualization of results generated by all three models. Run with:

```bash
jupyter notebook analysis.ipynb
```

---

## Epidemiological Constants (`bib/define.h`)

The C model uses the following SEIR/SEIHR compartment parameters:

| Parameter | Value | Description |
|---|---|---|
| `beta1` | 0.5 | Transmission rate (phase 1) |
| `beta2` | 0.41 | Transmission rate (phase 2) |
| `sigma` | 1/5.1 | Exposure progression rate |
| `gamma_A` | 1/7 | Recovery rate (asymptomatic) |
| `gamma_I` | 1/7 | Recovery rate (symptomatic) |
| `gamma_H` | 1/12.1 | Recovery rate (hospitalized) |
| `delta` | 1/13.7 | Mortality rate |
| `recupera` | 1/40 | Immunity loss rate |
| `dias` | 465 | Simulation duration (days) |

---

## Data

The raw **POLYMOD** study data (Mossong et al., 2008) is available in the `POLYMOD/` folder and corresponds to a human contact survey conducted across eight European countries. The **Political Blogs** network dataset (Adamic & Glance, 2005) is available through `graph-tool` and the `polblogs.zip` file.

---

## License

See [LICENSE.txt](LICENSE.txt) for details.

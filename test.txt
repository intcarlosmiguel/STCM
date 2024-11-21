#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <omp.h>

#include <igraph.h>
#include "mtwister.h"
#include "calc.h"
struct Graph{
    int **viz;
    int Nodes;
    double edges;
};

struct constant{
    double** prob_in;
    double** prob_out;
    double** distribution_in;
    double** distribution_out;
    double* p_faixas;
};

struct mce{
    bool is_undirect;
    int** degree_in;
    int** degree_out;
    int* n_faixas;
    int** site_per_faixas;
};

void free_CONSTANT(struct constant *CONSTANT,int categorias,bool is_undirect){
    int i;
    for (i = 0; i < categorias; i++){
        free(CONSTANT->prob_in[i]);
        free(CONSTANT->prob_out[i]);
        free(CONSTANT->distribution_in[i]);
        free(CONSTANT->distribution_out[i]);
    }
    free(CONSTANT->prob_in);
    free(CONSTANT->prob_out);
    free(CONSTANT->distribution_in);
    free(CONSTANT->distribution_out);
    free(CONSTANT->p_faixas);
}

void init_MCE(struct constant *CONSTANT,struct mce *MCE,struct Graph *G,int N,igraph_vector_t* faixas,bool is_undirect,int categorias){
    int i,j,k;
    double r;
    int* grau_in = (int*) calloc(N,sizeof(int));
    int* arr = (int*) calloc(N,sizeof(int));
    int* grau_out = (int*) calloc(N,sizeof(int));

    MCE->n_faixas = (int *) calloc(categorias, sizeof(int));
    MCE->site_per_faixas = (int**) malloc(categorias*sizeof(int*));
    MCE->degree_in = (int**) malloc(N*sizeof(int*));
    MCE->degree_out = (int**) malloc(N*sizeof(int*));
    for (i = 0; i < categorias; i++) MCE->site_per_faixas[i] = (int*) malloc(0*sizeof(int));
    for (i = 0; i < N; i++){
        j = 0;
        r = genrand64_real1();
        while(r > CONSTANT->p_faixas[j]) j++;
        VECTOR(*faixas)[i] = j;
        MCE->n_faixas[j]++;
        MCE->site_per_faixas[j] = realloc(MCE->site_per_faixas[j],MCE->n_faixas[j]*sizeof(int));
        MCE->site_per_faixas[j][MCE->n_faixas[j] - 1] = i;
        k = 0;

        MCE->degree_in[i] = (int*) calloc(categorias,sizeof(int));

        while(k == 0) k = empiric_distribution(CONSTANT->distribution_in[ (int) VECTOR(*faixas)[i]]);
        grau_in[i] = k;
        if(!is_undirect){
            MCE->degree_out[i] = (int*) calloc(categorias,sizeof(int));
            k = 0;
            while(k == 0)k = empiric_distribution(CONSTANT->distribution_out[ (int) VECTOR(*faixas)[i]]);
            grau_out[i] = k;
        }
        arr[i] = i;
    }
    int* temp_weights = (int*) calloc(N,sizeof(int));
    for (int i = 0; i < N; i++)temp_weights[i] = grau_in[i];
    mergeSort(temp_weights,arr,0,N-1);
    free(temp_weights);
    int site;
    for (i = 0; i < N; i++){

        site = arr[i];

        generate_multinomial(grau_in[site], categorias, CONSTANT->prob_in[(int) VECTOR(*faixas)[site]], MCE->degree_in[i]);

        G->viz[i] =(int*) malloc((grau_in[site]+1)*sizeof(int));
        G->viz[i][0] = 0;
        if(!is_undirect)generate_multinomial(grau_out[site], categorias, CONSTANT->prob_out[(int) VECTOR(*faixas)[site]], MCE->degree_out[i]);
    }
    free(grau_in);
}

void init_CONSTANT(struct constant *CONSTANT,int categorias,bool is_undirect){

    CONSTANT->p_faixas = (double *) calloc(categorias, sizeof(double));
    CONSTANT->distribution_in = (double**) malloc(categorias*sizeof(double*));
    CONSTANT->prob_in = (double**) malloc(categorias*sizeof(double*));
    int i,j;
    int N;
    char filename[200];
    FILE *file,*arquivo,*datei;

    if(is_undirect)arquivo = fopen("./input/multi_probability_in.txt","r");
    else arquivo = fopen("./input/polblogs/multi_probability_in.txt","r");
    if(is_undirect) datei = fopen("./input/polblogs/faixas.txt","r");
    else datei = fopen("./input/polblogs/faixas.txt","r");
    int a;
    for (i = 0; i < categorias; i++){
        if(is_undirect) sprintf(filename,"./input/distribution_%d_in.txt",i);
        else sprintf(filename,"./input/polblogs/distribution_%d_in.txt",i);
        file = fopen(filename,"r");

        N = size_txt(filename);
        CONSTANT->distribution_in[i] = (double*) calloc(N,sizeof(double));
        for (j = 0; j < N; j++) if(fscanf(file, "%d %lf\n",&a, &CONSTANT->distribution_in[i][j]));

        CONSTANT->prob_in[i] = (double*) malloc(sizeof(double)*categorias);

        for(j = 0; j < categorias; j++) if(fscanf(arquivo, "%lf ", &CONSTANT->prob_in[i][j]));
        if(fscanf(arquivo, "\n"));
        if(fscanf(datei, "%lf\n", &CONSTANT->p_faixas[i]));
        
        fclose(file);

    }
    for ( j = 0; j < categorias; j++) for (i = 1; i < categorias; i++) CONSTANT->prob_in[j][i] += CONSTANT->prob_in[j][i-1];

    for (i = 1; i < categorias; i++) CONSTANT->p_faixas[i] += CONSTANT->p_faixas[i-1];

    fclose(arquivo);
    if(!is_undirect){
        CONSTANT->distribution_out = (double**) malloc(categorias*sizeof(double*));
        CONSTANT->prob_out = (double**) malloc(categorias*sizeof(double*));
        if(is_undirect)arquivo = fopen("./input/multi_probability_out.txt","r");
        else arquivo = fopen("./input/polblogs/multi_probability_out.txt","r");
        for (i = 0; i < categorias; i++){
            if(is_undirect) sprintf(filename,"./input/distribution_%d_out.txt",i);
            else sprintf(filename,"./input/polblogs/distribution_%d_out.txt",i);
            file = fopen(filename,"r");

            N = size_txt(filename);
            CONSTANT->distribution_out[i] = (double*) calloc(N,sizeof(double));
            for (j = 0; j < N; j++) if(fscanf(file, "%d %lf\n",&a, &CONSTANT->distribution_out[i][j]));

            CONSTANT->prob_out[i] = (double*) malloc(sizeof(double)*categorias);

            for(j = 0; j < categorias; j++) if(fscanf(arquivo, "%lf ", &CONSTANT->prob_out[i][j]));
            if(fscanf(arquivo, "\n"));
            fclose(file);
        }
        for ( j = 0; j < categorias; j++) for (i = 1; i < categorias; i++) CONSTANT->prob_out[j][i] += CONSTANT->prob_out[j][i-1];
        fclose(arquivo);
    }

}

int somatorio(int** array,int elemento,int N){
    int soma = 0;
    for (int i = 0; i < N; i++) soma += array[elemento][i];
    return soma;
}

void append_vizinhos(struct Graph *G,int site,int vizinho,int plus){
    G->viz[site][0]++;
    G->viz[vizinho][0]++;
    G->viz[site] = (int*) realloc(G->viz[site],(G->viz[site][0]+1)*sizeof(int));
    G->viz[vizinho] = (int*) realloc(G->viz[vizinho],(G->viz[vizinho][0]+1)*sizeof(int));
    G->viz[site][G->viz[site][0]] = plus*vizinho;
    G->viz[vizinho][G->viz[vizinho][0]] = plus*site;
}


double generate_conections(struct Graph *G,int** degree, igraph_vector_t* faixas,double p,igraph_t* Grafo){
    int ligacoes_total = 0;
    int i,j;
    //int** matriz = (int**) malloc(5*sizeof(int*));
    //for (i = 0; i < 5; i++) matriz[i] = (int*) calloc(5,sizeof(int));
    for (i = 0; i < G->Nodes; i++){
        //print_vetor(degree[i],5,sizeof(int));
        //for (j = 0; j < 5; j++) matriz[(int) VECTOR(*faixas)[i]][j] += degree[i][j];
        ligacoes_total += somatorio(degree,i,5);
    }
    //for (i = 0; i < 5; i++) free(matriz[i]);
    //free(matriz);

    return ligacoes_total/G->edges;
}

struct Graph weighted_add_edge(struct Graph *G,double p,struct mce *MCE,int site,int faixa,igraph_vector_int_t* edges,igraph_vector_t* faixas,bool** liga){
    int i,vizinho,faixa1;
    faixa1 = (int) VECTOR(*faixas)[site];
    for ( i = 0; i < MCE->n_faixas[faixa]; i++){

        if(MCE->degree_in[site][faixa] == 0) break;

        vizinho = MCE->site_per_faixas[faixa][i];

        if(site == vizinho) continue;

        if(genrand64_real1() <= p){
            
            
            if(!MCE->is_undirect) if(MCE->degree_out[vizinho][faixa1] == 0) continue;
            else if(MCE->degree_in[vizinho][faixa1] == 0) continue;

            if(liga[site][vizinho]) continue;

            liga[site][vizinho] = true;
            if(MCE->is_undirect)liga[vizinho][site] = true;
            
            MCE->degree_in[site][faixa]--;

            if(!MCE->is_undirect) MCE->degree_out[vizinho][faixa1]--;
            else MCE->degree_in[vizinho][faixa1]--;
            
            igraph_vector_int_push_back(edges, site);
            igraph_vector_int_push_back(edges, vizinho);

            G->viz[site][0]++;
            if(MCE->is_undirect) G->viz[vizinho][0]++;
            G->viz[site][G->viz[site][0]] = vizinho;
            if(MCE->is_undirect) G->viz[vizinho][G->viz[vizinho][0]] = site;
            G->edges += 1;
        }
        else{
            
            if(liga[site][vizinho]) continue;
            liga[site][vizinho] = true;
            if(MCE->is_undirect) liga[vizinho][site] = true;
        }
    }
    return *G;
}

igraph_t local_configuration_model(int N,int categorias,bool is_undirect, double p,int seed,double* perca){

    struct Graph G;
    struct constant CONSTANT;
    struct mce MCE;
    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    int i = 0,j = 0,k = 0;
    
    G.edges = 0;

    igraph_vector_t faixas;
    igraph_vector_init(&faixas, G.Nodes);

    init_CONSTANT(&CONSTANT,categorias,is_undirect);
    init_MCE(&CONSTANT,&MCE,&G,N,&faixas,is_undirect,categorias);
    MCE.is_undirect = is_undirect;
    free_CONSTANT(&CONSTANT,categorias,is_undirect);
    //for ( i = 0; i < N; i++) printf("In: %d %d, Out: %d %d\n",MCE.degree_in[i][0],MCE.degree_in[i][1],MCE.degree_out[i][0],MCE.degree_out[i][1]);
    

    bool** liga = (bool**) malloc(N*sizeof(bool*));
    for (i = 0; i < N; i++) liga[i] = (bool*) calloc(N,sizeof(bool));
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    int site;
    uint8_t faixa;
    int vizinho;
    uint8_t faixa1;
    exit(0);
    for (site = 0; site < N; site++){

        for (faixa = 0; faixa < 5; faixa++){
            seed++;
            if(MCE.degree_in[site][faixa] == 0) continue;
            MCE.site_per_faixas[faixa] = randomize(MCE.site_per_faixas[faixa], MCE.n_faixas[faixa],seed);
            G = weighted_add_edge(&G,1,&MCE,site,faixa,&edges,&faixas,liga);
        }
        
        /*if(p > 0){
            int **vizinhos = (int**) malloc(categorias*sizeof(int*));
            int *q_vizinhos = (int*) calloc(categorias,sizeof(int));
            for (i = 0; i < 5; i++) vizinhos[i] = (int*) malloc(0*sizeof(int));
            int n = 0;
            for (i = 0; i < G.viz[site][0]; i++){
                
                vizinho = G.viz[site][i+1];
                if(vizinho < 0) continue;
                if(somatorio(degree,vizinho,categorias) !=0){
                    faixa1 = (int)VECTOR(faixas)[vizinho];
                    q_vizinhos[faixa1]++;
                    vizinhos[faixa1] = realloc(vizinhos[faixa1],q_vizinhos[faixa1]*sizeof(int));
                    vizinhos[faixa1][q_vizinhos[faixa1] - 1] = vizinho;
                    n++;
                }
            }
            if(n > 1){
                
                for (int viz = 0; viz < G.viz[site][0]; viz++){
                    vizinho = G.viz[site][viz+1];
                    
                    if(vizinho < 0) continue;
                    for (faixa = 0; faixa < 5; faixa++){
                        seed++;
                        if(degree[vizinho][faixa] == 0) continue;
                        vizinhos[faixa] = randomize(vizinhos[faixa], q_vizinhos[faixa],seed);
                        G = weighted_add_edge(&G,p,degree,vizinhos,q_vizinhos,vizinho,faixa,&edges,&faixas,liga);
                    }
                }
            }
            for(i = 0; i < 5; i++) free(vizinhos[i]);
            free(vizinhos);
            free(q_vizinhos);
        }*/

    }
    igraph_t Grafo;
    igraph_set_attribute_table(&igraph_cattribute_table);
    if(MCE.is_undirect) igraph_empty(&Grafo, G.Nodes, IGRAPH_UNDIRECTED);
    else igraph_empty(&Grafo, G.Nodes, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    igraph_cattribute_VAN_setv(&Grafo,"faixa",&faixas);
    //*perca = generate_conections(&G,degree,&faixas,p,&Grafo);
    for (i = 0; i < G.Nodes; i++){
        free(G.viz[i]);
        free(liga[i]);
    }
    free(liga);
    free(G.viz);
    for(i = 0;i < 5;i++) free(MCE.site_per_faixas[i]);
    free(MCE.site_per_faixas);
    free(MCE.n_faixas);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&faixas);

    
    return Grafo;
}

double calcular_correlacao(igraph_vector_t* vetor1, igraph_vector_t* vetor2) {
    double soma1 = 0, soma2 = 0, soma1Q = 0, soma2Q = 0, prodSoma = 0;
    int n = igraph_vector_size(vetor1);

    for (int i = 0; i < n; i++) {
        soma1 += VECTOR(*vetor1)[i];
        soma2 += VECTOR(*vetor2)[i];
        soma1Q += pow(VECTOR(*vetor1)[i], 2);
        soma2Q += pow(VECTOR(*vetor2)[i], 2);
        prodSoma += VECTOR(*vetor1)[i] * VECTOR(*vetor2)[i];
    }

    double numerador = prodSoma - (soma1 * soma2 / n);
    double denominador = sqrt((soma1Q - pow(soma1, 2) / n) * (soma2Q - pow(soma2, 2) / n));

    if (denominador == 0) {
        return 0;  // Evita divisão por zero
    } else {
        return numerador / denominador;
    }
}


void calcula_propriedades(igraph_t *Grafo,double p,int N, double *resultados) {

    double caminho_medio;
    double diametro;
    double agrupamento_medio;
    double agrupamento_total;
    double correlation;

    igraph_vector_int_t degrees;
    igraph_vector_int_init(&degrees, N);
    igraph_degree(Grafo, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

    igraph_vector_t clustering;
    igraph_vector_init(&clustering, N);
    igraph_transitivity_local_undirected(Grafo,&clustering,igraph_vss_all(),IGRAPH_TRANSITIVITY_ZERO);

    igraph_vector_t graus;
    igraph_vector_init(&graus, N);
    for (int i = 0; i < N; i++) VECTOR(graus)[i] = (double) VECTOR(degrees)[i];
    
    correlation = calcular_correlacao(&graus,&clustering);

    double k_medio = (double) igraph_vector_sum(&graus)/N;

    igraph_vector_add_constant(&graus,-k_medio);
    igraph_vector_mul(&graus,&graus);
    igraph_vector_scale(&graus,1/N);

    double std = 0 ;
    for (int i = 0; i < N; i++) std += pow((double)VECTOR(degrees)[i],2);
    std = std/N;
    // Mediana dos graus
    igraph_vector_sort(&graus);
    double mediana = (N%2 == 0)? VECTOR(degrees)[(int)(N-1)/2] : (VECTOR(degrees)[(int) N/2] + VECTOR(degrees)[(int) N/2+1])*0.5;

    // Calcula o menor caminho médio da Rede
    igraph_average_path_length(Grafo, &caminho_medio, NULL, IGRAPH_UNDIRECTED, 1);

    // Calcula o diâmetro da rede
    igraph_diameter(Grafo, &diametro, 0, 0, 0, 0, IGRAPH_UNDIRECTED, 1);

    // Calcula o agrupamento médio
    igraph_transitivity_avglocal_undirected(Grafo,&agrupamento_medio,IGRAPH_TRANSITIVITY_ZERO);

    // Calcula o agrupamento total
    igraph_transitivity_undirected(Grafo,&agrupamento_total,IGRAPH_TRANSITIVITY_NAN);
    resultados[0] = k_medio;
    resultados[1] = mediana;
    resultados[2] = std;
    resultados[3] = agrupamento_medio;
    resultados[4] = agrupamento_total;
    resultados[5] = correlation;
    resultados[6] = caminho_medio;
    resultados[7] = diametro;
    igraph_vector_destroy(&graus);
    igraph_vector_destroy(&clustering);
    igraph_vector_int_destroy(&degrees);

}


void generate_local_configuration_model(int N,int categorias,bool is_undirect,double p, int redes,int seed){

    int i;
    double a;
    double** resultados;
    if(redes > 1){
        resultados = (double **)malloc(redes * sizeof(double *));
        for(i = 0; i < redes; i++) resultados[i] = (double *)calloc(10 , sizeof(double ));  // Aponta para o início de cada linha no vetor resultados

    }
    clock_t inicio, fim;
    double perca;
    int count = 0;
    //omp_set_num_threads(7);
    //#pragma omp parallel for
    for (i = 0; i < redes; i++){

        //if(redes!=1) printf("\e[1;1H\e[2J");

        igraph_t G;
        inicio = clock();
        G = local_configuration_model(N,categorias,is_undirect,p,seed+i,&perca);
        if(redes!=1){
            calcula_propriedades(&G,p,N,resultados[i]);
            fim = clock();
            resultados[i][8] = ((double) (fim - inicio)) / CLOCKS_PER_SEC;
            resultados[i][9] = perca;
        }
        igraph_destroy(&G);
        count += 1;
    }


    if(redes > 1){
        FILE *file;
        int j;
        char filecheck[800];
        sprintf(filecheck,"./output/modelo/resultados_%d_%.2f.txt",N,p);
        file = fopen(filecheck,"w");
        char linha[10000];
        linha[0] = '\0';
        for ( i = 0; i < redes; i++){
            
            for(j = 0; j < 10 ;j++){
                char buffer[200]; // Buffer para armazenar a representação do número como string
                snprintf(buffer, sizeof(buffer), "%.2f ", resultados[i][j]); // Converte o número para string com duas casas decimais
                strcat(linha, buffer); // Adiciona o número à string final
            }
            fprintf(file,"%s\n",linha);
            linha[0] = '\0';
        }
        
        
        fclose(file);
        for(i = 0; i < redes; i++) free(resultados[i]);
        free(resultados);
    }
    //igraph_vector_int_destroy(&centralidade);
    printf("Terminou %d %.2f\n",N,p);
}
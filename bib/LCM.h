#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#include <unistd.h>

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
    int limit;
};

struct mce{
    bool is_undirect;
    int** degree_in;
    int** degree_out;
    int categorias;
};

void free_CONSTANT(struct constant *CONSTANT,int categorias,bool is_undirect){
    int i;
    for (i = 0; i < categorias; i++){
        if(!is_undirect)free(CONSTANT->prob_in[i]);
        free(CONSTANT->prob_out[i]);
        if(!is_undirect)free(CONSTANT->distribution_in[i]);
        free(CONSTANT->distribution_out[i]);
    }
    if(!is_undirect) free(CONSTANT->prob_in);
    free(CONSTANT->prob_out);
    if(!is_undirect) free(CONSTANT->distribution_in);
    free(CONSTANT->distribution_out);
    free(CONSTANT->p_faixas);
}

void init_MCE(struct constant *CONSTANT,struct mce *MCE,struct Graph *G,igraph_vector_t* faixas,bool is_undirect,int categorias){
    int i,j,k;
    double r;
    MCE->categorias = categorias;

    int* arr = (int*) calloc(G->Nodes,sizeof(int));

    int* grau_out = (int*) calloc(G->Nodes,sizeof(int));

    int* grau_in;

    int* faixas2 = (int*) calloc(G->Nodes,sizeof(int));

    MCE->degree_out = (int**) malloc(G->Nodes*sizeof(int*));
    if(!is_undirect){
        grau_in = (int*) calloc(G->Nodes,sizeof(int));
        MCE->degree_in = (int**) malloc(G->Nodes*sizeof(int*));
    }
    for (i = 0; i < G->Nodes; i++){
        j = 0;
        r = genrand64_real1();
        while(r > CONSTANT->p_faixas[j]) j++;
        faixas2[i] = j;
        
        k = 0;
        MCE->degree_out[i] = (int*) calloc(MCE->categorias,sizeof(int));
        if(is_undirect){
            while(k == 0) k = empiric_distribution(CONSTANT->distribution_out[faixas2[i]],CONSTANT->limit);
            grau_out[i] = k;
        }
        
        else{
            int k1 = 0;
            int k2 = 0;
            MCE->degree_in[i] = (int*) calloc(MCE->categorias,sizeof(int));
            while(k1 == 0){
                k = empiric_distribution(CONSTANT->distribution_out[faixas2[i]],CONSTANT->limit);
                k2 = empiric_distribution(CONSTANT->distribution_in[faixas2[i]],CONSTANT->limit);
                k1 = k+k2;
            }
            grau_out[i] = k;
            grau_in[i] = k2;
        }
        //printf("\n");
        arr[i] = i;
    }

    int* temp_weights = (int*) calloc(G->Nodes,sizeof(int));
    for (i = 0; i < G->Nodes; i++)temp_weights[i] = grau_out[i];
    mergeSort(temp_weights,arr,0,G->Nodes-1);
    free(temp_weights);

    int site;
    for (i = 0; i < G->Nodes; i++){
        site = arr[i];

        generate_multinomial(grau_out[site], MCE->categorias, CONSTANT->prob_out[faixas2[site]], MCE->degree_out[i]);

        G->viz[i] =(int*) malloc(1*sizeof(int));
        G->viz[i][0] = 0;
        if(!is_undirect)generate_multinomial(grau_in[site], MCE->categorias, CONSTANT->prob_in[faixas2[site]], MCE->degree_in[i]);
        VECTOR(*faixas)[i] = faixas2[site];
    }
    free(faixas2);
    free(grau_out);
    if(!is_undirect) free(grau_in);
    free(arr);
    
}

int* remove_element(int* vetor,int size,int i){
    int temp = vetor[i];
    vetor[i] = vetor[size-1];
    vetor[size-1] = temp;
    vetor = realloc(vetor,(size-1)*sizeof(int));
    return vetor;
}

void init_CONSTANT(struct constant *CONSTANT,bool is_undirect,int* category){
    int categorias;
    if(is_undirect) categorias = size_txt("./input/POLYMOD/faixas.txt");
    else categorias = size_txt("./input/polblogs/faixas.txt");
    CONSTANT->p_faixas = (double *) calloc(categorias, sizeof(double));
    CONSTANT->distribution_out = (double**) malloc(categorias*sizeof(double*));
    CONSTANT->prob_out = (double**) malloc(categorias*sizeof(double*));
    int i,j;
    int N;
    char filename[200];
    FILE *file,*arquivo,*datei;

    if(is_undirect)arquivo = fopen("./input/POLYMOD/multi_probability_out.txt","r");
    else arquivo = fopen("./input/polblogs/multi_probability_out.txt","r");

    if(is_undirect) datei = fopen("./input/POLYMOD/faixas.txt","r");
    else datei = fopen("./input/polblogs/faixas.txt","r");
    int a;
    for (i = 0; i < categorias; i++){

        if(is_undirect) sprintf(filename,"./input/POLYMOD/distribution_%d_out.txt",i);
        else sprintf(filename,"./input/polblogs/distribution_%d_out.txt",i);

        file = fopen(filename,"r");

        N = size_txt(filename);
        CONSTANT->limit = N;
        CONSTANT->distribution_out[i] = (double*) calloc(N,sizeof(double));
        for (j = 0; j < N; j++) if(fscanf(file, "%d %lf\n",&a, &CONSTANT->distribution_out[i][j]));
        
        CONSTANT->prob_out[i] = (double*) malloc(sizeof(double)*categorias);

        for(j = 0; j < categorias; j++) if(fscanf(arquivo, "%lf ", &CONSTANT->prob_out[i][j]));
        if(fscanf(arquivo, "\n"));
        if(fscanf(datei, "%lf\n", &CONSTANT->p_faixas[i]));
        
        fclose(file);

    }

    fclose(arquivo);
    fclose(datei);

    if(!is_undirect){
        arquivo = fopen("./input/polblogs/multi_probability_in.txt","r");
        CONSTANT->distribution_in = (double**) malloc(categorias*sizeof(double*));
        CONSTANT->prob_in = (double**) malloc(categorias*sizeof(double*));

        for (i = 0; i < categorias; i++){
            sprintf(filename,"./input/polblogs/distribution_%d_in.txt",i);
            file = fopen(filename,"r");

            N = size_txt(filename);

            CONSTANT->distribution_in[i] = (double*) calloc(N,sizeof(double));
            for (j = 0; j < N; j++) if(fscanf(file, "%d %lf\n",&a, &CONSTANT->distribution_in[i][j]));

            CONSTANT->prob_in[i] = (double*) malloc(sizeof(double)*categorias);

            for(j = 0; j < categorias; j++) if(fscanf(arquivo, "%lf ", &CONSTANT->prob_in[i][j]));

            if(fscanf(arquivo, "\n"));
            fclose(file);
        }
        for ( j = 0; j < categorias; j++) for (i = 1; i < categorias; i++) CONSTANT->prob_in[j][i] += CONSTANT->prob_in[j][i-1];
        fclose(arquivo);
    }
    for ( j = 0; j < categorias; j++) for (i = 1; i < categorias; i++) CONSTANT->prob_out[j][i] += CONSTANT->prob_out[j][i-1];
    for (i = 1; i < categorias; i++) CONSTANT->p_faixas[i] += CONSTANT->p_faixas[i-1];
    *category = categorias;
}

void print_lig(struct mce *MCE, igraph_vector_t* faixas,int N){
    int i,j;
    int** matriz_out = (int**) malloc(MCE->categorias*sizeof(int*));
    int** matriz_in = (int**) malloc(MCE->categorias*sizeof(int*));

    for (i = 0; i < MCE->categorias; i++) matriz_out[i] = (int*) calloc(MCE->categorias,sizeof(int));

    if(!MCE->is_undirect) for (i = 0; i < MCE->categorias; i++) matriz_in[i] = (int*) calloc(MCE->categorias,sizeof(int));
    for (i = 0; i < N; i++){
        //for (j = 0; j < 5; j++) matriz[(int) VECTOR(*faixas)[i]][j] += degree[i][j];
        for (j = 0; j < MCE->categorias; j++) matriz_out[(int) VECTOR(*faixas)[i]][j] += MCE->degree_out[i][j];
        if(!MCE->is_undirect){
            for (j = 0; j < MCE->categorias; j++) matriz_in[(int) VECTOR(*faixas)[i]][j] += MCE->degree_in[i][j];
        }
    }
    printf("\n");
    // Create a for to print the matrix matriz_out
    for (i = 0; i < MCE->categorias; i++){
        for (j = 0; j < MCE->categorias; j++) printf("%d ", matriz_out[i][j]);
        printf("\n");
    }
    printf("\n");
    if(!MCE->is_undirect){
        for (i = 0; i < MCE->categorias; i++){
            for (j = 0; j < MCE->categorias; j++) printf("%d ", matriz_in[i][j]);
            printf("\n");
        }
        printf("\n");
    }
    
    for (i = 0; i < MCE->categorias; i++) free(matriz_out[i]);
    free(matriz_out);

}

double generate_conections(struct Graph *G,struct mce *MCE, igraph_vector_t* faixas,double p,igraph_t* Grafo,double total_edges){
    double ligacoes_total = 0;
    int i,j;
    for (i = 0; i < G->Nodes; i++) ligacoes_total += somatorio(MCE->degree_out,i,MCE->categorias);
    if(MCE->is_undirect) (double)ligacoes_total/(2*total_edges);
    return (double)ligacoes_total/(total_edges);
}

bool check_ligacao(struct Graph *G,int site,int vizinho){
    int i;
    for( i = 0;i < G->viz[site][0]; i++) if(vizinho == G->viz[site][i+1]) return true;
    return false; 
}

struct Graph weighted_add_edge(struct Graph *G,double p,struct mce *MCE,int *n,int* site_per_faixas,int site,int faixa_vizinho,igraph_vector_t* faixas,igraph_vector_int_t* edges){

    int i,vizinho,faixa_site;
    faixa_site = (int) VECTOR(*faixas)[site];
    for ( i = 0; i < n[faixa_site]; i++){
        if(MCE->degree_out[site][faixa_vizinho] == 0) break;

        vizinho = site_per_faixas[i];
        if(vizinho < 0) continue;
        if(site == vizinho) continue;
        //printf("Vizinho = %d\n", vizinho);
        if(genrand64_real1() <= p){

            if(MCE->is_undirect) if(MCE->degree_out[vizinho][faixa_site] == 0) continue; 
            if(!MCE->is_undirect) if(MCE->degree_in[vizinho][faixa_site] == 0) continue;
            if(check_ligacao(G,site,vizinho)) continue;
            
            MCE->degree_out[site][faixa_vizinho]--;

            if(MCE->is_undirect){
                
                MCE->degree_out[vizinho][faixa_site]--;
                if(MCE->degree_out[vizinho][faixa_site] == 0){
                    site_per_faixas = remove_element(site_per_faixas,n[faixa_site],i);
                    i--;
                    n[faixa_site]--;
                }
            }
            if(!MCE->is_undirect){
                MCE->degree_in[vizinho][faixa_site]--;
                if(MCE->degree_in[vizinho][faixa_site] == 0){
                    site_per_faixas = remove_element(site_per_faixas,n[faixa_site],i);
                    i--;
                    n[faixa_site]--;
                }
            }
            
            igraph_vector_int_push_back(edges, site);
            igraph_vector_int_push_back(edges, vizinho);

            G->viz[site][0]++;
            G->viz[site] = (int*) realloc(G->viz[site],(G->viz[site][0]+1)*sizeof(int));

            G->viz[site][G->viz[site][0]] = vizinho;

            if(MCE->is_undirect){
                G->viz[vizinho][0]++;
                G->viz[vizinho] = realloc(G->viz[vizinho],(G->viz[vizinho][0]+1)*sizeof(int));
                G->viz[vizinho][G->viz[vizinho][0]] = site;
            }

            G->edges += 1;
        }
    }
    
    return *G;
}

void load_test(struct mce *MCE,igraph_vector_t* faixas,int N,struct Graph *G){

    FILE* arquivo = fopen("./input/test/in.txt", "r");
    FILE* file = fopen("./input/test/out.txt", "r");
    FILE* datei = fopen("./input/test/faixas.txt", "r");
    int a,b,i,idx;

    int n = size_txt("./input/test/faixas.txt");
    MCE->degree_out = (int**) malloc(N*n*sizeof(int*));
    MCE->degree_in = (int**) malloc(N*n*sizeof(int*));
    for (i = 0; i < N*n; i++) {
        MCE->degree_in[i] = (int*) calloc(2, sizeof(int));
        MCE->degree_out[i] = (int*) calloc(2,sizeof(int));
        if (i < n) {
            if (fscanf(arquivo, "%d %d\n", &MCE->degree_in[i][0], &MCE->degree_in[i][1]) != 2) {
                printf("Formato de arquivo in inválido\n");
                exit(1);
            }
            if (fscanf(file, "%d %d\n", &MCE->degree_out[i][0], &MCE->degree_out[i][1]) != 2){
                printf("Formato de arquivo out inválido\n");
                exit(1);
            }
            if(fscanf(datei, "%d\n",&a));
            VECTOR(*faixas)[i] = a;
        } 
        else {
            idx = i % n;
            MCE->degree_in[i][0] = MCE->degree_in[idx][0];
            MCE->degree_in[i][1] = MCE->degree_in[idx][1];
            MCE->degree_out[i][0] = MCE->degree_out[idx][0];
            MCE->degree_out[i][1] = MCE->degree_out[idx][1];
            VECTOR(*faixas)[i] = VECTOR(*faixas)[idx];
        }
        G->viz[i] = (int*) malloc(sizeof(int));
        G->viz[i][0] = 0;
    }
    fclose(arquivo);
    fclose(file);
    fclose(datei);


}

igraph_t local_configuration_model(int N,bool is_undirect, double p,int seed,double* perca){

    init_genrand64(seed);

    struct Graph G;
    double total_edges = 0;
    struct mce MCE;
    igraph_vector_t faixas;
    int i = 0,j = 0,k = 0,categorias = 0,f;
    G.edges = 0;
    FILE* in;
    FILE* out;
    MCE.is_undirect = is_undirect;
    char filename[200];
    if(!is_undirect) sprintf(filename,"./output/blogs/matrix/in_matrix_%d_%.2f.txt",N,p);
    if(!is_undirect) in = fopen(filename,"a");
    
    if(is_undirect) sprintf(filename,"./output/POLYMOD/matrix/out_matrix_%d_%.2f.txt",N,p);
    else sprintf(filename,"./output/blogs/matrix/out_matrix_%d_%.2f.txt",N,p);
    out = fopen(filename,"a");
    if(is_undirect){
        G.Nodes = N;
        G.viz = (int **)malloc(N*sizeof(int*));
        igraph_vector_init(&faixas, G.Nodes);
        struct constant CONSTANT;
        init_CONSTANT(&CONSTANT,is_undirect,&categorias);
        MCE.categorias = categorias;
        init_MCE(&CONSTANT,&MCE,&G,&faixas,is_undirect,categorias);
        free_CONSTANT(&CONSTANT,categorias,is_undirect);

    }
    else{
        categorias = 2;
        MCE.categorias = categorias;
        G.viz = (int **)malloc(N*1224*sizeof(int*));
        igraph_vector_init(&faixas, 1224*N);
        load_test(&MCE,&faixas,N,&G);
        
        
        N *= 793;
        G.Nodes = N;
    }
    
    int** n_faixas = calloc(categorias,sizeof(int*));
    int*** site_per_faixas = malloc(categorias*sizeof(int**));
    for ( i = 0; i < MCE.categorias; i++){
        n_faixas[i] = calloc(categorias,sizeof(int));
        site_per_faixas[i] = malloc(categorias*sizeof(int*));
        for ( j = 0; j < MCE.categorias; j++)site_per_faixas[i][j] = malloc(0*sizeof(int));
    }
    for(i = 0; i < G.Nodes; i++){
        f = VECTOR(faixas)[i];
        total_edges += somatorio(MCE.degree_out,i,MCE.categorias);
        for(j = 0;j < MCE.categorias;j++){
            if(is_undirect){
                if(MCE.degree_out[i][j] != 0){
                    n_faixas[f][j]++;
                    site_per_faixas[f][j] = realloc(site_per_faixas[f][j],n_faixas[f][j]*sizeof(int));
                    site_per_faixas[f][j][n_faixas[f][j] - 1] = i;
                }

            }
            if(!is_undirect){
                if(MCE.degree_in[i][j] != 0){
                    n_faixas[f][j]++;
                    site_per_faixas[f][j] = realloc(site_per_faixas[f][j],n_faixas[f][j]*sizeof(int));
                    site_per_faixas[f][j][n_faixas[f][j] - 1] = i;
                }
            }
        }
        
    }
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    int site;
    uint8_t faixa;
    int vizinho;
    uint8_t faixa1;
    for (site = 0; site < N; site++){
        f = VECTOR(faixas)[site];
        for (faixa = 0; faixa < categorias; faixa++){
            seed++;

            if(MCE.degree_out[site][faixa] == 0) continue;

            // Se têm nós na faixa etária "faixa" e podem receber ligações da faixa etária f
            if(n_faixas[faixa][f] == 0) continue;

            site_per_faixas[faixa][f] = randomize(site_per_faixas[faixa][f], n_faixas[faixa][f],seed);

            G = weighted_add_edge(&G,1,&MCE,n_faixas[faixa],site_per_faixas[faixa][f],site,faixa,&faixas,&edges);
            if(is_undirect){
                if(MCE.degree_out[site][faixa] == 0){
                    for (i = 0; i < n_faixas[f][faixa]; i++) if( site_per_faixas[f][faixa][i] == site) break;
                    site_per_faixas[f][faixa] = remove_element(site_per_faixas[f][faixa],n_faixas[f][faixa],i);
                    n_faixas[f][faixa]--;
                }

            }
        }
        if(p > 0){
            int** n_faixas_vizinhos = calloc(categorias,sizeof(int*));
            int*** site_per_faixas_vizinhos = malloc(categorias*sizeof(int**));

            for ( i = 0; i < MCE.categorias; i++){
                n_faixas_vizinhos[i] = calloc(categorias,sizeof(int));
                site_per_faixas_vizinhos[i] = malloc(categorias*sizeof(int*));
                for ( j = 0; j < MCE.categorias; j++)site_per_faixas_vizinhos[i][j] = malloc(0*sizeof(int));
                
            }

            for(i = 0; i <  G.viz[site][0]; i++){
                vizinho = G.viz[site][i+1];

                f = VECTOR(faixas)[vizinho];
                for(j = 0;j < MCE.categorias;j++){
                    if(is_undirect){
                        if(MCE.degree_out[vizinho][j] != 0){
                            n_faixas_vizinhos[f][j]++;
                            site_per_faixas_vizinhos[f][j] = realloc(site_per_faixas_vizinhos[f][j],n_faixas_vizinhos[f][j]*sizeof(int));
                            site_per_faixas_vizinhos[f][j][n_faixas_vizinhos[f][j] - 1] = vizinho;
                        }

                    }
                    if(!is_undirect){
                        if(MCE.degree_in[vizinho][j] != 0){
                            n_faixas_vizinhos[f][j]++;
                            site_per_faixas_vizinhos[f][j] = realloc(site_per_faixas_vizinhos[f][j],n_faixas_vizinhos[f][j]*sizeof(int));
                            site_per_faixas_vizinhos[f][j][n_faixas_vizinhos[f][j] - 1] = vizinho;
                        }
                    }
                }
                
            }
            for (i = 0; i < G.viz[site][0]; i++){
                vizinho = G.viz[site][i+1];
                f = VECTOR(faixas)[vizinho];
                for (faixa = 0; faixa < categorias; faixa++){
                    seed++;

                    if(MCE.degree_out[vizinho][faixa] == 0) continue;
                    if(n_faixas_vizinhos[faixa][f] == 0) continue;
                    site_per_faixas_vizinhos[faixa][f] = randomize(site_per_faixas_vizinhos[faixa][f], n_faixas_vizinhos[faixa][f],seed);
                    G = weighted_add_edge(&G,p,&MCE,n_faixas_vizinhos[faixa],site_per_faixas_vizinhos[faixa][f],vizinho,faixa,&faixas,&edges);
                    if(is_undirect){
                        if(MCE.degree_out[vizinho][faixa] == 0){
                            for (k = 0; k < n_faixas_vizinhos[f][faixa]; k++) if( site_per_faixas_vizinhos[f][faixa][k] == vizinho) break;
                            site_per_faixas_vizinhos[f][faixa] = remove_element(site_per_faixas_vizinhos[f][faixa],n_faixas_vizinhos[f][faixa],k);
                            n_faixas_vizinhos[f][faixa]--;
                        }

                    }
                    
                }
            }
            for (i = 0; i < categorias; i++) {
                for (j = 0; j < categorias; j++)if(n_faixas_vizinhos[i][j]!=0) free(site_per_faixas_vizinhos[i][j]);
                free(n_faixas_vizinhos[i]);
                free(site_per_faixas_vizinhos[i]);
            }
            free(n_faixas_vizinhos);
            free(site_per_faixas_vizinhos);
        } 

    }
    
    igraph_t Grafo;
    igraph_set_attribute_table(&igraph_cattribute_table);
    if(MCE.is_undirect) igraph_empty(&Grafo, G.Nodes, IGRAPH_UNDIRECTED);
    else igraph_empty(&Grafo, G.Nodes, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    
    *perca = generate_conections(&G,&MCE,&faixas,p,&Grafo,total_edges);
    for (i = 0; i < G.Nodes; i++){
        //get the sucessors of the node i
        if(!MCE.is_undirect){
            igraph_vector_int_t sucessors;
            igraph_vector_int_init(&sucessors, 0);
            igraph_neighbors(&Grafo, &sucessors, i, IGRAPH_OUT);

            int* connections = (int*) calloc(MCE.categorias,sizeof(int));
            for (j = 0; j < igraph_vector_int_size(&sucessors); j++) connections[(int) VECTOR(faixas)[VECTOR(sucessors)[j]]]++;
            fprintf(out,"%d %d %d\n",(int) VECTOR(faixas)[i],connections[0],connections[1]);
            free(connections);
            igraph_vector_int_destroy(&sucessors);

            // get the predecessors of the node i
            igraph_vector_int_t predecessors;
            igraph_vector_int_init(&predecessors, 0);
            igraph_neighbors(&Grafo, &predecessors, i, IGRAPH_IN);

            int* connections_in = (int*) calloc(MCE.categorias,sizeof(int));
            for (j = 0; j < igraph_vector_int_size(&predecessors); j++) connections_in[(int) VECTOR(faixas)[VECTOR(predecessors)[j]]]++;
            fprintf(in,"%d %d %d\n",(int) VECTOR(faixas)[i],connections_in[0],connections_in[1]);
            
            free(connections_in);
            igraph_vector_int_destroy(&predecessors);

            if(!MCE.is_undirect)free(MCE.degree_in[i]);
            free(MCE.degree_out[i]);
            free(G.viz[i]);

        }
        else{
            
            int* connections = (int*) calloc(MCE.categorias,sizeof(int));
            for (j = 0; j < G.viz[i][0]; j++) connections[(int) VECTOR(faixas)[G.viz[i][j+1]]]++;
            fprintf(out,"%d %d %d %d %d %d\n",(int) VECTOR(faixas)[i],connections[0],connections[1],connections[2],connections[3],connections[4]);
            
            /*for (j = 0; j < MCE.categorias; j++){
                if(j != MCE.categorias-1)fprintf(out,"%d ",connections[j]);
                else fprintf(out,"%d\n",connections[j]);
            }*/
            
            free(connections);
        }
    }
    if(!MCE.is_undirect)fclose(in);
    fclose(out);
    free(G.viz);
    if(!MCE.is_undirect) free(MCE.degree_in);
    free(MCE.degree_out);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&faixas);

    for (i = 0; i < categorias; i++) {
        for (j = 0; j < categorias; j++)if(n_faixas[i][j]!=0) free(site_per_faixas[i][j]);
        free(n_faixas[i]);
        free(site_per_faixas[i]);
    }
    free(n_faixas);
    free(site_per_faixas);



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

    double caminho_medio = 0;
    double diametro = 0;
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
    int n = 0;
    for (int i = 0; i < N; i++) std += pow((double)VECTOR(degrees)[i],2);
    std = std/N;
    // Mediana dos graus
    igraph_vector_sort(&graus);
    double mediana = (N%2 == 0)? VECTOR(degrees)[(int)(N-1)/2] : (VECTOR(degrees)[(int) N/2] + VECTOR(degrees)[(int) N/2+1])*0.5;


    // calculate the mean shortest path
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


void generate_local_configuration_model(int N,bool is_undirect,double p, int redes,int seed){

    char filename[200];
    if(is_undirect) {
        sprintf(filename, "./output/POLYMOD/matrix/out_matrix_%d_%.2f.txt", N, p);
        if (access(filename, F_OK) == 0) {
            remove(filename);
            printf("Removendo arquivo %s\n", filename);
        }

    }
    else{
        sprintf(filename, "./output/blogs/matrix/in_matrix_%d_%.2f.txt", N, p);
        if (access(filename, F_OK) == 0) {
            remove(filename);
            printf("Removendo arquivo %s\n", filename);
        }
        sprintf(filename, "./output/blogs/matrix/out_matrix_%d_%.2f.txt", N, p);
        if (access(filename, F_OK) == 0) {
            remove(filename);
            printf("Removendo arquivo %s\n", filename);
        }
    }

    int i;
    double a;
    double** resultados;
    if(redes > 1){
        resultados = (double **)malloc(redes * sizeof(double *));
        for(i = 0; i < redes; i++) resultados[i] = (double *)calloc(10 , sizeof(double ));  // Aponta para o início de cada linha no vetor resultados

    }
    
    double perca;
    int count = 0;
    //omp_set_num_threads(2);
    //#pragma omp parallel for schedule(dynamic)
    for (i = 0; i < redes; i++){
        double tempo = omp_get_wtime();
        //if(redes!=1) printf("\e[1;1H\e[2J");
        igraph_t G;
        G = local_configuration_model(N,is_undirect,p,seed+i,&perca);
        int Nodes = igraph_vcount(&G);
        if(redes!=1){
            //calcula_propriedades(&G,p,Nodes,resultados[i]);
            resultados[i][8] = omp_get_wtime() - tempo;
            resultados[i][9] = perca;
        }
        else{
            double agrupamento_total;
            igraph_transitivity_undirected(&G,&agrupamento_total,IGRAPH_TRANSITIVITY_NAN);
            printf("Agrupamento: %f\n",agrupamento_total);
            printf("Perda: %f\n",perca);
        }
        igraph_destroy(&G);
        count++;
        if(count%25 == 0) printf("%d/%d\n",count,redes);
    }

    /* if(redes > 1){
        FILE *file;
        int j;
        char filecheck[800];
        if(is_undirect){
            sprintf(filecheck,"./output/POLYMOD/metrics/resultados_POYLMOD_%d_%.2f.txt",N,p);
        }
        else sprintf(filecheck,"./output/blogs/metrics/resultados_blogs_%d_%.2f.txt",N,p);
        file = fopen(filecheck,"w");
        char linha[10000];
        linha[0] = '\0';
        for ( i = 0; i < redes; i++){
            
            for(j = 0; j < 10 ;j++){
                char buffer[200]; // Buffer para armazenar a representação do número como string
                snprintf(buffer, sizeof(buffer), "%.4f ", resultados[i][j]); // Converte o número para string com duas casas decimais
                strcat(linha, buffer); // Adiciona o número à string final
            }
            fprintf(file,"%s\n",linha);
            linha[0] = '\0';
        }
        
        
        fclose(file);
        for(i = 0; i < redes; i++) free(resultados[i]);
        free(resultados);
    } */
    printf("Terminou N: %d Probability:%.2f\n",N,p);
}
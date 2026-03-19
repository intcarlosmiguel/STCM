#pragma once

#include <stdlib.h>
#include <stdio.h>
#include "mtwister.h"
#include <string.h>
#include <math.h>
#include <igraph/igraph.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>

#define M_PI 3.14159265358979323846

static void swap(int *a, int *b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

static int* randomize (int* array, int n,int seed){
    init_genrand64(seed);
    for (int i = n-1; i > 0; i--){
        int j = genrand64_int64() % (i+1);
        swap(&array[i], &array[j]);
    }
    return array;
}
static void print_vetor(void* array,int N,int check){
    if(check == sizeof(int)){
        int* intArray = (int*)array;
        for (int i = 0; i < N; i++){
            if(i!=N-1) printf("%d ",intArray[i]);
            else printf("%d\n",intArray[i]);
        }
    }
    if(check == sizeof(double)){
        double* doubleArray = (double*)array;
        for (int i = 0; i < N; i++){
            if(i!=N-1) printf("%f ",doubleArray[i]);
            else printf("%f\n",doubleArray[i]);
        }
    }
}


static int size_txt(char *str){
    FILE* f;
    int L = 0;
    char c;
    f = fopen(str,"r");
    for (c = getc(f); c != EOF; c = getc(f)) if (c == '\n') L = L + 1;
    return L;
}

static void print_matrix(void** mat,int N,int n,int check){
    if(check == sizeof(int)){
        int** intMat = (int**)mat;
        printf("========================================================\n");
        for (int i = 0; i < N; i++) print_vetor(intMat[i],n,sizeof(intMat[i][0]));
        printf("========================================================\n");
    }
    if(check == sizeof(double)){
        double** doubleMat = (double**)mat;
        printf("========================================================\n");
        for (int i = 0; i < N; i++) print_vetor(doubleMat[i],n,sizeof(doubleMat[i][0]));
        printf("========================================================\n");
    }
    
}
static void generate_multinomial(int n, int k, double *probabilities, int *outcomes) {
    int i,j;

    // Gera os k-1 primeiros valores
    for (i = n; i > 0; i-- ) {
        double r = genrand64_real1();
        for(j = 0; j < k; j++) if(r < probabilities[j]) break;
        outcomes[j]++;
    }

    // O último valor é determinado pelo resto
}

static int empiric_distribution(double* distribution,int limit){
    double r = genrand64_real1();
    int n = 0;
    while(r > distribution[n]){
        n++;
        if(n > limit) return 0;
    }
    
    return n;
}

static void merge(int *temp_weights, int *labels, int l, int m, int r) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    // Cria vetores temporários
    int *L_weights = (int *)malloc(n1 * sizeof(int));
    int *R_weights = (int *)malloc(n2 * sizeof(int));
    int *L_labels = (int *)malloc(n1 * sizeof(int));
    int *R_labels = (int *)malloc(n2 * sizeof(int));

    if (L_weights == NULL || R_weights == NULL || L_labels == NULL || R_labels == NULL) {
        fprintf(stderr, "Erro ao alocar memória\n");
        exit(EXIT_FAILURE);
    }

    // Copia os dados para os vetores temporários L_weights[] e R_weights[], L_labels[] e R_labels[]
    for (i = 0; i < n1; i++) {
        L_weights[i] = temp_weights[l + i];
        L_labels[i] = labels[l + i];
    }
    for (j = 0; j < n2; j++) {
        R_weights[j] = temp_weights[m + 1 + j];
        R_labels[j] = labels[m + 1 + j];
    }

    // Mescla os vetores temporários de volta para temp_weights[l..r] e labels[l..r]
    i = 0; // Índice inicial do primeiro subarray
    j = 0; // Índice inicial do segundo subarray
    k = l; // Índice inicial do subarray mesclado

    while (i < n1 && j < n2) {
        if (L_weights[i] >= R_weights[j]) {  // Ordem decrescente
            temp_weights[k] = L_weights[i];  // Note que isso atualiza apenas o vetor temporário
            labels[k] = L_labels[i];
            i++;
        } else {
            temp_weights[k] = R_weights[j];
            labels[k] = R_labels[j];
            j++;
        }
        k++;
    }

    // Copia os elementos restantes de L_weights[] e L_labels[], se houver
    while (i < n1) {
        temp_weights[k] = L_weights[i];
        labels[k] = L_labels[i];
        i++;
        k++;
    }

    // Copia os elementos restantes de R_weights[] e R_labels[], se houver
    while (j < n2) {
        temp_weights[k] = R_weights[j];
        labels[k] = R_labels[j];
        j++;
        k++;
    }

    // Libera a memória alocada
    free(L_weights);
    free(R_weights);
    free(L_labels);
    free(R_labels);
}

// Função principal que implementa o MergeSort
static void mergeSort(int *temp_weights, int *labels, int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;

        // Ordena a primeira e a segunda metade
        mergeSort(temp_weights, labels, l, m);
        mergeSort(temp_weights, labels, m + 1, r);

        // Combina as duas metades
        merge(temp_weights, labels, l, m, r);
    }
}
static int somatorio(int** array,int elemento,int N){
    int soma = 0;
    for (int i = 0; i < N; i++) soma += array[elemento][i];
    return soma;
}
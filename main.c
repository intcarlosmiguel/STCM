#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <omp.h>

#include "bib/mtwister.h"
#include "bib/calc.h"
#include "bib/LCM.h"

int main(int argc, char *argv[]) {
    // Verificar se todos os argumentos necessários estão presentes
    if (argc < 5) {
        fprintf(stderr, "Erro: Número insuficiente de argumentos.\n");
        fprintf(stderr, "Uso: %s <seed> <N> <is_undirect> <redes> <prob>\n", argv[0]);
        return 1;
    }
    
    int seed = atoi(argv[1]);
    int N = atoi(argv[2]);
    bool is_undirect = atoi(argv[3]);
    int redes = atoi(argv[4]);
    double prob = atof(argv[5]) / 100;
    printf("Seed: %d, N: %d, is_undirect: %d, redes: %d, prob: %.2f\n", seed, N, is_undirect, redes, prob);
    // Chamar a função de geração do modelo com os parâmetros validados
    generate_local_configuration_model(N, is_undirect, prob, redes, seed);
    
    return 0;
}
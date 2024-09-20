#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <omp.h>

#include "bib/mtwister.h"
#include "bib/calc.h"
#include "bib/LCM.h"
//#include "bib/define.h"



int main(int argc,char *argv[ ]){
    int seed = atoi(argv[1]);
    int N = atoi(argv[2]);
    int categorias = atoi(argv[3]);
    bool is_undirect = atoi(argv[4]);
    int redes = atoi(argv[5]);
    double prob = atof(argv[6])/100;
    generate_local_configuration_model(N,categorias,is_undirect,prob ,redes,seed);

    
}

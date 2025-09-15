gcc -I/home/carlos/igraph/build/include -I/home/carlos/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=102

./main 6437 30000 1 100 0
#./main 752611724138 20 0 100 0
#time ./main 2835 8 0 2 100
gcc -g -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++

gdb ./main
#  235 1000 1 500 0
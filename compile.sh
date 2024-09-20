gcc -I/home/carlos/igraph/build/include -I/home/carlos/igraph/include main.c bib/mtwister.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
./main 100 100 2 0 1 0
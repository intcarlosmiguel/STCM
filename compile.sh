gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=843
for ((i = 0; i <= 100; i += 1)); do
    ./main $seed 10000 1 500 $i
done
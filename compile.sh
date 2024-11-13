gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=183
for ((i = 26; i <= 100; i += 1)); do
    ./main $seed 100000 1 100 $i
    ((seed += 500))
done
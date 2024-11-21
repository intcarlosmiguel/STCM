gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=118
for ((i = 0; i <= 100; i += 1)); do
    ./main $seed 50000 1 100 $i
    ((seed += 100))
done
for ((i = 0; i <= 100; i += 1)); do
    ./main $seed 75000 1 100 $i
    ((seed += 100))
done

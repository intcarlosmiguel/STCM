gcc -I/home/carlos/igraph/build/include -I/home/carlos/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=218
for ((i = 54; i <= 100; i += 1)); do
    ./main $seed 20 0 100 $i
    ((seed += 100))
done

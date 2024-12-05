gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=20
for j in {1..100}
do
    ./main $seed 100 0 100 $j
    ((seed+=100))
done

#time ./main 2835 8 0 2 100
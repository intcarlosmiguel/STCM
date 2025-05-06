gcc -I/home/carlos/igraph/build/include -I/home/carlos/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=102

for j in 0 25 50 75 100
do
    ./main $seed 20 0 100 $j
    ((seed+=100))
done

#time ./main 2835 8 0 2 100
gcc -I/home/carlos/igraph/build/include -I/home/carlos/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=20
for i in {1..6}
do
    for j in {1..100}
    do
        ./main $seed $i 0 100 $j
        ((seed+=100))
    done
done
#time ./main 2835 8 0 2 100
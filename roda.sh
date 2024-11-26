gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=2835
for i in {1..100..2}
do
    for j in {0..100}
    do
        ./main $seed $i 0 100 $j
        ((seed+=100))
    done
done

#time ./main 2835 8 0 2 100
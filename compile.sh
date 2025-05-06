gcc -I/home/carlos/igraph/build/include -I/home/carlos/igraph/include main.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=208
for ((i = 17000; i < 20000; i += 1000)); do


    for ((j = 0; j <= 100; j += 25));
    do
        ./main $seed $i 1 50 $j
        ((seed += 100))
    done
    
done

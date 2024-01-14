#set -x
NITER=1
for((iter=0; iter<$NITER; iter++))
do
    counter=1;
    for((j=1; j<=20; j++))
    do
        ./ms2geno 10000000 1000 1 slim.output_s${counter}_it${iter}.p1.txt > slim.output_s${counter}_it${iter}.p1.txt.geno.txt
        ./ms2geno 10000000 1000 1 slim.output_s${counter}_it${iter}.p2.txt > slim.output_s${counter}_it${iter}.p2.txt.geno.txt
            counter=$((counter+1));
    done
done


  #set -x
NITER=1
for((iter=0; iter<$NITER; iter++))
do
    counter=1;
    for((j=1; j<=20; j++))
    do
        perl ./merge_geno_tables.pl -i1 slim.output_s${counter}_it${iter}.p1.txt.geno.txt -i2 slim.output_s${counter}_it${iter}.p2.txt.geno.txt -t1 slim.output_s${counter}_it${iter}.p1_table.txt -t2 slim.output_s${counter}_it${iter}.p2_table.txt -n1 500 -n2 500 -o1 slim.output_s${counter}_it${iter}.p1p2.txt.geno.txt -o2 slim.output_s${counter}_it${iter}.p1p2_table.txt
        counter=$((counter+1));
    done
done

for((iter=0; iter<$NITER; iter++))
do
    counter=1;
    for((j=1; j<=20; j++))
    do
        nr=`awk 'END{print NR}' slim.output_s${counter}_it${iter}.p1p2.txt.geno.txt`
        nrows=$(($nr-2))
        #nrows=$(($nrows/100))
        ./Rsbki slim.output_s${counter}_it${iter}.p1p2.txt.geno.txt $nrows 1000 0.1 0.05 $RANDOM 2 2 500 500 POP1 POP2
        ./Rsbki slim.output_s${counter}_it${iter}.p1p2.txt.geno.txt $nrows 1000 0.1 0.05 $RANDOM 2 1 500 500 POP1 POP2
        counter=$((counter+1));
    done
done


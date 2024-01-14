NITER=1
PB2_VALUES=( 0.0010 0.0010 0.0010 0.0010 0.0010 0.01000 0.01000 0.01000 0.01000 0.01000 );
MB2_VALUES=( 0.0125 0.0125 0.0125 0.0125 0.0125 0.00125 0.00125 0.00125 0.00125 0.00125 );
NSWEEPS=( 0 0 0 0 0 5 5 5 5 5 );

DISP_VALUES=(  0.00 0.05 0.20 0.20 0.20 0.00 0.05 0.20 0.20 0.20 );
TAU_VALUES=(   0.00 0.10 0.25 0.50 0.80 0.00 0.10 0.25 0.50 0.80 );
SIGMA_VALUES=( 0.60 0.50 0.40 0.30 0.10 0.60 0.50 0.40 0.30 0.10 );

for((iter=0; iter<$NITER; iter++))
do
    counter=11;
    for((j=0; j<10; j++))
    do
        RN=$RANDOM
        sleep 1
echo slim -t -m -d \"seed=$RN\" -d \"nloci=1000\" -d \"nchrom=1\" -d \"nblank_chr=0\" -d \"size_chrom=10000000\" -d \"size_gene=1500\" -d \"Ne_pop=500\" -d \"t0=5\*Ne_pop\" -d \"tend=asInteger\(0.10\*Ne_pop\)\" -d \"maprec=\'flat\'\" -d \"mutation_rate=1.25e-7\" -d \"max_recrate=1.25e-7\" -d \"min_recrate=1.25e-7\" -d \"s_mean_beneficial_p1=0.0125\" -d \"s_mean_beneficial_p2=${MB2_VALUES[$j]}\" -d \"s_mean_deleterious=-1e-6\" -d \"shape_deleterious=0.2\" -d \"shape_beneficial=1\" -d \"prop_beneficial_p1=0.001\" -d \"prop_beneficial_p2=${PB2_VALUES[$j]}\" -d \"tau=${TAU_VALUES[$j]}\" -d \"sigma=${SIGMA_VALUES[$j]}\" -d \"disp=${DISP_VALUES[$j]}\" -d \"h_mean=0.36\" -d \"refn=${counter}\" -d \"nsweeps=${NSWEEPS[$j]}\" -d \"s_beneficial_sweep=0.1\" -d \"file_input1=\'EW_function_neg_plus_posB.R\'\" -d \"file_output1=\'slim.output_s${counter}_it${iter}\'\" ./EHH-GWAS-sims4_blankChr.slim
slim -t -m -d "seed=$RN" -d "nloci=1000" -d "nchrom=1" -d "nblank_chr=0" -d "size_chrom=10000000" -d "size_gene=1500" -d "Ne_pop=500" -d "t0=5*Ne_pop" -d "tend=asInteger(0.10*Ne_pop)" -d "maprec='flat'" -d "mutation_rate=1.25e-7" -d "max_recrate=1.25e-7" -d "min_recrate=1.25e-7" -d "s_mean_beneficial_p1=0.0125" -d "s_mean_beneficial_p2=${MB2_VALUES[$j]}" -d "s_mean_deleterious=-1e-6" -d "shape_deleterious=0.2" -d "shape_beneficial=1" -d "prop_beneficial_p1=0.001" -d "prop_beneficial_p2=${PB2_VALUES[$j]}" -d "tau=${TAU_VALUES[$j]}" -d "sigma=${SIGMA_VALUES[$j]}" -d "disp=${DISP_VALUES[$j]}" -d "h_mean=0.36" -d "refn=${counter}" -d "nsweeps=${NSWEEPS[$j]}" -d "s_beneficial_sweep=0.1" -d "file_input1='EW_function_neg_plus_posB.R'" -d "file_output1='slim.output_s${counter}_it${iter}'" ./EHH-GWAS-sims4_blankChr.slim
        counter=$((counter+1));
    done
done


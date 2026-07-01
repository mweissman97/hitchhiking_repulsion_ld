#!/bin/bash
# An example pipeline for how to run and process simulations using hitchhiking_basic_sim.txt
# The paramter values used here are those used in Fig. 1B,E and Fig. 3

# Load modules
. /u/local/Modules/default/init/modules.sh
module load python
module load slim/4.0.1
module load conda
module load R

# Create new output directory
OUTPUT_DIR="hitchhiking_basic_example"
mkdir -p $OUTPUT_DIR

# Define parameter values
NREP=100 # Number of replicates. Throughout the paper, we used 1000 per simulation scenario, which will take a long time
S_BEN=0.05
S_DEL=0.001
RHO_EXP=6
MU_EXP=6
THETA=0.01
H_DEL=0.5
H_BEN=0.5

# Run simulations NREP times
for i in $(seq 1 $NREP); do 
    FILE_PATH="$OUTPUT_DIR/replicate_$i"
    # 1, run slim simulation
    slim -d S_BEN=$S_BEN -d S_DEL=$S_DEL -d RHO_EXP=$RHO_EXP -d MU_EXP=$MU_EXP -d THETA=$THETA -d H_DEL=$H_DEL -d H_BEN=$H_BEN -d "FILE_PATH='$FILE_PATH'" simulations/hitchhiking_basic_sim.txt

    # 2, process SLiM's outputs to create a haplotype csv for downstream analyses
    python sim_processing/make_ms_haplotypes_slim4.py --file_path "$FILE_PATH"

    # 3, calculate pairwise linkage disequilibrium
    Rscript --vanilla sim_processing/pairwise_ld.R $FILE_PATH
done

# Across replicates, summarize pairwise LD
# Concatenate all individual replicate pairwise files into 1 file 
cd $OUTPUT_DIR
all_pair="pairwise_all.csv"
header=$(head -n 1 replicate_1_pairwise.csv)
echo "$header" > "$all_pair"
for file in replicate_*_pairwise.csv; do
    tail -n +2 "$file" >> "$all_pair"
done
cd ..

# Now, calculate mean LD across replicates.
# This will create 3 csv files that are all named something like "hitchhiking_basic_example/hitchhiking_basic_example_mean_r2_global.csv"
Rscript --vanilla sim_processing/mean_r2.R "$OUTPUT_DIR/$all_pair" $S_BEN $S_DEL $RHO_EXP $MU_EXP $THETA $H_DEL $H_BEN 0 1 "$OUTPUT_DIR/$OUTPUT_DIR"

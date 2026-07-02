# Overview
This repository contains code for the manuscript "The counteracting effects of hitchhiking and repulsion of deleterious alleles during a selective sweep" (Weissman, Lohmueller, and Garud 2026). It is split into 4 subdirectories:
```simulations```: contains the [SLiM 4]([https://messerlab.org/slim/](https://www.journals.uchicago.edu/doi/10.1086/723601?__cf_chl_f_tk=f1IcuQK8F2UzTzjpkZifidwsr5MgAmy3uZmmbfgvJHA-1782863354-1.0.1.1-yryJ4SmSgxPEfLUbhizaPaTbJFCCU265agQZjlOsF9M)) simulation scripts.
* ```sim_processing```: contains R and Python scripts used to convert SLiM 4 outputs into haplotype csv files and then measure linkage disequilibrium.
* ```make_sim_figures```: contains the csv files created through ```simulations``` and ```sim_processing``` and associated R files to recreate all simulation based figures in the paper
* ```c_diff```: contains the *Clostridium difficile* haplotype csv, results of the [iLDS](https://github.com/garudlab/iLDS) scan, and R scripts needed to create figures related to the *Clostridium difficile* analysis (Fig. 4, Supplemental Fig. 4)

In the main directory, the file ```example_hitchhiking_basic_pipeline.sh``` demonstrates the steps to run and process simulations.

# File Directory
## simulations/
Simulations were written using [SLiM 4]([https://messerlab.org/slim/](https://www.journals.uchicago.edu/doi/10.1086/723601?__cf_chl_f_tk=f1IcuQK8F2UzTzjpkZifidwsr5MgAmy3uZmmbfgvJHA-1782863354-1.0.1.1-yryJ4SmSgxPEfLUbhizaPaTbJFCCU265agQZjlOsF9M)). Simulations were run on [hoffman2](https://www.hoffman2.idre.ucla.edu/About/System-overview.html), a Linux compute cluster.
* [hitchhiking_basic_sim.txt](https://github.com/mweissman97/hitchhiking_repulsion_ld/blob/7b8e5472f4625886472a2ce7886f03fdc6c7846b/simulations/hitchhiking_basic_sim.txt): Primary, constant Ne model that simulates a partial sweep of a single, central adaptive locus. Used to generate data for figures: 1B-E, 2A-F, 3A-B, S1, S2, S3A-B, S4A-D.
  * Arguments:
     * S_BEN: selection coefficient of beneficial variant, float. 
     * S_DEL: absolute value of selection coefficient of background deleterious variants, float
     * RHO_EXP: log10 recombination rate, can be a float or an integer
     * MU_EXP: log10 mutation rate for background delterious and neutral variants, can be a float or an integer
     * THETA: population scaled mutation rate for focal beneficial mutation, float
     * H_DEL: dominance coefficient for deleterious variants, float
     * H_BEN: dominance coefficient for beneficial variant, float
     * FILE_PATH: a string that specifies what to call the output files, ie "path/replicate_1.1". This should not include the suffix, as SLiM will create multiple outputs with the same root name
  * Outputs:
     * FILE_PATH.burnin.txt: the results of SLiM's outputFull() after the burn-in, but before the sweep has begun
     * FILE_PATH.burnin.ms: the results of SLiM's outputMSSample() after the burn-in, but before the sweep has begun
     * FILE_PATH.txt: the results of SLiM's outputFull() at the conclusion of the simulation, once the mutation has partially swept
     * FILE_PATH.ms: the results of SLiM's outputMSSample() at the conclusion of the simulation, once the mutation has partially swept
* [epistasis_introgression_sim.txt](https://github.com/mweissman97/hitchhiking_repulsion_ld/blob/7b8e5472f4625886472a2ce7886f03fdc6c7846b/simulations/epistasis_introgression_sim.txt): Simulates a selective sweep of a multi-locus adaptation adaptation that is introduced via introgression / horizontal gene transer. Used to generate data for figure: 5.
  * Arguments:
     * MODEL: an integer 1-4 that corresponds to the fitness interactions of mutations in the adaptive fragment
     * NUM_SITES: an integer that represents the number of interacting loci in the adaptive fragment. The paper uses 1 and 10.
     * FILE_PATH: a string that specifies what to call the output files, ie "path/replicate_1.1"
  * Outputs:
     * FILE_PATH.txt: the results of SLiM's outputFull() at the conclusion of the simulation, once the mutation has partially swept
     * FILE_PATH.ms: the results of SLiM's outputMSSample() at the conclusion of the simulation, once the mutation has partially swept
* [longshallow.txt](https://github.com/mweissman97/hitchhiking_repulsion_ld/blob/7b8e5472f4625886472a2ce7886f03fdc6c7846b/simulations/longshallow.txt): Simulates a long, shallow demographic contraction where the population size is halved and the contraction proceeds for 5000 generations. Used to generate data for figure: S1
  * Arguments:
     * FILE_PATH: a string that specifies what to call the output files, ie "path/replicate_1.1"
  * Outputs:
     * FILE_PATH.txt: the results of SLiM's outputFull() at the conclusion of the simulation, once the mutation has partially swept
     * FILE_PATH.ms: the results of SLiM's outputMSSample() at the conclusion of the simulation, once the mutation has partially swept
* [sharpshort.txt](https://github.com/mweissman97/hitchhiking_repulsion_ld/blob/7b8e5472f4625886472a2ce7886f03fdc6c7846b/simulations/sharpshort.txt): Simulates a short, sharp demographic contraction where the population size is decreased to 0.1*N and the contraction proceeds for 20 generations. Used to generate data for figure: S1
  * Arguments:
     * FILE_PATH: a string that specifies what to call the output files, ie "path/replicate_1.1"
  * Outputs:
     * FILE_PATH.txt: the results of SLiM's outputFull() at the conclusion of the simulation, once the mutation has partially swept
     * FILE_PATH.ms: the results of SLiM's outputMSSample() at the conclusion of the simulation, once the mutation has partially swept

## sim_processing/
Python and R files used to process simulation outputs.
* [make_ms_haplotypes_slim4.py](https://github.com/mweissman97/hitchhiking_repulsion_ld/blob/419bdbbd97f59bf80411f43b21a6ad6e0b04f02c/sim_processing/make_ms_haplotypes_slim4.py): Creates a haplotype csv to calculate LD with information about whether SNPs or synonymous or nonsynonymous.
  * Arguments:
     * --file_path: a string that specifies what to call the output files, ie "path/replicate_1.1". It will read both the FILE_PATH.txt and FILE_PATH.ms created in the SLiM simulations.
  * Outputs:
     * FILE_PATH_raw.csv: a haplotype csv where each row corresponds to a SNP, the first column gives the position of the site, the second gives whether the SNP is syn or nonsyn, and the remaining columns (X0:X99) are states for a given individual where 0 is ancestral and 1 is derived.
* [pairwise_ld.R](https://github.com/mweissman97/hitchhiking_repulsion_ld/blob/ec5793f9b887d11ef165cadd7d2bfcca7fd24bba/sim_processing/pairwise_ld.R): uses a haplotype csv to calculate pairwise LD for all SNP pairs.
  * Arguments:
     * FILE_PATH: path to haplotype csv, minus the suffix "_raw.csv"
  * Outputs:
     * FILE_PATH_pairwise.csv: a csv dataframe with the columns:
       * snp1: site position of first variant (A)
       * snp2: site position of second variant (B)
       * dist: distance, in base pairs, between snp1 and snp2
       * D: LD, calculated as D = P<sub>AB</sub> - P<sub>A</sub> * P<sub>B</sub>
       * r2: LD, calculated as D<sup>2</sup> / (P<sub>A</sub> * (1-P<sub>A</sub>) * P<sub>B</sub> * (1-P<sub>B</sub>))
       * D_prime: LD, calculated as D / min(P<sub>A</sub> * (1-P<sub>A</sub>), P<sub>B</sub> * (1-P<sub>B</sub>))
       * site_type: whether both snp1 and snp2 are "syn" or "nonsyn"
       * afreq_bin: coarse allele frequency of both snp1 and snp2, where "common" corresponds to MAF >= 0.2, "rare" corresponds to a MAF <= 0.05, and "int" corresponds to intermediate allele frequencies
       * snp1_freq: absolute allele frequency of snp1
       * snp2_freq: absolute allele frequency of snp2
* [mean_r2.R](https://github.com/mweissman97/hitchhiking_repulsion_ld/blob/ec5793f9b887d11ef165cadd7d2bfcca7fd24bba/sim_processing/mean_r2.R): Calculates mean r2 from pairwise LD file. Arguments include parameters used in simulation to more easily concatenate files downstream.
  * Arguments:
     * FILE_NAME: full path to the FILE_PATH_pairwise.csv file created by pairwise_ld.R
     * S_BEN: selection coefficient of beneficial variant, float. 
     * S_DEL: absolute value of selection coefficient of background deleterious variants, float
     * RHO_EXP: log10 recombination rate, can be a float or an integer
     * MU_EXP: log10 mutation rate for background delterious and neutral variants, can be a float or an integer
     * THETA: population scaled mutation rate for focal beneficial mutation, float
     * H_DEL: dominance coefficient for deleterious variants, float
     * H_BEN: dominance coefficient for beneficial variant, float
     * MODEL: for introgression epistasis simulations: an integer 1-4 that corresponds to the fitness interactions of mutations in the adaptive fragment. For single locus, additive simulations, set a dummy 0 value
     * NUM_SITES: an integer that represents the number of interacting loci in the adaptive fragment
     * OUT_NAME: desired prefix for output filename. Output file name should not include suffix (ie .csv) because we will make several output files that all start with the same root.
  * Outputs: All output files also include columns corresponding to above parameters (ie s_ben and model)
     * OUT_NAME_mean_r2_global.csv: mean r2 by mutation type, allele frequency, and distance bin. Includes columns:
       * site_type: whether both snp1 and snp2 are "syn" or "nonsyn"
       * allele_freq: coarse allele frequency of both snp1 and snp2, where "common" corresponds to MAF >= 0.2, "rare" corresponds to a MAF <= 0.05, and "int" corresponds to intermediate allele frequencies
       * dist_bin: maximum distance in fine-scaled, non-overlapping distance bins
       * mean_D: mean D across all pairs
       * mean_r2: mean r2 across all pairs
       * mean_D_prime: mean D_prime across all pairs
       * num_pairs: number of pairs
     * OUT_NAME_paired_delta_r2: mean r2 for nonsyn and syn pairs that are precisely matched for distance apart and allele frequency. Includes columns:
       * dist: precise distance between variants
       * afreq_bin: coarse allele frequency of both snp1 and snp2, where "common" corresponds to MAF >= 0.2, "rare" corresponds to a MAF <= 0.05, and "int" corresponds to intermediate allele frequencies
       * snp1_freq: precise frequency of the first SNP
       * snp2_freq: precise frequency of the second SNP
       * rN: mean r2 for all nonsynonymous pairs that are precisely dist base pairs apart and have the same snp1_freq and snp2_freq
       * rS: mean r2 for all synonymous pairs that are precisely dist base pairs apart and have the same snp1_freq and snp2_freq
       * delta_r2: rN-rS for that combination of dist, snp1_freq, and snp2_freq
     * OUT_NAME_spatial_ld.csv: mean r2 for common (MAF > 0.2) deleterious-deleterious nonsyn and syn pairs as a function of both distance apart and distance from the focal beneficial variant. Includes columns:
       * site_type: whether both snp1 and snp2 are "syn" or "nonsyn"
       * dist_bin: maximum distance between variants in fine-scaled, non-overlapping distance bins
       * ben_dist_breaks: cut-off for distance from beneficial variant to be considered "near" vs. "far"
       * ben_proximity: whether pairs are "near" (ie both SNPs in pair are < ben_dist_breaks away from beneficial variant) are "far" (ie both SNPs in pair are > ben_dist_breaks away from beneficial variant)
       * mean_D: mean D across all pairs
       * mean_r2: mean r2 across all pairs
       * mean_D_prime: mean D_prime across all pairs
       * num_pairs: number of pairs

## make_sim_figures/
Csv files generated by simulations and R scripts to recreate all figures in manuscript.

## c_diff/
Data and post processing files for running iLDS and plotting LD decay for *Clostridium difficile*.
* c_diff.txt: Haplotype file for *Clostridium difficile* genomes from the Unified Human Gastrointenstinal Genome (UHGG) project (alias MGYG-HGUT-02369)
* [pi_calcs_simple.py](https://github.com/mweissman97/hitchhiking_repulsion_ld/blob/a406167c12051687c7b18b85c3d8a82110c7c975/c_diff/pi_calcs_simple.py): Calculates average nucleotide diversity (pi) in overlapping windows of size 1000 bp with a step size of 500 bp between windows. Uses c_diff.txt.
* [C_diff_pi_windows.tsv](https://github.com/mweissman97/hitchhiking_repulsion_ld/blob/a406167c12051687c7b18b85c3d8a82110c7c975/c_diff/C_diff_pi_windows.tsv): result of pi_calcs_simple.py
* [C_diff_full_scan.txt](https://github.com/mweissman97/hitchhiking_repulsion_ld/blob/a406167c12051687c7b18b85c3d8a82110c7c975/c_diff/C_diff_full_scan.txt): result of running [iLDS](https://github.com/garudlab/iLDS) pipeline on c_diff.txt
* [figures_cdiff.R](https://github.com/mweissman97/hitchhiking_repulsion_ld/blob/a406167c12051687c7b18b85c3d8a82110c7c975/c_diff/figures_cdiff.R): R files used to create all figures related to *C. difficile*; Figure 4, Supplemental Figure 4. This includes plotting iLDS results, pi in sliding windows, rN - rS graphs for all sweeps, and haplotype visualizations.

## make_figures/
R scripts used to generate all figures for paper.

# License
This project is covered under the MIT License.

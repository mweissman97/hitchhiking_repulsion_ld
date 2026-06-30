# This file takes the SLiM outputs of the full output and the MS style to create a csv where each snp is annotated as synonymous or nonsynonymous
# Import necessary libraries
import pandas as pd  # For data manipulation and analysis
import numpy as np  # For numerical operations
import os

# Function to read haplotype data from simulation files
def read_haplotype(file_path, G):
    """
    Reads haplotype data from a simulation file, processes it, and returns a DataFrame.
    """
    
    # Construct the file path for the simulation file
    file_loc = f"{file_path}"
    filename = os.path.abspath(f"{file_path}.ms")
    if not os.path.exists(filename):
        raise FileNotFoundError(f"ERROR: MS file not found: {filename}")

    # Read the simulation file line by line
    with open(filename) as file:
        lines = [line.rstrip() for line in file]
        
    # Extract positions of mutations and scale them by genome length
    pos = [int(np.round(float(i) * G)) for i in lines[2].split(" ")[1:]]
    #print("Parsed positions:", pos)
    #print("Number of positions:", len(pos))

    # Extract haplotype data from the file
    H = []
    for l in lines[3:]:
        H.append(list(l))
        
    # Transpose the haplotype matrix and convert it to a DataFrame
    Hs = pd.DataFrame(zip(*H))

    # Convert haplotype data to integers
    Hs = Hs.astype(int)
    
    # Set the mutation positions as the index
    Hs.index = pos
    Hs.index.name = "site_pos"
    
    # Group by site position and sum the haplotypes
    Hs = Hs.groupby("site_pos").sum()
    print("Hs shape after transpose:", Hs.shape)
    print(Hs.head())

    # Ensure all values are binary (0 or 1)
    Hs = Hs.where(Hs <= 1, 1)
    #print("Hs shape after binary:", Hs.shape)
    #print(Hs.head())

    # Filter out invariant sites (those with all 0s or all 1s)
    #F = (Hs.T.mean() > 0) & (Hs.T.mean() < 1)
    #Hs = Hs.loc[F]
    Hs = Hs.loc[Hs.apply(lambda row: row.nunique() > 1, axis=1)]

    #print("Hs shape after invariant site filtering:", Hs.shape)

    # Print the last 20 rows of the processed haplotype data
    #print(Hs.tail(20))
    
    # Read the mutation information from the corresponding `.txt` file
    mutfile = os.path.abspath(f"{file_loc}.txt")
    if not os.path.exists(mutfile):
        raise FileNotFoundError(f"ERROR: TXT mutation file not found: {mutfile}")    
    with open(mutfile) as file:
        lines = [line.rstrip() for line in file]
        
    # Extract mutation data from the file
    #begin_M = next(i for i, l in enumerate(lines) if "Mutations:" in l) + 1
    #end_M = next(i for i, l in enumerate(lines) if "Haplosomes:" in l)
    begin_M = next(i for i, l in enumerate(lines) if l.strip() == "Mutations:") + 1
    end_M = next(i for i, l in enumerate(lines) if l.strip() == "Individuals:")
    #print("Mutation file first line: ", begin_M)
    L_dict = {}
    print("Begin line idx: ", begin_M)
    print("End line idx: ", end_M)
    for l in lines[begin_M:end_M]:
        #print("Mutfile line: ", l, "/n")
        parts = l.split()
        if len(parts) >= 4:
            mut_pos = int(parts[3])
            mut_type = parts[2]
            L_dict[mut_pos] = mut_type

    L = pd.Series(L_dict, dtype=str)  # values are strings (mutation types)
    L.index = L.index.astype(int)     # ensure positions are integers

    intr = Hs.index.intersection(L.index)
    print("Hs.index types:", [type(x) for x in Hs.index[:10]])
    print("Hs.index values:", list(Hs.index[:10]))

    print("L.index types:", [type(x) for x in L.index[:10]])
    print("L.index values:", list(L.index[:10]))

    print("Intersection size:", len(intr))
    print("Intersection values:", list(intr))
    
    #begin_M = lines.index('Mutations:') + 1
    #end_M = lines.index("Individuals:")
    #L = [l.split(" ") for l in lines[begin_M:end_M]]
    #L = {int(l[3]): l[2] for l in L}
    #L = pd.Series(L)
    
    # Set the mutation information as the index of the haplotype data
    # Find the intersection of mutation positions with the haplotype data
    #matched_positions = []
    #for pos in Hs.index:
        # find closest mutation position
     #   closest = L.index[np.argmin(np.abs(L.index - pos))]
        # optionally: only accept if difference is <= 1 bp
      #  if abs(closest - pos) <= 1:
       #     matched_positions.append(pos)    

    Hs_filtered = Hs.loc[intr]
    # subset mutation types accordingly
    L_filtered = L.loc[[L.index[np.argmin(np.abs(L.index - p))] for p in intr]]


    # Map mutation types to site positions
    M = L_filtered.reset_index()
    M.columns = ["site_pos", "site_type"]
    M["site_type"] = M["site_type"].map({"m1": "syn", "m2": "nonsyn", "m3": "nonsyn", "m4": "nonsyn"})
    M = pd.MultiIndex.from_frame(M)

    # Set the mutation information as the index of the haplotype data
    Hs_filtered.index = M    
    # Return the processed haplotype data

    # Print the shape of the haplotype and mutation data
    print(Hs_filtered.shape)
    #print(M.shape)
    #print(M[-10:])

    return Hs_filtered

def compute_ld_matrix(geno_df):
    """
    Computes the pairwise LD (r^2) matrix and returns it in long format (upper triangle only).
    
    Parameters:
    - geno_df: DataFrame where rows are samples and columns are SNPs (genotypes).
    
    Returns:
    - cor_long_upper: DataFrame with columns ['snp1', 'snp2', 'r2'] for the upper triangle.
    """
    # Compute the Pearson correlation matrix and square it to get r^2 values
    cor_matrix = geno_df.corr(method="pearson") ** 2

    # Convert the correlation matrix to a long DataFrame
    cor_df = cor_matrix.reset_index()
    cor_df = cor_df.rename(columns={"index": "snp1_idx"})  # Rename the index column to 'snp1'
    cor_long = cor_df.melt(id_vars=["snp1_idx"], var_name="snp2_idx", value_name="r2")

    # Filter to include only the upper triangle (snp1 < snp2)
    cor_long_upper = cor_long[cor_long["snp1_idx"] < cor_long["snp2_idx"]]
    #convert snp indexes to actual snp values

    return cor_long_upper

def add_meta_to_ld(cor_long_upper, df_meta, rep):
    """
    Adds metadata to the LD matrix.
    
    Parameters:
    - cor_long_upper: DataFrame with columns ['snp1', 'snp2', 'r2'] (output of compute_ld_matrix).
    - df_meta: DataFrame with metadata for SNPs, including 'site_pos', 'site_type', and 'site_freq'.
    - rep: Replicate identifier (e.g., an integer or string).
    
    Returns:
    - ld_matrix_full: DataFrame with additional metadata columns.
    """
   
    # Merge to map "snp1_idx" to "site_pos" for "snp1_pos"
    cor_long_upper = cor_long_upper.merge(
        df_meta.rename(columns={"site_pos": "snp1_pos",
                           "site_type": "site_type_1",
                           "site_freq": "site_freq_1"}), 
        left_on="snp1_idx", 
        right_index=True, 
        how="left"
    )

    # Merge to map "snp2_idx" to "site_pos" for "snp2_pos"
    cor_long_upper = cor_long_upper.merge(
        df_meta.rename(columns={"site_pos": "snp2_pos",
                               "site_type": "site_type_2",
                               "site_freq": "site_freq_2"}), 
        left_on="snp2_idx", 
        right_index=True, 
        how="left"
    )

    # Calculate distance between SNPs
    cor_long_upper["dist"] = abs(cor_long_upper["snp1_pos"] - cor_long_upper["snp2_pos"])

    # Determine SNP type (match or mismatch)
    cor_long_upper["snp_type"] = np.where(
        cor_long_upper["site_type_1"] == cor_long_upper["site_type_2"],
        cor_long_upper["site_type_2"],
        "mismatch"
    )

    # Add replicate identifier
    cor_long_upper["rep"] = rep

    return cor_long_upper

# Main script execution
if __name__ == "__main__":
    # Import argparse for command-line argument parsing
    import argparse
    
    # Set up argument parser
    parser = argparse.ArgumentParser()
    
    # Define command-line arguments
    parser.add_argument('--file_path', help="", type=str)
    parser.add_argument('--repNum', help="", type=str, default="1")
    parser.add_argument('--G', help="", type=int, default=1e4)    
    
    # Parse command-line arguments
    args = parser.parse_args()
        
    file_path = os.path.expandvars(args.file_path)  # expand any $VAR
    file_path = os.path.expanduser(file_path)      # expand ~
    file_path = os.path.abspath(file_path)         # absolute path
    repNum = args.repNum
    G = args.G
    
    # Read and process haplotype data
    df = read_haplotype(file_path, G)
    #print(df.head())
    df.to_csv(f"{file_path}_raw.csv")

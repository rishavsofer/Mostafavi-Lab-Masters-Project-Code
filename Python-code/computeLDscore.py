"""
Written by Rishav Eliyahu Sofer Dasgupta
March 4th, 2025

"""
import pandas as pd
import numpy as np
import os
import gzip
import gc
from pandarallel import pandarallel
from compress_and_uncompress_fns import read_compressed_file

# =============================================================================
# STEP 0: SET GLOBAL PARAMETERS AND READ IN DATA FILES
# =============================================================================

# Define file paths
root_fp = "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/"
ldsc_fp = root_fp + "ldsc/"

# Data files
gwas_fp = "/gpfs/data/mostafavilab/GWAS/50_irnt_standing_height_both_sexes/filtered_height_gwas_data.tsv.gz"
gene_features_fp = root_fp + "gene_features.bed"

# Load GWAS data
gwas_df = pd.read_csv(gwas_fp, sep="\t", compression="gzip")

# Load gene features data
gene_df = pd.read_csv(gene_features_fp, sep="\t", header=None, names=["chr", "start", "end", "gene"], comment="#")
gene_df["chr"] = gene_df["chr"].apply(lambda x: int(x.replace("chr", "")))

# Chromosome range
chromosomes = list(range(1, 23))  # Autosomal chromosomes only

# =============================================================================
# STEP 1: GENERATE ANNOTATION FILE
# =============================================================================

def generate_annotation_file(chrom, gwas_df, gene_df, cm_fp):
    """Generates annotation file for LD Score regression."""
    gwas_chr = gwas_df[gwas_df.CHR == chrom]
    gene_chr = gene_df[gene_df.chr == chrom]
    
    # Read CM values
    cm_df = pd.read_csv(cm_fp, sep="\t", header=None, names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    
    # Find 5 nearest genes
    def nearest_gene_finder(m):
        dists = np.abs(gene_chr["start"] - m.POS)
        top_5 = dists.nsmallest(5)
        out = dict(m)
        for i, (idx, val) in enumerate(top_5.items()):
            out[f"GeneSymbol.f{i+1}"] = gene_chr.iloc[idx]["gene"]
            out[f"dist.f{i+1}"] = max(val / 1000, 1/1000)  # Convert to kb
        return pd.Series(out)
    
    pandarallel.initialize(progress_bar=True, nb_workers=6)
    nearest_genes_chr = gwas_chr.parallel_apply(nearest_gene_finder, axis=1)
    gc.collect()
    
    # Prepare annotation dataframe
    annot_df_prep = nearest_genes_chr.copy(deep=True)
    annot_df_prep = annot_df_prep.rename(columns={"POS": "BP", "variant": "SNP"})
    annot_df_prep = annot_df_prep.merge(cm_df[["SNP", "CM"]], on="SNP", how="left").fillna({"CM": 1.0})
    
    return annot_df_prep

# =============================================================================
# STEP 2: FILTER & ORDER ANNOTATION FILE
# =============================================================================

def filter_annotation_file(chrom, root_fp, ldsc_fp, annot_df_prep):
    """Filters and orders the annotation file to match SNP ordering."""
    filtered_bim_fp = f"{root_fp}/filtered_plink_triplets/{chrom}_filt.bim.gz"
    filtered_bim = read_compressed_file(filtered_bim_fp)
    
    annot_df = annot_df_prep.set_index("SNP").loc[filtered_bim.SNP.values].reset_index()
    
    annot_fp = f"{ldsc_fp}/chr{chrom}/s_het_w_chr{chrom}.annot"
    annot_df.to_csv(annot_fp, sep="\t", index=False)
    print(f"Chromosome {chrom} annotation file written.")

# =============================================================================
# STEP 3: RUN FOR ALL CHROMOSOMES
# =============================================================================

for chrom in chromosomes:
    branch_fp = f"{ldsc_fp}/chr{chrom}/"
    cm_fp = f"{root_fp}/ld_score_outputs/CM_values/chr{chrom}.bim"
    
    annot_df_prep = generate_annotation_file(chrom, gwas_df, gene_df, cm_fp)
    filter_annotation_file(chrom, root_fp, ldsc_fp, annot_df_prep)

# =============================================================================
# STEP 4: PREPARE LD SCORE COMPUTATION
# =============================================================================

def prep_ld_score_compute(chrom, root_fp, ldsc_fp):
    """Prepares command for LD Score computation."""
    out_fp = f"{root_fp}/ld_score_outputs/extracted_plink_triplet/{chrom}_filt"
    annot_fp = f"{ldsc_fp}/chr{chrom}/s_het_w_chr{chrom}.annot"
    ld_score_fp = f"{ldsc_fp}/chr{chrom}/ld_scores/s_het_w_chr{chrom}"
    
    ld_score_cmd = f"python {ldsc_fp}/ldsc.py --l2 --bfile {out_fp} --ld-wind-cm 1 --annot {annot_fp} --out {ld_score_fp}"
    print(ld_score_cmd)

# Generate LD Score commands for all chromosomes
for chrom in chromosomes:
    prep_ld_score_compute(chrom, root_fp, ldsc_fp)


# Install required packages
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
library(biomaRt)

# Connect to Ensembl database
ensembl <- useEnsembl(biomart="ensembl", 
                        dataset="hsapiens_gene_ensembl", 
                        host="https://grch37.ensembl.org")

# Get the transcript for ABL1
transcript_info <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name"),
  filters = "external_gene_name",
  values = "ABL1",
  mart = ensembl
)

# Select a canonical transcript (first one for simplicity)
canonical_transcript <- transcript_info$ensembl_transcript_id[1]

# Query Ensembl for coding sequence and peptide sequence for ABL1
coding_info <- getBM(
  attributes = "coding",
  filters = "external_gene_name",
  values = "ABL1",
  mart = ensembl
)

peptide_info <- getBM(
  attributes = "peptide",
  filters = "external_gene_name",
  values = "ABL1",
  mart = ensembl
)

# Extract the coding sequence and peptide
coding_sequence <- coding_info$coding[1]
peptide_sequence <- peptide_info$peptide[1]

# Split coding sequence into codons
codons <- strsplit(coding_sequence, "")[[1]]
codons <- paste0(codons[c(TRUE, FALSE, FALSE)], codons[c(FALSE, TRUE, FALSE)], codons[c(FALSE, FALSE, TRUE)])

# Map codons to their respective amino acids
amino_acids <- unlist(strsplit(peptide_sequence, ""))

# Check SNP feasibility for each amino acid mutation
library(dplyr)

# Codon to amino acid map
codon_table <- list(
  "A" = c("GCT", "GCC", "GCA", "GCG"),
  "R" = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
  "N" = c("AAT", "AAC"),
  "D" = c("GAT", "GAC"),
  "C" = c("TGT", "TGC"),
  "Q" = c("CAA", "CAG"),
  "E" = c("GAA", "GAG"),
  "G" = c("GGT", "GGC", "GGA", "GGG"),
  "H" = c("CAT", "CAC"),
  "I" = c("ATT", "ATC", "ATA"),
  "L" = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
  "K" = c("AAA", "AAG"),
  "M" = c("ATG"),
  "F" = c("TTT", "TTC"),
  "P" = c("CCT", "CCC", "CCA", "CCG"),
  "S" = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
  "T" = c("ACT", "ACC", "ACA", "ACG"),
  "W" = c("TGG"),
  "Y" = c("TAT", "TAC"),
  "V" = c("GTT", "GTC", "GTA", "GTG")
)

# Function to check SNP feasibility
is_snp_feasible <- function(original_codon, target_codon) {
  sum(unlist(strsplit(original_codon, "")) != unlist(strsplit(target_codon, ""))) == 1
}

# Iterate over each codon and amino acid
results <- data.frame()

for (i in seq_along(codons)) {
  original_codon <- codons[i]
  original_aa <- amino_acids[i]
  
  for (target_aa in names(codon_table)) {
    if (target_aa != original_aa) {
      for (target_codon in codon_table[[target_aa]]) {
        if (is_snp_feasible(original_codon, target_codon)) {
          results <- rbind(
            results,
            data.frame(
              Position = i,
              OriginalCodon = original_codon,
              OriginalAA = original_aa,
              TargetCodon = target_codon,
              TargetAA = target_aa,
              SNPFeasible = TRUE
            )
          )
          break
        }
      }
    }
  }
}

# Display results
print(results)
df.out = results %>% 
  mutate(mut = paste0(OriginalAA,Position,TargetAA))

# Save results to a file
write.csv(df.out, "./code/Refs/ABL1_SNP_Accessibility.csv", row.names = FALSE)

# # Append snp accessibility to Twist Mut dataframe
# df.twist = read.csv('code/Refs/Twist_ABL_All_Data.csv')
# 
# df.merge = df.out %>% select(mut,SNPFeasible)
# 
# df.updated = merge(df.twist, df.merge, by = 'mut', all.x = T) %>% 
#   mutate(snp = if_else(is.na(SNPFeasible), FALSE, SNPFeasible)) %>% 
#   select(-SNPFeasible)
# 
# write.csv(df.updated, "./code/Refs/Twist_ABL_All_Data.csv", row.names = FALSE)

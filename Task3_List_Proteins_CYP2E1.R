# ============================================================
# TASK 3: List All Proteins for CYP2E1 (ENSG00000130649)
# Gene: CYP2E1 | Tissue: Liver | Breaks down toxins & alcohol
# ============================================================

# ---- Install Bioconductor packages if not already installed ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "ensembldb",
  "EnsDb.Hsapiens.v86",
  "AnnotationFilter"
), ask = FALSE)

# ---- Load Libraries ----
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(AnnotationFilter)

# ---- Set up the Ensembl database ----
edb <- EnsDb.Hsapiens.v86

# ---- Define CYP2E1 gene ID ----
gene_id <- "ENSG00000130649"

# ============================================================
# LIST ALL PROTEINS FOR CYP2E1
# ============================================================

# Method 1: proteins() function — returns protein info directly
proteins_cyp2e1 <- proteins(
  edb,
  filter  = GeneIdFilter(gene_id),
  columns = c(
    "protein_id",          # UniProt / Ensembl protein ID
    "tx_id",               # Associated transcript
    "gene_id",             # Gene Ensembl ID
    "gene_name",           # Gene symbol
    "protein_sequence"     # Amino acid sequence (FASTA)
  )
)

# ---- Display results ----
cat("==============================================\n")
cat("  CYP2E1 (ENSG00000130649) - All Proteins\n")
cat("==============================================\n\n")

cat("Total number of proteins:", length(proteins_cyp2e1), "\n\n")

# Convert to data frame
prot_df <- as.data.frame(proteins_cyp2e1)
print(prot_df)

# ---- Show protein IDs and their transcript links ----
cat("\n----------------------------------------------\n")
cat("Protein ID <-> Transcript ID Mapping:\n")
cat("----------------------------------------------\n")
mapping <- prot_df[, c("protein_id", "tx_id", "gene_name")]
print(mapping)

# ---- Protein sequence lengths ----
cat("\n----------------------------------------------\n")
cat("Protein Sequence Lengths (amino acids):\n")
cat("----------------------------------------------\n")
if ("protein_sequence" %in% colnames(prot_df)) {
  prot_df$aa_length <- nchar(as.character(prot_df$protein_sequence))
  print(prot_df[, c("protein_id", "tx_id", "aa_length")])
}

# ---- Method 2: Using select() ----
cat("\n----------------------------------------------\n")
cat("Tabular View via select():\n")
cat("----------------------------------------------\n")
prot_table <- select(
  edb,
  keys    = gene_id,
  keytype = "GENEID",
  columns = c("GENEID", "GENENAME", "TXID", "PROTEINID")
)
# Keep only rows that have a protein (removes non-coding transcripts)
prot_table <- prot_table[!is.na(prot_table$PROTEINID), ]
print(prot_table)

# ---- Save to CSV ----
prot_df_out <- prot_df
if ("protein_sequence" %in% colnames(prot_df_out)) {
  # Truncate sequence for CSV readability; full seq in FASTA below
  prot_df_out$protein_sequence <- substr(as.character(prot_df_out$protein_sequence), 1, 50)
  prot_df_out$protein_sequence <- paste0(prot_df_out$protein_sequence, "...")
}
write.csv(prot_df_out, file = "CYP2E1_proteins.csv", row.names = FALSE)
cat("\nResults saved to: CYP2E1_proteins.csv\n")

# ---- Write FASTA file for each protein sequence ----
if ("protein_sequence" %in% colnames(prot_df)) {
  fasta_file <- "CYP2E1_protein_sequences.fasta"
  con <- file(fasta_file, open = "w")
  for (i in seq_len(nrow(prot_df))) {
    header <- paste0(">", prot_df$protein_id[i],
                     " | tx:", prot_df$tx_id[i],
                     " | gene:", prot_df$gene_name[i])
    writeLines(header, con)
    writeLines(as.character(prot_df$protein_sequence[i]), con)
  }
  close(con)
  cat("Protein sequences saved to:", fasta_file, "\n")
}

# ============================================================
# EXPECTED OUTPUT:
# CYP2E1 has 2 protein-coding transcripts, so 2 proteins:
#   ENSP00000349420 — linked to ENST00000360332 (493 aa, canonical)
#   ENSP00000374303 — linked to ENST00000389728 (shorter isoform)
#
# Non-coding transcripts (retained_intron, processed_transcript)
# do NOT produce proteins and will be absent from this list.
# ============================================================

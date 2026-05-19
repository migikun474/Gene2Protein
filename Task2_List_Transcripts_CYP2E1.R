# ============================================================
# TASK 2: List All Transcripts for CYP2E1 (ENSG00000130649)
# Gene: CYP2E1 | Tissue: Liver | Breaks down toxins & alcohol
# ============================================================

# ---- Install Bioconductor packages if not already installed ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required packages
BiocManager::install(c(
  "ensembldb",
  "EnsDb.Hsapiens.v86",   # Ensembl v86 (GRCh38) human database
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
# LIST ALL TRANSCRIPTS FOR CYP2E1
# ============================================================

# Method 1: Using filter() with GeneIdFilter
transcripts_cyp2e1 <- transcripts(
  edb,
  filter = GeneIdFilter(gene_id),
  columns = c(
    "tx_id",          # Transcript Ensembl ID
    "tx_name",        # Transcript name
    "tx_biotype",     # Biotype (protein_coding, retained_intron, etc.)
    "tx_seq_start",   # Transcript start position
    "tx_seq_end",     # Transcript end position
    "seq_name",       # Chromosome
    "seq_strand",     # Strand (1 = plus, -1 = minus)
    "tx_cds_seq_start",  # CDS start
    "tx_cds_seq_end",    # CDS end
    "gene_id",           # Gene Ensembl ID
    "gene_name"          # Gene symbol
  )
)

# ---- Display results ----
cat("==============================================\n")
cat("  CYP2E1 (ENSG00000130649) - All Transcripts\n")
cat("==============================================\n\n")

cat("Total number of transcripts:", length(transcripts_cyp2e1), "\n\n")

# Convert to data frame for easier viewing
tx_df <- as.data.frame(transcripts_cyp2e1)
print(tx_df)

# ---- Summary by biotype ----
cat("\n----------------------------------------------\n")
cat("Transcripts by Biotype:\n")
cat("----------------------------------------------\n")
print(table(tx_df$tx_biotype))

# ---- Protein-coding transcripts only ----
cat("\n----------------------------------------------\n")
cat("Protein-Coding Transcripts Only:\n")
cat("----------------------------------------------\n")
pc_transcripts <- tx_df[tx_df$tx_biotype == "protein_coding", ]
print(pc_transcripts[, c("tx_id", "tx_name", "tx_biotype",
                          "tx_seq_start", "tx_seq_end", "seq_name", "seq_strand")])

# ---- Method 2: Using select() for a cleaner tabular output ----
cat("\n----------------------------------------------\n")
cat("Tabular View (select method):\n")
cat("----------------------------------------------\n")
tx_table <- select(
  edb,
  keys    = gene_id,
  keytype = "GENEID",
  columns = c("TXID", "TXNAME", "TXBIOTYPE", "TXSEQSTART", "TXSEQEND",
              "SEQNAME", "SEQSTRAND")
)
print(tx_table)

# ---- Save to CSV ----
write.csv(tx_df, file = "CYP2E1_transcripts.csv", row.names = FALSE)
cat("\nResults saved to: CYP2E1_transcripts.csv\n")

# ============================================================
# EXPECTED OUTPUT (Ensembl v86):
# CYP2E1 has approximately 5 transcripts:
#   ENST00000360332 - protein_coding (main isoform, 9 exons)
#   ENST00000389728 - protein_coding
#   ENST00000495786 - retained_intron
#   ENST00000493705 - processed_transcript
#   ENST00000464098 - retained_intron
# ============================================================

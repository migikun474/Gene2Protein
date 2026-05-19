# ============================================================
# TASK 4: Protein-to-Genome Map for Each CYP2E1 Protein
# Gene: CYP2E1 | ENSG00000130649 | Liver | Toxin metabolism
#
# Generates a coordinate map: each amino acid position in the
# protein → its corresponding genomic coordinates (GRCh38)
# ============================================================

# ---- Install packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "ensembldb",
  "EnsDb.Hsapiens.v86",
  "AnnotationFilter",
  "IRanges",
  "GenomicRanges"
), ask = FALSE)

# ---- Load libraries ----
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(AnnotationFilter)
library(IRanges)
library(GenomicRanges)

# ---- Database and gene setup ----
edb     <- EnsDb.Hsapiens.v86
gene_id <- "ENSG00000130649"

# ============================================================
# STEP 1: Get all proteins for CYP2E1
# ============================================================
prot_df <- as.data.frame(
  proteins(edb,
           filter  = GeneIdFilter(gene_id),
           columns = c("protein_id", "tx_id", "gene_name", "protein_sequence"))
)

cat("==============================================\n")
cat("  CYP2E1 Protein → Genome Coordinate Mapping\n")
cat("==============================================\n")
cat("Proteins found:", nrow(prot_df), "\n")
cat("Protein IDs:", paste(prot_df$protein_id, collapse = ", "), "\n\n")

# ============================================================
# STEP 2: Map protein coordinates to genome for each protein
# using proteinToGenome() from ensembldb
# ============================================================

# Helper function: map one protein fully and save to file
map_protein_to_genome <- function(edb, protein_id, tx_id, gene_name,
                                  protein_sequence) {

  cat("----------------------------------------------\n")
  cat("Processing protein:", protein_id, "\n")
  cat("Linked transcript:", tx_id, "\n")

  # Protein length
  prot_len <- nchar(as.character(protein_sequence))
  cat("Protein length:", prot_len, "amino acids\n")

  # Define IRanges covering the FULL protein (aa position 1 to end)
  # Each position = 1 amino acid
  prot_range <- IRanges(start = 1, end = prot_len)
  names(prot_range) <- protein_id

  # ---- Map to genome ----
  # proteinToGenome() converts protein aa coordinates → genomic coordinates
  genome_map <- proteinToGenome(prot_range, edb)

  cat("Mapping complete. Genomic ranges returned:", length(genome_map[[1]]), "\n\n")

  # Convert to data frame
  map_df <- as.data.frame(genome_map[[1]])

  # Add metadata columns
  map_df$protein_id       <- protein_id
  map_df$transcript_id    <- tx_id
  map_df$gene_name        <- gene_name
  map_df$protein_aa_start <- 1
  map_df$protein_aa_end   <- prot_len

  # Reorder columns for clarity
  cols_order <- c("protein_id", "transcript_id", "gene_name",
                  "protein_aa_start", "protein_aa_end",
                  "seqnames", "start", "end", "width", "strand",
                  setdiff(colnames(map_df),
                          c("protein_id", "transcript_id", "gene_name",
                            "protein_aa_start", "protein_aa_end",
                            "seqnames", "start", "end", "width", "strand")))
  map_df <- map_df[, intersect(cols_order, colnames(map_df))]

  return(map_df)
}

# ============================================================
# STEP 3: Loop over all proteins and save individual map files
# ============================================================

all_maps <- list()

for (i in seq_len(nrow(prot_df))) {

  pid   <- prot_df$protein_id[i]
  txid  <- prot_df$tx_id[i]
  gname <- prot_df$gene_name[i]
  pseq  <- prot_df$protein_sequence[i]

  # Generate the map
  map_i <- map_protein_to_genome(edb, pid, txid, gname, pseq)

  # Print summary
  cat("Genomic mapping for", pid, ":\n")
  print(map_i[, c("protein_id", "seqnames", "start", "end", "strand")])
  cat("\n")

  # Save individual file per protein
  out_file <- paste0("CYP2E1_protein2genome_", pid, ".csv")
  write.csv(map_i, file = out_file, row.names = FALSE)
  cat("Saved:", out_file, "\n\n")

  all_maps[[pid]] <- map_i
}

# ============================================================
# STEP 4: Save a combined map file (all proteins together)
# ============================================================
combined_map <- do.call(rbind, all_maps)
write.csv(combined_map,
          file = "CYP2E1_protein2genome_ALL_proteins.csv",
          row.names = FALSE)

cat("==============================================\n")
cat("Combined map saved: CYP2E1_protein2genome_ALL_proteins.csv\n")
cat("==============================================\n\n")

# ============================================================
# STEP 5: Visualise the exon-level mapping (optional summary)
# ============================================================
cat("----------------------------------------------\n")
cat("Exon-level genomic blocks per protein:\n")
cat("----------------------------------------------\n")
for (pid in names(all_maps)) {
  m <- all_maps[[pid]]
  cat("\nProtein:", pid, "\n")
  cat(sprintf("  Chromosome: %s | Strand: %s\n",
              unique(m$seqnames), unique(m$strand)))
  cat(sprintf("  Genomic span: %s - %s\n",
              format(min(m$start), big.mark=","),
              format(max(m$end),   big.mark=",")))
  cat(sprintf("  Exonic blocks: %d\n", nrow(m)))
  cat(sprintf("  Total coding bases: %d\n", sum(m$width)))
}

# ============================================================
# EXPECTED OUTPUT FILES:
#   CYP2E1_protein2genome_ENSP00000349420.csv  — canonical protein
#   CYP2E1_protein2genome_ENSP00000374303.csv  — shorter isoform
#   CYP2E1_protein2genome_ALL_proteins.csv     — combined
#
# EXPECTED CONTENT (per file):
#   protein_id, transcript_id, gene_name,
#   protein_aa_start, protein_aa_end,
#   seqnames (chr10), start, end, width, strand (-)
#
# NOTE:
#   CYP2E1 is on the MINUS strand of chromosome 10.
#   proteinToGenome() accounts for this and returns
#   correct genomic coordinates for each exonic segment.
# ============================================================

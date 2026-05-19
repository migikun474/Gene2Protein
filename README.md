README
Gene: CYP2E1
Ensembl ID: ENSG00000130649
Tissue: Liver
Function: Breaks down toxins and alcohol
Chromosome: 10 (minus strand) | Genome build: GRCh38 / hg38

File Overview
FileTaskDescriptionTask1_IGV_CYP2E1_Instructions.txtTask 1IGV setup & navigation guideTask2_List_Transcripts_CYP2E1.RTask 2R script — list all transcriptsTask3_List_Proteins_CYP2E1.RTask 3R script — list all proteinsTask4_Protein2Genome_Map_CYP2E1.RTask 4R script — protein-to-genome maps

Requirements
Software

R (version ≥ 4.0)
IGV — download from https://software.broadinstitute.org/software/igv/download (install locally, NOT on Google Drive)

R Packages
All scripts install their own dependencies automatically on first run. The core packages used are:
rBiocManager::install(c(
  "ensembldb",
  "EnsDb.Hsapiens.v86",   # Ensembl v86, GRCh38 human annotation
  "AnnotationFilter",
  "IRanges",
  "GenomicRanges"
))

How to Run
Run the scripts in order (Tasks 2 → 3 → 4). Each script is self-contained and can also be run independently.
bashRscript Task2_List_Transcripts_CYP2E1.R
Rscript Task3_List_Proteins_CYP2E1.R
Rscript Task4_Protein2Genome_Map_CYP2E1.R
Or open each .R file in RStudio and run line by line.

Task Details
Task 1 — IGV View
Open Task1_IGV_CYP2E1_Instructions.txt and follow the step-by-step guide to:

Install IGV on your local machine
Load the hg38 human genome
Navigate to CYP2E1 using gene name or coordinates (chr10:133,520,000–133,561,000)
Take a screenshot for your report


Task 2 — List All Transcripts
Script: Task2_List_Transcripts_CYP2E1.R
Uses transcripts() from ensembldb with a GeneIdFilter to retrieve all transcripts for ENSG00000130649.
Output printed to console:

Total transcript count
Per-transcript: ID, name, biotype, start/end, chromosome, strand
Breakdown by biotype (protein_coding, retained_intron, etc.)

Output file generated:

CYP2E1_transcripts.csv

Expected transcripts (~5 total):
Transcript IDBiotypeENST00000360332protein_codingENST00000389728protein_codingENST00000495786retained_intronENST00000493705processed_transcriptENST00000464098retained_intron

Task 3 — List All Proteins
Script: Task3_List_Proteins_CYP2E1.R
Uses proteins() from ensembldb to retrieve only the protein products (non-coding transcripts are automatically excluded).
Output printed to console:

Total protein count
Protein ID ↔ Transcript ID mapping
Amino acid length of each protein

Output files generated:

CYP2E1_proteins.csv
CYP2E1_protein_sequences.fasta

Expected proteins (~2 total):
Protein IDTranscript IDLengthENSP00000349420ENST00000360332493 aa (canonical)ENSP00000374303ENST00000389728shorter isoform

Task 4 — Protein-to-Genome Map
Script: Task4_Protein2Genome_Map_CYP2E1.R
Uses proteinToGenome() from ensembldb to map each amino acid position of every CYP2E1 protein back to its genomic coordinates on GRCh38. CYP2E1 is on the minus strand of chromosome 10 — this is handled automatically.
Output files generated (one per protein + one combined):

CYP2E1_protein2genome_ENSP00000349420.csv
CYP2E1_protein2genome_ENSP00000374303.csv
CYP2E1_protein2genome_ALL_proteins.csv

Columns in each map file:
ColumnDescriptionprotein_idEnsembl protein IDtranscript_idLinked transcript IDgene_nameCYP2E1protein_aa_startAmino acid start positionprotein_aa_endAmino acid end positionseqnamesChromosome (chr10)startGenomic start coordinate (GRCh38)endGenomic end coordinate (GRCh38)widthLength of genomic block in basesstrandStrand (− for CYP2E1)

Key Notes

CYP2E1 sits on the minus (−) strand of chromosome 10. IGV will display transcription running right to left.
Only protein_coding transcripts produce proteins. Retained intron and processed transcript biotypes are non-coding and will not appear in Task 3 or Task 4 output.
All scripts use Ensembl v86 (GRCh38) via EnsDb.Hsapiens.v86. If your course uses a different Ensembl version, swap in the appropriate EnsDb package.
proteinToGenome() maps each exonic coding block separately — the number of rows in the map file equals the number of coding exons in that transcript.


Reference

ensembldb vignette (coordinate mapping): https://www.bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/coordinate-mapping.html
Ensembl gene page: https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000130649
IGV download: https://software.broadinstitute.org/software/igv/download

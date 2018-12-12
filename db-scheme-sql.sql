CREATE TABLE cell_line
(
id INT COMMENT 'Primary key.',
cl_id CHAR UNIQUE  COMMENT 'Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
cell_line CHAR COMMENT 'Human readable cell line identifier (could be a catalog number or a different standard cell line term).',
disease CHAR COMMENT 'Disease associated with cell line sample. Set to NORMAL if no disease presented in sample. Disease terms not derived from controlled vocabulary at present.',
tissue CHAR COMMENT 'Source tissue of cell line sample.',
PRIMARY KEY (id)
);

CREATE TABLE exome_mutation_annotation
(
id INT,
exome_id INT COMMENT 'Foreign key. Unique key for each exome sequencing experiment (per cell line per condition).',
HUGO_symbol VARCHAR(255) COMMENT 'HUGO symbol for the gene (HUGO symbols are always in all caps). Unknown is used for regions that do not correspond to a gene.',
Chromosome VARCHAR(255) COMMENT 'The affected chromosome (chr1).',
Start_position INT COMMENT 'Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate.',
End_position INT COMMENT 'Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate.',
Reference_Allele VARCHAR(255) COMMENT 'The plus strand reference allele at this position. Includes the deleted sequence for a deletion or - for an insertion.',
Variant_Classification VARCHAR(255) COMMENT 'Translational effect of variant allele.',
Variant_Type CHAR COMMENT 'Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more (SNP, DNP, TNP, ON',
Mutation_Status CHAR COMMENT 'An assessment of the mutation as somatic, germline, LOH, post transcriptional modification, unknown, or none. The values allowed in this field are constrained by the value in the Validation_Status field.',
Exon_Number INT COMMENT 'The exon number (out of total number).',
Transcript_ID INT COMMENT 'Ensembl ID of the transcript affected by the variant.',
Protein_Position CHAR COMMENT 'Relative position of affected amino acid in protein. A - symbol is displayed as the numerator if the variant does not appear in coding sequence.',
Codons CHAR COMMENT 'The alternative codons with the variant base in upper case',
BIOTYPE CHAR COMMENT 'Biotype of transcript.',
PRIMARY KEY (id)
);

CREATE TABLE exome_sequencing_files
(
id INT COMMENT 'Primary key',
exome_id INT COMMENT 'Foreign key. Unique key for each exome sequencing experiment (per cell line per condition).',
bam_file CHAR COMMENT 'Exome mapping results file.',
bai_file CHAR COMMENT 'Exome mapping index file.',
raw_fastqc_files CHAR COMMENT 'Raw reads and FastQC reports file.',
trim_fastqc_files CHAR COMMENT 'Trimmed reads and FastQC reports file.',
nt_sequences_fa_file CHAR COMMENT 'Nucleotide sequences file.',
coding_sequences_fa_file CHAR COMMENT 'Coding sequences file.',
peptide_translation_aa_file CHAR COMMENT 'Peptide translations file.',
gene_accession_coordinates_bed_file CHAR COMMENT 'Gene accession and coordinates file containing read level coverage and exon accessions and coordinates.',
variant_analysis_calls_vcf CHAR COMMENT 'Variant analysis results file (per cell line).',
copy_number_variation_vcf CHAR COMMENT 'Copy number analysis results file (per cell line).',
PRIMARY KEY (id)
);

CREATE TABLE exome_structural_variants
(
id INT,
exome_id INT COMMENT 'Foreign key. Unique key for each exome sequencing experiment (per cell line per condition).',
cluster_id CHAR,
left_chromosome INT,
left_breakpoint CHAR,
right_chromosome INT,
right_breakpoint CHAR,
num_prs INT,
localization FLOAT,
type VARCHAR(255),
PRIMARY KEY (id)
);

CREATE TABLE exome_variant_analysis_annotation
(
id INT,
exome_id INT COMMENT 'Foreign key. Unique key for each exome sequencing experiment (per cell line per condition).',
dbSNP VARCHAR(255) COMMENT 'The rs-IDs from the   dbSNP database, \"novel\" if not found in any database used, or null if there is no dbSNP record, but it is found in other databases',
chromosome VARCHAR(255) COMMENT 'The affected chromosome (chr1).',
position INT,
ref_base VARCHAR(255),
sample_genotype VARCHAR(255),
sample_alleles CHAR,
alleles_dbSNP CHAR,
accession INT,
function_GVS INT,
function_dbSNP INT,
rs_ID INT,
protein_position INT,
cDNA_position CHAR,
grantham_score CHAR,
score_phast_cons CHAR,
cons_score_GERP CHAR,
dbSNP_validation CHAR,
repeat_masker CHAR,
tandem_repeat CHAR,
kegg_pathway CHAR,
cpg_islands CHAR,
amino_acids CHAR,
poly_phen FLOAT,
african_hap_map_freq CHAR,
european_hap_map_freq CHAR,
asian_hap_map_freq CHAR,
has_genotype CHAR,
distance_to_splice CHAR,
micro_RNAs CHAR,
tfbs CHAR,
genome_esp CHAR,
ppi CHAR,
protein_sequence CHAR,
PRIMARY KEY (id)
);

CREATE TABLE miRNA_mapping_annotation
(
id INT,
miRNA_mapping_id INT,
miRDeep2_score FLOAT,
provisional_id CHAR,
true_positive_prob CHAR,
rfam_alert CHAR,
total_read_count INT,
mature_read_count INT,
loop_read_count INT,
star_read_count INT,
significant_randfold_p_value FLOAT,
miRBase_miRNA VARCHAR(255),
same_seed_example_miRBase_miRNA VARCHAR(255),
UCSC_browser VARCHAR(255),
NCBI_blastn VARCHAR(255),
consensus_mature_seq CHAR,
consensus_star_seq CHAR,
consensus_precursor_seq CHAR,
precursor_coordinate CHAR,
novel BOOLEAN,
mature BOOLEAN,
PRIMARY KEY (id)
);

CREATE TABLE miRNA_mapping_summary
(
id INT,
miRNA_mapping_id INT,
miRDeep2_score FLOAT,
novel_miRNAs_reported INT,
novel_miRNAs_false_positive_est CHAR,
novel_miRNAs_true_positive_est CHAR,
known_miRNAs_total INT,
known_miRNAs_data INT,
known_miRNA_detected INT,
SNR_est FLOAT,
excision_gearing INT,
PRIMARY KEY (id)
);

CREATE TABLE mRNA_gene_transcript_annotation
(
id INT,
mRNA_gene_transcript_map_id INT COMMENT 'Foreign key. Unique ID mapping mRNA abundance analysis results to raw files from sequence experiments per cell line and condition.',
gene_id VARCHAR(255) COMMENT 'Ensemble gene ID.',
transcript_ids CHAR COMMENT 'Expressed transcript ensemble IDs.',
length FLOAT,
effective_length FLOAT,
expected_count INT,
TPM FLOAT COMMENT 'Normalized transcript abundance - Transcripts Per Kilobase Million.',
FPKM FLOAT COMMENT 'Normalized transcript abundance - Fragments Per Kilobase Million.',
PRIMARY KEY (id)
);

CREATE TABLE mRNAseq_variant_call_annotation
(
id INT,
mRNAseq_id INT COMMENT 'Foreign key. Unique key for each mRNA sequencing experiment (per cell line per condition).',
HUGO_symbol VARCHAR(255) COMMENT 'HUGO symbol for the gene (HUGO symbols are always in all caps). Unknown is used for regions that do not correspond to a gene.',
Chromosome VARCHAR(255) COMMENT 'The affected chromosome (chr1).',
Start_position INT COMMENT 'Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate.',
End_position INT COMMENT 'Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate.',
Reference_Allele VARCHAR(255) COMMENT 'The plus strand reference allele at this position. Includes the deleted sequence for a deletion or - for an insertion.',
Variant_Classification VARCHAR(255) COMMENT 'Translational effect of variant allele.',
Variant_Type CHAR COMMENT 'Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more (SNP, DNP, TNP, ON',
Mutation_Status CHAR COMMENT 'An assessment of the mutation as somatic, germline, LOH, post transcriptional modification, unknown, or none. The values allowed in this field are constrained by the value in the Validation_Status field.',
Exon_Number INT COMMENT 'The exon number (out of total number).',
Transcript_ID INT COMMENT 'Ensembl ID of the transcript affected by the variant.',
Protein_Position CHAR COMMENT 'Relative position of affected amino acid in protein. A - symbol is displayed as the numerator if the variant does not appear in coding sequence.',
Codons CHAR COMMENT 'The alternative codons with the variant base in upper case',
BIOTYPE CHAR COMMENT 'Biotype of transcript.',
PRIMARY KEY (id)
);

CREATE TABLE proliferation_measurement_timing
(
id INT,
cl_id CHAR UNIQUE  COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
p1t1 FLOAT COMMENT 'Plate 1 - sample at time 1 [hours]',
p1t2 FLOAT COMMENT 'Plate 1 - sample at time 2  [hours]',
p1t3 FLOAT COMMENT 'Plate 1 - sample at time 3 [hours]',
p2t1 FLOAT COMMENT 'Plate 2 - sample at time 1 [hours]',
p2t2 FLOAT COMMENT 'Plate 2 - sample at time 2 [hours]',
p2t3 FLOAT COMMENT 'Plate 2 - sample at time 3 [hours]',
PRIMARY KEY (id)
);

CREATE TABLE proteomics_tmt
(
id INT COMMENT 'Primary key.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
accessions CHAR COMMENT 'Indicates that this preplicate had NOT been spiked with UPS2.',
score FLOAT COMMENT 'Replicate number. Technical replicate identifier. Also called preplicate for preparation replicate. Options: 1; 2; N/A. I.e. two internal replicates were prepared to make the TMT-10-plex mixture. 1 - indicates that this is the first of 3 preplicates; 2 - ',
sequence CHAR COMMENT 'The name of the sequence from the protein database.',
modifications CHAR COMMENT 'The score assigned by PeptideProphet. This score represents the probability that the peptide identification is correct. A higher score indicates a better match.',
x_corr FLOAT COMMENT 'The sequence of the peptide match. The previous and next amino acids in the database sequence are printed before/after the identified peptide, separated by periods.',
description TEXT COMMENT 'Short description of the protein nature and function.',
modifications_all_sites CHAR COMMENT 'The score assigned by PeptideProphet. This score represents the probability that the peptide identification is correct. A higher score indicates a better match.',
sample INT COMMENT 'Replicate number. Technical replicate identifier. Also called preplicate for preparation replicate. Options: 1; 2; N/A. I.e. two internal replicates were prepared to make the TMT-10-plex mixture. 1 - indicates that this is the first of 3 preplicates; 2 - ',
abundance_ratio_F1_127N_F1_126 FLOAT,
abundance_ratio_F1_127C_F1_126 FLOAT,
abundance_ratio_F1_128N_F1_126 FLOAT,
abundance_ratio_F1_128C_F1_126 FLOAT,
abundance_ratio_F1_129N_F1_126 FLOAT,
abundance_ratio_F1_129C_F1_126 FLOAT,
abundance_ratio_F1_130N_F1_126 FLOAT,
abundance_ratio_F1_130C_F1_126 FLOAT,
abundance_ratio_F1_131_F1_126 FLOAT,
abundance_ratio_log2_F1_127N_F1_126 FLOAT,
abundance_ratio_log2_F1_127C_F1_126 FLOAT,
abundance_ratio_log2_F1_128N_F1_126 FLOAT,
abundance_ratio_log2_F1_128C_F1_126 FLOAT,
abundance_ratio_log2_F1_129N_F1_126 FLOAT,
abundance_ratio_log2_F1_129C_F1_126 FLOAT,
abundance_ratio_log2_F1_130N_F1_126 FLOAT,
abundance_ratio_log2_F1_130C_F1_126 FLOAT,
abundance_ratio_log2_F1_131_F1_126 FLOAT,
PRIMARY KEY (id)
);

CREATE TABLE substrate
(
id INT COMMENT 'Primary key.',
sub_id CHAR UNIQUE  COMMENT 'Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
substrate CHAR COMMENT 'Human readable string identifying each gel condition, as used across experiment data labels.',
stiffness CHAR COMMENT 'Gel matrix stiffness measured in Pascals.',
integrin_ligand CHAR COMMENT 'Integrin ligand type, e.g. Colagen I (col), Fibronectin (Fn).',
matrix CHAR COMMENT 'Hydrogel matrix type, e.g. polyacrylamide (PA), hyaluronic acid (HA)',
PRIMARY KEY (id)
);

CREATE TABLE stiffness
(
id INT COMMENT 'Primary key.',
atomic_force_measurement_id CHAR COMMENT 'Concatenation of each cell condition, cell line and trial combination - unique identifier for each atomic force measurement measurement.',
calc_ramp_ex_nm FLOAT COMMENT 'Movement of the AFM piezo during extension phase (nanometers).',
calc_ramp_rt_nm FLOAT COMMENT 'Movement of the AFM piezo during retraction phase (nanometers).',
defl_nm_ex FLOAT COMMENT 'Deflection of the cantilever during extension phase (nanometers).',
defl_nm_rt FLOAT COMMENT 'Deflection of the cantilever during retraction phase (nanometers).',
time_s_ex FLOAT COMMENT 'Duration of extension phase (seconds).',
time_s_rt FLOAT COMMENT 'Duration of retraction phase (seconds).',
PRIMARY KEY (id)
);

CREATE TABLE miRNA_differential_expression
(
id INT,
sub_id CHAR,
tissue_1 CHAR COMMENT 'Source tissue 1 (out of 2)  for differential expression analysis (note all cell lines from that tissue are included in the analysis).',
tissue_2 CHAR COMMENT 'Source tissue 2 (out of 2)  for differential expression analysis (note all cell lines from that tissue are included in the analysis).',
miRNA_accession_id CHAR COMMENT 'miRNA miRBase accession ID.',
base_mean FLOAT COMMENT 'A measure of read abundance, i.e. expression level (analogous to  \'mean of normalized counts\' in DESeq2, and \'concentration\' or \'counts-per-million\' in edgeR).',
log2_fold_change FLOAT COMMENT 'A log2 fold change (in expression level) for each transcript.',
lfcSE FLOAT COMMENT 'Log fold change standard error.',
stat FLOAT,
p_value FLOAT COMMENT 'p-value (not adjusted).',
p_adj FLOAT COMMENT 'Adjusted p-value.',
NCBI_genome_build CHAR COMMENT 'Reference human genome build.',
PRIMARY KEY (id)
);

CREATE TABLE proteomics_phospho
(
id INT COMMENT 'Primary key.',
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
accession CHAR COMMENT 'Indicates that this preplicate had NOT been spiked with UPS2.',
score FLOAT COMMENT 'Replicate number. Technical replicate identifier. Also called preplicate for preparation replicate. Options: 1; 2; N/A. I.e. two internal replicates were prepared to make the TMT-10-plex mixture. 1 - indicates that this is the first of 3 preplicates; 2 - ',
coverage FLOAT COMMENT 'Probability score assigned to the protein group by ProteinProphet.',
sequence CHAR COMMENT 'The name of the sequence from the protein database.',
num_peptides INT COMMENT 'The number of filtered peptides in the run that were matched to this sequence.',
unique_peptides INT COMMENT 'The number of unique filtered peptides in the run that were matched to this sequence.',
modifications CHAR COMMENT 'The score assigned by PeptideProphet. This score represents the probability that the peptide identification is correct. A higher score indicates a better match.',
x_corr FLOAT COMMENT 'The sequence of the peptide match. The previous and next amino acids in the database sequence are printed before/after the identified peptide, separated by periods.',
description TEXT COMMENT 'Short description of the protein nature and function.',
replicate INT COMMENT 'Replicate number. Technical replicate identifier. Also called preplicate for preparation replicate. Options: 1; 2; N/A. I.e. two internal replicates were prepared to make the TMT-10-plex mixture. 1 - indicates that this is the first of 3 preplicates; 2 - ',
PRIMARY KEY (id)
);

CREATE TABLE contractility
(
id INT COMMENT 'Primary key.',
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
traction_force FLOAT COMMENT 'Traction force is used as a proxy for cell contractility, which is measured as total maximum force (nN) using the traction force microscopy (TFM). All measreuments performed by a Laser scanning confocal microscope (Leica TCS SP5, Wetzlar, Germany).',
nuclear_volume FLOAT COMMENT 'Nuclear volume derived from 3D-image stacks recorder by 633 nm laser lines.',
cell_volume FLOAT COMMENT 'Cell volume derived from 3D-image stacks recorded by 488nm laser lines.',
cell_area FLOAT COMMENT 'Cell area based on identified cell boundary from cell volume calculations.',
PRIMARY KEY (id)
);

CREATE TABLE exome_sequencing_summary
(
id INT,
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
num_exons INT COMMENT 'Number of exons per cell line in the given condition.',
total_size_exons INT COMMENT 'Total size of exons per sequenced cell line in the given condition.',
longest_exon INT COMMENT 'Longest exon length per sequenced cell line in the given condition.',
shortest_exon INT COMMENT 'Shortest exon length per sequenced cell line in the given condition.',
mean_exon_size INT COMMENT 'Average exon size per sequenced cell line in the given condition.',
median_exon_size INT COMMENT 'Median exon size per sequenced cell line in the given condition.',
N50_exon_length INT COMMENT 'N50 contig length (weighted mean of lengths - greater length to longer reads) per sequenced cell line in the given condition.',
L50_exon_count INT COMMENT 'L50 number of contigs whose summed length is N50, per sequenced cell line in the given condition.',
exon_percent_A FLOAT COMMENT 'Percentage of nucleotide A across exons.',
exon_percent_C FLOAT COMMENT 'Percentage of nucleotide C across exons.',
exon_percent_G FLOAT COMMENT 'Percentage of nucleotide G across exons.',
exon_percent_T FLOAT COMMENT 'Percentage of nucleotide T across exons.',
exon_percent_N FLOAT COMMENT 'Percentage of nucleotide N across exons.',
exon_percent_GC FLOAT COMMENT 'Percentage of nucleotides that are either G or C across exons.',
exon_percent_other_nt FLOAT COMMENT 'Percentage of unidentified/other nucleotides.',
exon_num_other_nt FLOAT COMMENT 'Number of exon non-ACGTN nucleotides.',
mean_exon_completeness FLOAT COMMENT 'Average exon completeness per sequenced cell line in the given condition.',
mean_exon_coverage_per_base FLOAT COMMENT 'Average exon coverage per base, per sequenced cell line in the given condition.',
NCBI_genome_build CHAR COMMENT 'Reference human genome build.',
PRIMARY KEY (id)
);

CREATE TABLE miRNA_differential_expression_files
(
id INT,
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR,
miRNA_differential_expression_id INT COMMENT 'Foreign key. Unique id of mRNA differential expression analysis. Unique for each tissue and condition pair (a tissue may encompass multiple cell lines).',
raw_fastq_files CHAR COMMENT 'The raw Illumina sequence reads converted to .fastq files (individual files are comma-separated).',
trim_reads_fq_files CHAR COMMENT 'Trimmed reads .fq files (individual files are comma-separated).',
sorted_25nt_bam_file CHAR COMMENT 'Sequence reads mapped to the reference genome (HG19/37) using Bowtie2 (bam file).',
sorted_25nt_bai_file CHAR COMMENT 'Sequence reads mapped to the reference genome (HG19/37) using Bowtie2 index file (bai).',
unmapped_25nt_fq_file CHAR COMMENT 'Unmapped reads .fq file.',
miRanalyzer_files CHAR COMMENT 'Outputs of miRanalyzer (individual files are comma-separated).',
NCBI_genome_build CHAR,
PRIMARY KEY (id)
);

CREATE TABLE miRNA_mapping
(
id INT,
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
raw_fastq_file CHAR,
NCBI_genome_build CHAR COMMENT 'Reference human genome build.',
PRIMARY KEY (id)
);

CREATE TABLE morphology
(
id INT COMMENT 'Primary key.',
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
aspect_ratio FLOAT COMMENT 'Single cell aspect ratio measurements using ImageJ photo stacks.',
circularity FLOAT COMMENT 'Single cell circularity measurements using ImageJ photo stacks.',
area FLOAT COMMENT 'Single cell area measurements using ImageJ photo stacks (micrometers squared).',
PRIMARY KEY (id)
);

CREATE TABLE motility
(
id INT COMMENT 'Primary key',
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
avg_end_to_end_dist FLOAT COMMENT 'Average straight line (end-to-end) distance traversed by a cell of a given cell line in a given culture condition. Measured in micrometers.',
avg_total_dist FLOAT COMMENT 'Average total distance traversed by a cell of a given cell line in a given culture condition. Measured in micrometers.',
avg_speed FLOAT COMMENT 'Average velocity of a given cell line in a given culture condition. Measured in micrometers/min.',
end_to_end_dist_st_dev FLOAT COMMENT 'Standard deviation of straight line (end-to-end) distance distribution for a given cell line in a given culture condition.',
total_dist_st_dev FLOAT COMMENT 'Standard deviation of total distance distribution for a given cell line in a given culture condition.',
speed_st_dev FLOAT COMMENT 'Standard deviation of velocity distribution of a given cell line in a given culture condition. ',
end_to_end_dist_st_err FLOAT COMMENT 'Standard error of straight line (end-to-end) distance distribution for a given cell line in a given culture condition.',
total_dist_st_err FLOAT COMMENT 'Standard error of total distance distribution for a given cell line in a given culture condition.',
speed_st_err FLOAT COMMENT 'Standard error of velocity distribution of a given cell line in a given culture condition. ',
PRIMARY KEY (id)
);

CREATE TABLE motility_measurements
(
id INT COMMENT 'Primary key
',
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
distance FLOAT COMMENT 'Distance value reported by ImageJ - not used in downstream summary motility calculations for now. Multiple measurements per cell line and culture condition.',
pixel_value FLOAT COMMENT 'Pixel value reported by ImageJ - not used in downstream summary motility calculations for now. Multiple measurements per cell line and culture condition. Please refer to ImageJ\'s manual for this value\'s description.',
slice INT COMMENT 'Frame number for the position indicated by pixel (x, y) coordinates; time between frames is 5 minutes for all but cell line DU145, plate 1 = 6.5 min; and cell line Caov3 = 1 hour.',
track INT COMMENT 'Unique ID of a cell being tracked (i.e. indicates which cell is tracked in this measurement).',
velocity FLOAT COMMENT 'Velocity value reported by ImageJ - not used in downstream summary motility calculations for now. Multiple measurements per cell line and culture condition. Please refer to ImageJ\'s manual for this value\'s description.',
x FLOAT COMMENT 'x-coordinate (in pixels) of this cell position in this frame; each pixel corresponds to 1.5 micrometers.',
y FLOAT COMMENT 'y-coordinate (in pixels) of this cell position in this frame; each pixel corresponds to 1.5 micrometers.',
PRIMARY KEY (id)
);

CREATE TABLE mRNA_gene_transcript_map
(
id INT,
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
raw_fastq_files CHAR COMMENT 'Reads and FastQC reports files (individual files are comma-separated).',
NCBI_genome_build CHAR COMMENT 'Reference human genome build.',
sorted_STAR_bam_file CHAR COMMENT 'Mapped reads (sorted STAR) bam file.',
sorted_STAR_bai_file CHAR COMMENT 'Mapped reads (sorted STAR) bai file.',
PRIMARY KEY (id)
);

CREATE TABLE mRNA_sequencing_files
(
id INT,
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
raw_fastqc_files CHAR,
trim_reads_fq_files CHAR,
accepted_hits_bam_file CHAR,
accepted_hits_bai_file CHAR,
unmapped_bam_file CHAR,
deletions_bed_file CHAR,
junctions_bed_file CHAR,
insertions_bed_file CHAR,
transcripts_transdecoder_files CHAR,
cufflinks_files CHAR,
NCBI_genome_build CHAR,
variant_analysis_calls_vcf_file CHAR,
transcript_sequences_file CHAR,
PRIMARY KEY (id)
);

CREATE TABLE proliferation
(
id INT COMMENT 'Primary key.',
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
plate INT COMMENT 'Replicate plate number. All measurements were simultaneously performed in two plates replicating the experiment.',
trial INT COMMENT 'Time point in which measurement was obtained (i.e. three possible timepoints). Trials at timepoints 2 and 3 only recorded total number of cells in culture (w/o distinction between touching and single cells). Please refer to proliferation plate timing for ',
single INT COMMENT 'The number of single cells (cells not touching other cells) were manually counted for each frame that had clear enough cells to count for each condition, replicate, and cell line for time point 1 (see proliferation plate timing for time point description)',
touching INT COMMENT 'The number of touching cells (cells not isolated from other cells) were manually counted for each frame that had clear enough cells to count for each condition, replicate, and cell line for time point 1 (see proliferation plate timing for time point descr',
frame_num INT COMMENT 'Number of ImageJ frames counted for this cell line and gel culture condition.',
PRIMARY KEY (id)
);

CREATE TABLE proteomics_iBAQ
(
id INT COMMENT 'Primary key.',
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
UPS2 BOOLEAN COMMENT 'Indicates that this preplicate had NOT been spiked with UPS2.',
replicate INT COMMENT 'Replicate number. Technical replicate identifier. Also called preplicate for preparation replicate. Options: 1; 2; N/A. I.e. two internal replicates were prepared to make the TMT-10-plex mixture. 1 - indicates that this is the first of 3 preplicates; 2 - ',
prob FLOAT COMMENT 'Probability score assigned to the protein group by ProteinProphet.',
protein CHAR COMMENT 'The name of the sequence from the protein database.',
total_filtered_peptides INT COMMENT 'The number of filtered peptides in the run that were matched to this sequence.',
unique_filtered_peptides INT COMMENT 'The number of unique filtered peptides in the run that were matched to this sequence.',
pep_prophet FLOAT COMMENT 'The score assigned by PeptideProphet. This score represents the probability that the peptide identification is correct. A higher score indicates a better match.',
peptide CHAR COMMENT 'The sequence of the peptide match. The previous and next amino acids in the database sequence are printed before/after the identified peptide, separated by periods.',
description TEXT COMMENT 'Short description of the protein nature and function.',
PRIMARY KEY (id)
);

CREATE TABLE stiffness_measurements
(
id INT COMMENT 'Primary key.',
sub_id CHAR COMMENT 'Foreign key. Normalized string unique for each culture condition; concatenation of gel matrix type, stiffness (in Pascals), integrin ligands type; uses glass for uncoated glass culture condition.',
cl_id CHAR COMMENT 'Foreign key. Cell line id: normalized string uniquely identifying each cell line; human readable but not equivalent to standard catalog cell line names.',
atomic_force_measurement_id CHAR UNIQUE  COMMENT 'Concatenation of each cell condition, cell line and trial combination - unique identifier for each atomic force measurement measurement.',
date DATE COMMENT 'Date the measurement was performed.',
young_modulus FLOAT COMMENT 'Elastic modulus calculation for cell membrane stiffness (i.e. Young modulus) based on the movement of the AFM piezo and deflection of the cantilever during extension phase. Please refer to https://www.synapse.org/#!Synapse:syn7248585 for more details o',
spring_constant FLOAT COMMENT 'Spring constant calculation of the cantilever used in Young modulus measurement. The resonant frequency of the cantilever was first measured by applying a range of vibrational frequencies to the cantilever(10-30,000 Hz) and determining the frequency wh',
PRIMARY KEY (id)
);
library('dplyr')
library('tidyr')
library('biomaRt')

BASE = '/u/project/zarlab/nlapier2/geuvadis_mirror/'
#ensembl = useMart("ensembl")
#listEnsembl(GRCh=37)
grch37 = useEnsembl(biomart="ensembl",GRCh=37)
#listDatasets(grch37)
ds = useDataset('hsapiens_gene_ensembl', mart = grch37)
#listAttributes(ds)

source('./fastqtl_utils.R')

# method = 'fastqtl'
method = 'qtltools'

# get these fields from the biomart used above
gene_tss = getBM(
  c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position',
    'strand'),
  mart = ds)
# a few renaming and transforming operations on the dataframe
gene_tss <- dplyr::rename(gene_tss, ens_gene = ensembl_gene_id, chr = chromosome_name)
gene_tss <- dplyr::mutate(gene_tss, tss = ifelse(strand == 1, start_position, end_position))
gene_tss = dplyr::mutate(gene_tss, strand = ifelse(strand == 1, '+', '-'))

if (method == 'qtltools') {
  tss_info = dplyr::filter(gene_tss, grepl('^[0-9]+', chr))
  tss_info = dplyr::mutate(tss_info, pid = ens_gene, gid = ens_gene, start_position = start_position - 1)
  tss_info = dplyr::select(tss_info, Chr = chr, start = start_position,
    end = end_position, pid, gid, strand)
  tss_info = dplyr::arrange(tss_info, Chr, start, end)
} else if (method == 'fastqtl') {
  tss_info <- dplyr::filter(gene_tss, grepl('^[0-9]+$', chr))
  tss_info <- mutate(tss_info, start = tss - 1, end = tss)
  tss_info <- dplyr::select(tss_info, Chr = chr, start, end, ID = ens_gene)
  tss_info = dplyr::arrange(tss_info, Chr, start, end)
} else {
  stop('unrecognized method')
}



# # put in a format for fastqtl
# # NOTE: 0 based indexing
# tss_info <- dplyr::filter(gene_tss, grepl('^[0-9]+$', chr))
# tss_info <- mutate(tss_info, start = tss - 1, end = tss, pid = ens_gene, gid = ens_gene)
# tss_info <- dplyr::select(tss_info, Chr = chr, start, end, pid, gid, strand)
# tss_info = dplyr::arrange(tss_info, Chr, start, end)

########################################################################
# first process kallisto
########################################################################

# counts_TPM.rds scaledTPM object
# counts.rds the raw counts object
ti_scaled = readRDS(paste0(BASE, 'counts_TPM.rds'))
ti_raw = readRDS(paste0(BASE, 'counts.rds'))


# conclusion: ti_scaled$counts is different, while abundance is not
all.equal(ti_scaled$counts, ti_raw$counts)
all.equal(ti_scaled$abundance, ti_raw$abundance)

all.equal(rownames(ti_scaled$counts), rownames(ti_raw$counts))

raw_counts = as.matrix(ti_raw$counts)
tmp = sub('\\..*$', '', rownames(raw_counts))
rownames(raw_counts) = tmp
scaled_counts = as.matrix(ti_scaled$counts)
rownames(scaled_counts) = tmp

# now filter on the raw counts
mean_minimum = 0.50
min_counts = 5

pass_count_filter = apply(raw_counts, 1,
  function(x) mean(x >= min_counts) >= mean_minimum)
pass_count_filter = names(pass_count_filter[pass_count_filter])

# now normalize the scaled counts
scaled_filter = scaled_counts[pass_count_filter, ]
sf = DESeq2::estimateSizeFactorsForMatrix(scaled_filter)
norm_counts = t(t(scaled_filter) / sf)

filter_low_variance <- function(m) {
  q_drop <- 0.05
  all_sd <- apply(m, 1, sd)
  which_include <- all_sd >= quantile(all_sd, probs = q_drop)
  m[which_include, ]
}

tmp = expression_to_bed(norm_counts,
  tss_info,
  paste0('beds/gene_level_norm/', method, '/kallisto_scaled_tpm'),
  n_pc = 20,
  normalize = TRUE,
  tool = method)

tmp = expression_to_bed(norm_counts,
  tss_info,
  paste0('beds/gene_level_tpm/', method, '/kallisto_scaled_tpm'),
  n_pc = 20,
  normalize = FALSE,
  tool = method)

# now compute the shared covariates
metadata = read.table(paste0(BASE, 'metadata/clean.tsv'), header = FALSE)
colnames(metadata) <- c('sample', 'left', 'right', 'ena_sample', 'assay_name', 'population',
  'laboratory')
metadata <- filter(metadata, population == 'YRI')
metadata = mutate(metadata, laboratory = paste0('lab', laboratory))
covariates = dplyr::select(metadata, sample, laboratory)
covariates = pivot_wider(covariates, names_from = sample, values_from = laboratory)
covariates = data.frame(id = 'laboratory', covariates)
dir.create('beds/gene_level_norm/qtl_shared', showWarnings = FALSE)
dir.create('beds/gene_level_tpm/qtl_shared', showWarnings = FALSE)
write.table(covariates, file = 'beds/gene_level_norm/qtl_shared/covariates.txt', quote = FALSE,
  sep = '\t', row.names = FALSE)
write.table(covariates, file = 'beds/gene_level_tpm/qtl_shared/covariates.txt', quote = FALSE,
  sep = '\t', row.names = FALSE)




# now process transcript-level kallisto counts
# get these fields from the biomart used above
trans_tss = getBM(
  c('ensembl_transcript_id_version', 'chromosome_name', 'transcript_start', 'transcript_end', 'strand'),
  mart = ds)
# a few renaming and transforming operations on the dataframe
trans_tss <- dplyr::rename(trans_tss, ens_trans = ensembl_transcript_id_version, chr = chromosome_name)
trans_tss <- dplyr::mutate(trans_tss, tss = ifelse(strand == 1, transcript_start, transcript_end))
trans_tss = dplyr::mutate(trans_tss, strand = ifelse(strand == 1, '+', '-'))

if (method == 'qtltools') {
  tss_info = dplyr::filter(trans_tss, grepl('^[0-9]+', chr))
  tss_info = dplyr::mutate(tss_info, pid = ens_trans, gid = ens_trans, transcript_start = transcript_start - 1)
  tss_info = dplyr::select(tss_info, Chr = chr, start = transcript_start,
    end = transcript_end, pid, gid, strand)
  tss_info = dplyr::arrange(tss_info, Chr, start, end)
} else if (method == 'fastqtl') {
  tss_info <- dplyr::filter(trans_tss, grepl('^[0-9]+$', chr))
  tss_info <- mutate(tss_info, start = tss - 1, end = tss)
  tss_info <- dplyr::select(tss_info, Chr = chr, start, end, ID = ens_trans)
  tss_info = dplyr::arrange(tss_info, Chr, start, end)
} else {
  stop('unrecognized method')
}

# transcript_counts_TPM.rds scaledTPM object
# transcript_counts.rds the raw counts object
ti_scaled = readRDS(paste0(BASE, 'transcript_counts_TPM.rds'))
ti_raw = readRDS(paste0(BASE, 'transcript_counts.rds'))


# conclusion: ti_scaled$counts is different, while abundance is not
all.equal(ti_scaled$counts, ti_raw$counts)
all.equal(ti_scaled$abundance, ti_raw$abundance)

all.equal(rownames(ti_scaled$counts), rownames(ti_raw$counts))

raw_counts = as.matrix(ti_raw$counts)
#tmp = sub('\\..*$', '', rownames(raw_counts))
#rownames(raw_counts) = tmp
scaled_counts = as.matrix(ti_scaled$counts)
#rownames(scaled_counts) = tmp

# now filter on the raw counts
mean_minimum = 0.50
min_counts = 5

pass_count_filter = apply(raw_counts, 1,
  function(x) mean(x >= min_counts) >= mean_minimum)
pass_count_filter = names(pass_count_filter[pass_count_filter])

# now normalize the scaled counts
scaled_filter = scaled_counts[pass_count_filter, ]
sf = DESeq2::estimateSizeFactorsForMatrix(scaled_filter)
norm_counts = t(t(scaled_filter) / sf)

filter_low_variance <- function(m) {
  q_drop <- 0.05
  all_sd <- apply(m, 1, sd)
  which_include <- all_sd >= quantile(all_sd, probs = q_drop)
  m[which_include, ]
}

tmp = expression_to_bed(norm_counts,
  tss_info,
  paste0('beds/transcript_level_norm/', method, '/kallisto_scaled_tpm'),
  n_pc = 20,
  normalize = TRUE,
  tool = method)

tmp = expression_to_bed(norm_counts,
  tss_info,
  paste0('beds/transcript_level_tpm/', method, '/kallisto_scaled_tpm'),
  n_pc = 20,
  normalize = FALSE,
  tool = method)

# now compute the shared covariates
metadata = read.table(paste0(BASE, 'metadata/clean.tsv'), header = FALSE)
colnames(metadata) <- c('sample', 'left', 'right', 'ena_sample', 'assay_name', 'population',
  'laboratory')
metadata <- filter(metadata, population == 'YRI')
metadata = mutate(metadata, laboratory = paste0('lab', laboratory))
covariates = dplyr::select(metadata, sample, laboratory)
covariates = pivot_wider(covariates, names_from = sample, values_from = laboratory)
covariates = data.frame(id = 'laboratory', covariates)
dir.create('beds/transcript_level_norm/qtl_shared', showWarnings = FALSE)
dir.create('beds/transcript_level_tpm/qtl_shared', showWarnings = FALSE)
write.table(covariates, file = 'beds/transcript_level_norm/qtl_shared/covariates.txt', quote = FALSE,
  sep = '\t', row.names = FALSE)
write.table(covariates, file = 'beds/transcript_level_tpm/qtl_shared/covariates.txt', quote = FALSE,
  sep = '\t', row.names = FALSE)





# now process featureCounts
#fc = read.table('./featureCounts.txt', header = TRUE, row.names = 1)
#fc = as.matrix(fc)

#pass_count_filter_fc = apply(fc, 1,
#  function(x) mean(x >= min_counts) >= mean_minimum)
#pass_count_filter_fc = names(pass_count_filter_fc[pass_count_filter_fc])

#fc_filtered = fc[pass_count_filter_fc, ]

#fc_sf = DESeq2::estimateSizeFactorsForMatrix(fc_filtered)
#norm_fc = t(t(fc_filtered) / fc_sf)

#tmp = expression_to_bed(norm_fc,
#  tss_info,
#  paste0(method, '/featureCounts'),
#  n_pc = 20,
#  normalize = TRUE,
#  tool = method)


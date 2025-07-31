library(Seurat)
options(warn = -1)

# get args
args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
species <- args[2]
sample <- args[3]
script_dir <- args[4]

# def function
process_downsampling <- function(AllCounts, mtgenelist, downsampling_type, indir, output_file) {
  results <- data.frame(
    downsample = character(),
    median_gene = numeric(),
    median_umis = numeric(),
    stringsAsFactors = FALSE
  )
  
  downsampling_names <- names(AllCounts[[downsampling_type]]$inex$downsampling)
  
  for (ds_name in downsampling_names) {
    umis <- as.matrix(AllCounts[[downsampling_type]]$inex$downsampling[[ds_name]])
    seu <- CreateSeuratObject(counts = umis)
    mt_genes <- intersect(rownames(seu), as.character(mtgenelist[, 1]))
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mt_genes)
    
    median_gene <- median(seu@meta.data$nFeature_RNA)
    median_umis <- median(seu@meta.data$nCount_RNA)
    
    ds_name_clean <- sub('downsampled_', '', ds_name)
    
    results <- rbind(results, data.frame(
      downsample = ds_name_clean,
      median_gene = median_gene,
      median_umis = median_umis
    ))
  }
  
  write.table(results, output_file, sep='\t', quote=FALSE,row.names=FALSE)
}

# get data
AllCounts <- readRDS(paste0(indir, "/zUMIs_output/expression/", sample, ".dgecounts.rds"))

# get MTgent list
if (species == 'Human') {
  mtgenelist <- read.table(file.path(script_dir, 'human_MTgeneList'))
} else if (species == 'Mouse') {
  mtgenelist <- read.table(file.path(script_dir, 'mouse_MTgeneList'))
}

# process data
for (data_type in c("umicount", "readcount")) {
  umis <- as.matrix(AllCounts[[data_type]]$inex$all)
  seu <- CreateSeuratObject(counts = umis)
  mt_genes <- intersect(rownames(seu), as.character(mtgenelist[, 1]))
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mt_genes)
  
  output_file <- paste0(indir, "/summary/", data_type, "_all.xls")
  write.table(seu@meta.data, output_file, sep='\t', quote=FALSE,col.names=NA)
}

# process downmsampling data
process_downsampling(AllCounts, mtgenelist, "readcount", indir, paste0(indir, "/summary/readcount_saturation.xls"))
process_downsampling(AllCounts, mtgenelist, "umicount", indir, paste0(indir, "/summary/umicount_saturation.xls"))

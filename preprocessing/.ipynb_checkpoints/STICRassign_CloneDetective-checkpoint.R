library(CloneDetective)
library(DropletUtils)
library(scater)
library(data.table)
library(Seurat)

dir <- '/wynton/group/pollen/jding/brainchromatin/HM2D/cellranger-7.2.0/'
fig_dir='/wynton/home/pollenlab/jding/BrainChromatin/macaque/figures/'
dir.create(fig_dir, showWarnings = FALSE)

#for (f in c('HM2D_L1','HM2D_L2','HM2D_L3')){
for (f in c('HM2D_2nd_L1','HM2D_2nd_L2')){
  sce <- read10xCounts(file.path(dir, f,  'outs/filtered_feature_bc_matrix'))
  cell_clone_bcode_dt <- fread(file.path(dir,f,'outs/nextclone/clone_barcodes.csv'))
  cell_clone_bcode_dt$UID = paste0(cell_clone_bcode_dt$CellBarcode, cell_clone_bcode_dt$UMI)
  wl = unique(names(table(cell_clone_bcode_dt$UID) > 1))
  cell_clone_bcode_dt = cell_clone_bcode_dt[cell_clone_bcode_dt$UID %in% wl,]
  cell_by_clone_mat = generate_cell_clone_barcode_matrix(
    cell_clone_bcode_dt = cell_clone_bcode_dt,
    cell_bcode_col = "CellBarcode",
    clone_bcode_col = "CloneBarcode",
    umi_col = "UMI",
    umi_clone_consensus_threshold = 0.7
  )
  cell_by_clone_mat = cell_by_clone_mat[cell_by_clone_mat$n_reads > 1,]
  plt <- draw_treemap(
    cell_by_clone_matrix = cell_by_clone_mat,
    valid_cells_bcodes = colData(sce)$Barcode
    )
  pdf(paste(fig_dir, f,'_Clones_1.pdf',sep=''),width=6, height=6)
  print(plt)
  dev.off()
  cell_by_clone_mat = cell_by_clone_mat[cell_by_clone_mat$n_reads > 2,]
  plt <- draw_treemap(
    cell_by_clone_matrix = cell_by_clone_mat,
    valid_cells_bcodes = colData(sce)$Barcode
    )
  pdf(paste(fig_dir, f,'_Clones_2.pdf',sep=''),width=6, height=6)
  print(plt)
  dev.off()
  cell_by_clone_mat$CellClone = paste0(cell_by_clone_mat$CellBarcode,cell_by_clone_mat$CloneBarcode)
  cell_clone_bcode_dt$CellClone = paste0(cell_clone_bcode_dt$CellBarcode,cell_clone_bcode_dt$CloneBarcode)
  cell_clone_bcode_dt = cell_clone_bcode_dt[cell_clone_bcode_dt$CellClone %in% cell_by_clone_mat$CellClone,]
  sce_with_clone <- assign_and_embed_clones(
    cell_by_gene_mat = sce,
    cell_clone_reads_dt = cell_clone_bcode_dt)
  write.csv(colData(sce_with_clone),file.path(dir, f, 'outs/nextclone/assignments.csv'))
}




.libPaths("/mnt/beegfs/mccarthy/backed_up/general/rlyu/Software/Rlibs/4.1.2")
library(scater)
library(ggplot2)
library(mclust)
library(slingshot)
library(tradeSeq)
library(grDevices)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(Seurat)
library(dplyr)
library(BiocParallel)
library(mgcv) 

args <- (commandArgs(trailingOnly = TRUE))

for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

print(celltype)
print(outrds)
print(out_knots_png)
print(filter_by)
print(use_major_ct)
filter_by <- as.numeric(filter_by)
print(min_lumen)
min_lumen <- as.numeric(min_lumen)

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 4 # use 4 cores
BPPARAM$stop.on.error <- FALSE
lumen_celltype_sce <- readRDS(file="output/savedRDS/revisions2024/lumen_cell_types_sce.all.ct.rds")

lumen_nuclei_RNA_sce <- readRDS(file = "output/savedRDS/revisions2024/lumen_nuclei_RNA_sce_t4_ordered.rds")

lumen_nuclei_RNA_sce_per_ct <- readRDS("output/savedRDS/revisions2024/lumen_nuclei_RNA_per_celltype_sce.rds")
lumen_nuclei_RNA_sce_per_ct$tma <-
  lumen_nuclei_RNA_sce$tma[match(lumen_nuclei_RNA_sce_per_ct$sample_name,
                                 lumen_nuclei_RNA_sce$sample)]

table(lumen_nuclei_RNA_sce_per_ct$final_CT)
lumen_nuclei_RNA_sce_per_ct$Lineage1 <- lumen_nuclei_RNA_sce[,lumen_nuclei_RNA_sce_per_ct$lumen_id]$Lineage1

assays(lumen_celltype_sce)$final_CT_prop[1:4,1:5]
final_CT_prop <- t(assays(lumen_celltype_sce)$final_CT_prop)
final_CT_count <- t(assays(lumen_celltype_sce)$final_CT_count)

final_CT_prop_filter <- final_CT_prop[,colnames(final_CT_prop) %in%
                                          names(colSums(final_CT_prop !=0) >= 100)[colSums(final_CT_prop !=0) >= 100]]
colnames(final_CT_prop_filter)
final_CT_prop <- data.frame(final_CT_prop)
colnames(final_CT_prop)

final_CT_prop_filter <- data.frame(final_CT_prop_filter)
dim(final_CT_prop_filter)
## >=3 in 300 lumens

final_CT_prop_major_celltypes <- colnames(final_CT_count)[colSums(final_CT_count >= 3) >= 300]
final_CT_prop_major <- data.frame(final_CT_prop)[,make.names(final_CT_prop_major_celltypes)]
final_CT_prop_major <- data.frame(final_CT_prop_major)
dim(final_CT_prop_major)
colnames(final_CT_prop_major)
if(use_major_ct == "true"){
  final_CT_prop_filter = final_CT_prop_major
}
dim(final_CT_prop_filter)
assignCells_test <- function(cellWeights) {
  if (is.null(dim(cellWeights))) {
    if (any(cellWeights == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      return(matrix(1, nrow = length(cellWeights), ncol = 1))
    }
  } else {
    if (any(rowSums(cellWeights) == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      # normalize weights
      normWeights <- sweep(cellWeights, 1,
                           FUN = "/",
                           STATS = apply(cellWeights, 1, sum)
      )
      # sample weights
      wSamp <- apply(normWeights, 1, function(prob) {
        stats::rmultinom(n = 1, prob = prob, size = 1)
      })
      # If there is only one lineage, wSamp is a vector so we need to adjust for that
      if (is.null(dim(wSamp))) {
        wSamp <- matrix(wSamp, ncol = 1)
      } else {
        wSamp <- t(wSamp)
      }
      return(wSamp)
    }
  }
}
findKnots_test <- function(nknots, pseudotime, wSamp) {
  # Easier to recreate them all here than to pass them on
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("t",ii), pseudotime[,ii])
  }
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("l",ii),1*(wSamp[,ii] == 1))
  }
  
  # Get the times for the knots
  tAll <- c()
  for (ii in seq_len(nrow(pseudotime))) {
    tAll[ii] <- pseudotime[ii, which(as.logical(wSamp[ii,]))]
  }
  
  knotLocs <- stats::quantile(tAll, probs = (0:(nknots - 1)) / (nknots - 1))
  if (any(duplicated(knotLocs))) {
    # fix pathological case where cells can be squeezed on one pseudotime value.
    # take knots solely based on longest lineage
    knotLocs <- stats::quantile(t1[l1 == 1],
                                probs = (0:(nknots - 1)) / (nknots - 1))
    # if duplication still occurs, get average btw 2 points for dups.
    if (any(duplicated(knotLocs))) {
      dupId <- duplicated(knotLocs)
      # if it's the last knot, get duplicates from end and replace by mean
      if (max(which(dupId)) == length(knotLocs)) {
        dupId <- duplicated(knotLocs, fromLast = TRUE)
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                  knotLocs[which(dupId) + 1]))
      } else {
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                  knotLocs[which(dupId) + 1]))
      }
    }
    # if this doesn't fix it, get evenly spaced knots with warning
    if (any(duplicated(knotLocs))) {
      knotLocs <- seq(min(tAll), max(tAll), length = nknots)
    }
  }
  
  maxT <- max(pseudotime[,1])
  if (ncol(pseudotime) > 1) {
    maxT <- c()
    # note that first lineage should correspond to the longest, hence the
    # 100% quantile end point is captured.
    for (jj in 2:ncol(pseudotime)) {
      maxT[jj - 1] <- max(get(paste0("t", jj))[get(paste0("l",jj)) == 1])
    }
  }
  # if max is already a knot we can remove that
  if (all(maxT %in% knotLocs)) {
    knots <- knotLocs
  } else {
    maxT <- maxT[!maxT %in% knotLocs]
    replaceId <- vapply(maxT, function(ll){
      which.min(abs(ll - knotLocs))
    }, FUN.VALUE = 1)
    knotLocs[replaceId] <- maxT
    if (!all(maxT %in% knotLocs)) {
      # if not all end points are at knots, return a warning, but keep
      # quantile spaced knots.
      warning(paste0("Impossible to place a knot at all endpoints.",
                     "Increase the number of knots to avoid this issue."))
    }
    knots <- knotLocs
  }
  # guarantees that first knot is 0 and last knot is maximum pseudotime.
  knots[1] <- min(tAll)
  knots[nknots] <- max(tAll)
  knotList <- lapply(seq_len(ncol(pseudotime)), function(i){
    knots
  })
  names(knotList) <- paste0("t", seq_len(ncol(pseudotime)))
  return(knotList)
}
#celltype = "Basal"
dim(lumen_celltype_sce)
lumen_nuclei_RNA_sce_per_ct <- lumen_nuclei_RNA_sce_per_ct[,make.names(lumen_nuclei_RNA_sce_per_ct$final_CT)==celltype]
lumen_nuclei_RNA_sce_per_ct$celltype_count <- as.numeric(final_CT_count[lumen_nuclei_RNA_sce_per_ct$lumen_id,
                                                                         make.names(colnames(final_CT_count))==celltype])
dim(lumen_nuclei_RNA_sce_per_ct)
## filter lumens, keep lumen with >=3 such cell types
lumen_nuclei_RNA_sce_per_ct <- lumen_nuclei_RNA_sce_per_ct[,lumen_nuclei_RNA_sce_per_ct$celltype_count>=3]

ncoefficients <- 4*ncol(final_CT_prop_filter)+4+1+4

if(min_lumen < ncoefficients){
  min_lumen = ncoefficients
}

min_lumen

stopifnot(ncol(lumen_nuclei_RNA_sce_per_ct) >= min_lumen)

dim(lumen_nuclei_RNA_sce_per_ct)
cellWeights <- matrix(rep(1,ncol(lumen_nuclei_RNA_sce_per_ct)),ncol = 1)
wSamp  <- assignCells_test(cellWeights)
# pseudotim knots 5
knotList <- findKnots_test(nknots=5,
                           pseudotime = matrix(lumen_nuclei_RNA_sce_per_ct$Lineage1,
                                               nrow=length(lumen_nuclei_RNA_sce_per_ct$Lineage1)),
                           wSamp = wSamp)
## accounting for 15 cell types proportions
final_ct_prop_filter_ct <- final_CT_prop_filter[lumen_nuclei_RNA_sce_per_ct$lumen_id,]
## manually set the knotlist, and use the same knotlist for each gene
## celltype proportion nknots =5  

knot_list_per_ct <- sapply(colnames(final_ct_prop_filter_ct),
                           function(ct_name){
  each_ct <- final_ct_prop_filter_ct[,ct_name]
  each_ct <- as.matrix(each_ct,nrow=length(each_ct))
  each_res <- findKnots_test(5,each_ct,wSamp)
  #names(each_res) <- ct_name
  each_res
})

names(knot_list_per_ct) <- gsub(".t1","",names(knot_list_per_ct))
knot_list_per_ct
knot_list_t1_per_ct <-  c(knotList,knot_list_per_ct)

png(out_knots_png,
    width = 1400,height = 1200,res = 100,units = "px")
par(mfrow=c(5,7))
for(x in names(knot_list_t1_per_ct)) {
  plot(x = knot_list_t1_per_ct[[x]],
       y=rep(1,length(knot_list_t1_per_ct[[x]])),main = x,ylab = "")
}
dev.off()

gam_models_list = list()
set.seed(100)
sterms <- paste0("s(t1,bs = 'cr',k=5, by = l,id=1)") 

smooth_terms <- paste0(sterms,"+",paste0(vapply(colnames(final_CT_prop_filter),
                                                function(ct){
                                                  paste0("s(",ct,",bs = 'cr',k=5)") 
                                                  },
                                                FUN.VALUE = "formula"),collapse = "+"),
                       " + offset(log(total_tx_count))")

as.formula(paste0("gene_count ~ -1 + U +",smooth_terms))

## only test genes that have >=3 counts in 50% of the lumens

# test_genes <- rowSums(lumen_nuclei_RNA_sce_per_ct@assays@data$nuclei_trans_count_perCT >= 3)/ncol(lumen_nuclei_RNA_sce_per_ct) 
# names(test_genes)[test_genes >= filter_by]

## only test genes that were sig in bulk results
bulk_sig_genes <- read.csv(
  "./analysis/revision_2024/lumen_order_t4/offset_total_trans/gene_specc/gene_specc_c4_257genes.csv")
bulk_sig_genes <- bulk_sig_genes$gene
length(bulk_sig_genes)
test_genes_c <- rowSums(lumen_nuclei_RNA_sce_per_ct@assays@data$nuclei_trans_count_perCT >= 3)/ncol(lumen_nuclei_RNA_sce_per_ct) 
test_genes_c <- names(test_genes_c)[test_genes_c >= filter_by]
length(test_genes_c)
test_genes <- intersect(bulk_sig_genes,test_genes_c)

aCap_gam_models_list = bplapply(test_genes,
                                function(gene_name){
  gene_df <- data.frame(gene_count = lumen_nuclei_RNA_sce_per_ct@assays@data$nuclei_trans_count_perCT[gene_name,],
                        #cell_count = lumen_nuclei_RNA_sce_per_ct$celltype_count,
                        total_tx_count = as.numeric(
                          colSums(lumen_nuclei_RNA_sce_per_ct@assays@data$nuclei_trans_count_perCT)),
                        t1 = lumen_nuclei_RNA_sce_per_ct$Lineage1,
                        final_ct_prop_filter_ct,
                        #U = model.matrix(~0+tma, colData(lumen_nuclei_RNA_sce_per_ct)),
                        l = factor(1)) 
  U = model.matrix(~0+tma,colData(lumen_nuclei_RNA_sce_per_ct))
  model_res <- gam(formula =as.formula(paste0("gene_count ~ -1 + U + ",smooth_terms)),
                   knotList = c(knotList,knot_list_per_ct),
                   control = mgcv::gam.control(),
                   family = "nb",data=gene_df)
  model_res[["gene"]] <- gene_name
  model_res
},BPPARAM = BPPARAM)

saveRDS(aCap_gam_models_list,outrds)
write.table(test_genes, file = paste0(outrds,".txt"))


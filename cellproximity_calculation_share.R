library(sf)
library(dplyr)

filter <- dplyr::filter
select <- dplyr::select

## load functions 
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

## calculate distance function
calc_d <- function(pt1, pt2) {
  d <- sqrt( (( pt2[2] - pt1[2] )^2 ) + (( pt2[1] - pt1[1] )^2 ) )
  return(d)
}

getDir <- function(origin, pt2){
  pt2 <- pt2 - origin
  deg <- abs(atan(pt2[2]/pt2[1]) * (180/pi))
  if(pt2[1] > 0 & pt2[2] > 0){deg <- deg
  } else if(pt2[1] < 0 & pt2[2] > 0){ deg <- 180-deg
  } else if(pt2[1] < 0 & pt2[2] < 0){ deg <- 180 + deg
  } else if(pt2[1] > 0 & pt2[2] < 0){ deg <- 360 - deg}
  return(deg)
}

assign.angle <- function(degrees, range){
  require(data.table)
  a=data.table(degrees=degrees)
  a[,merge:=degrees]
  
  b=data.table(range=range)
  b[,merge:=range]
  
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  Merge_a_b=b[a,roll='nearest']
  return(Merge_a_b)
}

anchorPoint <- function(degrees, anchor){
  object <- (degrees + (360 - anchor)) 
  object <- ifelse(object > 360, object - 360, object)
  return(object)
}


# prepare data for runs
# x <-  readRDS('/scratch/nhadad/lung_datasets/prox_analysis_av/xenium_final_CT_plus_niches_062724.rds')
## generate niche list
# tmp <- x$CNiche %>% unique() %>% as.data.frame()
# write.table(tmp, file = file.path(inDir, "avannan_niche_list.txt"),
#             quote = F, sep = "\t", col.names = F, row.names = F)
# ## generate tniche list 
# tmp <- x$TNiche %>% unique() %>% as.data.frame()
# write.table(tmp, file = file.path(inDir, "avannan_tniche_list.txt"),
#             quote = F, sep = "\t", col.names = F, row.names = F)
# ## generate celltype list 
# tmp <- x$final_CT %>% unique() %>% as.data.frame()
# write.table(tmp, file = file.path(inDir, "avannan_celltype_list.txt"),
#             quote = F, sep = "\t", col.names = F, row.names = F)
# 
# meta <- x@meta.data
# meta$cid <- rownames(meta)
# saveRDS(meta, file.path(inDir, "xenium_final_CT_plus_niches_062724_metadata.rds"))


args <- commandArgs(trailingOnly = TRUE)
cell_type_id <- if (length(args) > 0) args[1] else ''
cell_type_id_out <- gsub("\\/", "\\_", cell_type_id)
print(cell_type_id)

## load data and define code paramters
inDir <- './input/'
x <- read.csv(file.path(inDir, "xenium_final_CT_plus_niches_062724_metadata.csv"), row.names = 1) # xenium object meta.data 
# code loops on samples
sample_ids <- x %>% 
  pull(sample) %>% as.character %>% unique() 

nCores <- 24
library(doParallel)
registerDoParallel(cores=nCores)
proximal_cell_pop <- foreach (i=1:length(sample_ids), .combine = rbind) %dopar% {
  
  require("sf"); require("dplyr")
  sid <- sample_ids[i]; print(sid)
  obj.sid <- subset(x, sample == sid)
  
  # 1. get cell coordinates
  meta <- obj.sid
  cells_coord <- meta %>% 
    filter(final_CT == cell_type_id)
  if(nrow(cells_coord) == 0){return(NULL)} ## skip if sample does not contain celltype
  
  # 2. define radius
  cell_counts <- c()
  r <- 60 ## set radius ~30 is 3 cells 
  print(r)
  for(j in 1:nrow(cells_coord)){
    print(j)
    cell_id <- rownames(cells_coord)[j]
    split_cell.id <- strsplit(cell_id, ";") %>% unlist()
    cell_x <- cells_coord$x_centroid[j]
    cell_y <- cells_coord$y_centroid[j]
    
    circ_coord <- circleFun(c(cell_x, cell_y), r, 1000)
    circ_coord <- circ_coord[!duplicated(circ_coord),]
    
    # 3. get cell composition within circle
    # 3a. subset meta data to include all points with ring boundary 
    all_cell_coord <- data.frame(x = meta$x_centroid, y = meta$y_centroid)
    tmp <- circ_coord[,c(1:2)]
    io <- splancs::inout(all_cell_coord, tmp)
    
    if(mean(io) == 0){next}
    extract_cells <- as.data.frame(meta)[which(io == TRUE),] ## get cell within radius
    if(nrow(extract_cells) == 1){next}
    
    # calculate distance and angle between cells
    tmp <- extract_cells[,c('x_centroid', 'y_centroid')]
    tmp <- rbind(tmp[which(rownames(tmp) %in% cell_id),], tmp)
    tmp <- tmp[!duplicated(tmp),]
    
    toCheck <- combn(rownames(tmp), 2, simplify = FALSE)
    names(toCheck) <-
      sapply(toCheck, paste, collapse = " - ")
    toCheck <- toCheck[grep(pattern = cell_id, toCheck)]
    
    ## distance
    tmp.d <- sapply(toCheck, function(j){
      calc_d(tmp[j[1],c(1,2)], tmp[j[2],c(1,2)]) })
    names(tmp.d) <- gsub('\\.y_centroid', '', names(tmp.d))
    tmp.d <- as.data.frame(do.call(rbind, tmp.d))
    tmp.d$cell_id <- rownames(tmp.d)
    tmp.d <- tidyr::separate(tmp.d, cell_id, sep = "\\ - ", c('cell.a', 'cell.b'))
    rownames(tmp.d) <- NULL
    
    ## angle
    tmp.degree <- sapply(toCheck, function(j){
      getDir(c(cell_x, cell_y), tmp[j[2],c(1,2)]) })
    names(tmp.degree) <- gsub('\\.y_centroid', '', names(tmp.degree))
    tmp.degree <- as.data.frame(do.call(rbind, tmp.degree))
    tmp.degree$cell_id <- rownames(tmp.degree)
    tmp.degree <- tidyr::separate(tmp.degree, cell_id, sep = "\\ - ", c('cell.a', 'cell.b'))
    rownames(tmp.degree) <- NULL
    
    ## add degrees to distance table
    tmp.d$degree <- tmp.degree$V1[match(tmp.d$cell.b, tmp.degree$cell.b)]
    tmp.d$celltypeA <- extract_cells$final_CT[match(tmp.d$cell.a, rownames(extract_cells))]
    tmp.d$celltypeB <- extract_cells$final_CT[match(tmp.d$cell.b, rownames(extract_cells))]
    
    
    ## assign each cell to a range category - this will be calculated in the logistic regression code
    range <- seq(30, 360, 30) - 15
    Merge_a_b <- assign.angle(degrees = tmp.d$degree, range = range)
    tmp.d$angle <- Merge_a_b$range[match(tmp.d$degree, Merge_a_b$degrees)]
    tmp.d$cell.a <- cell_id
    tmp.d <- tmp.d[!duplicated(tmp.d),]
    # append
    cell_counts <- rbind(cell_counts, tmp.d)
  }
  cell_counts$sid <- sid
  return(cell_counts)
}

stopImplicitCluster() 
out_path <- './output/prox_files'
out_fname <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
saveRDS(proximal_cell_pop, file.path(out_path, out_fname))  



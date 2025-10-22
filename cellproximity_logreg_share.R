library(data.table)
library(dplyr)
`%notin%` <- Negate(`%in%`)

# quick or calculation
getOR <- function(log_odds, se_log_odds) {
  OR <- exp(log_odds)
  ci_low <- exp(log_odds - 1.96 * se_log_odds)
  ci_hi <- exp(log_odds + 1.96 * se_log_odds)
  return(cbind(OR = OR, ci_low = ci_low, ci_high = ci_hi))
}

# Direct namespace imports
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
group_by <- dplyr::group_by
summarise <- dplyr::summarise
arrange <- dplyr::arrange

# Load data and define parameters
inDir <- './input/'
outDir <- './output/prox_enrichment'
meta <- read.csv(file.path(inDir, "xenium_final_CT_plus_niches_062724_metadata.csv"), row.names = 1)
meta <- meta[, !grepl('^leiden', colnames(meta))]

# define contrasts
options(contrasts = c("contr.treatment", "contr.poly"))

# Parameters
comparisons <- c('all', 'LF', 'MF', 'Unaffected') # other comps c('all', 'disease', 'LF', 'MF', 'ILD', 'Unaffected')
cell_type_niche_proximity <- list()  

for (comp in comparisons) {
  # Get unique cell types 
  cell_types <- names(table(meta$final_CT))
  
  celltype_logit_res <- data.frame()
  
  for (celltype_id in cell_types) {
    cat("Processing:", celltype_id, "\n")
    
    # Load proximity data
    celltype_id_out <- gsub("\\/", "\\_", celltype_id)
    inpath <- '/scratch/nhadad/lung_datasets/prox_analysis_av/output/prox_files'
    f <- paste0(celltype_id_out, '_cell_dist_adjdegree.rds')
    df <- readRDS(file.path(inpath, f))
    df <- df[order(df$V1), ]  # Sort by distance 
    
    # Match metadata
    df$NicheA <- meta$CNiche[match(df$cell.a, rownames(meta))]
    df$NicheB <- meta$CNiche[match(df$cell.b, rownames(meta))]
    df$celltypeA <- meta$final_CT[match(df$cell.a, rownames(meta))]
    df$celltypeB <- meta$final_CT[match(df$cell.b, rownames(meta))]
    
    # Compute global threshold for closest neighbors (across all cells)
    df_ranked <- df %>%
      group_by(cell.a) %>%
      mutate(rank_overall = row_number()) %>%  # Rank neighbors by distance
      ungroup()
    d_threshold <- df_ranked %>%
      filter(rank_overall == 1) %>%
      summarise(d_threshold = mean(V1) + 2 * sd(V1), .groups = 'drop') %>%
      pull(d_threshold)
    # cat("Distance threshold:", d_threshold, "\n")
    
    # Bin angles and rank within bins
    degrees <- df$degree
    range <- seq(30, 360, 30) - 15
    a <- data.table(degrees = degrees)
    a[, merge := degrees]
    b <- data.table(range = range)
    b[, merge := range]
    setkeyv(a, 'merge')
    setkeyv(b, 'merge')
    Merge_a_b <- b[a, roll = 'nearest']
    df$angle <- Merge_a_b$range[match(df$degree, Merge_a_b$degrees)]
    
    # Rank within each cell.a and angle
    df <- df %>%
      group_by(cell.a, angle) %>%
      arrange(V1) %>%  # within-group ordering
      mutate(idx = row_number()) %>%
      ungroup()
    
    # Keep first neighbor per angle and within distance threshold 
    tmp0 <- df %>%
      filter(V1 < d_threshold & idx == 1)
    if (nrow(tmp0) == 0) next
    
    tmp0 <- tmp0 %>%
      select(cell.b, celltypeB, sid) %>%
      distinct()
    
    # Subset metadata based on comparison
    meta_tmp <- switch(comp,
                       'all' = meta,
                       'disease' = filter(meta, sample_type != 'Unaffected'),
                       'LF' = filter(meta, sample_type == 'LF'),
                       'MF' = filter(meta, sample_type == 'MF'),
                       'ILD' = filter(meta, sample_type == 'ILD'),
                       'Unaffected' = filter(meta, sample_type == 'Unaffected')
    )
    
    # exclude self from analysis 
    meta_tmp <- filter(meta_tmp, final_CT != celltype_id)
    
    # skip if n celltypes is one or less
    levels <- names(table(meta_tmp$final_CT))
    if (length(levels) <= 1) next
    
    # binary proximity
    meta_tmp$proximal <- ifelse(rownames(meta_tmp) %in% tmp0$cell.b, 1, 0)
    
    # exclude celltypes with 0 proximity (avoid model convergence)
    prox_table <- table(meta_tmp$proximal, meta_tmp$final_CT)
    no_prox_ct <- names(which(prox_table[2, ] == 0))
    if (length(no_prox_ct) > 0) {
      meta_tmp <- filter(meta_tmp, final_CT %notin% no_prox_ct)
      levels <- names(table(meta_tmp$final_CT))
      if (length(levels) <= 1) next
    }
    
    # set levels as factors
    meta_tmp$final_CT <- factor(meta_tmp$final_CT, levels = levels)
    
    # log reg model
    mylogit <- glm(proximal ~ final_CT, data = meta_tmp, family = "binomial")
    
    # predict log-odds 
    newdata1 <- data.frame(final_CT = factor(levels, levels = levels))
    pred_link <- predict(mylogit, newdata = newdata1, type = "link", se.fit = TRUE)
    lp_fit <- pred_link$fit
    se_lp <- pred_link$se.fit
    n_levels <- length(levels)
    avg_lp <- mean(lp_fit)
    
    # get deviation relative to mean
    delta <- lp_fit - avg_lp
    
    # SE
    se_delta <- numeric(n_levels)
    for (k in seq_len(n_levels)) {
      var_lp_k <- se_lp[k]^2
      var_others <- sum(se_lp[-k]^2)
      se_delta[k] <- sqrt( ((n_levels - 1) / n_levels)^2 * var_lp_k + (1 / n_levels)^2 * var_others )
    }
    
    # Z-values and p-values for null
    z_val <- delta / se_delta
    p_val <- 2 * pnorm(-abs(z_val))
    
    # Odds ratios and CIs relative to mean
    or_results <- getOR(delta, se_delta) # for smaller n or more accurate OR use confint()
    
    # compile results
    model_coefficients <- data.frame(
      Estimate = delta,
      `Std. Error` = se_delta,
      `z value` = z_val,
      `Pr(>|z|)` = p_val,
      row.names = levels
    )
    colnames(or_results) <- c("OR", "ci5", "ci95")
    
    # get probabilities
    newdata1$rankP <- predict(mylogit, newdata = newdata1, type = "response")
    
    # Combine
    results <- cbind(newdata1, model_coefficients, or_results)
    results$Celltype <- celltype_id
    rownames(results) <- NULL
    
    # Append
    celltype_logit_res <- rbind(celltype_logit_res, results)
    
    cat("Cell type complete\n")
  }
  
  #  save results per comparison
  cell_type_niche_proximity[[comp]] <- celltype_logit_res
  # write.csv(celltype_logit_res, file.path(outDir, paste0(comp, '_LogReg_enrichment.csv')), row.names = FALSE)
}
gc()

# table for plot
colnames(celltype_logit_res) <- c('target', 'prob', 'coefficient', 'std.err', 'z.value', 'p.value', 'OR', 'ci05', 'ci95', 'source')
# FDR per celltype
celltype_logit_res <- celltype_logit_res %>% 
  group_by(source) %>% 
  mutate(fdr = p.adjust(p.value, method = 'fdr'))
celltype_logit_res$logOR <- log(celltype_logit_res$OR)

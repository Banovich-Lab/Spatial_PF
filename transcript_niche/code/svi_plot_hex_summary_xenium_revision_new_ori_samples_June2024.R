### generate hex plot for xenium sample

#.libPaths("/data/gpfs/projects/punim0741/rlyu/Software/Rlib/4.2.0")
.libPaths("/mnt/beegfs/mccarthy/backed_up/general/rlyu/Software/Rlibs/4.1.2")

library(hexbin)
library(ggplot2)
library(tibble)
library(cowplot)
#library(ggtree)
library(Seurat)
library(dplyr)
library(Polychrome)

args <- (commandArgs(trailingOnly = TRUE))

for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

print(sample)
print(node_meta_file)
print(cid_results)

print(outpng)
print(bin_thresh)
print(bin_width)
bin_width = as.numeric(bin_width)
bin_thresh = as.numeric(bin_thresh)


# s_n <- paste0(sample,"_embedded_",tma,".csv")
#node_meta_file <- paste0("./output/xenium/fullPanel/graphs/",s_n)
node_meta <- readr::read_csv(node_meta_file,
                             num_threads = 4,
                             show_col_types = FALSE)
node_meta <- node_meta[,c("transcript_id" ,"x_location", "y_location", "z_location",
                          "feature_name")]


gmm_id <- read.table(file =cid_results)

node_meta$gmm <- as.character(gmm_id$V1)
kmeans13_colors <- c(
  "#2f4f4f","#000000","#008000","#4b0082","#ff0000","#ffd700", '#b5bd61',"#00ffff",
  "#0000ff","#ff69b4","#1e90ff","#ffdab9","#ff00ff")

names(kmeans13_colors) <- 0:12

k_15_colors <- c("#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919",
                 "#005C31", "#2BCE48", 
                 "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", 
                 "#C20088", "#003380", "#FFA405" )
names(k_15_colors) <- paste0("c",0:14)
gmm8_colors <- k_15_colors[paste0("c",0:7)]
gmm8_k <- 0:7
names(gmm8_colors) <- gmm8_k


gmm9_colors <- k_15_colors[c(paste0("c",0:6),"c8","c7")]
gmm9_k <- c(4,5,1,6,2,3,8,0,7)

names(gmm9_colors) <- (gmm9_k)

gmm10_colors <- k_15_colors[c(paste0("c",0:4),"c9","c5","c10","c6","c7")]
gmm10_k <- c(3,7,6,0,1,4,2,5,8,9)
names(gmm10_colors) <- gmm10_k


gmm11_colors <- k_15_colors[c("c0","c12","c2","c3","c4","c9","c5","c10","c6","c7","c11")]
gmm11_k <- c(1,8,6,4,3,7,2,9,10,0,5)
names(gmm11_colors) <- gmm11_k

gmm12_colors <- k_15_colors[c("c0","c12","c2","c3","c4","c9","c5","c1","c10","c6","c7","c11")]
gmm12_k <- c(8,9,5,10,1,0,6,7,2,4,3,11)
names(gmm12_colors) <- gmm12_k


gmm13_colors <- k_15_colors[c("c0","c12","c2","c3","c4","c9","c5","c13","c1","c10","c6","c7","c11")]
gmm13_k <- c(9,4,1,3,7,6,5,8,10,12,2,0,11)
names(gmm13_colors) <- gmm13_k


matched_color_list <- list("gmm8_colors" = gmm8_colors,
                           "gmm9_colors" = gmm9_colors,
                           "gmm10_colors" = gmm10_colors,
                           "gmm11_colors" = gmm11_colors,
                           "gmm12_colors" = gmm12_colors,
                           "gmm13_colors" = gmm13_colors )

if(paste0("gmm",length(unique(node_meta$gmm)),"_colors") %in% names(matched_color_list)){
  k_colors <- matched_color_list[[paste0("gmm",length(unique(node_meta$gmm)),"_colors")]]
} else {
  k_colors = DiscretePalette(n = length(unique(node_meta$gmm)))
}

#k_colors[]
k_colors_new = k_colors
k_colors_new["3"] = k_colors['0']
k_colors_new["5"] = k_colors['1']
k_colors_new["1"] = k_colors['10']
k_colors_new["7"] = k_colors['11']
k_colors_new["0"] = k_colors['2']
k_colors_new["2"] = k_colors['3']
k_colors_new["10"] = k_colors['4']
k_colors_new["6"] = k_colors['5']
k_colors_new["8"] = k_colors['6']
k_colors_new["4"] = k_colors['7']
k_colors_new["11"] = k_colors['8']


count_nodes <- function(x,threshold=bin_thresh){
  ifelse( length(x)>threshold,  length(x),NA)
}
get_major_cluster <- function(x,threshold=bin_thresh){
  ifelse( length(x)>threshold, names(table(x))[which.max(table(x))],NA)
  
}
get_major_cluster <- function(x,threshold=bin_thresh){
  ifelse( length(x)>threshold, names(table(x))[which.max(table(x))],NA)
  
}
node_meta$gmm <- as.numeric(node_meta$gmm)

dim(node_meta)

p1 = node_meta %>% ggplot() + 
  stat_summary_hex( aes(x = x_location, y = y_location, z = gmm),
                    size = 0.05, fun = count_nodes, binwidth = bin_width)+
  theme_bw()+xlab(sample)

color_labels <- paste0("T",1:length(k_colors_new))
names(color_labels) <- as.character(0:(length(k_colors_new) - 1))

p2 = node_meta %>% ggplot() + 
  stat_summary_hex( aes(x = x_location, y = y_location, z = gmm),
                    size = 0.05, 
                    fun = get_major_cluster, binwidth = bin_width)+
  theme_bw(base_size=18)+theme(panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank(), 
                               panel.border = element_blank(),
                               axis.line = element_line(color="black"),
                               axis.ticks = element_line(color="black"),
                               axis.text = element_text(color = "black"))+
  scale_fill_manual(values = k_colors_new, 
                    breaks = as.character(0:(length(k_colors_new) - 1)),
                    labels =color_labels) + 
  guides(fill=guide_legend(title="Niche"))+
  scale_y_reverse()+coord_equal()+
  xlab(sample)
## legend Niche, remove labels.
## 
p3 = node_meta %>% ggplot() + 
  stat_summary_hex( aes(x = x_location, y = y_location, z = gmm),
                    size = 0.1, fun = get_major_cluster, binwidth = bin_width)+
  theme_bw()+scale_fill_manual(values = k_colors_new,
                               breaks = as.character(0:(length(k_colors_new) - 1)),
                               labels =color_labels)+coord_equal()+
  xlab(sample)+
  facet_wrap(.~gmm,ncol = 5)


p = p1 + p2

ggsave(plot = p, filename = outpng,width = 10,height = 4,dpi = 150)
ggsave(p2, filename = paste0(outpng,".clusteronly.png"),width = 9,height = 7,dpi = 300)
ggsave(p2+theme(legend.position="None",
                axis.title.y=element_blank(),
                axis.text=element_blank(),
                axis.ticks=element_blank(),
                axis.line = element_blank()), filename = paste0(outpng,".clusteronly_nolegend.png"),
       width = 8,height = 8,dpi = 300)
ggsave(p2+theme(legend.position="None",
                axis.title.y=element_blank(),
                axis.text=element_blank(),
                axis.ticks=element_blank(),
                axis.line = element_blank()), filename = paste0(outpng,".clusteronly_nolegend.pdf"),
       width = 8,height = 8,dpi = 300,device="pdf")
ggsave(p3, filename = paste0(outpng,".split_gmm.png"),width = 20,height = 10,dpi = 150)
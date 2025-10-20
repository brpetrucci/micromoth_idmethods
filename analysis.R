###
# load packages

# dplyr
library(dplyr)

# tidyr
library(tidyr)

# ggplot
library(ggplot2)

# vegan
library(vegan)

# treeio
library(treeio)

# ape
library(ape)

###
# load files

# base directory
base_dir <- "/Users/petrucci/Documents/research/micromoths_idmethods/"

# ids
ids_raw <- read.delim(paste0(base_dir, "final_ids.csv"), sep = ",")
ids_raw <- ids_raw[!(ids_raw$coi_id %in% c("too short", "too many ambiguous bases")), ]

# separate first column into 4 columns
ids <- ids_raw %>%
  extract(
    ACC.,
    into = c("Site", "City", "Room", "N"),
    regex = "^([A-Z]{2})([A-Z])\\-([A-Z])\\-([0-9]+)$"
  ) %>%
  mutate(N = as.integer(N)) %>%
  select(Site, City, Room, N, morphological_id, genitalia_id, coi_id)

# alter city names
city_names <- c(C = "Cambridge", N = "Newton", W = "Waltham")
ids$City <- city_names[ids$City]

# add ACC
ids$ACC <- ids_raw$ACC.

# genitalia measurements table
gen_measures <- na.omit(read.delim(paste0(base_dir, 
                                          "genitalia_measurements.csv"),
                           sep = ","))

###
# rename alignment

# read alignment
coi <- read.FASTA(paste0(base_dir, "mafft_coi.fasta"))

# match names
matches <- unlist(lapply(ids$coi_id[match(names(coi), ids$ACC)],
                         function(x) paste(strsplit(x, " ")[[1]], 
                                           collapse = "_")))

# get new name
names_align <- paste0(matches, "_",
                      ave(seq_along(matches), matches, FUN = seq_along))

# change name of alignment
named_coi <- coi
names(named_coi) <- names_align

# save named alignment
write.FASTA(named_coi, paste0(base_dir, "named_coi.fasta"))

# save names
write.table(data.frame(ACC = names(coi), names = names_align),
            file = paste0(base_dir, "names_align.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

###
# calculate correctness percentages

# correctness columns
ids$morphological_correct <- ids$morphological_id == ids$coi_id
ids$genitalia_correct <- ids$genitalia_id == ids$coi_id

# total percent correct
percent_correct <- c(morpho = sum(ids$morphological_correct) / nrow(ids),
                     genitalia = sum(ids$genitalia_correct) / nrow(ids))

# percent correct per species table
correct_species <- ids %>%
  select(morphological_correct, genitalia_correct, coi_id) %>%
  group_by(coi_id) %>%
  summarise(morphological_percentage = mean(morphological_correct),
            genitalia_percentage = mean(genitalia_correct),
            .groups = "drop") 

###
# Shannon's diversity index

# function to calculate shannon's div
shannon_div <- function(div) {
  # make vector of results
  res <- c()
  
  # iterate through sites
  for (i in 1:nrow(div)) {
    # vector of prop * ln
    vals <- c()
    
    # iterate through species
    for (j in 1:ncol(div)) {
      # calculate proportion
      prop <- div[i, j] / sum(div[i, ])
      
      # get val
      vals <- c(vals, ifelse(prop == 0, 0, prop * log(prop)))
    }
    
    # add to results
    res <- c(res, -1 * sum(vals))
  }
  
  # name res
  names(res) <- rownames(div)
  
  # get final value
  return(res)
}

# function to calculate species richness
n_sp <- function(div) {
  # vector to hold richness
  n_sp <- unlist(lapply(1:nrow(div), function(i) sum(div[i, ] > 0)))
  
  # name vector
  names(n_sp) <- rownames(div)
  
  # return it
  return(n_sp)
}

# function to calculate pielou's evenness
pielous_evenness <- function(div) {
  # get species number per pop
  richness <- n_sp(div)
  
  # get shannon's div per pop
  shannon_div_pop <- shannon_div(div)
  
  # divide one by the other
  return(shannon_div_pop / log(richness))
}

# function to get all 3 into a data frame for a given column in ids
div_df <- function(ids, col) {
  # matrix of diversity
  pop_div <- table(ids[, c(col, "coi_id")])

  # species richness
  pop_richness <- n_sp(pop_div)

  # get site and city shannon's
  pop_shannon_div <- shannon_div(pop_div)

  # get evenness
  pop_evenness <- pielous_evenness(pop_div)
  
  # make data frame
  df <- na.omit(data.frame(Diversity = pop_shannon_div,
                           Evenness = pop_evenness,
                           Richness = pop_richness))
  
  # add column name and remove rownames
  df[[col]] <- rownames(df)
  rownames(df) <- 1:nrow(df)
  
  # return data frame
  return(df)
}

# make a data frame for site, city and room
site_df <- div_df(ids, "Site")
city_df <- div_df(ids, "City")
room_df <- div_df(ids, "Room")

# add city to site_df
site_df$City <- unlist(lapply(site_df$Site, function(x) ids$City[ids$Site == x][1]))

###
# make figures

# site and city palettes
site_palette <- colorRampPalette(palette.colors(8, palette = "R4"))(length(unique(ids$Site)))
names(site_palette) <- unique(ids$Site)
city_palette <- palette.colors(3, palette = "Okabe-Ito")
names(city_palette) <- unique(ids$City)

# function to make even-rich plot and div barplot
div_plots <- function(df, col, shape, palette) {
  # make a richness by evenness plot
  even_rich <- ggplot(df, aes(x = Richness, y = Evenness,
                              colour = .data[[col]], shape = .data[[shape]])) +
    scale_colour_manual(values = palette) +
    geom_point(size = 5) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill = "transparent"),
          panel.background = element_blank(),
          )
  
  # make diversity barplot
  div <- ggplot(df, aes(y = Diversity, x = .data[[col]], colour = .data[[col]],
                        fill = .data[[col]])) +
    scale_colour_manual(values = palette) +
    scale_fill_manual(values = palette) +
    geom_bar(stat = "identity") +
    theme_bw()
  
  # return plots
  return(list(EVENRICH = even_rich, DIV = div))
}

# function to save a plot as a pdf
pdf_plot <- function(fig, filename, width, height) {
  # start pdf
  pdf(filename, width = width, height = height)
  
  # plot
  print(fig)
  
  # close connection
  dev.off()
}

# get plots for each pop
site_plots <- div_plots(site_df, "Site", "City", site_palette)
city_plots <- div_plots(city_df, "City", "City", city_palette)

# save plots separately
pdf_plot(site_plots[[1]], paste0(base_dir, "site_evenness.pdf"),
         width = 10, height = 10)
pdf_plot(site_plots[[2]], paste0(base_dir, "site_diversity.pdf"),
         width = 10, height = 10)
pdf_plot(city_plots[[1]], paste0(base_dir, "city_evenness.pdf"),
         width = 10, height = 10)
pdf_plot(city_plots[[2]], paste0(base_dir, "city_diversity.pdf"),
         width = 10, height = 10)


###
# make PCAs for genitalia measurement

# get principal components
pp <- capscale(gen_measures[, 5:9] ~ 1, distance = "euclidean")

# get loadings
loadings <- as.data.frame(scores(pp, display = "species", scaling = 1)) * 0.7
loadings$Species <- rownames(loadings)

# get sample scores
pca_s <- as.data.frame(scores(pp, display = "sites", scaling = 1))

# get eigenvalues
eigs <- as.data.frame(pp$CA$eig)

# make rownames a column
eigs$MDS <- rownames(eigs)

# get variance explained by each PCA vector
eig_var <- eigs[, 1] / sum(eigs[, 1])

# # get colors
# colors_pop <- colors_vec

# get axes names
x_axis <- colnames(pca_s[1])
y_axis <- colnames(pca_s[2])

# pca plot
pca_plot <- ggplot(pca_s, aes_string(x = x_axis, y = y_axis)) +
  theme_bw() +
  geom_point(aes(colour = gen_measures[, "Species"], 
                 shape = gen_measures[, "City"]), size = 2, stroke = 1)  +
  stat_ellipse(type = "t",
               aes(color = gen_measures[, "Species"]), show.legend = NA, lwd = 1) + 
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = MDS1, yend = MDS2),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black"
  ) +
  geom_text(
    data = loadings,
    aes(x = MDS1, y = MDS2 + 2, label = Species),
    color = "black"
  ) +
  scale_colour_brewer(palette = "Set2") +
  xlab(paste0(rownames(eigs)[1], " (", 
              round(eig_var[1] * 100, 2), 
              "% explained variance)")) +
  ylab(paste0(rownames(eigs)[2], " (", 
              round(eig_var[2] * 100, 2), 
              "% explained variance)")) +
  guides(color = guide_legend(title = "Species"),
         shape = guide_legend(title = "City")) + 
  ggtitle(paste0("Genitalia measurements PCA colored by city"))

# print it
print(pca_plot)

# save pdf
pdf_plot(pca_plot, paste0(base_dir, "pca.pdf"), 10, 7)

###
# get a neighbor-joining tree

# read fasta
named_coi <- read.FASTA(paste0(base_dir, "named_coi.fasta"))

# make a nj tree
tree <- nj(dist.dna(named_coi, model = "F84"))
tree <- root(tree, c("Ephestia_elutella_1"))

# write tree as is
write.nexus(tree, file = paste0(base_dir, "nj_tree.nex"))

# get names for the alignment
names_align <- read.table(paste0(base_dir, "names_align.tsv"), header = TRUE)

# open pdf
pdf(file = paste0(base_dir, "tree_plot.pdf"),
    width = 10, height = 7)

# plot tree without tip labels
plot(tree, show.tip.label = FALSE,
     x.lim = c(-0.2 * max(node.depth.edgelength(tree)), 
               max(node.depth.edgelength(tree))))

# factor for site
site_factor <- factor(ids$Site[match(names_align$ACC, ids$ACC)])

# make colors based on site
tip.cols <- adjustcolor(site_palette[site_factor], alpha = 0.5)

# get tip coordinates
coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# add colored circles at tips
points(coords$xx[1:coords$Ntip], coords$yy[1:coords$Ntip],
       pch = 21,                      
       bg  = tip.cols,                
       col = "black",                 
       cex = 1.5)

# find unique factor levels in the order you want
levs <- levels(site_factor)

# plot legend
par(xpd = TRUE)
legend("topleft",                     
       legend = levs,               
       pt.bg  = site_palette[levs], 
       pch    = 21,                 
       pt.cex = 1.5,                
       bty    = "n")     
par(xpd = FALSE)
dev.off()

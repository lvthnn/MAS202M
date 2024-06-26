############################################################
# PART 1: Data cleaning
############################################################

library(dplyr)
library(scales)

source("pair_plot.R")

set.seed(1)

islrblue <- "#3333b2"
islrred <- "#990000"
islrgreen <- "#009900"

########################
# READ IN DATA
########################

synthetic_dat <- read.table(
  "MAS202M_s24_synthetic_data_20240321.csv",
  sep = ",",
  head = TRUE
)

pair_plot(synthetic_dat[, -c(1, 2)], col = islrblue)

########################
# DATA CLEANING
########################

synthetic_dat <- synthetic_dat[, -1]

# data split across col3 and col6 -- mutate col3 and delete col6
synthetic_dat <- synthetic_dat |>
  mutate(col3 = ifelse(is.na(col3), col6, col3)) |>
  select(-col6)

# order the samples by increasing index
n <- nrow(synthetic_dat)
ord <- gsub("sample", "", synthetic_dat$id) |>
  as.numeric() |>
  order()
synthetic_dat <- synthetic_dat[ord, ]
rownames(synthetic_dat) <- 1:n

dup <- duplicated(synthetic_dat)
sum(dup) # 97 identical observations

# remove identical observations
synthetic <- synthetic_dat |> distinct(.keep_all = TRUE)

dup_ids <- which(duplicated(synthetic$id))
dup_ids_tmp <- synthetic[c(dup_ids, dup_ids - 1), ] |>
  arrange(id) |>
  filter(id != "")

dup_ids_drop <- which((rowSums(apply(dup_ids_tmp, 2, is.na)) == TRUE)) |>
  names() |>
  as.numeric()

synthetic <- synthetic[-dup_ids_drop, ] # 2 dropped (NA-identical)
synthetic <- synthetic[-1669, ] # drop identical obs. with no ID

# rename samples
n <- nrow(synthetic)
rownames(synthetic) <- 1:n
synthetic$id <- paste0("sample", 1:n)

########################
# IMPUTE MISSING VALUES
########################

synthetic <- synthetic |> mutate(
  across(col1:col5, function(c) ifelse(is.na(c), mean(c, na.rm = TRUE), c))
)

########################
# SUMMARY STATICS & PLOTS
########################

summary(synthetic)

pair_plot(synthetic[, -1])

############################################################
# PART 2: Clustering methods
############################################################

# pca for decision boundary?

pr_out <- prcomp(synthetic[, 2:6], scale = TRUE)

pr_var <- pr_out$sdev^2
pve <- pr_var / sum(pr_var)

biplot(
  pr_out,
  scale = 0,
  col = c(alpha(islrblue, 0.25), islrblue),
  cex = c(0.5, 1.25)
)
abline(h = 0, lty = 2, col = "#cccccc")
abline(v = 0, lty = 2, col = "#cccccc")

par(mfrow = c(1, 2))
plot(
  pve,
  type = "b",
  xlab = "Principal Component",
  ylab = "PVE",
  ylim = c(0, 1)
)

plot(
  cumsum(pve),
  type = "b",
  xlab = "Princap Component",
  ylab = "Cumulative PVE",
  ylim = c(0, 1)
)

# Hierarchical clustering - Out of the box (Average, Complete)
hclust <- hclust(dist(synthetic[, 2:6]), method = "average")
plot(hclust)

hclust_out <- cutree(hclust, 2)

pair_plot(synthetic[, 2:6], col = hclust_out + 2, cex = 0.2)

# Hierarchical clustering -- Scaled (Complete)
hclust_sc <- hclust(dist(scale(synthetic[, 2:6])), method = "complete")
plot(hclust_sc)

hclust_sc_out <- cutree(hclust_sc, 2)
hclust_sc_out <- ifelse(hclust_sc_out == 1, 2, 1)

pair_plot(synthetic[, 2:6], col = hclust_sc_out + 2, cex = 0.2)

# Hierarchical clustering -- PC1 & PC2 (Average)
hclust_pca <- hclust(dist(pr_out$x[, 1:2]), method = "complete")
plot(hclust_pca)

hclust_pca_out <- cutree(hclust_pca, 2)

pair_plot(synthetic[, 2:6], col = hclust_pca_out + 2)

# K Means -- Out of the box performance
kmeans <- kmeans(synthetic[, 2:6], centers = 2, nstart = 50)
kmeans$cluster <- ifelse(kmeans$cluster == 1, 2, 1)
pair_plot(synthetic[, 2:6], col = kmeans$cluster + 2)

# K Means -- Scaled data
kmeans_sc <- kmeans(scale(synthetic[, 2:6]), centers = 2, nstart = 50)
pair_plot(synthetic[, 2:6], col = kmeans_sc$cluster + 2)

# K Means -- PC1 & PC2
kmeans_pr <- kmeans(pr_out$x[, 1:2], centers = 2, nstart = 50)
pair_plot(synthetic[, 2:6], col = kmeans_pr$cluster + 2)

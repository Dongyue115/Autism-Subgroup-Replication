rm(list = ls())

library(NbClust)
library(proxy)
library(R.matlab)

cca_data <- readMat("para_com_out.mat")

all_data <- cca_data$para.com

cluster_index <- list()

for (i in 1:length(all_data)){
  replicate_data <- all_data[[i]][[1]][[4]]
  replicate_data <- as.data.frame(replicate_data)
  cca_rs_data <- scale(replicate_data)
  cca_rs_data_dis <- dist(cca_rs_data, method = "cosine")
  cluster_index[[i]] <- list(
    ch = NbClust(cca_rs_data, diss = cca_rs_data_dis, distance = NULL, method="average" ,min.nc = 2, max.nc = 7, index = "ch"),
    db = NbClust(cca_rs_data, diss = cca_rs_data_dis, distance = NULL, method="average",min.nc = 2, max.nc = 7, index = "db"))
}

ch_all <- matrix(nrow = 1000, ncol = 6)
db_all <- matrix(nrow = 1000, ncol = 6)

for (i in 1:length(cluster_index)){
  a <- cluster_index[[i]]$ch$All.index
  ch_all[i,] <- a
  b <- cluster_index[[i]]$db$All.index
  db_all[i,] <- b
  }

save(ch_all, file = "ch_all_zfc_1000.RData")
save(db_all, file = "db_all_zfc_1000.RData")

######plot########

load("ch_all_zfc_1000.RData")
load("db_all_zfc_1000.RData")
png("cluster/ch_all.png", width = 400*(300/72), height = 600*(300/72), res = 300)
par(mgp = c(2.5, 1, 0))
boxplot(ch_all,notch = TRUE,
        main = "CH index in 1000 replicates",
        xlab = "Number of clusters",
        ylab = "Calinski Harabasz Criterion",
        names = 2:7,
        cex.main=2.2,
        cex.lab=2.2,
        cex.axis=1.8)

dev.off()

png("cluster/db_all.png", width = 400*(300/72), height = 600*(300/72), res = 300)
par(mgp = c(2.5, 1, 0))
boxplot(db_all, notch = TRUE,
        main = "DB index in 1000 replicates",
        xlab = "Number of clusters",
        ylab = "Davies Bouldin Criterion",
        names = 2:7,
        cex.main=2.2,
        cex.lab=2.2,
        cex.axis=1.8)
dev.off()


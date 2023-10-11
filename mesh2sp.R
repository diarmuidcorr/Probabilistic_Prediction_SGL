## Libraries and Settings: 
library(pROC)
library(terra)
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(INLA)
inla.setOption(inla.mode="experimental")
library(geoR)
library(ggplot2)
library(gridExtra)
library(inlabru)
#> Loading required package: fmesher
library(fmesher)
library(ggplot2)

wd <- ("/mnt/Cryo-Chem/diarmuid/Paper3_Dev/")
setwd(wd)

load('Rdata/inla.spde.3conditions.boundary.mesh.RData')
bound = "Defined"

i <- 1

boundary_path <- paste0("/data/Boundary_Data/", bound ,"_Boundary/Data")
label_path <- paste0("/data/Boundary_Data/", bound ,"_Boundary/Labels/")
files_path <- paste0("/data/Boundary_Data/", bound ,"_Boundary/Files/")
data_rasters <- Sys.glob(file.path(paste0(wd, boundary_path), "*.tif"))


data_ras <- data_rasters[i]
label_ras <- file.path(paste0(wd, label_path, substr(data_rasters[i], 89, 160)))
data_ras <- rast(paste0(data_ras))
label_ras <- rast(paste0(label_ras))

names(data_ras) <- c('Band_01', 'Band_02', 'Band_03', 'Band_04', 'Band_05',
                     'Band_06', 'Band_07', 'Band_08', 'Band_09', 'Band_10',
                     'Band_11', 'Band_12', 'Band_13', 'Band_14', 'Band_15',
                     'Band_16', 'Band_17', 'Band_18', 'Band_19', 'Band_20',
                     'Band_21', 'Band_22', 'Band_23', 'Band_24', 'Band_25')

names(label_ras) <- c('Labels')

point_data <- as.data.frame(data_ras, 'SpatialPointsDataFrame')
point_label <- as.data.frame(label_ras, 'SpatialPointsDataFrame')

locations_ <- crds(label_ras)

locations_[,1] <- locations_[,1] - min(locations_[,1]) + 15
locations_[,2] <- locations_[,2] - min(locations_[,2]) + 15

r <- point_data$Band_16
g <- point_data$Band_15
b <- point_data$Band_14

r <- (r - min(r)) / (max(r) - min(r))
g <- (g - min(g)) / (max(g) - min(g))
b <- (b - min(b)) / (max(b) - min(b))

#r <- matrix(r, nrow = 192, ncol = 192)
#g <- matrix(g, nrow = 192, ncol = 192)
#b <- matrix(b, nrow = 192, ncol = 192)

dat.df <- data.frame(lon <- locations_[,1],
                     lat <- locations_[,2],
                     r <- r,
                     g <- g,
                     b <- b)


ggpl1 <- ggplot(data = dat.df, aes(x = lon, y =lat)) +
geom_raster(fill=rgb(r = r,
                     g = g,
                     b = b),
            #width=10, height=10,
            show.legend = FALSE) +
scale_fill_identity() +
gg(mesh_general,
   alpha = 0.5,
   edge.color = "azure4",
   edge.linewidth = 0.0000001,
   int.color = "black",
   ext.color = "black",
   ext.linewidth = 0.025) +
theme_void() +
scale_fill_identity() + 
coord_equal() +
#ggtitle("RGB with INLA Mesh") +
theme(plot.title = element_text(hjust = 0.5)) +
xlab("Longitude") +
ylab("Latitude")

ggsave(plot=ggpl1, filename="Rdata/mesh.pdf")

#pdf("Rdata/mesh.pdf") 
#plot(raster::as.raster(my_img))
#plot(mesh_general)

# Close the pdf file
#dev.off()
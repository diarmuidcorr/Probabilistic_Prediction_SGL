## This code is written in R and allows spatial data 
## processing and modeling using the INLA package in R.
## The script process geospatial data for three different 
## lake boundary conditions: "Defined," "Blurred," and "NoClear."
## Creates a consistent mesh determined for input data 192*192 pixels.
## It defines an SPDE model by constructinf a Matérn covariance model 
## on the SPDE mesh.
## Then INLA is used to create a separate model for each condition,
## with different spatial parameters for each.
## The range and standard deviation parameters for each condition's 
## model are extracted from the fitted models.
## New data is loaded and processed in a similar way to the training data,
## and the range and standard deviation are calculated.
## The distances between the new data's parameters and the reference 
## models' parameters are calculated, with the smallest distance used
## to identify the most appropriate model. 
## The selected model is then used to predict over the new data. 
## This code was compiled by Diarmuid Corr, Lancaster University.
## Contact dcorr103@gmail.com for any further information.

## Libraries and Settings: 
library(pROC)
library(terra)
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(ggplot2)
library(lattice)
library(caret)  # For data splitting and model evaluation
library(nnet)   # For multinomial logistic regression
library(moments)
library(abind)
library(geoR)
library(gridExtra)
library(entropy)
library(fractaldim)



wd <- ("/mnt/Cryo-Chem/diarmuid/Paper3_Dev/")
setwd(wd)


formula <- boundary_type ~ mean_reflectance + sd_reflectance + skewness + kurtosis + smoothness + contrast


formula <- as.formula(formula)

conditions <- list(
    Condition_A = list(label = "Defined"),
    Condition_B = list(label = "Blurred"),
    Condition_C = list(label = "NoClear"),
    Condition_D = list(label = "Testing")
)

train_data <-  list()
test_data <-  list()

est.LOGR.Boundary.Labels <- function(){
    
    a <- 1
    b <- 1
    for (condition  in names(conditions)) {
        bound <- conditions[[condition]]$label
        boundary_path <- paste0("/data/Boundary_Data/", bound ,"_Boundary/Data")
        label_path <- paste0("/data/Boundary_Data/", bound ,"_Boundary/Labels/")
        files_path <- paste0("/data/Boundary_Data/", bound ,"_Boundary/Files/")
        data_rasters <- Sys.glob(file.path(paste0(wd, boundary_path), "*.tif"))
        
        length.ras <- length(data_rasters)
        
        for(i in 1:length.ras) { 

            data_ras <- data_rasters[i]
            label_ras <- file.path(paste0(wd, label_path, 
                                          substr(data_rasters[i], 89, 160)))
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

            locations_[,1] <- ((locations_[,1] - min(locations_[,1]))) + 15
            locations_[,2] <- ((locations_[,2] - min(locations_[,2]))) + 15

            dat <- data.frame(x.coord = locations_[,1], 
                              y.coord = locations_[,2],
                              reflectance_01 = point_data$Band_01,                     
                              reflectance_02 = point_data$Band_02,                     
                              reflectance_03 = point_data$Band_03,                     
                              reflectance_04 = point_data$Band_04,                     
                              reflectance_05 = point_data$Band_05, 
                              reflectance_06 = point_data$Band_06, 
                              reflectance_07 = point_data$Band_07, 
                              reflectance_08 = point_data$Band_08, 
                              reflectance_09 = point_data$Band_09, 
                              reflectance_10 = point_data$Band_10, 
                              reflectance_11 = point_data$Band_11, 
                              reflectance_12 = point_data$Band_12, 
                              reflectance_13 = point_data$Band_13, 
                              reflectance_14 = point_data$Band_14,
                              reflectance_15 = point_data$Band_15, 
                              reflectance_16 = point_data$Band_16, 
                              reflectance_17 = point_data$Band_17,
                              reflectance_18 = point_data$Band_18, 
                              reflectance_19 = point_data$Band_19, 
                              reflectance_20 = point_data$Band_20,
                              reflectance_21 = point_data$Band_21, 
                              reflectance_22 = point_data$Band_22, 
                              reflectance_23 = point_data$Band_23, 
                              reflectance_24 = point_data$Band_24,
                              reflectance_25 = point_data$Band_25
                             )

            og_dat <- data.frame(x.coord = locations_[,1],
                                 y.coord = locations_[,2],
                                 Labels = point_label$Labels)
            
            shannon_entopy <- 0
            mean_reflectance <- 0
            sd_reflectance <- 0
            fractal_dim <- 0
            skewness <- 0
            kurtosis <- 0
            smoothness <- 0
            contrast <- 0
            
            
            for (j in 3:length(dat)){
                
                ref = dat[,j]
                shannon_entopy <- shannon_entopy + entropy.empirical(ref, unit="log2")
                mean_reflectance <- mean_reflectance + mean(ref)
                sd_reflectance <- sd_reflectance + sd(ref)
                fractal_dim <- fractal_dim + fd.estim.variogram(ref)$fd
                skewness <- skewness + skewness(ref)
                kurtosis <- kurtosis + kurtosis(ref)
                smoothness <- smoothness + (1 - 1/(1 + sd_reflectance**2))
                contrast <- contrast + ((sd_reflectance**2)/(kurtosis(ref)**0.25))
                
                }
            
            if(bound == 'Testing'){
                
                bound.dat <- data.frame(
                    
                    boundary_type = NA,
                    shannon_entopy = shannon_entopy,
                    mean_reflectance = mean_reflectance,
                    sd_reflectance = sd_reflectance,
                    fractal_dim = fractal_dim,
                    skewness = skewness,
                    kurtosis = kurtosis,
                    smoothness = smoothness,
                    contrast = contrast
                    
                    )
                test_data[[a]] <- bound.dat
                a <- a + 1
                
            } else {
                
                bound.dat <- data.frame(
                    
                    boundary_type = bound,
                    shannon_entopy = shannon_entopy,
                    mean_reflectance = mean_reflectance,
                    sd_reflectance = sd_reflectance,
                    fractal_dim = fractal_dim,
                    skewness = skewness,
                    kurtosis = kurtosis,
                    smoothness = smoothness,
                    contrast = contrast
                    
                    )
                
                
                train_data[[b]] <- bound.dat
                b <- b + 1
                
            }
                        
        }
        

    }
    test_data.df <- do.call("rbind",test_data) #combine all vectors into a matrix
    train_data.df <- do.call("rbind",train_data) #combine all vectors into a matrix
    
    # Setting the reference
    train_data.df$boundary_type <- relevel(factor(train_data.df$boundary_type), 
                                           ref = "NoClear")
    
    # Specify all levels of the response variable
    all_levels <- unique(train_data.df$boundary_type)
    
    
    # Train multinomial logistic regression model
    model <- multinom(boundary_type ~ ., data = train_data.df, levels = all_levels)
    #print(summary(model))
    
    predicted_probs <- predict(model, newdata = test_data.df, type = "probs")

    #print(predicted_probs)

    # Determine most likely boundary type for each image
    most_likely_boundary <- colnames(predicted_probs)[apply(predicted_probs, 1, which.max)]

    #print(most_likely_boundary)
    return(most_likely_boundary)
    

}

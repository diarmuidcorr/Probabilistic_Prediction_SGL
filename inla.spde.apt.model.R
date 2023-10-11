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
library(INLA)
inla.setOption(inla.mode="experimental")
library(geoR)
library(ggplot2)
library(gridExtra)


wd <- ("/mnt/Cryo-Chem/diarmuid/Paper3_Dev/")
setwd(wd)
source("scripts/multinomial.regression.boundary.model.R", 
       encoding = 'UTF-8')

#functions to define the most likely boundary using multinomial regression
most_likely_boundary <- est.LOGR.Boundary.Labels()
print(most_likely_boundary)
source("scripts/est.REG.INLA.R", encoding = 'UTF-8')

formula_Defined <- Labels ~ -1 + intercept_Defined + reflectance_01 + reflectance_02 + reflectance_03 + reflectance_04 + reflectance_05 + reflectance_06 + reflectance_07 + reflectance_08 + reflectance_09 + reflectance_10 + reflectance_11 + reflectance_12 + reflectance_13 + reflectance_14 + reflectance_15 + reflectance_16 + reflectance_17 + reflectance_18 + reflectance_19 + reflectance_20 + reflectance_21 + reflectance_22 + reflectance_23 + reflectance_24 + reflectance_25 + f(i_Defined, model = spde)

formula_Blurred <- Labels ~ -1 + intercept_Blurred + reflectance_01 + reflectance_02 + reflectance_03 + reflectance_04 + reflectance_05 + reflectance_06 + reflectance_07 + reflectance_08 + reflectance_09 + reflectance_10 + reflectance_11 + reflectance_12 + reflectance_13 + reflectance_14 + reflectance_15 + reflectance_16 + reflectance_17 + reflectance_18 + reflectance_19 + reflectance_20 + reflectance_21 + reflectance_22 + reflectance_23 + reflectance_24 + reflectance_25 + f(i_Blurred, model = spde)

formula_NoClear <- Labels ~ -1 + intercept_NoClear + reflectance_01 + reflectance_02 + reflectance_03 + reflectance_04 + reflectance_05 + reflectance_06 + reflectance_07 + reflectance_08 + reflectance_09 + reflectance_10 + reflectance_11 + reflectance_12 + reflectance_13 + reflectance_14 + reflectance_15 + reflectance_16 + reflectance_17 + reflectance_18 + reflectance_19 + reflectance_20 + reflectance_21 + reflectance_22 + reflectance_23 + reflectance_24 + reflectance_25 + f(i_NoClear, model = spde)

conditions <- list(
    Condition_A = list(label = "Defined", formula = formula_Defined),
    Condition_B = list(label = "Blurred", formula = formula_Blurred),
    Condition_C = list(label = "NoClear", formula = formula_NoClear),
    Condition_D = list(label = "Testing", formula = NA)
)

## Create mesh, SPDE, A matrix, and INLA stack for each model if needed:
if(FALSE){
    i=1

    bound = "Defined"

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


    step_size <- 6
    label_ras_gen <- aggregate(label_ras, step_size)
    locations_nodes <- crds(label_ras_gen)

    locations_nodes[,1] <- locations_nodes[,1] - min(locations_nodes[,1]) + 10
    locations_nodes[,2] <- locations_nodes[,2] - min(locations_nodes[,2]) + 10

    data_x <- data.frame(x = seq(10, 1930, by=step_size*10), 
                         y = 1930)
    #data_x <- rbind(data_x, data.frame(x = 192, y = 192))

    data_y <- data.frame(x = 1930, 
                         y = seq(10, 1920, by=step_size*10))

    locations_nodes <- rbind(locations_nodes, data_x)  
    locations_nodes <- rbind(locations_nodes, data_y)

    mesh_general <- inla.mesh.2d(loc = locations_nodes, #determine nodes
                                 loc.domain = locations_, #determine boundary/domain
                                 #offset = c(10, 20), #determines how strict domain is 
                                          #(in dom, out dom)
                                 #max.n=c(1000, 100), #Safeguard against large meshes.
                                 #max.n.strict=c(128000, 128000), #Don't build a huge mesh!
                                 max.edge=1500, #determines max edge length
                                 #cutoff = 10,
                                )

    #save(mesh_general, file = 'Rdata/inla.spde.3conditions.boundary.mesh.RData')
    #pdf("Rdata/mesh.pdf") 
    #plot(mesh_general)

    # Close the pdf file
    #dev.off()


    spde <- inla.spde2.pcmatern(
        mesh = mesh_general, ## mesh
        alpha = 2,   ## smoothness parameter
        prior.range=c(500, 0.8), ## P(practic.range<500)=0.5
        prior.sigma=c(10, 0.01))  ## P(sigma>4)=0.01

    #save(spde, file = 'Rdata/inla.spde.3conditions.boundary.spde.RData')

    #save(mesh_general, file = 'Rdata/mesh_general_8step.RData')
    #save(spde, file = 'Rdata/spde_8step.RData')

    mesh <- spde$mesh
    A_ <- inla.spde.make.A(mesh= mesh,
                           loc=locations_)
    #save(A_, file = 'Rdata/inla.spde.3conditions.boundary.A_mat.RData')
    

    models <- list()
    
    
    test.i <- 1
    for (condition  in names(conditions)) {
        bound <- conditions[[condition]]$label
        boundary_path <- paste0("/data/Boundary_Data/", bound ,"_Boundary/Data")
        label_path <- paste0("/data/Boundary_Data/", bound ,"_Boundary/Labels/")
        files_path <- paste0("/data/Boundary_Data/", bound ,"_Boundary/Files/")
        data_rasters <- Sys.glob(file.path(paste0(wd, boundary_path), "*.tif"))
        length.ras <- length(data_rasters)
     
        a <- 1
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
            
            #functions to create the INLA stacks
            if (bound == 'Defined') {
                
                stack.def <- est.LOGR.INLA(dat, og_dat, bound,
                                           spde, A_, i
                                          )
                
                if (a == 1) {
                    stack.def_ <- stack.def
                    a <- 2
                } else {
                    stack.def_ <- inla.stack(stack.def_, stack.def)
                }
            
            } else if (bound == 'Blurred') {
            
                stack.blur <- est.LOGR.INLA(dat, og_dat, bound,
                                            spde, A_, i
                                           )
                
                if (a == 1) {
                    stack.blur_ <- stack.blur
                    a <- 2
                } else {
                    stack.blur_ <- inla.stack(stack.blur_, stack.blur)
                } 
            
            } else if (bound == 'NoClear') {
                            
                stack.noc <- est.LOGR.INLA(dat, og_dat, bound,
                                           spde, A_, i
                                          )
                
                if (a == 1) {
                    stack.noc_ <- stack.noc
                    a <- 2
                } else {
                    stack.noc_ <- inla.stack(stack.noc_, stack.noc)
                }
            
            }
            
            if(bound == 'Testing'){
                bound_ <- most_likely_boundary[test.i]
                if (bound_ == 'Defined__'){
                
                    stack.def <- est.LOGR.INLA(dat, og_dat, bound_,
                                           spde, A_, i,
                                           labels_=FALSE
                                           )
                    stack.def_ <- inla.stack(stack.def_, stack.def)
                
                } else if (bound_ == 'Blurred__'){
                
                    stack.blur <- est.LOGR.INLA(dat, og_dat, bound_,
                                           spde, A_, i,
                                           labels_=FALSE
                                           )
                    stack.blur_ <- inla.stack(stack.blur_, stack.blur)
                
                } else if (bound_ == 'NoClear__'){
                
                    stack.noc <- est.LOGR.INLA(dat, og_dat, bound_,
                                           spde, A_, i,
                                           labels_=FALSE
                                           )
                    stack.noc_ <- inla.stack(stack.noc_, stack.noc)
                }
                
                stack.def <- est.LOGR.INLA(dat, og_dat, 'Defined',
                                           spde, A_, i,
                                           labels_=FALSE
                                           )
                stack.def_ <- inla.stack(stack.def_, stack.def)
                
                stack.blur <- est.LOGR.INLA(dat, og_dat, 'Blurred',
                                           spde, A_, i,
                                           labels_=FALSE
                                           )
                stack.blur_ <- inla.stack(stack.blur_, stack.blur)
                
                stack.noc <- est.LOGR.INLA(dat, og_dat, 'NoClear',
                                           spde, A_, i,
                                           labels_=FALSE
                                           )
                stack.noc_ <- inla.stack(stack.noc_, stack.noc)

                
                test.i <- test.i + 1
            
            }
        }
    }
    save(stack.def_, 
         file = paste0('Rdata/inla.spde.3conditions.boundary.DefinedData_.RData'))
    save(stack.blur_, 
         file = paste0('Rdata/inla.spde.3conditions.boundary.BlurredData_.RData'))
    save(stack.noc_, 
         file = paste0('Rdata/inla.spde.3conditions.boundary.NoClearData_.RData'))
    
    def.model <- inla(formula_Defined, family="binomial",
                      data = inla.stack.data(stack.def_, spde=spde),
                      control.predictor = list(A = inla.stack.A(stack.def_), 
                                               compute = TRUE, link=1),
                      control.compute=list(cpo=TRUE, config=TRUE, dic=TRUE,
                                           return.marginals.predictor=TRUE),
                      control.inla=list(strategy='gaussian',
                                        int.strategy='eb'),
                      num.threads=12)
    
    blur.model <- inla(formula_Blurred, family="binomial",
                       data = inla.stack.data(stack.blur_, spde=spde),
                       control.predictor = list(A = inla.stack.A(stack.blur_), 
                                                compute = TRUE, link=1),
                       control.compute=list(cpo=TRUE, config=TRUE, dic=TRUE,
                                            return.marginals.predictor=TRUE),
                       control.inla=list(strategy='gaussian',
                                         int.strategy='eb'),
                       num.threads=12)
    
    noc.model <- inla(formula_NoClear, family="binomial",
                      data = inla.stack.data(stack.noc_, spde=spde),
                      control.predictor = list(A = inla.stack.A(stack.noc_), 
                                               compute = TRUE, link=1),
                      control.compute=list(cpo=TRUE, config=TRUE, dic=TRUE,
                                           return.marginals.predictor=TRUE),
                      control.inla=list(strategy='gaussian',
                                        int.strategy='eb'),
                      num.threads=12)
    
    models[["Defined"]] <- def.model
    models[["Blurred"]] <- blur.model
    models[["NoClear"]] <- noc.model
    
    #save(models, file = 'Rdata/inla.spde.3boundary.conditions.models.RData')

} else {
    print('Model training already done, loading from memory instead')
    #load('Rdata/inla.spde.3boundary.conditions.models.RData')
    load('Rdata/inla.spde.3conditions.boundary.spde.RData')
    load('Rdata/inla.spde.3conditions.boundary.A_mat.RData')
    
    load('Rdata/inla.spde.3conditions.boundary.DefinedData_.RData')
    load('Rdata/inla.spde.3conditions.boundary.BlurredData_.RData')
    load('Rdata/inla.spde.3conditions.boundary.NoClearData_.RData')
    
    print('SPDEs loaded')
    
    def.model <- inla(formula_Defined, family="binomial",
                      data = inla.stack.data(stack.def_, spde=spde),
                      control.predictor = list(A = inla.stack.A(stack.def_), 
                                               compute = TRUE, link=1),
                      control.compute=list(cpo=TRUE, config=TRUE, dic=TRUE,
                                           return.marginals.predictor=TRUE),
                      control.inla=list(strategy='gaussian',
                                        int.strategy='eb'),
                      num.threads=12)
    
    print('def.model loaded')
    
    blur.model <- inla(formula_Blurred, family="binomial",
                       data = inla.stack.data(stack.blur_, spde=spde),
                       control.predictor = list(A = inla.stack.A(stack.blur_), 
                                                compute = TRUE, link=1),
                       control.compute=list(cpo=TRUE, config=TRUE, dic=TRUE,
                                            return.marginals.predictor=TRUE),
                       control.inla=list(strategy='gaussian',
                                         int.strategy='eb'),
                       num.threads=12)
    
    print('blur.model loaded')
    
    noc.model <- inla(formula_NoClear, family="binomial",
                      data = inla.stack.data(stack.noc_, spde=spde),
                      control.predictor = list(A = inla.stack.A(stack.noc_), 
                                               compute = TRUE, link=1),
                      control.compute=list(cpo=TRUE, config=TRUE, dic=TRUE,
                                           return.marginals.predictor=TRUE),
                      control.inla=list(strategy='gaussian',
                                        int.strategy='eb'),
                      num.threads=12)
    
    print('noc.model loaded')
    
}


boundary_path <- paste0("/data/Boundary_Data/Testing_Boundary/Data")
label_path <- paste0("/data/Boundary_Data/Testing_Boundary/Labels/")
files_path <- paste0("/data/Boundary_Data/Testing_Boundary/Files/")

data_rasters <- Sys.glob(file.path(paste0(wd, boundary_path), "*.tif"))

test_i <- 1
for (bound in most_likely_boundary) {
    
    data_ras <- data_rasters[test_i]
    label_ras <- file.path(paste0(wd, label_path,
                                  substr(data_rasters[test_i], 89, 160)))
    
    label_ras <- file.path(paste0(wd, label_path,
                                  substr(data_rasters[test_i], 89, 160)))

    data_ras <- rast(paste0(data_ras))
    label_ras <- rast(paste0(label_ras))
    source.srs <- crs(data_ras)
    source.ext <- ext(data_ras)

    names(data_ras) <- c('Band_01', 'Band_02', 'Band_03', 'Band_04', 'Band_05',
                         'Band_06', 'Band_07', 'Band_08', 'Band_09', 'Band_10',
                         'Band_11', 'Band_12', 'Band_13', 'Band_14', 'Band_15',
                         'Band_16', 'Band_17', 'Band_18', 'Band_19', 'Band_20',
                         'Band_21', 'Band_22', 'Band_23', 'Band_24', 'Band_25')

    names(label_ras) <- c('Labels')

    locations_new <- crds(label_ras)

    locations_new[,1] <- locations_new[,1] - min(locations_new[,1]) + 15
    locations_new[,2] <- locations_new[,2] - min(locations_new[,2]) + 15

    
    if (bound == 'Defined__') {
    
        id_allef <- inla.stack.index(stack.def_,
                                     paste0(bound, "pred_allef", test_i))$data
    
    } else if (bound == 'Blurred__') {
    
        id_allef <- inla.stack.index(stack.blur_, 
                                     paste0(bound, "pred_allef", test_i))$data
        
    } else if (bound == 'NoClear__') {
        
        id_allef <- inla.stack.index(stack.noc_, 
                                     paste0(bound, "pred_allef", test_i))$data
        
    }
    
    id_allef.def_ <- inla.stack.index(stack.def_,
                                      paste0('Defined', "pred_allef", test_i))$data
    id_allef.blur_ <- inla.stack.index(stack.blur_, 
                                       paste0('Blurred', "pred_allef", test_i))$data
    id_allef.noc_ <- inla.stack.index(stack.noc_,
                                      paste0('NoClear', "pred_allef", test_i))$data

    
    newdat <- as.data.frame(xyFromCell(data_ras, cell = 1:ncell(data_ras)))

    newdat$pred.def_ <- def.model$summary.fitted.values[id_allef.def_, "mean"]
    
    newdat$pred.blur_ <- blur.model$summary.fitted.values[id_allef.blur_, "mean"]
    
    newdat$pred.noc_ <- noc.model$summary.fitted.values[id_allef.noc_, "mean"]

    pred.df <- data.frame(lon = locations_new[,1],
                          lat = locations_new[,2],
                          x = crds(label_ras)[,1],
                          y = crds(label_ras)[,2],
                          pred = newdat$pred.def_)
    
    out_path <- file.path(paste0(wd, files_path, substr(data_rasters[test_i], 89, 156)))


    #convert to matrix
    mat.pred = as.matrix(pred.df)

    #define non-uniform spatial extent
    ext.both = terra::ext(source.ext)
    #create raster
    r.both = terra::rast(ext.both, nrow=192, ncol=192,
                         resolution=10, crs=source.srs)

    #populate raster
    x.pred = terra::rasterize(mat.pred[,3:4], r.both, mat.pred[,5])


    out_path.pred <- file.path(paste0(out_path, bound, 
                                      '_3Models_Prediction_ID_def_', test_i, '.tif'))

    writeRaster(x.pred, filename=out_path.pred, overwrite=TRUE)
    
    
    
    
    
    
    
    pred.df <- data.frame(lon = locations_new[,1],
                          lat = locations_new[,2],
                          x = crds(label_ras)[,1],
                          y = crds(label_ras)[,2],
                          pred = newdat$pred.blur_)


    #convert to matrix
    mat.pred = as.matrix(pred.df)

    #define non-uniform spatial extent
    ext.both = terra::ext(source.ext)
    #create raster
    r.both = terra::rast(ext.both, nrow=192, ncol=192,
                         resolution=10, crs=source.srs)

    #populate raster
    x.pred = terra::rasterize(mat.pred[,3:4], r.both, mat.pred[,5])


    out_path.pred <- file.path(paste0(out_path, bound, 
                                      '_3Models_Prediction_ID_blur_', test_i, '.tif'))

    writeRaster(x.pred, filename=out_path.pred, overwrite=TRUE)
    
    
    
    
    
    
    pred.df <- data.frame(lon = locations_new[,1],
                          lat = locations_new[,2],
                          x = crds(label_ras)[,1],
                          y = crds(label_ras)[,2],
                          pred = newdat$pred.noc_)


    #convert to matrix
    mat.pred = as.matrix(pred.df)

    #define non-uniform spatial extent
    ext.both = terra::ext(source.ext)
    #create raster
    r.both = terra::rast(ext.both, nrow=192, ncol=192,
                         resolution=10, crs=source.srs)

    #populate raster
    x.pred = terra::rasterize(mat.pred[,3:4], r.both, mat.pred[,5])


    out_path.pred <- file.path(paste0(out_path, bound, 
                                      '_3Models_Prediction_ID_noc_', test_i, '.tif'))

    writeRaster(x.pred, filename=out_path.pred, overwrite=TRUE)
    
    test_i <- test_i + 1
        
}

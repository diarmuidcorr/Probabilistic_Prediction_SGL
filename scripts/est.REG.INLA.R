library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(INLA)

est.LOGR.INLA <- function(dat, og_dat, bound,
                          spde, A_, i,
                          labels_=TRUE
                         ){
    
    #mesh.index = inla.spde.make.index(name="field",
     #                                 n.spde=mesh$n)
    
        
    if(labels_){
        y_factor <- og_dat$Labels
        tag.label <- paste0("est",bound,i)

        if(bound=='Defined'){
            print(tag.label)
            stack.est <- inla.stack(data=list(Labels = y_factor), #og_dat, 
                                    A=list(A_,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1
                                           ),
                                    effects=list(i_Defined = 1:spde$n.spde,
                                                 #i = 1:spde$n.spde,
                                                 intercept_Defined = rep(1, nrow(dat)),
                                                 reflectance_01 = dat$reflectance_01,
                                                 reflectance_02 = dat$reflectance_02,
                                                 reflectance_03 = dat$reflectance_03,
                                                 reflectance_04 = dat$reflectance_04,
                                                 reflectance_05 = dat$reflectance_05,
                                                 reflectance_06 = dat$reflectance_06,
                                                 reflectance_07 = dat$reflectance_07,
                                                 reflectance_08 = dat$reflectance_08,
                                                 reflectance_09 = dat$reflectance_09,
                                                 reflectance_10 = dat$reflectance_10,
                                                 reflectance_11 = dat$reflectance_11,
                                                 reflectance_12 = dat$reflectance_12,
                                                 reflectance_13 = dat$reflectance_13,
                                                 reflectance_14 = dat$reflectance_14,
                                                 reflectance_15 = dat$reflectance_15,
                                                 reflectance_16 = dat$reflectance_16,
                                                 reflectance_17 = dat$reflectance_17,
                                                 reflectance_18 = dat$reflectance_18,
                                                 reflectance_19 = dat$reflectance_19,
                                                 reflectance_20 = dat$reflectance_20,
                                                 reflectance_21 = dat$reflectance_21,
                                                 reflectance_22 = dat$reflectance_22,
                                                 reflectance_23 = dat$reflectance_23,
                                                 reflectance_24 = dat$reflectance_24,
                                                 reflectance_25 = dat$reflectance_25
                                                ),
                                    tag=tag.label)
        }

        if(bound=='Blurred'){
            print(tag.label)
            stack.est <- inla.stack(data=list(Labels = y_factor), #og_dat, 
                                    A=list(A_,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1
                                          ),
                                    effects=list(i_Blurred = 1:spde$n.spde,
                                                 #i_Defined_1 = 1:spde$n.spde,
                                                 #i = 1:spde$n.spde,
                                                 intercept_Blurred = rep(1, nrow(dat)),
                                                 reflectance_01 = dat$reflectance_01,
                                                 reflectance_02 = dat$reflectance_02,
                                                 reflectance_03 = dat$reflectance_03,
                                                 reflectance_04 = dat$reflectance_04,
                                                 reflectance_05 = dat$reflectance_05,
                                                 reflectance_06 = dat$reflectance_06,
                                                 reflectance_07 = dat$reflectance_07,
                                                 reflectance_08 = dat$reflectance_08,
                                                 reflectance_09 = dat$reflectance_09,
                                                 reflectance_10 = dat$reflectance_10,
                                                 reflectance_11 = dat$reflectance_11,
                                                 reflectance_12 = dat$reflectance_12,
                                                 reflectance_13 = dat$reflectance_13,
                                                 reflectance_14 = dat$reflectance_14,
                                                 reflectance_15 = dat$reflectance_15,
                                                 reflectance_16 = dat$reflectance_16,
                                                 reflectance_17 = dat$reflectance_17,
                                                 reflectance_18 = dat$reflectance_18,
                                                 reflectance_19 = dat$reflectance_19,
                                                 reflectance_20 = dat$reflectance_20,
                                                 reflectance_21 = dat$reflectance_21,
                                                 reflectance_22 = dat$reflectance_22,
                                                 reflectance_23 = dat$reflectance_23,
                                                 reflectance_24 = dat$reflectance_24,
                                                 reflectance_25 = dat$reflectance_25
                                                ),
                                    tag=tag.label)
            }

        if(bound=='NoClear'){
            print(tag.label)
            stack.est <- inla.stack(data=list(Labels = y_factor), #og_dat, 
                                    A=list(A_,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1
                                          ),
                                    effects=list(i_NoClear = 1:spde$n.spde,
                                                 #i_Blurred_2 = 1:spde$n.spde,
                                                 #i_NoClear_3 = 1:spde$n.spde,
                                                 #i = 1:spde$n.spde,
                                                 intercept_NoClear = rep(1, nrow(dat)),
                                                 reflectance_01 = dat$reflectance_01,
                                                 reflectance_02 = dat$reflectance_02,
                                                 reflectance_03 = dat$reflectance_03,
                                                 reflectance_04 = dat$reflectance_04,
                                                 reflectance_05 = dat$reflectance_05,
                                                 reflectance_06 = dat$reflectance_06,
                                                 reflectance_07 = dat$reflectance_07,
                                                 reflectance_08 = dat$reflectance_08,
                                                 reflectance_09 = dat$reflectance_09,
                                                 reflectance_10 = dat$reflectance_10,
                                                 reflectance_11 = dat$reflectance_11,
                                                 reflectance_12 = dat$reflectance_12,
                                                 reflectance_13 = dat$reflectance_13,
                                                 reflectance_14 = dat$reflectance_14,
                                                 reflectance_15 = dat$reflectance_15,
                                                 reflectance_16 = dat$reflectance_16,
                                                 reflectance_17 = dat$reflectance_17,
                                                 reflectance_18 = dat$reflectance_18,
                                                 reflectance_19 = dat$reflectance_19,
                                                 reflectance_20 = dat$reflectance_20,
                                                 reflectance_21 = dat$reflectance_21,
                                                 reflectance_22 = dat$reflectance_22,
                                                 reflectance_23 = dat$reflectance_23,
                                                 reflectance_24 = dat$reflectance_24,
                                                 reflectance_25 = dat$reflectance_25
                                                ),
                                    tag=tag.label)
        }
    } else {
        if(bound=='Defined'){
            y_factor <- NA
            tag.label <- paste0(bound, "pred_allef",i)
            print(tag.label)
            stack.est <- inla.stack(data=list(Labels = y_factor), #og_dat, 
                                    A=list(A_,
                                           1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1,
                                           1,1,1,1,1
                                          ),
                                    effects=list(i_Defined = 1:spde$n.spde,
                                                 intercept_Defined = rep(1, nrow(dat)),
                                                 reflectance_01 = dat$reflectance_01,
                                                 reflectance_02 = dat$reflectance_02,
                                                 reflectance_03 = dat$reflectance_03,
                                                 reflectance_04 = dat$reflectance_04,
                                                 reflectance_05 = dat$reflectance_05,
                                                 reflectance_06 = dat$reflectance_06,
                                                 reflectance_07 = dat$reflectance_07,
                                                 reflectance_08 = dat$reflectance_08,
                                                 reflectance_09 = dat$reflectance_09,
                                                 reflectance_10 = dat$reflectance_10,
                                                 reflectance_11 = dat$reflectance_11,
                                                 reflectance_12 = dat$reflectance_12,
                                                 reflectance_13 = dat$reflectance_13,
                                                 reflectance_14 = dat$reflectance_14,
                                                 reflectance_15 = dat$reflectance_15,
                                                 reflectance_16 = dat$reflectance_16,
                                                 reflectance_17 = dat$reflectance_17,
                                                 reflectance_18 = dat$reflectance_18,
                                                 reflectance_19 = dat$reflectance_19,
                                                 reflectance_20 = dat$reflectance_20,
                                                 reflectance_21 = dat$reflectance_21,
                                                 reflectance_22 = dat$reflectance_22,
                                                 reflectance_23 = dat$reflectance_23,
                                                 reflectance_24 = dat$reflectance_24,
                                                 reflectance_25 = dat$reflectance_25
                                                ),
                                    tag=tag.label)
        }


        if(bound=='Blurred'){
        y_factor <- NA
        tag.label <- paste0(bound, "pred_allef",i)
        print(tag.label)
        stack.est <- inla.stack(data=list(Labels = y_factor), #og_dat, 
                                A=list(A_,
                                       1,
                                       1,1,1,1,1,
                                       1,1,1,1,1,
                                       1,1,1,1,1,
                                       1,1,1,1,1,
                                       1,1,1,1,1
                                      ),
                                effects=list(i_Blurred = 1:spde$n.spde,
                                             intercept_Blurred = rep(1, nrow(dat)),
                                             reflectance_01 = dat$reflectance_01,
                                             reflectance_02 = dat$reflectance_02,
                                             reflectance_03 = dat$reflectance_03,
                                             reflectance_04 = dat$reflectance_04,
                                             reflectance_05 = dat$reflectance_05,
                                             reflectance_06 = dat$reflectance_06,
                                             reflectance_07 = dat$reflectance_07,
                                             reflectance_08 = dat$reflectance_08,
                                             reflectance_09 = dat$reflectance_09,
                                             reflectance_10 = dat$reflectance_10,
                                             reflectance_11 = dat$reflectance_11,
                                             reflectance_12 = dat$reflectance_12,
                                             reflectance_13 = dat$reflectance_13,
                                             reflectance_14 = dat$reflectance_14,
                                             reflectance_15 = dat$reflectance_15,
                                             reflectance_16 = dat$reflectance_16,
                                             reflectance_17 = dat$reflectance_17,
                                             reflectance_18 = dat$reflectance_18,
                                             reflectance_19 = dat$reflectance_19,
                                             reflectance_20 = dat$reflectance_20,
                                             reflectance_21 = dat$reflectance_21,
                                             reflectance_22 = dat$reflectance_22,
                                             reflectance_23 = dat$reflectance_23,
                                             reflectance_24 = dat$reflectance_24,
                                             reflectance_25 = dat$reflectance_25
                                            ),
                                tag=tag.label)
        }


        if(bound=='NoClear'){
        y_factor <- NA
        tag.label <- paste0(bound, "pred_allef",i)
        print(tag.label)
        stack.est <- inla.stack(data=list(Labels = y_factor), #og_dat, 
                                A=list(A_,
                                       1,
                                       1,1,1,1,1,
                                       1,1,1,1,1,
                                       1,1,1,1,1,
                                       1,1,1,1,1,
                                       1,1,1,1,1
                                      ),
                                effects=list(i_NoClear = 1:spde$n.spde,
                                             intercept_NoClear = rep(1, nrow(dat)),
                                             reflectance_01 = dat$reflectance_01,
                                             reflectance_02 = dat$reflectance_02,
                                             reflectance_03 = dat$reflectance_03,
                                             reflectance_04 = dat$reflectance_04,
                                             reflectance_05 = dat$reflectance_05,
                                             reflectance_06 = dat$reflectance_06,
                                             reflectance_07 = dat$reflectance_07,
                                             reflectance_08 = dat$reflectance_08,
                                             reflectance_09 = dat$reflectance_09,
                                             reflectance_10 = dat$reflectance_10,
                                             reflectance_11 = dat$reflectance_11,
                                             reflectance_12 = dat$reflectance_12,
                                             reflectance_13 = dat$reflectance_13,
                                             reflectance_14 = dat$reflectance_14,
                                             reflectance_15 = dat$reflectance_15,
                                             reflectance_16 = dat$reflectance_16,
                                             reflectance_17 = dat$reflectance_17,
                                             reflectance_18 = dat$reflectance_18,
                                             reflectance_19 = dat$reflectance_19,
                                             reflectance_20 = dat$reflectance_20,
                                             reflectance_21 = dat$reflectance_21,
                                             reflectance_22 = dat$reflectance_22,
                                             reflectance_23 = dat$reflectance_23,
                                             reflectance_24 = dat$reflectance_24,
                                             reflectance_25 = dat$reflectance_25
                                            ),
                                tag=tag.label)
        }
    }
    
    return(stack.est)
}


    
    
    

# Probabilistic_Prediction_SGL
A probabilistic prediction of supraglacial lakes on the southwest Greenland Ice Sheet.

This code is written in R and allows spatial data processing and modeling using the INLA package in R. 
The script process geospatial data for three different lake boundary conditions: "Defined," "Blurred," and "NoClear." 
Creates a consistent mesh determined for input image data 192*192 pixels.
It defines an SPDE model by constructinf a Matérn covariance model on the SPDE mesh.
Then INLA is used to create a separate model for each condition, with different spatial parameters for each.
A multinomial logistic regression is used to discern one of three lake border conditions (well defined, blurred, and without a clear border) within the images.
The border condition predicted by the multinomial regression on the new data is used to identify the most appropriate model. 
The selected model is then used to predict over the new data. 
This code was compiled by Diarmuid Corr, Lancaster University.
Contact dcorr103@gmail.com for any further information.

The main script, inla.spde.apt.model.R, calls on the multinomial.regression.boundary.model.R script to make a lake border type prediction based on the outcome of multinomial regression, and the est.REG.INLA.R script to create the stacks for the INLA-SPDE model. 

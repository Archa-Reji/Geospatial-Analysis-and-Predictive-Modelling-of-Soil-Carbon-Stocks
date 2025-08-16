#set wd
# let us start with a legacy data
# Soil C density to a given depth
library(ithir)
data(oneProfile)
str(oneProfile)

#We can fit a spline to the max soil depth or any depth less than that
# Fitting a ea_Spline
#ea_spline function will predict a continuous function 
#from the top of the soil profile to the maximum soil depth,
#such that it will interpolate values both within the observed depths and between the depths where there is no observation
#parametrs: lam (smoothness) and d
eaFit <- ea_spline(oneProfile, var.name = "C.kg.m3.",
                   d = t(c(0, 5, 15, 30, 60, 100, 200)), lam = 0.1, vlow = 0,
                   show.progress = FALSE)
str(eaFit)

#Plotting the outputs of ea spline
par(mar = c(2, 2, 2, 2))
par(mfrow = c(3, 1))
for (i in 1:3) {
  plot_ea_spline(splineOuts = eaFit, d = t(c(0, 5, 15, 30, 60,
                                             100, 200)), maxd = 200, type = i, plot.which = 1,
                 label = "carbon density") }




#Environmental covariates
data(edgeroi_splineCarbon)
str(edgeroi_splineCarbon)

#Accessing environmental covariates
#elevation,terrain wetness index (twi), gamma radiometric potassium (radK), and Landsat 7 spectral bands 3 and 4 (landsat_b3 and landsat_b4).
# Common CRS, dimension and resolution (Ideal for DSM)
# When variable CRS, dim and resolution use "projectRaster" and "resample".
data(edgeroiCovariates)
library(raster)
elevation
twi
radK
landsat_b3
landsat_b4

#plot raster

plot(elevation, main = "Edgeroi elvation map with overlayed
point locations")

coordinates(edgeroi_splineCarbon) <- ~east + north


## plot points


plot(edgeroi_splineCarbon, add = T)


#stacking covariates with common CRS,resolution and dimension

covStack <- stack(elevation, twi, radK, landsat_b3, landsat_b4)
covStack

#intersecting between the soil observations and covariate layers

DSM_data <- extract(covStack, edgeroi_splineCarbon, sp = 1,
                    method = "simple")


#export the soil and covariate data intersect object to file for later use.
#First we convert the spatial object to a dataframe, then export as .csv file

DSM_data <- as.data.frame(DSM_data)
write.table(DSM_data, "edgeroiSoilCovariates_C.TXT",
            col.names = T, row.names = FALSE, sep = ",")





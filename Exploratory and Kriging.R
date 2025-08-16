
#Some Exploratory Data Analysis
# work with the same wd
#
edge.dat <- read.table("edgeroiSoilCovariates_C.TXT",
                       sep = ",", header = T)
str(edge.dat)

#summary statistics of SOC
round(summary(edge.dat$X0.5.cm), 1)


install.packages("fBasics")
install.packages("nortest")
library(fBasics)
library(nortest)

# skewness
sampleSKEW(edge.dat$X0.5.cm)

# kurtosis
sampleKURT(edge.dat$X0.5.cm)

#Normality test
ad.test(edge.dat$X0.5.cm)

#plot untransformed data
par(mfrow = c(1, 2))
hist(edge.dat$X0.5.cm)
qqnorm(edge.dat$X0.5.cm, plot.it = TRUE, pch = 4, cex = 0.7)
qqline(edge.dat$X0.5.cm, col = "red", lwd = 2)


#data transformation
sampleSKEW(log(edge.dat$X0.5.cm))
sampleKURT(log(edge.dat$X0.5.cm))
ad.test(log(edge.dat$X0.5.cm))


#plot transformed data
par(mfrow = c(1, 2))
hist(log(edge.dat$X0.5.cm))
qqnorm(log(edge.dat$X0.5.cm), plot.it = TRUE, pch = 4, cex = 0.7)
qqline(log(edge.dat$X0.5.cm), col = "red", lwd = 2)


#Spatial interpolation
# First step: prepare a grid of points upon which the interpolators will be used.
# can be done by extracting the coordinates from either of the 90m resolution rasters we have for the Edgeroi
install.packages("gstat")
library (raster)
library(ithir)

data(edgeroiCovariates)

tempD <- data.frame(cellNos = seq(1:ncell(elevation)))
tempD$vals <- getValues(elevation) # get the pixels with value associated (discards NA)
tempD <- tempD[complete.cases(tempD), ]
cellNos <- c(tempD$cellNos)
gXY <- data.frame(xyFromCell(elevation, cellNos, spatial = FALSE))

#IDI/ IDW interpolation

library(gstat)
names(edge.dat)[2:3] <- c("x", "y")  #specify the observed data with thier locations
IDW.pred <- idw(log(edge.dat$X0.5.cm) ~ 1, locations = ~x + y,
                data = edge.dat, newdata = gXY, idp = 2) #specify the spatial points we want to interpolate onto

#Remember: idp=inverse distance weighting parameter, default value is 2, but you can play with it

#Plot
IDW.raster.p <- rasterFromXYZ(as.data.frame(IDW.pred[, 1:3]))
plot(IDW.raster.p)


#kriging

vgm1 <- variogram(log(X0.5.cm) ~ 1, ~x + y, edge.dat, width = 400,
                  cressie = TRUE, cutoff = 10000)
mod <- vgm(psill = var(log(edge.dat$X0.5.cm)),
           "Exp", range = 5000, nugget = 0)
model_1 <- fit.variogram(vgm1, mod)
model_1

#plot variogram
plot(vgm1, model = model_1)

krig.pred <- krige(log(edge.dat$X0.5.cm) ~ 1, locations = ~x + y,
                   data = edge.dat, newdata = gXY, model = model_1)

par(mfrow = c(2, 1))
krig.raster.p <- rasterFromXYZ(as.data.frame(krig.pred[, 1:3]))
krig.raster.var <- rasterFromXYZ(as.data.frame
                                 (krig.pred[, c(1:2, 4)]))
plot(krig.raster.p, main = "ordinary kriging predictions")
plot(krig.raster.var, main = "ordinary kriging variance")

# If we want to know the correlation between target property and covariates

edge.dat$logC_0_5 <- log(edge.dat$X0.5.cm)
names(edge.dat)

cor(edge.dat[, c("elevation", "twi", "radK",
                 "landsat_b3", "landsat_b4")], edge.dat[, "logC_0_5"])

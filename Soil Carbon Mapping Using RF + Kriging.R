library(ithir)
library(MASS)
library(Cubist)
library(gstat)
library(ithir)
library(raster)
library(randomForest)

# point data
data(edgeroi_splineCarbon)
names(edgeroi_splineCarbon)[2:3] <- c("x", "y")
# natural log transform
edgeroi_splineCarbon$log_cStock <- log(edgeroi_splineCarbon
                                          $X15.30.cm)
# grids
data(edgeroiCovariates)
coordinates(edgeroi_splineCarbon) <- ~x + y
# stack the rasters
covStack <- stack(elevation, twi, radK, landsat_b3, landsat_b4)
# extract
DSM_data <- extract(covStack, edgeroi_splineCarbon, sp = 1,
                    method = "simple")
DSM_data <- as.data.frame(DSM_data)
set.seed(123)
training <- sample(nrow(DSM_data), 0.7 * nrow(DSM_data))
mDat <- DSM_data[training, ]
# fit the model
edge.RF.Exp <- randomForest(log_cStock ~ elevation + twi + radK +
                              landsat_b3 + landsat_b4, data = mDat,
                            importance = TRUE, ntree = 1000)
print(edge.RF.Exp)

mDat$residual <- mDat$log_cStock - predict(edge.RF.Exp,
                                              newdata = mDat)
mean(mDat$residual)


coordinates(mDat) <- ~x + y
crs(mDat) <- "+proj=utm +zone=55 +south +ellps=WGS84 +datum=WGS84
+units=m +no_defs"
vgm1 <- variogram(residual ~ 1, mDat, width = 250, cressie = TRUE,
                  cutoff = 10000)
mod <- vgm(psill = var(mDat$residual), "Sph", range = 5000,
           nugget = 0)
model_1 <- fit.variogram(vgm1, mod)
model_1

# Residual kriging model
gRK <- gstat(NULL, "RKresidual", residual ~ 1, mDat,
             model = model_1)
crs(mDat)
#External Validation

RandomForest.pred.V <- predict(edge.RF.Exp, newdata = DSM_data[training, ])
# RF model with residual variogram
vDat <- DSM_data[-training, ]
coordinates(vDat) <- ~x + y
crs(vDat) <- "+proj=utm +zone=55 +south +ellps=WGS84 +datum=WGS84
+units=m +no_defs"
# make the residual predictions
RK.preds.V <- as.data.frame(krige(residual ~ 1, mDat, model = model_1,
                                  newdata = vDat))
## [using ordinary kriging]
# Sum the two components together
RK.preds.fin <- RandomForest.pred.V + RK.preds.V[, 3]
# validation 
goof(observed = DSM_data$log_cStock[-training],
     predicted = RandomForest.pred.V)
# validation regression kriging with RF model
goof(observed = DSM_data$log_cStock[-training],
     predicted = RK.preds.fin)
par(mfrow = c(3, 1))
crs(covStack)= crs(mDat)
map.RK1 <- predict(covStack, edge.RF.Exp,
                   filename = "cStock_15_30_RF_RK.tif",
                   format = "GTiff", datatype = "FLT4S", overwrite = TRUE)
plot(map.RK1, main = "Random Forest model predicted 15-30cm log
carbon stocks (15-30cm)")
map.RK2 <- interpolate(covStack, gRK, xyOnly = TRUE, index = 1,
                       filename = "cStock_15_30_residualRK.tif", format = "GTiff",
                       datatype = "FLT4S", overwrite = TRUE)

plot(map.RK2, main = "Kriged residual")
pred.stack <- stack(map.RK1, map.RK2)
map.RK3 <- calc(pred.stack, fun = sum,
                filename = "cStock_15_30_finalPredRK.tif",
                format = "GTiff", progress = "text", overwrite = T)
plot(map.RK3,
     main = "Regression kriging predicted 15-30cm log
carbon stocks (15-30cm)")




#depth 30- 60 cm

# point data
data(edgeroi_splineCarbon)
names(edgeroi_splineCarbon)[2:3] <- c("x", "y")
# natural log transform
edgeroi_splineCarbon$log_cStock <- log(edgeroi_splineCarbon
                                          $X30.60.cm)
# grids
data(edgeroiCovariates)
coordinates(edgeroi_splineCarbon) <- ~x + y
# stack the rasters
covStack <- stack(elevation, twi, radK, landsat_b3, landsat_b4)
# extract
DSM_data <- extract(covStack, edgeroi_splineCarbon, sp = 1,
                    method = "simple")
DSM_data <- as.data.frame(DSM_data)
set.seed(123)
training <- sample(nrow(DSM_data), 0.7 * nrow(DSM_data))
mDat <- DSM_data[training, ]
# fit the model

edge.RF.Exp <- randomForest(log_cStock ~ elevation + twi + radK +
                              landsat_b3 + landsat_b4, data = mDat,
                            importance = TRUE, ntree = 1000)
print(edge.RF.Exp)

mDat$residual <- mDat$log_cStock - predict(edge.RF.Exp,
                                           newdata = mDat)
mean(mDat$residual)


coordinates(mDat) <- ~x + y
crs(mDat) <- "+proj=utm +zone=55 +south +ellps=WGS84 +datum=WGS84
+units=m +no_defs"
vgm1 <- variogram(residual ~ 1, mDat, width = 250, cressie = TRUE,
                  cutoff = 10000)
mod <- vgm(psill = var(mDat$residual), "Sph", range = 5000,
           nugget = 0)
model_1 <- fit.variogram(vgm1, mod)
model_1
# Residual kriging model
gRK <- gstat(NULL, "RKresidual", residual ~ 1, mDat,
             model = model_1)
crs(mDat)
#External Validation

RandomForest.pred.V <- predict(edge.RF.Exp, newdata = DSM_data[-training, ])
#RF model with residual variogram
vDat <- DSM_data[-training, ]
coordinates(vDat) <- ~x + y
crs(vDat) <- "+proj=utm +zone=55 +south +ellps=WGS84 +datum=WGS84
+units=m +no_defs"
# make the residual predictions
RK.preds.V <- as.data.frame(krige(residual ~ 1, mDat, model = model_1,
                                  newdata = vDat))
## [using ordinary kriging]
# Sum the two components together
RK.preds.fin <- RandomForest.pred.V + RK.preds.V[, 3]
# validation
goof(observed = DSM_data$log_cStock[-training],
     predicted = RandomForest.pred.V)
# validation regression kriging with RF model
goof(observed = DSM_data$log_cStock[-training],
     predicted = RK.preds.fin)
par(mfrow = c(3, 1))
crs(covStack)= crs(mDat)
map.RK1 <- predict(covStack, edge.RF.Exp,
                   filename = "cStock_30_60_RF_RK.tif",
                   format = "GTiff", datatype = "FLT4S", overwrite = TRUE)
plot(map.RK1, main = "Random Forest model predicted 30-60cm log
carbon stocks (30-60cm)")
map.RK2 <- interpolate(covStack, gRK, xyOnly = TRUE, index = 1,
                       filename = "cStock_30_60_residualRK.tif", format = "GTiff",
                       datatype = "FLT4S", overwrite = TRUE)

plot(map.RK2, main = "Kriged residual")
pred.stack <- stack(map.RK1, map.RK2)
map.RK3 <- calc(pred.stack, fun = sum,
                filename = "cStock_30_60_finalPredRK.tif",
                format = "GTiff", progress = "text", overwrite = T)
plot(map.RK3,
     main = "Regression kriging predicted 30-60cm log
carbon stocks (30-60cm)")




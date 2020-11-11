library(ncdf4)
library(viridis)
library(raster)
library(rgdal)
library(sf)

# Load ERA5 monthly averaged data on single levels:
# Total precipitation, total downward sw radiation at the surface, and direct downward sw
precip_brick <- brick("C:/Russell/Projects/Amazon_GOME-2/Data/ERA5_Monthly_Mean_Precip_SWRad_2007-2017.nc", varname = "tp")
sw_brick <- brick("C:/Russell/Projects/Amazon_GOME-2/Data/ERA5_Monthly_Mean_Precip_SWRad_2007-2017.nc", varname = "msdwswrf")
direct_sw_brick <- brick("C:/Russell/Projects/Amazon_GOME-2/Data/ERA5_Monthly_Mean_Precip_SWRad_2007-2017.nc", varname = "msdrswrf")

# Amazon shapefile
amazon <- shapefile("C:/Russell/Projects/Amazon_GOME-2/Data/Amazon_poly.shp")

# Moist Forest (2000+ MAP)
moist_forest <- raster("C:/Russell/Projects/Amazon_GOME-2/Data/Amazon_Consistent_Forest80_2000+.tif")
moist_forest <- projectRaster(moist_forest, crs = "+proj=longlat +datum=WGS84 +no_defs")
moist_forest_poly <- rasterToPolygons(moist_forest)
moist_forest_poly <- crop(moist_forest_poly, amazon)

# ENSO Index MEI v2 data
mei_data <- read.csv("C:/Russell/Projects/Amazon_GOME-2/Data/ENSO_Monthly_2007_MEIv2.csv", header = TRUE) # Load the file

# Clip to Amazon
precip_brick_a <- crop(precip_brick, amazon)
precip_brick_a <- mask(precip_brick_a, amazon)
sw_brick_a <- crop(sw_brick, amazon)
sw_brick_a <- mask(sw_brick_a, amazon)
direct_sw_brick_a <- crop(direct_sw_brick, amazon)
direct_sw_brick_a <- mask(direct_sw_brick_a, amazon)
# Clip to Moist Forest
precip_brick_mf <- crop(precip_brick, moist_forest_poly)
precip_brick_mf <- mask(precip_brick_mf, moist_forest_poly)
sw_brick_mf <- crop(sw_brick, moist_forest_poly)
sw_brick_mf <- mask(sw_brick_mf, moist_forest_poly)
direct_sw_brick_mf <- crop(direct_sw_brick, moist_forest_poly)
direct_sw_brick_mf <- mask(direct_sw_brick_mf, moist_forest_poly)

# Make diffuse brick
diffuse_sw_brick_a <- sw_brick_a - direct_sw_brick_a
diffuse_sw_brick_mf <- sw_brick_mf - direct_sw_brick_mf

# Calc means
precip_mean_a <- calc(precip_brick_a, mean)
sw_mean_a <- calc(sw_brick_a, mean)
direct_sw_mean_a <- calc(direct_sw_brick_a, mean)
diffuse_sw_mean_a <- calc(diffuse_sw_brick_a, mean)
# Moist Forest
precip_mean_mf <- calc(precip_brick_mf, mean)
sw_mean_mf <- calc(sw_brick_mf, mean)
direct_sw_mean_mf <- calc(direct_sw_brick_mf, mean)
diffuse_sw_mean_mf <- calc(diffuse_sw_brick_mf, mean)

# Calc sd
precip_sd_a <- calc(precip_brick_a, sd)
sw_sd_a <- calc(sw_brick_a, sd)
direct_sw_sd_a <- calc(direct_sw_brick_a, sd)
diffuse_sw_sd_a <- calc(diffuse_sw_brick_a, sd)
# Moist Forest
precip_sd_mf <- calc(precip_brick_mf, sd)
sw_sd_mf <- calc(sw_brick_mf, sd)
direct_sw_sd_mf <- calc(direct_sw_brick_mf, sd)
diffuse_sw_sd_mf <- calc(diffuse_sw_brick_mf, sd)

# Monthly Standardized Anomalies
precip_anom_a <- (precip_brick_a - precip_mean_a) / precip_sd_a
sw_anom_a <- (sw_brick_a - sw_mean_a) / sw_sd_a
direct_sw_anom_a <- (direct_sw_brick_a - direct_sw_mean_a) / direct_sw_sd_a
diffuse_sw_anom_a <- (diffuse_sw_brick_a - diffuse_sw_mean_a) / diffuse_sw_sd_a
# Moist Forest
precip_anom_mf <- (precip_brick_mf - precip_mean_mf) / precip_sd_mf
sw_anom_mf <- (sw_brick_mf - sw_mean_mf) / sw_sd_mf
direct_sw_anom_mf <- (direct_sw_brick_mf - direct_sw_mean_mf) / direct_sw_sd_mf
diffuse_sw_anom_mf <- (diffuse_sw_brick_mf - diffuse_sw_mean_mf) / diffuse_sw_sd_mf

# Basin scale timeseries of anomalies
precip_basin_anom_timeseries <- vector()
sw_basin_anom_timeseries <- vector()
direct_sw_basin_anom_timeseries <- vector()
diffuse_sw_basin_anom_timeseries <- vector()
for (i in 1:(nlayers(precip_anom_a))){
    precip_basin_anom_timeseries <- c(precip_basin_anom_timeseries, cellStats(precip_anom_a[[i]], stat = 'mean', na.rm = TRUE))
    sw_basin_anom_timeseries <- c(sw_basin_anom_timeseries, cellStats(sw_anom_a[[i]], stat = 'mean', na.rm = TRUE))
    direct_sw_basin_anom_timeseries <- c(direct_sw_basin_anom_timeseries, cellStats(direct_sw_anom_a[[i]], stat = 'mean', na.rm = TRUE))
    diffuse_sw_basin_anom_timeseries <- c(diffuse_sw_basin_anom_timeseries, cellStats(diffuse_sw_anom_a[[i]], stat = 'mean', na.rm = TRUE))
}
# Moist Forest timeseries of anomalies
precip_mf_anom_timeseries <- vector()
sw_mf_anom_timeseries <- vector()
direct_sw_mf_anom_timeseries <- vector()
diffuse_sw_mf_anom_timeseries <- vector()
for (i in 1:(nlayers(precip_anom_a))){
    precip_mf_anom_timeseries <- c(precip_mf_anom_timeseries, cellStats(precip_anom_mf[[i]], stat = 'mean', na.rm = TRUE))
    sw_mf_anom_timeseries <- c(sw_mf_anom_timeseries, cellStats(sw_anom_mf[[i]], stat = 'mean', na.rm = TRUE))
    direct_sw_mf_anom_timeseries <- c(direct_sw_mf_anom_timeseries, cellStats(direct_sw_anom_mf[[i]], stat = 'mean', na.rm = TRUE))
    diffuse_sw_mf_anom_timeseries <- c(diffuse_sw_mf_anom_timeseries, cellStats(diffuse_sw_anom_mf[[i]], stat = 'mean', na.rm = TRUE))
}


# Plot basin Scale
pdf("C:/Russell/Projects/Amazon_GOME-2/Figures/Precip_PARtoc_Anomalies.pdf", width=7.5, height=5, compress=FALSE)

par(mfrow=c(1,1),oma=c(2.5,2.5,0.5,2.75))
xaxis <- c(1:132)
names <- c("2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017")
ENSO <- mei_data$MEI
ENSO <- head(ENSO, -12)

# Axes margin
op <- par(mar = c(0, 0, 0.85, 0.5))

# Plot MEI
barplot(ENSO, col = ifelse(ENSO > 0.5, rgb(255,191,191,max=255),
                           ifelse(ENSO < -0.5, rgb(191,191,225,max=255), rgb(235,235,235,max=255))), ylim=c(-3,3), axes=F, tck=F,
        xpd=F, mgp=c(3,0.3,0), ann=FALSE, xaxs = "i", border = NA, space = 0)
abline(h = 0)

# Add anomaly lines
par(new=T)
bar <- barplot(precip_basin_anom_timeseries, axes = FALSE, ylim = c(-2.5, 2.5), col = NA, border = NA, ann = FALSE, xaxs = "i", space = 0)
lines(bar, precip_basin_anom_timeseries, type = "l", col = "#045a8d", lwd = 1.5, ylim = c(-2.5, 2.5), ann = FALSE)
lines(bar, sw_basin_anom_timeseries, type="l", col="#b30000", lwd=1.5, ylim=c(-2.5, 2.5), ann=FALSE)
# lines(bar, direct_sw_basin_anom_timeseries, type="l", lty=2, col="#ef6548", lwd=1.5, ylim=c(-2.5, 2.5), ann=FALSE)
# lines(bar, diffuse_sw_basin_anom_timeseries, type="l", lty=2, col="#3690c0", lwd=1.5, ylim=c(-2.5, 2.5), ann=FALSE)

# x-axis labels, ticks, and title
axis(side=1, labels=names, tck=0.03, cex.axis=0.85, mgp=c(3,0.3,0), at=c(seq(from=0,to=120,by=12)))
# y-axis labels, ticks, and title
axis(side=4, mgp=c(3,0.3,0), las=1, tck=0.03, at=c(seq(from=-2.5,to=2.5,by=1)))
axis(side=4, mgp=c(3,0.3,0), las=1, tck=0.02, at=c(seq(from=-2,to=2,by=1)), labels = FALSE)
axis(side=2, tck = 0.03, las=1, mgp=c(3,0.3,0), at=c(seq(from=-3,to=3,by=1)))
axis(side=2, tck = 0.02, labels=FALSE, c(seq(from=-2.5,to=2.5,by=1)))
# Build the legend
# legend("topleft", lty=c(1,1,2,2), lwd = c(1.5, 1.5, 1.5, 1.5), horiz=FALSE, col=c("#045a8d", "#b30000", "#ef6548", "#3690c0"), c("Precipitation", expression(paste("PAR"['TOC']*"")),
#        expression(paste("Direct PAR"['TOC']*"")), expression(paste("Diffuse PAR"['TOC']*""))), bty="n", cex=0.85, pt.cex=0.85)
# legend("topleft", lty=c(1), lwd = c(1.5), horiz=FALSE, col=c("#045a8d"), c("Precipitation"), bty="n", cex=0.85, pt.cex=0.85)
legend("topright", lty=c(1,1), lwd = c(1.5,1.5), horiz=FALSE, col=c("#045a8d", "#b30000"), c("Precipitation", expression(paste("PAR"['TOC']*""))), bty="n", cex=0.85, pt.cex=0.85)
legend("topleft", fill = c(rgb(255,191,191,max=255), rgb(235,235,235,max=255), rgb(191,191,225,max=255)), c("Warm Phase", "Neutral", "Cool Phase"), bty="n", cex=0.85, pt.cex=0.85)
box()

mtext("Standardized Anomaly", side=4, line=2)
mtext("Multivariate ENSO Index Version 2", side=2, line=1.25)
mtext(expression(paste("Year")), 1, 1.5, outer=TRUE)
mtext(expression(paste("Amazon Basin")), 3, -0.5, outer=TRUE)

dev.off()


# Plot Moist Forest
pdf("C:/Russell/Projects/Amazon_GOME-2/Figures/Precip_PARtoc_Anomalies_MF.pdf", width=7.5, height=5, compress=FALSE)

par(mfrow=c(1,1),oma=c(2.5,2.5,0,2.5))
xaxis <- c(1:132)
names <- c("2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017")
ENSO <- mei_data$MEI
ENSO <- head(ENSO, -12)

# Axes margin
op <- par(mar = c(0, 0, 0.85, 0.5))

# Plot MEI
barplot(ENSO, col = ifelse(ENSO > 0.5, rgb(255,191,191,max=255),
                           ifelse(ENSO < -0.5, rgb(191,191,225,max=255), rgb(235,235,235,max=255))), ylim=c(-3,3), axes=F, tck=F,
        xpd=F, mgp=c(3,0.3,0), ann=FALSE, xaxs = "i", border = NA, space = 0)
abline(h = 0)

# Add anomaly lines
par(new=T)
bar <- barplot(precip_mf_anom_timeseries, axes = FALSE, ylim = c(-2.5, 2.5), col = NA, border = NA, ann = FALSE, xaxs = "i", space = 0)
lines(bar, precip_mf_anom_timeseries, type = "l", col = "#045a8d", lwd = 1.5, ylim = c(-2.5, 2.5), ann = FALSE)
lines(bar, sw_mf_anom_timeseries, type="l", col="#b30000", lwd=1.5, ylim=c(-2.5, 2.5), ann=FALSE)
# lines(bar, direct_sw_mf_anom_timeseries, type="l", lty=2, col="#ef6548", lwd=1.5, ylim=c(-2.5, 2.5), ann=FALSE)
# lines(bar, diffuse_sw_mf_anom_timeseries, type="l", lty=2, col="#3690c0", lwd=1.5, ylim=c(-2.5, 2.5), ann=FALSE)

# x-axis labels, ticks, and title
axis(side=1, labels=names, tck=0.03, cex.axis=0.85, mgp=c(3,0.3,0), at=c(seq(from=0,to=120,by=12)))
# y-axis labels, ticks, and title
axis(side=4, mgp=c(3,0.3,0), las=1, tck=0.03, at=c(seq(from=-2.5,to=2.5,by=1)))
axis(side=4, mgp=c(3,0.3,0), las=1, tck=0.02, at=c(seq(from=-2,to=2,by=1)), labels = FALSE)
axis(side=2, tck = 0.03, las=1, mgp=c(3,0.3,0), at=c(seq(from=-3,to=3,by=1)))
axis(side=2, tck = 0.02, labels=FALSE, c(seq(from=-2.5,to=2.5,by=1)))
# Build the legend
# legend("topleft", lty=c(1,1,2,2), lwd = c(1.5, 1.5, 1.5, 1.5), horiz=FALSE, col=c("#045a8d", "#b30000", "#ef6548", "#3690c0"), c("Precipitation", expression(paste("PAR"['TOC']*"")),
#        expression(paste("Direct PAR"['TOC']*"")), expression(paste("Diffuse PAR"['TOC']*""))), bty="n", cex=0.85, pt.cex=0.85)
# legend("topleft", lty=c(1), lwd = c(1.5), horiz=FALSE, col=c("#045a8d"), c("Precipitation"), bty="n", cex=0.85, pt.cex=0.85)
legend("topleft", lty=c(1,1), lwd = c(1.5,1.5), horiz=FALSE, col=c("#045a8d", "#b30000"), c("Precipitation", expression(paste("PAR"['TOC']*""))), bty="n", cex=0.85, pt.cex=0.85)
legend("topright", fill = c(rgb(255,191,191,max=255), rgb(235,235,235,max=255), rgb(191,191,225,max=255)), c("Warm Phase", "Neutral", "Cool Phase"), bty="n", cex=0.85, pt.cex=0.85)
box()

mtext("Standardized Anomaly", side=4, line=1.75)
mtext("Multivariate ENSO Index Version 2", side=2, line=1.25)
mtext(expression(paste("Year")), 1, 1.5, outer=TRUE)

dev.off()
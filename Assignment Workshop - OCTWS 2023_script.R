## Assignment Workshop - OCTWS 2023
## Jackson Kusack
## March 24, 2023

rm(list=ls())
set.seed(725)
dev.off()

################################
## PACKAGES
if (!require('assignR')) install.packages('assignR'); library('assignR')
if (!require('raster')) install.packages('raster'); library('raster')
if (!require('terra')) install.packages('terra'); library('terra')
if (!require('sf')) install.packages('sf'); library('sf')
if (!require('sp')) install.packages('sp'); library('sp')
if (!require('rasterVis')) install.packages('rasterVis'); library('rasterVis')
if (!require('rnaturalearth')) install.packages('rnaturalearth'); library('rnaturalearth')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('plotly')) install.packages('plotly'); library('plotly')
if (!require('rdryad')) install.packages('rdryad'); library('rdryad')

custom.palette <- colorRampPalette(brewer.pal(9, "YlGnBu")) # Color palette


##############################
## ISOTOPE DATA

# Download data directly from Dryad
duck.d2h <- read.csv(rdryad::dryad_files_download(1697805)[[1]]) 
glimpse(duck.d2h) # Check data structure

# Select juvenile black ducks harvested in Ontario
duck.d2h <- filter(duck.d2h, age == "I") %>% 
  filter(state.prov == "ON") %>% 
  select(id.lab, VSMOW, long.dec, lat.dec)

# Create point object from the black duck data
duck.points <- vect(duck.d2h, geom = c("long.dec", "lat.dec")) # create vect (points) object

nrow(duck.d2h) # sample size


##############################
## SHAPEFILES

northamerica <- ne_countries(continent = "North America", scale = 50, returnclass = "sf") %>% # Countries
  st_transform(st_crs('EPSG:4326')) # Project to WGS84

northamerica.states <- ne_states(country =  c("Canada","United States of America"), returnclass = "sf") %>% # States and provinces
  st_transform(st_crs('EPSG:4326')) # Project to WGS84

breeding.range <- st_read("/vsicurl/https://github.com/JacksonKusack/Assignment-Workshop/raw/main/Shapefiles/ABDU_baldassarre_lowdens.shp") %>% # Load a shapefile from GitHub
  st_transform(st_crs('EPSG:4326')) # Project to WGS84

plot(st_geometry(breeding.range), col = brewer.pal(9, "YlGnBu")[5], border = NA)
plot(st_geometry(northamerica), add = T)

plot(duck.points, col = 'red', pch = 17)
plot(st_geometry(northamerica.states), add = T)


##############################
## ISOCAPES

# Download isoscape directly
gsd <- getIsoscapes(isoType = "GlobalPrecipGS", timeout = 1200) %>% 
  projectRaster(crs = CRS(SRS_string = 'EPSG:4326'))

# Pull out the d2H surfaces (1 and 2)
gsd <- gsd[[1:2]] 

plot(gsd, xlab="Longitude", ylab="Latitude", col = custom.palette(16))


##############################
## CALIBRATION

# Load known origin data
data("knownOrig")

knownOrig$sources[1:2] # list the dataset names and ids

unique(knownOrig$samples$Taxon) # list the taxa

# Extract the data from dataset 12 (van Dijk et al. 2014 - Mallards)
cal.mall <- subOrigData(marker = 'd2H', dataset = 12, ref_scale = NULL) # Extract data, keeping the ref_scale the same
cal.mall$data <- spTransform(cal.mall$data, CRS(SRS_string = 'EPSG:4326')) # reproject the data
str(cal.mall)

# Test the transformed vs untransformed d2H data
test.transformed <- subOrigData(marker = 'd2H', dataset = 12, ref_scale = "VSMOW_H")
mean(cal.mall$data$d2H) # Untransformed mean
mean(test.transformed$data$d2H) # Transformed mean

# Calibrate the mean and sd values from our isoscape
r <- calRaster(known = cal.mall, isoscape = gsd, interpMethod = 1, verboseLM = F, genplot = F) 
r.model <- r$lm.model # Extract model results
summary(r.model) 

ggplot(data = r$lm.data, aes(y = tissue.iso, x = isoscape.iso)) + 
  geom_point()  + 
  stat_smooth(method = "lm", formula = 'y ~ x') + 
  theme_classic()

plot(r$isoscape.rescale, xlab="Longitude", ylab="Latitude", col = custom.palette(16))

r$isoscape.rescale <- mask(rast(r$isoscape.rescale), breeding.range) %>% # Mask to breeding range
  crop(extent(breeding.range)) %>% # Crop to breeding range
  raster::stack() # Convert back to raster object

plot(r$isoscape.rescale, xlab="Longitude", ylab="Latitude", col = custom.palette(16))


####################################
## ASSIGNMENT

# Assignment model
origins <- pdRaster(r, data.frame(duck.d2h)[1:2], genplot = F)
cellStats(origins[[10]], 'sum') # Posterior probabilities across a single raster layer should sum to 1

plot(origins[[3]], col = custom.palette(8))
plot(st_geometry(breeding.range), add = T)
plot(st_geometry(northamerica), add = T)

#########################
## Priors

# Make a new prior, with arbitrary probabilities for different states
prior <- northamerica.states
prior$prob[prior$name == "Ontario"] <- 0.1
prior$prob[prior$name == "Québec"] <- 0.4
prior$prob[prior$name %in% c("New Brunswick","Nova Scotia","Prince Edward Island")] <- 0.2
prior$prob[prior$name %in% c("Newfoundland and Labrador")] <- 0.3
prior$prob[prior$admin == "United States of America"] <- 0

prior <- rasterize(prior, r$isoscape.rescale[[1]], field = "prob") %>% # rasterize the polygons to match the calibrated isoscape
  mask(breeding.range) # mask

plot(prior, col = custom.palette(16))

# Assignment model, with a prior
origins.prior <- pdRaster(r, data.frame(duck.d2h)[1:2], prior = prior, genplot = F) # assignment with the prior

# Compare the two assigments (with and without prior)
par(mfrow = c(1,2))
plot(origins[[12]], col = custom.palette(8))
plot(st_geometry(breeding.range), add = T)
plot(st_geometry(northamerica), add = T)

plot(origins.prior[[12]], col = custom.palette(8))
plot(st_geometry(breeding.range), add = T)
plot(st_geometry(northamerica), add = T)

############################
## Binary surfaces

# Determine odds ratio
odds <- 2/3 # select upper 66% of cells

# Create binary regions
binary.origins <- qtlRaster(origins, threshold = odds, thresholdType = "prob", genplot = F)

# Compare the probability surface and the binary surface 
par(mfrow = c(1,2))
plot(origins[[3]], col = custom.palette(8)) 
plot(st_geometry(northamerica), add = T)

plot(binary.origins[[3]], col = custom.palette(8))
plot(st_geometry(northamerica), add = T)

# Sum the binary surfaces
origins <- calc(binary.origins, sum)

# Save the raster
writeRaster(origins, "ABDU_origins", format = "ascii", overwrite = TRUE) 

# Plot the results
(p1 <- levelplot(origins, col.regions = custom.palette(16), margin=FALSE, 
                 legend=list(top=list(fun=grid::textGrob("Count", y=0.3, x=1.09)))) +
    latticeExtra::layer(sp.polygons(as(northamerica, "Spatial"), size = 0.5)) +
    latticeExtra::layer(sp.polygons(as(breeding.range, "Spatial"), fill = NA, alpha = 1, col = "black", lwd = 1.5, lty = 3)))

# Save the plot
png(filename = "ABDU_origins.png", units = "in", width = 6, height = 5, res=1200)  
p1
invisible(dev.off())

sessionInfo()
# This script generates the map of Waimea and the sampling sites for the study.

# Load Libraries-----------------------------------------------------------------------------------------------------------------------------------------
spatial_libs <- c("rgdal",
                  "FedData",
                  "raster",
                  "sp",
                  "RStoolbox",
                  "automap",
                  "geoviz",
                  "sf",
                  "ggspatial"
                  )

# general
general_libs <- c("data.table",
                  "dplyr",
                  "tidyr",
                  "ggplot2",
                  "viridisLite",
                  "grid",
                  "RColorBrewer",
                  "cowplot",
                  "ggnewscale")

# load libraries invisibly
invisible( lapply( c(spatial_libs,general_libs),
                   library, character.only = T) )
## Set Theme -------------------------------------------------------------------------------------------------------------------------------------------------------
theme_set(theme_bw())

# Load Data -------------------------------------------------------------------------------------------------------------------------------------------------------
# read phyloseq object
meta_file <- "data/processed/table_exports/all_marine_metabolite_sample_flat_table.csv"
meta<- fread(meta_file)
meta <- meta[ , c("lat","long"):= list(as.numeric(lat), as.numeric(long))]

# extract sites
point_data <- meta[ , .N , by = .(site_name, lat,long)]

## Set Boundaries -------------------------------------------------------------------------------------------------------------------------------------------------------
# generate the bounding box for the map in several formats for use with various mapping packages

# Check geographic range of sampling points
limits <- c(
  min(meta$long),
  min(meta$lat), 
  max(meta$long),
  max(meta$lat) 
)

# define a bounding box with a small cushion around the minimum and maximum
bbox <- list(
  p1 = list(long = limits[1] - 0.03, lat = limits[2] - 0.04) ,
  p2 = list(long = limits[3] + 0.03, lat = limits[4] + 0.04)
)


# genreate a spatial polygon template of the bounding box for use with FedData
bbox_extent <-
  polygon_from_extent(raster::extent(bbox$p1$long, bbox$p2$long, bbox$p1$lat, bbox$p2$lat),
                      proj4string = "+proj=longlat +datum=WGS84 +no_defs")


std_proj <- crs(bbox_extent)

# set up plot limits
plot_lims <- list(x=c(bbox$p1$long, bbox$p2$long), y = c(bbox$p1$lat, bbox$p2$lat))

##  Base Map -------------------------------------------------------------------------------------------------------------------------------------------------------

# hawaii coastline (state GIS, projection utm-4, nad83)
hi_coast <- st_read("data/raw/spatial/hawaii_coastline/coast_n83.shp",
                    crs = "epsg:26904")

# project to std proj
hi_coast <-st_transform(hi_coast,std_proj)

## Streams -----------------------------------------------------------------------------------------------------------------------------------
# check stream metadata to identify stream name
streams      <- st_read('data/raw/spatial/darstreams.kml')
waimea_river <- streams[streams$STREAM_NAM == "Waimea R",]

# crop river to plot area
waimea_river <- st_crop(waimea_river, c(xmin = plot_lims$x[1],
                                        xmax = plot_lims$x[2],
                                        ymin = plot_lims$y[1],
                                        ymax = plot_lims$y[2]))

test <- waimea_river[10,]

## Elevation Raster  -------------------------------------------------------------------------------------------------------------------------------------------------------

# get extent (federal 1 arc second elevation, lat long GRS80)
elev_extent <- spTransform(bbox_extent, CRSobj = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs ")

# pull elevation data in raster format
elev_ras_full <-
  get_ned(
    template = bbox_extent,
    label = "elev",
    res = "1",
    force.redo = F,
    raw.dir = "/data/raw/spatial",
    extraction.dir = "../data/processed/interim/"
  )

crs(elev_ras_full)


# project into standard format
elev_ras_full <- projectRaster(elev_ras_full, crs=std_proj)

crs(elev_ras_full)


## Bathymetry Rasters -------------------------------------------------------------------------------------------------------------------------------------------------------
# Bathymetry data from 
## Bathymetry
bath_ras_50 <- raster("data/raw/spatial/mhi_mbsyn_bathyonly_50m_v21.nc")
bath_ras_50 <- crop(bath_ras_50, bbox_extent)

# set projection and resolution to match elevation
bath_ras_50 <- projectRaster(bath_ras_50, elev_ras_full)
bath_dt_50 <- as(raster::as.data.frame(bath_ras_50, xy=T), "data.table")

## Mask Bathymetry and Elevation ----------------------------------------------------------------------------------------------------------------------------

# Mask Bathymetry
bath_ras <- mask(bath_ras_50, elev_ras_full)

# convert raster to data.table
# set a name for the 3rd dimension of the raster (i.e. depth)
bath_dt <- as(raster::as.data.frame(bath_ras, xy=T), "data.table")
setnames(bath_dt, 3, "depth")
bath_dt[ , alpha := ifelse(is.na(depth), 0, 1)]

# as matrix
bath_mat <- as.matrix(bath_ras) %>% t()

# test map 
bath_map <- ggplot() +
  geom_raster(data= bath_dt, mapping =aes(x = x, y = y, fill = depth))
bath_map

# Mask Elevation
elev_ras <- mask(elev_ras_full, bath_ras, inverse = T)

# and convert raster to data table and matrix 
elev_dt  <- raster::as.data.frame(elev_ras, xy=T) %>% as.data.table()
setnames(elev_dt, 3, "elev")
elev_dt[ , alpha := ifelse(is.na(elev), 0, 1)]
elev_mat <- as.matrix(elev_ras) %>% t()

# test map
elev_map <- ggplot() +
  geom_raster(data= elev_dt, mapping =aes(x = x, y = y, fill = elev))
elev_map

# try merging
bath_elev_merge <- merge(bath_ras_50, elev_ras)
bath_elev_dt <- raster::as.data.frame(bath_elev_merge , xy=T) %>% as.data.table()
bath_elev_mat <- as.matrix(bath_elev_merge)

 ggplot()+
  geom_raster(data = bath_elev_dt, mapping  =aes(x = x, y = y, fill = layer))

## Satelite Image -------------------------------------------------------------------------------------------------------------------------------------------------------

# function crop_sat will read in RGB landsat raster, crop it to map extent, write out as png

# use this fuction to define image size for using PNG images as map overlays

define_image_size <- function(match_to, major_dim = 400) {
  
  # calculate aspect ration (width/height) from lat/long bounding box
  lims <- extent(match_to)
  
  aspect_ratio <- abs((lims[1]- lims[2]) / (lims[3] - lims[4]))
  
  # define dimensions
  img_width <- ifelse(aspect_ratio > 1, major_dim, major_dim*aspect_ratio) %>% round()
  img_height <- ifelse(aspect_ratio < 1, major_dim, major_dim/aspect_ratio) %>% round()
  
  size_str <- paste(img_width, img_height, sep = ",")
  
  list(height = img_height, width = img_width, size = size_str)
  }


# calculate image dimensions based on size of elevation raster
img_size <- define_image_size(match_to = bath_elev_merge,
                              major_dim = max(dim(bath_elev_merge)))


## Plot 2D Map -------------------------------------------------------------------------------------------------------------------------------

# Plot on raster background
p1 <- ggplot() + 
  # elevation raster
  geom_raster(data= elev_dt,
              mapping = aes(x = x, y = y, fill = elev),
              alpha = elev_dt$alpha,
              interpolate = T) +
  scale_fill_gradient( low = "khaki", high = "darkgreen") +
  
  # new scale for bathymetry
  new_scale_fill()+
  geom_raster(data = bath_dt,
              mapping = aes(x = x, y = y, fill = depth), 
              alpha = bath_dt$alpha,
              interpolate = T) +
  scale_fill_gradient(low = "darkblue",
                      high = "white")+

  # river line
  geom_sf(data = waimea_river,
          color = alpha("blue",0.8)) +
  
  # points
  geom_point(aes(x = long, y = lat), 
             data = point_data,
             color = "white",
             shape = 16,
             size = 3) +
  geom_point(aes(x = long, y = lat, color = site_name), 
             data = point_data,
             shape = 17,
             size = 1.5) +

  scale_color_manual(values = viridis(n = 5)) +

  # annotations
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(
    height = unit(0.5, "cm"),
    width = unit(0.5, "cm"),
    pad_y = unit(0.75, "cm"),
    style = north_arrow_minimal("text_size = 8")
  ) +
  
  # map settings
  coord_sf(xlim = plot_lims$x,
           ylim = plot_lims$y,
           clip = "on") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        #axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,0,0)) +
  labs(fill = "Depth (m)", color = "Site")

p1

# With Oahu inset
oahu_limits <- c(-158.305366, 21.212964, -157.624682, 21.766562)

p2 <-ggplot()+
  # coastline
  geom_sf(data = hi_coast) +
  scale_color_manual(values = c("orange","red","blue","yellow","black")) +
  # rectangle
  geom_rect(aes(xmin = plot_lims$x[1]- 0.003,
                xmax = plot_lims$x[2] + 0.003,
                ymin = plot_lims$y[1] - 0.003,
                ymax = plot_lims$y[2] + 0.003), color = "gray10", fill = NA)+
  # map settings
  theme_void()+
  theme(legend.position = "none")+
  coord_sf(xlim = c(oahu_limits[c(1,3)]), ylim= c(oahu_limits[c(2,4)])) 

p2

p1.5 <- p1 + annotation_custom(
                              grob = ggplotGrob(p2),
                              ymin = plot_lims$y[1],
                              ymax = plot_lims$y[1] + .02,
                              xmin = plot_lims$x[2] - 0.02,
                              xmax = plot_lims$x[2]
)



p3 <-  ggplot() + 
  
  # Satellite Image
  annotation_custom(g, xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf)+
  
  # river line
  geom_sf(data = waimea_river,
          color = alpha("blue",0.8)) +
  
  # points
  geom_point(aes(x = long, y = lat), 
             data = point_data,
             color = "white",
             shape = 16,
             size = 3) +
  geom_point(aes(x = long, y = lat, color = site_name), 
             data = point_data,
             shape = 17,
             size = 1.5)+
  scale_color_manual(values = viridis(n = 5)) +
  
  # annotations
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(
    height = unit(0.5, "cm"),
    width = unit(0.5, "cm"),
    pad_y = unit(0.75, "cm"),
    style = north_arrow_minimal("text_size = 8")
  ) +
  
  # map settings
  coord_sf(xlim = plot_lims$x,
           ylim = plot_lims$y,
           clip = "on") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        legend.position = "none")


plot_grid(p3,p1.5, labels = c("A","B"), label_size = 12, axis = "tblr",align = "v" )

ggsave("output/map/ST_marine_map.pdf", 
       width = 15,
       height = 15,
       units = "in")

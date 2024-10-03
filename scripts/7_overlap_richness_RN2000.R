# Proyectar y ver solapamiento mapas
library(terra)
library(sf)
library(tidyr)
library(dplyr)
options("sf_use_s2" = FALSE)

riqueza <- rast("../data/riqueza_esteparias_reproyectada.tif")

#riqueza_proj <- terra::project(riqueza, "EPSG:4326")

riqueza_points <- as.data.frame(riqueza, xy = TRUE, cell=TRUE)

extents <- ext(riqueza)
plot(riqueza)
#IBAS_2010 <- st_read("scripts/prueba_proyeccion/Ibas_2010_p.shp")
#ggplot(data = IBAS_2010) +
#  geom_sf()

#red_natura <- st_read("scripts/prueba_proyeccion/rn2000/PS_NATURA2000_2023_PB.shp")
red_natura <- st_read("../data/threats/red_natura_2000/Red_Natura_2000_de_Espa%C3%B1a.shp")
red_natura = st_make_valid(red_natura)
red_natura_proyect <- st_transform(x = red_natura, crs=25830)

class(red_natura_proyect)

# Overlap betweem richness and RN2000
rn_sp <- vect(red_natura_proyect)
#sf:::as_Spatial(red_natura_proyect$geom) # This works
class(rn_sp)
riqueza_rn <- terra::extract(x = riqueza, y = rn_sp, cells= TRUE)
rn_df <- red_natura_proyect %>% 
  st_drop_geometry() %>% 
  mutate(ID := seq_len(nrow(.)))

red_NATURA_richness <- 
  #--- back to sf ---#
  st_as_sf(rn_sp) %>% 
  #--- define ID ---#
  mutate(ID := seq_len(nrow(.))) %>% 
  #--- merge by ID ---#
  dplyr::left_join(., riqueza_rn, by = "ID")

str(riqueza_points)
str(riqueza_rn)

richness_red_NATURA <- riqueza_points %>% 
  dplyr::left_join(., riqueza_rn, by = "cell") %>% 
  dplyr::left_join(., rn_df, by = "ID") %>% 
  distinct(cell, .keep_all = TRUE)


# Calculate percentage of pixel from X richness value overlapped with areas Red Natura 20000

richness_red_NATURA_overlapped <- richness_red_NATURA

richness_values <- unique(richness_red_NATURA_overlapped$riqueza_esteparias_reproyectada.x)
overlaped_values <- data.frame(richness= as.integer(), percentage_overlaped= as.integer())

for(x in 1:length(richness_values)){
  rvalue <- richness_values[x]
  filter_r_value <- richness_red_NATURA_overlapped %>% 
    filter(riqueza_esteparias_reproyectada.x== rvalue)
  filter_red_natura <- filter_r_value %>% 
    drop_na(ID)
  percentage_r <- (nrow(filter_red_natura)/nrow(filter_r_value))*100
  s_value <- data.frame(richness= rvalue, percentage_overlaped= percentage_r)
  overlaped_values <- rbind(overlaped_values, s_value)
  
}

# only the 25 percentage of the areas with 6 grasslands species are protected y the RN
overlaped_values


library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

spain <- ne_countries(scale = "medium", returnclass = "sf", country = "spain")
#class(world)

spain_crop <- rmapshaper::ms_filter_islands(spain,
                                            min_area = 100000000000,
                                            drop_null_geometries=TRUE)

riqueza_binaria6 <- riqueza 
riqueza_binaria6[riqueza_binaria6 < 6] <- NA
riqueza_binaria6[riqueza_binaria6 >= 6] <- 6
  
plot(riqueza_binaria6)
riqueza_binaria6_points <- as.data.frame(riqueza_binaria6, xy = TRUE, cell=TRUE)



ggplot() +
  geom_raster(aes(x=x, y=y, fill = riqueza_esteparias_reproyectada), data = riqueza_binaria6_points) +
  ylab("") +
  xlab("") +
 # geom_sf(data=spain_crop, fill="grey70", color=NA, alpha= 0.1) +
  geom_sf(data=red_natura_proyect, size=0.2, fill=alpha("green4",0.5), color=alpha("black",0.2)) +
  coord_sf(xlim = c(-63390.7932, 1066474.6778), ylim = c(3985951.8493, 4863441.8715)) 
  
  
  
  
#-------------------------------------------------------------------------------
#           Modelamiento de la Distribucion Espacial de B. cockerelli
# authors: Felipe López-Hernández y Paula Bermeo-Fúquene
#
# affiliation: Corporación colombiana de investigación agropecuaria AGROSAVIA
#
# date: Mayo, 2025
# ------------------------------------------------------------------------------


#Cargar librerías/paquetes -------------------------------------------------------
pacman::p_load(BIEN, corrplot, cowplot, dplyr, dismo, ecospat, ellipsenm, ENMeval, envirem, geodata, geosphere, ggplot2, ggspatial, ggnewscale, ggpattern, grid, gridExtra, kuenm, maxnet, pROC, raster, readr, readxl, rgbif, rnaturalearth, rnaturalearthdata, RPostgreSQL, rworldxtra, sf, sp, spThin, stars, terra, tidyverse, tiff, tmap, viridis)

# Establecer el directorio de trabajo donde se guardarán y leerán archivos
#setwd("D:/AGROSAVIA_ 2020 _BUHO/PAPERS_2020/SDM BC") #PC Felipe
setwd("D:/AGROSAVIA - CORPORACION COLOMBIANA DE INVESTIGACION AGROPECUARIA/Mapas de riesgo B. cockerelli CC - Documentos/Mapas de riesgo B. cockerelli  CC/Meta24_2025/Modelos ecologicos de nicho/") #PC Paula



# Descarga, limpieza y filtraje de ocurrencias ----------------------------


# Seleccionar columnas de interés: longitud y latitud en la base de datos de presencia

name_all <- read_delim('Datos presencia/Presencia_total.csv', delim = ';', col_names = T) |> as.data.frame() |> na.omit() #Genera el dataframe y elimina filas con valores faltantes (NAs)
 

#Reducir 1 punto por cada km2
# Run spatial thinning, using 1 km distance
thinnned <-
  thin(loc.data = name_all,
       lat.col = "Latitude", long.col = "Longitude",
       spec.col = "Species",
       thin.par = 1, reps = 100,
       locs.thinned.list.return = TRUE,
       out.base = 'name_all',
       write.files = TRUE, max.files=1, out.dir="Datos presencia/",
       write.log.file = FALSE)
# Have a look at the first thinned data set
View(thinnned[[1]])

# Plot the first thinned data set over the full data set to see thinning
points(thinnned[[1]]$Longitude, thinnned[[1]]$Latitude, col = "red", pch = 20)

dev.off()
par(mfrow = c(2, 1), cex = 0.5, mar = rep(0.6, 4))
points(thinnned[[1]]$Longitude, thinnned[[1]]$Latitude, col = "red", pch = 20)
points(name_all[, 5:6], col = 'gray')

# Convertir los datos a un objeto espacial `sf`
puntos_sf <- st_as_sf(name_all, coords = c("Longitude", "Latitude"), crs = 4326) 
# Convierte las coordenadas en un objeto sf con sistema de referencia CRS WGS 84.

# Descargar el polígono de Colombia usando datos de Natural Earth
# Obtener los polígonos de primer nivel administrativo (departamentos) de Colombia
departamentos <- ne_states(country = "Colombia", returnclass = "sf")

# Obtener los departamentos de Colombia
departamentos <- ne_states(country = "Colombia", returnclass = "sf")

# Seleccionar departamento de Nariño
pais <- departamentos %>%
  filter(name %in% c("Nariño"))

# Supongamos que 'puntos_sf' ya está cargado como un objeto sf con geometría de puntos

# Filtrar los puntos que están dentro de alguno de los polígonos
puntos_filtrados <- puntos_sf %>%
  filter(rowSums(st_within(geometry, pais, sparse = FALSE)) > 0)

# Selecciona los puntos que están dentro de los límites del polígono de Colombia.

# Convertir los datos filtrados a formato tabular (opcional)
puntos_filtrados_tabla <- puntos_filtrados %>%
  mutate(
    longitude = st_coordinates(geometry)[, 1], # Extrae la longitud de las geometrías.
    latitude = st_coordinates(geometry)[, 2]  # Extrae la latitud de las geometrías.
  )

# Convertir a data.frame y eliminar la columna de geometría
puntos_filtrados_tabla <- as.data.frame(puntos_filtrados_tabla) # Convierte el objeto a data.frame.
data_name <- puntos_filtrados_tabla[, -1] # Elimina la columna de geometría si existe.


# Definición del área de calibración --------------------------------------

#Cargar datos de presencia filtrados
name_all <- read_delim("Datos presencia/name_all_thin1.csv", delim = ',', col_names = T) |> as.data.frame() |> na.omit()


#Descargar raster de elevación para Colombia
#elevation_30s(country="COL", path=("Meta24_2025/Modelos ecológicos de nicho/Rasters/"))

#elevacion <- list.files("D:/AGROSAVIA_ 2020 _BUHO/PAPERS_2020/AvocadoClust/Scripts_resultados_definitivos/Clustering_datos_Paula/2021_V3/Dataset 4/", pattern = 'elev.tif', full.names = TRUE) #PC Felipe
elevacion <- list.files("Rasters/elevation/", pattern = 'elev.tif', full.names = TRUE)

 elev1 <- stack(elevacion)
# # recortamos y enmascaramos usando un polígono
 elev <- elev1
 elev <- crop(elev, allac)
 elev <- mask(elev, allac)
 elev1 <- crop(elev1, allac)
# # creamos una máscara terrestre / marina estableciendo valores <0 a NA, y valores> = 0 a 1
# elev[elev < 0] <- NA
# elev[!is.na(elev)] <- 1
# 
# elev <- mask(elev1, elev)

# Generando Buffer

occs.sf <- sf::st_as_sf(recordsRomeron, coords = c("Longitude","Latitude"), crs = raster::crs(elev))

# Now, we project our point data to an equal-area projection, which converts our 
# degrees to meters, which is ideal for buffering (the next step). 
# We use the typical Eckert IV projection.
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
occs.sf <- sf::st_transform(occs.sf, crs = eckertIV)

# Buffer all occurrences by 500 m, union the polygons together 
occs.buf <- sf::st_buffer(occs.sf, dist = 500) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(elev))

plot(occs.buf)

## Opc 1. Concave area (Including buffer 500 Km)

occs.buf <- concave_area(name_all, longitude = "Longitude", latitude = "Latitude", length_threshold = 5,
                       buffer_distance = 0.5) |> vect() 

plot(occs.buf); points(name_all[, 2:3], col = 'red'); legend("topleft", legend = "Concave area", bty = "n"); points(name_all[, 2:3], col = 'black')

## Opc 2. Convex area
occs.buf1<- convex_area(name_all, longitude = "Longitude", latitude = "Latitude",
               buffer_distance = 0.5) |> vect()
plot(occs.buf1); points(name_all[, 5:6], col = 'red'); legend("topleft", legend = "Convex area", bty = "n")

#Opc 3. ecoregions analysis

#Load shape WWF
ecor <- "Vectores/wwf_terr_ecos.shp"
ecor = vect(ecor)
ecor <- as(ecor, "Spatial")

# Areas by selecting ecoregions (Including buffer 30 Km)
M_ecorreg <- polygon_selection(name_all, longitude = "Longitude", latitude = "Latitude",
                               polygons = ecor, buffer_distance = 0.5) |> vect() 

plot(M_ecorreg); points(name_all[, 5:6], col = 'red'); legend("topleft", legend = "Ecoregions", bty = "n")

#Plot compuesto
dev.off()
par(mfrow = c(2, 2), cex = 0.5, mar = rep(0.6, 4))
plot(occs.buf); points(name_all[, 2:3], col = 'red'); legend("top", legend = "Concave area", bty = "n"); points(name_all[, 2:3], col = 'black')
plot(occs.buf1); points(name_all[, 2:3], col = 'red'); legend("top", legend = "Convex area", bty = "n")
plot(M_ecorreg); points(name_all[, 2:3], col = 'red'); legend("top", legend = "Ecoregions", bty = "n")

#Conclusión: con el fin de ampliar el rango de background, teniendo en cuenta la larga extensión de las ecorregiones, se decide utilizar la opción 2. Convex area

#Save shapefile
writeVector(occs.buf1, "Vectores/calibrat_area.shp", overwrite=T)


# Descarga y selección de variables predictorias (OMITIR SÍ LAS BIOVARIABLES YA ESTÁN DESCARGADAS) --------------------------

#Descargar climate data Worldclim V2.1 

#worldclim_country(country="COL", var="bio", res=0.5, path="Rasters/")



# Selección de variables predictoras --------------------------------------


#Cargue del área de calibración/Nariño
shape <- 'Vectores/Servicios_Públicos_-_Departamentos.shp'
ca = vect(shape)
ca <- ca[ca$DEPTO == "NARIÑO"]
crs(ca) <- "EPSG:3857"

#plot(ca)

#Cargue puntos de presencia
name_all <- read_delim("Datos presencia/name_all_thin1.csv", delim = ',', col_names = T) |> as.data.frame() |> na.omit()

#List Worldclim (wc) rasters
wc <- rast('Rasters/bio/climate/wc2.1_country/COL_wc2.1_30s_bio.tif')
crs(wc) <- "epsg:4326"

#Change raster names
names(wc) <- c(paste0("Bio_", seq(1, 19, by=1)))

#Project rasters
ca <- project(ca, wc)


#Mask variables to the ca (Calibration Area)
wc.ca <- mask(crop(wc, ca), ca)
plot(wc.ca[[1]])
res(wc.ca)


#Extracting values from rasters

#Eliminar primera columna datos de presencia
name_all <- name_all[ ,-1]

Pointrast <- raster::extract(wc.ca, name_all)
PointRast <- Pointrast[,-1]
sdmdataP <- data.frame(PointRast)


#Table cleaning
valores_soma <- rowSums(sdmdataP)
valores_soma_validos <- 1:nrow(sdmdataP)
valores_soma_validos <- ifelse(is.na(valores_soma), NA, valores_soma_validos)
valores_soma_validos <- subset(valores_soma_validos, valores_soma_validos >0)

sdmdata_validosP <- sdmdataP[valores_soma_validos, ]


# Correlation 

cor(sdmdata_validosP) 
round(cor(sdmdata_validosP),2)
par(mfrow=c(1,1))
hist(cor(sdmdata_validosP), main="Correlation matrix", col="gray90")

write.table(round(cor(sdmdata_validosP), 2), 'Rasters/Analisis_predictores/correlac.xls', row.names = T, sep = '\t')
cor_matrix <- cor(sdmdata_validosP)
write.table(ifelse(cor_matrix>= 0.7 | cor_matrix  <= -0.7, 'Yes', 'No'), 'Rasters/Analisis_predictores/correlac_Y-N_BIO.xls', row.names = T, 
            sep = '\t')

tiff('Rasters/Analisis_predictores/grafica_correlac.tif', width = 20, height = 20, units = 'cm', res = 300)#, compression = 'lzw')
corrplot(cor(sdmdata_validosP), type = 'lower', diag = F, tl.srt = 45, mar = c(3, 0.5, 2, 1),
         title = 'Environmental Variables')
dev.off()

# figure export
tiff('Rasters/Analisis_predictores/grafica_correlac.tiff', width = 8, height = 9, units = 'in', res = 900, compression = 'lzw')
corrplot(cor(sdmdata_validosP), type = "lower", diag = F, title = 'Correlation Between Environmental Variables', 
       mar = c(3, 0.5, 2, 1), tl.srt = 45)
dev.off()


# 2. PCA 

pca <- prcomp(sdmdata_validosP, scale = T)

# eigenvalues - autovalores
summary(pca)

# Gráfico de barras con las contribuciones
screeplot(pca, main = 'Autovalores')
abline(h = 1, col = 'red', lty = 2)

tiff('Rasters/Analisis_predictores/screeplotPCA.tif', width = 20, height = 20, units = 'cm', res = 300)
screeplot(pca, main = 'Autovalores')
abline(h = 1, col = 'red', lty = 2)
dev.off()

# valores de cada eje (eigenvectors - autovectores - scores)
pca$x
plotPCA <- biplot(prcomp(sdmdata_validosP, scale. = T, ))
biplot(pca, col = c('darkblue', 'red'),
       scale = 0, xlabs = rep("*", 89)
)

# Relación de variables en cada eje (loadings - cargas)
pca$rotation[, 1:4]
abs(pca$rotation[, 1:4])

# exportar tabla con cada contribución
write.table(abs(pca$rotation[, 1:4]), 'Rasters/Analisis_predictores/PCA_BIO.xls', row.names = T, sep = '\t')

# plot
biplot(pca)
dev.off()

# 3.  VIF 
  vif <- usdm::vifstep(sdmdata_validosP, th=10, keep = c("Bio_1", "Bio_12"))

#save.output
capture.output(vif, file='Rasters/Analisis_predictores/VIFdata.txt')

#vif1 <- as.data.frame(vif)


#Save final variables

pred <- c(wc.ca$Bio_1, wc.ca$Bio_2, wc.ca$Bio_3, wc.ca$Bio_12, wc.ca$Bio_15)
raster::writeRaster(pred, "Rasters/wccurrent_Narino.tiff", overwrite=TRUE)

#pred <- stack(pred)
#wc.fv <- stack("Other_files/Rasters_wc/wccurrent_ca.tiff")

 raster::writeRaster(wc.ca$Bio_1, "Rasters/wccurrent_bio1.tiff", overwrite=TRUE) #T media
 raster::writeRaster(wc.ca$Bio_2, "Rasters/wccurrent_bio2.tif", overwrite=TRUE) #Rango T
 raster::writeRaster(wc.ca$Bio_3, "Rasters/wccurrent_bio3.tif", overwrite=TRUE) #Isotermalidad
 raster::writeRaster(wc.ca$Bio_12, "Rasters/wccurrent_bio12.tif", overwrite=TRUE) #Pp anual
 raster::writeRaster(wc.ca$Bio_15, "Rasters/wccurrent_bio15.tif", overwrite=TRUE) #SD Pp


 # Selección del mejor modelo MAxEnt usando KUENM ----------------------------------------
 pacman::p_load(kuenm, readr, raster, terra)
 
 # Establecer el directorio de trabajo donde se guardarán y leerán archivos
 setwd("D:\\AGROSAVIA - CORPORACION COLOMBIANA DE INVESTIGACION AGROPECUARIA\\Mapas de riesgo B. cockerelli CC - Documentos\\Mapas de riesgo B. cockerelli  CC\\Meta24_2025\\Modelos ecologicos de nicho")
 
 #Cargar datos de presencia filtrados
 # name_all <- read_delim("Datos presencia/name_all_thin1.csv", delim = ',', col_names = T) |> as.data.frame() |> na.omit()
 # 
 # colnames(name_all) <- c("species", "longitude", "latitude") 
 # 
 # # Set a random seed in order to be able to reproduce this analysis.
 # set.seed(2025)
 # 
 # split <- kuenm_occsplit(occ = name_all, train.proportion = 0.7, method = "random", save = TRUE, name = "occ")
 # 
 # # #Variables predictoras raster stack
 # pred <- terra::rast(raster::stack("Rasters/wccurrent_Narino.tiff"))
 # ## preparing sets of variables (complete the code)
 # help(kuenm_varcomb)
 # 
 # #Convert predictors in ascii format
 # 
 # #Variables predictoras raster stack
 # pred <- terra::rast(stack("Rasters/wccurrent_Narino.tiff"))
 # 
 # raster::writeRaster(pred$Bio_1, "Variables/Bio_1.asc", overwrite=TRUE, filetype = "AAIGrid", NAflag = -9999)
 # raster::writeRaster(pred$Bio_2, "Variables/Bio_2.asc", overwrite=TRUE, filetype = "AAIGrid", NAflag = -9999)
 # raster::writeRaster(pred$Bio_3, "Variables/Bio_3.asc", overwrite=TRUE, filetype = "AAIGrid", NAflag = -9999)
 # raster::writeRaster(pred$Bio_12, "Variables/Bio_12.asc", overwrite=TRUE, filetype = "AAIGrid", NAflag = -9999)
 # raster::writeRaster(pred$Bio_15, "Variables/Bio_15.asc", overwrite=TRUE, filetype = "AAIGrid", NAflag = -9999)
 # 
 # vs <- kuenm_varcomb(var.dir = "Variables", out.dir = "M_var", min.number = 4,  in.format = "ascii", out.format = "ascii")
 # 
 # #Markdown
 # kuenm_start(file.name= 'Bactericera_KUENM')
 
 # Calibration process
 
 set.seed(2025)
 occ_joint <- "occ_joint.csv"
 occ_tra <- "occ_train.csv"
 M_var_dir <- "M_var"
 batch_cal <- "Candidate_models"
 out_dir <- "Candidate_Models"
 reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1))
 f_clas <- c("l", "lq","lqp") # 
 maxent_path <-  "D:/OneDrive - AGROSAVIA - CORPORACION COLOMBIANA DE INVESTIGACION AGROPECUARIA/Agrosavia/Apps/maxent"
 wait <- F
 run <- T
 max.memory = 4000
 args = "maximumbackground=2500"
 
 kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
           out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
           maxent.path = maxent_path, wait = wait, run = run)
 
 
 
 #Models evaluation
 set.seed(2025)
 occ_test <- "occ_test.csv"
 out_eval <- "Calibration_results"
 threshold <- 5
 rand_percent <- 50
 iterations <- 500
 kept <- T
 selection <- "OR_AICc"
 paral_proc <- FALSE 
 
 cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                         out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                         kept = kept, selection = selection, parallel.proc = paral_proc) 
 

# Modelamiento MAxEnt usando Maxnet -----------------------------------------------------------


#Cargar datos de presencia filtrados
name_all <- read_delim("Datos presencia/name_all_thin1.csv", delim = ',', col_names = T) |> as.data.frame() |> na.omit()


#Variables predictoras raster stack
pred <- terra::rast(stack("Rasters/wccurrent_Narino.tiff"))

# Set a random seed in order to be able to reproduce this analysis.
set.seed(2025)

# Selección del número de background points
Number_background_points = 2500
os <- list(validation.bg = "partition") #Métricas a partir de los valores de testing

#Función para generar la métrica pROC test 10.1016/J.ECOLMODEL.2007.11.008 
pROC <- function(vars) {
  pROC <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
  out <- data.frame(pROC_auc_ratio = pROC$pROC_summary[1], 
                    pROC_pval = pROC$pROC_summary[2], row.names = NULL)
  return(out)
}


# Run ENMevaluate
Results <- ENMevaluate(occs =  name_all[,-1 ], envs = pred, n.bg = Number_background_points,
                       algorithm = 'maxnet', partitions = "block", 
                       tune.args = list(fc = c("LQ"), rm = 0.1, other.settings=os), user.eval=pROC 
                       ) 
                                        #rm = c(seq(0.1, 1, 0.1), 2:6)))
    #Revisar sí se puede optimizar con las proyecciones de las curvas de los predictores "Free"

# Modeling results
Results@results


## Best Model Prediction
Models <- Results@results
Models$ID <- 1:nrow(Models)
Models <- Models %>% arrange("pROC_pval.avg")
BestModels <- Results@models[[Models$ID[1]]] #|> as.data.frame() HERE!!!!!!!

saveRDS(Models, file = "Salidas MaxEnt/Bcokerelli_models_present_2500_backg_points")



#Prediction <- predict(predicted, BestModels, type = "cloglog") #Error: [predict] the number of values returned by 'fun' (model predict function) does not match the input. Try na.rm=TRUE?#Proyección al área de estudio

Prediction <- predict(pred, BestModels, type = "cloglog", na.rm=T) 

#Responses of the variables (predictors)

plot(BestModels,  type = "cloglog")


# Graficando el nicho potencial actual usando tmap

#Realized niche
RN <- Prediction$lyr1
palette <- viridis_pal(begin = 0, end = 1, option = 'D', direction = 1)(5)
labs   <- c("<0.2","0.2–0.4","0.4–0.6","0.6–0.8",">0.8")

#Load Occurrences data
Occ <- read_delim('Datos presencia/name_all_thin1.csv', delim = ',', col_names = T) |> as.data.frame()
wgs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
WGS84 <- sp::CRS(wgs)
occ <- sp::SpatialPointsDataFrame(Occ[, c('Longitude', 'Latitude')], Occ, proj4string = WGS84)

#Load Narino municipalities
NAR <- st_read('Vectores/Servicios_Públicos_-_Municipios_2005.shp')
NARIN <- NAR[NAR$DEPTO == "NARIÑO", ]
target_crs <- crs(RN, proj = TRUE) 
NARIN <- st_transform(NARIN, target_crs)

# change to tmap mode
tmap::tmap_mode(mode = "plot")
tmap_options(check.and.fix = TRUE)


#Niche Map

Nicho <-  
  tm_shape(RN, bbox = c(-80, 0.35, -76.5, 2.7)) +
  tm_raster(style = "fixed", legend.reverse= T,
            #legend.format = list(text.align = 'right'),
            palette = palette, labels=labs, title = 'Idoneidad', breaks=seq(0, 1, 0.2), midpoint=NA) +
  tm_shape(occ) +
  tm_symbols(size= 0.08, col= 'red', alpha= 0.8) +
  tm_shape(NARIN) +
  tm_borders(col = 'white', lwd = 0.05, lty=2) +
  #tm_text('ISO_SUB', size = 1, col = 'black', fontfamily = 'sans', fontface='bold') +  
  tm_compass(position = c(0.85, 0.7), size= 1.5) + 
  tm_scale_bar(text.size = 0.6, position = c(0.2, 0), width = 0.15) + #scale bar
  tm_graticules(lines = F, labels.rot = c(0, 90), labels.size = 0.55) +
  tm_add_legend(type = "symbol", 
                labels = "Ocurrencias de B. cockerelli", 
                col = "red", 
                lwd = 1,
                size = 0.5
  ) +
  tm_layout(main.title= expression(paste("Distribución potencial de ", italic("B. cockerelli"), " actual en Nariño")),
            legend.position = c("left", "top"),
            legend.hist.height = 0,
            inner.margins = c(0, 0, 0, 0),
            legend.title.size = 1,
            legend.text.size = 0.9,
            legend.width = 0.7
  )

Nicho

tmap_save(
  Nicho,
  Nicho,
  filename = "Mapas/Nicho_B_cockerelli.png",
  width = 7, height = 6, units = "in", dpi = 600
)



# Crear las carpetas para la lectura de los archivos en el futuro ---------


wcF <- rast(list.files(path='D:/SDM/Futur_data', pattern = ".tif$", full.names = TRUE))
crs(wcF) <- "epsg:4326"

#Shapefile G area
shape <- "Servicios_Públicos_-_Departamentos.shp"
G = vect(shape)
G <- G[G$DEPTO == "NARIÑO"]
crs(G) <- "EPSG:3857"
plot(G)

#Project rasters
G <- project(G, wcF)

## Future distribution with worldclim Data 
wcF.ca <- mask(crop(wcF, G), G)
plot(wcF.ca[[1]])
res(wcF.ca)

#Change names
wcF1 <- wcF.ca

names(wcF1) <- gsub('wc2.1_30s_bioc_','',names(wcF1))
names(wcF1) <- gsub('MIROC6', 'MI',names(wcF1))
names(wcF1) <- gsub('MPI-ESM1-2-HR', 'MP',names(wcF1))
names(wcF1) <- gsub('MRI-ESM2-0', 'MR',names(wcF1))
names(wcF1) <- gsub('ssp245', '2',names(wcF1))
names(wcF1) <- gsub('ssp370', '3',names(wcF1))
names(wcF1) <- gsub("2041-2060",'50',names(wcF1))
names(wcF1) <- gsub("2061-2080",'70',names(wcF1))



ssps <- c('2', "3")
years <- c('50', '70')
bands <- c(1L, 2L, 3L, 12L, 15L) 

  for (rc in 1: length(ssps)) {
    for (ye in 1: length(years)) {
      for (band in 1: length(bands)) {
        
        MI <- wcF1[[paste0("MI_",ssps[rc],'_',years[ye],'_',bands[band])]]
        MP  <- wcF1[[paste0("MP_",ssps[rc],'_',years[ye],'_',bands[band])]]
        MR <- wcF1[[paste0("MR_",ssps[rc],'_',years[ye],'_',bands[band])]]
        
        bio <- mean(MI, MP, MR)
        
        assign(paste0('bio_', bands[band]), bio)
      }
      
      Result <- stack(c(bio_1, bio_2, bio_3, bio_12, bio_15))
      
      names(Result) <- c('Bio_1', 'Bio_2', "Bio_3", 'Bio_12', "Bio_15")
      
      dir.create(paste0(ssps[rc],'_',years[ye]))
      
      dir <- (paste0(ssps[rc],'_',years[ye]))
      
      raster::writeRaster(Result, filename= file.path(dir, names(Result)), bylayer=TRUE, format="ascii", overwrite=T)
      
    }
  }


# Random Forest -----------------------------------------------------------

p_load::pacman(caret, ranger, randomForest)


# HIPERPAR?METROS, N?MERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICI?N

particiones  <- 5
repeticiones <- 10
hiperparametros <- expand.grid(C = c(1))

# Hiperparametros
hiperparametros <- expand.grid(mtry = c(3, 4, 5, 7),
                               min.node.size = c(2, 3, 4, 5, 10, 15, 22), #Ajustar de acuerdo con nuestros pará
                               splitrule = "gini")

set.seed(123)
seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)
for (i in 1:(particiones * repeticiones)) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros)) 
}
seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)


# DEFINICI?N DEL ENTRENAMIENTO

control_train <- trainControl(method = "repeatedcv", number = particiones,
                              repeats = repeticiones, seeds = seeds,
                              returnResamp = "final", verboseIter = FALSE,
                              classProbs = TRUE, allowParallel = TRUE)

#AJUSTE DEL MODELO

modelo_rf <- train(pb ~ ., data = datos_train_prep,
                   method = "ranger",
                   tuneGrid = hiperparametros,
                   metric = "Accuracy", #métrica probabilística
                   trControl = control_train,
                   # N?mero de ?rboles ajustados
                   num.trees = 700)


modelo_rf$finalModel

# REPRESENTACI?N GR?FICA

# pdf("RF_EXP4_V6_training.pdf")
ggplot(modelo_rf, highlight = TRUE) +
  scale_x_continuous(breaks = 1:30) +
  labs(title = "Evoluci?n del accuracy del modelo Random Forest") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

modelo_rf$finalModel
# dev.off()

# TRAINING : Error train

predicciones <- predict(object = modelo_rf,
                        newdata = datos_train_prep,
                        type = "prob")

AUC_RF_train <- auc(response = datos_train_prep$pb, 
                    predictor = predicciones$pres) 
AUC_RF_train

# TESTING : Error test

predicciones_raw <- predict(modelo_rf, newdata = datos_test_prep,
                            type = "raw")
confusionMatrix <- confusionMatrix(data = predicciones_raw, reference = datos_test_prep$pb,
                                   positive = "pres")

confusionMatrix

precision(data = predicciones_raw, reference = datos_test_prep$pb,
          positive = "pres")

error_test <- mean(predicciones_raw != datos_test_prep$pb)
paste("Model Error Test", round(error_test*100, 2), "%")


# ROC test ----------------------------------------------------------------

# ROC
# C?lculo de la curva
predicciones <- predict(object = modelo_rf,
                        newdata = datos_test_prep,
                        type = "prob")

curva_roc <- roc(response = datos_test_prep$pb, 
                 predictor = predicciones$pres) 

# tr <- threshold(curva_roc, 'spec_sens')
# Gr?fico de la curva
# plot(curva_roc)

d= cbind(curva_roc$specificities, curva_roc$sensitivities)
colnames(d) <- c("specificities","sensitivities")
d <- as.data.frame(d)
# pdf("ROC_v5.pdf")
p2 <- ggplot(d,aes(specificities, sensitivities), )+
  geom_line()+
  geom_point() + scale_x_reverse()
p2
# dev.off()

AUC_RF_test <- auc(response = datos_test_prep$pb, 
                   predictor = predicciones$pres) 
AUC_RF_test

RFpm  <-  raster::predict(predicted,modelo_rf, progress="text",type = "prob", 
                          args=c("extrapolate=T", "doclamp=TRUE"))

mapRF_correct <- 1 - RFpm 

# write.csv(RFpm,"RFpm.csv", row.names = FALSE)

### graficando ggplot2
mapRF.ggplot2 <- mapRF_correct %>% as("SpatialPixelsDataFrame") %>% 
  as.data.frame()
a2 <- ggplot() + geom_tile(data = mapRF.ggplot2, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis(option="turbo") + xlab("Longitud") + ylab("Latitud")+ 
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_nautical) + 
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "gray45"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "gray45"),
        panel.background = element_rect(fill = "gray54",colour = "gray85",size = 0.5, linetype = "solid"))+ labs(fill = "Probabilidad") + 
  labs(subtitle= "Modelamiento de la distribución de"~~italic("Rubus spp.")~~"",title="Random Forest") +
  # geom_sf(data = Romeron, fill="black", color="black", size=1) +
  scale_fill_viridis(option="turbo")
a2




# K-Nearest Neighbor (kNN) ------------------------------------------------



# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN

particiones  <- 5
repeticiones <- 10

# Hiperparámetros
hiperparametros <- data.frame(k = c(1, 2, 5, 10, 15, 20, 30, 50))

set.seed(123)
seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)
for (i in 1:(particiones * repeticiones)) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros)) 
}
seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO

control_train <- trainControl(method = "repeatedcv", number = particiones,
                              repeats = repeticiones, seeds = seeds,
                              returnResamp = "final", verboseIter = FALSE,
                              allowParallel = TRUE)

# AJUSTE DEL MODELO

set.seed(342)
modelo_knn <- train(pb ~ ., data = datos_train_prep,
                    method = "knn",
                    tuneGrid = hiperparametros,
                    metric = "Accuracy",
                    trControl = control_train)
modelo_knn



# REPRESENTACIÓN GRÁFICA

ggplot(modelo_knn, highlight = TRUE) +
  scale_x_continuous(breaks = hiperparametros$k) +
  labs(title = "Evolución del accuracy del modelo KNN", x = "K") +
  theme_bw()





# TRAINING : Error train

predicciones <- predict(object = modelo_knn,
                        newdata = datos_train_prep,
                        type = "prob")

AUC_KNN_train <- auc(response = datos_train_prep$pb, 
                     predictor = predicciones$pres) 
AUC_KNN_train




# TESTING : Error test

predicciones_raw <- predict(modelo_knn, newdata = datos_test_prep,
                            type = "raw")
confusionMatrix <- confusionMatrix(data = predicciones_raw, reference = datos_test_prep$pb,
                                   positive = "pres")

confusionMatrix



precision(data = predicciones_raw, reference = datos_test_prep$pb,
          positive = "pres")



# Error de test
error_test <- mean(predicciones_raw != datos_test_prep$pb)
paste("Model Error Test", round(error_test*100, 2), "%")



# ROC
# C?lculo de la curva
predicciones <- predict(object = modelo_knn,
                        newdata = datos_test_prep,
                        type = "prob")

curva_roc <- roc(response = datos_test_prep$pb, 
                 predictor = predicciones$pres) 

# tr <- threshold(curva_roc, 'spec_sens')
# Gr?fico de la curva
# plot(curva_roc)

d= cbind(curva_roc$specificities, curva_roc$sensitivities)
colnames(d) <- c("specificities","sensitivities")
d <- as.data.frame(d)
# pdf("ROC_v5.pdf")
p2 <- ggplot(d,aes(specificities, sensitivities), )+
  geom_line()+
  geom_point() + scale_x_reverse()
p2
# dev.off()

AUC_KNN_test <- auc(response = datos_test_prep$pb, 
                    predictor = predicciones$pres) 
AUC_KNN_test






knn_prec  <-  raster::predict(predicted,modelo_knn, progress="text",type = "prob", 
                              args=c("extrapolate=T", "doclamp=TRUE"))

knn_prec_corrected <- 1 - knn_prec


{r fig.asp = 1, fig.width = 9}

# write.csv(RFpm,"RFpm.csv", row.names = FALSE)

### graficando ggplot2
mapRF.ggplot2 <- knn_prec_corrected %>% as("SpatialPixelsDataFrame") %>% 
  as.data.frame()
b2 <- ggplot() + geom_tile(data = mapRF.ggplot2, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis(option="turbo") + xlab("Longitud") + ylab("Latitud")+ 
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_nautical) + 
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "gray45"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "gray45"),
        panel.background = element_rect(fill = "gray54",colour = "gray85",size = 0.5, linetype = "solid"))+ labs(fill = "Probabilidad") + 
  labs(subtitle= "Modelamiento de la distribución de"~~italic("Rubus sp.")~~"",title="K-Nearest Neighbor (kNN)")
b2



#Naive Bayes


# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN

particiones  <- 5
repeticiones <- 10

# Hiperparámetros
hiperparametros <- data.frame(usekernel = FALSE, fL = 0 , adjust = 0)

set.seed(123)
seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)
for (i in 1:(particiones * repeticiones)) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO

control_train <- trainControl(method = "repeatedcv", number = particiones,
                              repeats = repeticiones, seeds = seeds,
                              returnResamp = "final", verboseIter = FALSE,
                              allowParallel = TRUE)

# AJUSTE DEL MODELO

set.seed(342)
modelo_nb <- train(pb ~ ., data = datos_train_prep,
                   method = "nb",
                   tuneGrid = hiperparametros,
                   metric = "Accuracy",
                   trControl = control_train)
modelo_nb




# TRAINING : Error train

predicciones <- predict(object = modelo_nb,
                        newdata = datos_train_prep,
                        type = "prob")

AUC_NB_train <- auc(response = datos_train_prep$pb, 
                    predictor = predicciones$pres) 
AUC_NB_train




# TESTING : Error test

predicciones_raw <- predict(modelo_nb, newdata = datos_test_prep,
                            type = "raw")
confusionMatrix <- confusionMatrix(data = predicciones_raw, reference = datos_test_prep$pb,
                                   positive = "pres")

confusionMatrix



precision(data = predicciones_raw, reference = datos_test_prep$pb,
          positive = "pres")



# Error de test
error_test <- mean(predicciones_raw != datos_test_prep$pb)
paste("Model Error Test", round(error_test*100, 2), "%")




# ROC
# C?lculo de la curva
predicciones <- predict(object = modelo_nb,
                        newdata = datos_test_prep,
                        type = "prob")

curva_roc <- roc(response = datos_test_prep$pb, 
                 predictor = predicciones$pres) 

# tr <- threshold(curva_roc, 'spec_sens')
# Gr?fico de la curva
# plot(curva_roc)

d= cbind(curva_roc$specificities, curva_roc$sensitivities)
colnames(d) <- c("specificities","sensitivities")
d <- as.data.frame(d)
# pdf("ROC_v5.pdf")
p2 <- ggplot(d,aes(specificities, sensitivities), )+
  geom_line()+
  geom_point() + scale_x_reverse()
p2
# dev.off()

AUC_NB_test <- auc(response = datos_test_prep$pb, 
                   predictor = predicciones$pres) 
AUC_NB_test







nb_prec  <-  raster::predict(predicted,modelo_nb, progress="text",type = "prob", 
                             args=c("extrapolate=T", "doclamp=TRUE"))

nb_prec_corrected <- 1 - nb_prec



# write.csv(RFpm,"RFpm.csv", row.names = FALSE)

### graficando ggplot2
mapRF.ggplot2 <- nb_prec_corrected %>% as("SpatialPixelsDataFrame") %>% 
  as.data.frame()
c2 <- ggplot() + geom_tile(data = mapRF.ggplot2, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis(option="turbo") + xlab("Longitud") + ylab("Latitud")+ 
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_nautical) + 
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "gray45"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "gray45"),
        panel.background = element_rect(fill = "gray54",colour = "gray85",size = 0.5, linetype = "solid"))+ labs(fill = "Probabilidad") + 
  labs(subtitle= "Modelamiento de la distribución de"~~italic("Rubus sp.")~~"",title="Naive Bayes")
c2


#Regresión logística


# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN

particiones  <- 5
repeticiones <- 10

# Hiperparámetros
hiperparametros <- data.frame(parameter = "none")

set.seed(123)
seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)
for (i in 1:(particiones * repeticiones)) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO

control_train <- trainControl(method = "repeatedcv", number = particiones,
                              repeats = repeticiones, seeds = seeds,
                              returnResamp = "final", verboseIter = FALSE,
                              allowParallel = TRUE)

# AJUSTE DEL MODELO

set.seed(342)
modelo_logistic <- train(pb ~ ., data = datos_train_prep,
                         method = "glm",
                         tuneGrid = hiperparametros,
                         metric = "Accuracy",
                         trControl = control_train,
                         family = "binomial")
modelo_logistic




# TRAINING : Error train

predicciones <- predict(object = modelo_logistic,
                        newdata = datos_train_prep,
                        type = "prob")

AUC_Log_train <- auc(response = datos_train_prep$pb, 
                     predictor = predicciones$pres) 
AUC_Log_train




# TESTING : Error test

predicciones_raw <- predict(modelo_logistic, newdata = datos_test_prep,
                            type = "raw")
confusionMatrix <- confusionMatrix(data = predicciones_raw, reference = datos_test_prep$pb,
                                   positive = "pres")

confusionMatrix



precision(data = predicciones_raw, reference = datos_test_prep$pb,
          positive = "pres")



# Error de test
error_test <- mean(predicciones_raw != datos_test_prep$pb)
paste("Model Error Test", round(error_test*100, 2), "%")



# ROC
# C?lculo de la curva
predicciones <- predict(object = modelo_logistic,
                        newdata = datos_test_prep,
                        type = "prob")

curva_roc <- roc(response = datos_test_prep$pb, 
                 predictor = predicciones$pres) 

# tr <- threshold(curva_roc, 'spec_sens')
# Gr?fico de la curva
# plot(curva_roc)

d= cbind(curva_roc$specificities, curva_roc$sensitivities)
colnames(d) <- c("specificities","sensitivities")
d <- as.data.frame(d)
# pdf("ROC_v5.pdf")
p2 <- ggplot(d,aes(specificities, sensitivities), )+
  geom_line()+
  geom_point() + scale_x_reverse()
p2
# dev.off()

AUC_Log_test <- auc(response = datos_test_prep$pb, 
                    predictor = predicciones$pres) 
AUC_Log_test






log_prec  <-  raster::predict(predicted,modelo_logistic, progress="text",type = "prob", 
                              args=c("extrapolate=T", "doclamp=TRUE"))

log_prec_corrected <- 1 - log_prec




# write.csv(RFpm,"RFpm.csv", row.names = FALSE)

### graficando ggplot2
mapRF.ggplot2 <- log_prec_corrected %>% as("SpatialPixelsDataFrame") %>% 
  as.data.frame()
d2 <- ggplot() + geom_tile(data = mapRF.ggplot2, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis(option="turbo") + xlab("Longitud") + ylab("Latitud")+ 
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_nautical) + 
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "gray45"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "gray45"),
        panel.background = element_rect(fill = "gray54",colour = "gray85",size = 0.5, linetype = "solid"))+ labs(fill = "Probabilidad") + 
  labs(subtitle= "Modelamiento de la distribución de"~~italic("Rubus sp.")~~"",title="Regresión logística")
d2




LDA


# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN

particiones  <- 5
repeticiones <- 10

# Hiperparámetros
hiperparametros <- data.frame(parameter = "none")

set.seed(123)
seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)
for (i in 1:(particiones * repeticiones)) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO

control_train <- trainControl(method = "repeatedcv", number = particiones,
                              repeats = repeticiones, seeds = seeds,
                              returnResamp = "final", verboseIter = FALSE,
                              allowParallel = TRUE)

# AJUSTE DEL MODELO

set.seed(342)
modelo_lda <- train(pb ~ ., data = datos_train_prep,
                    method = "lda",
                    tuneGrid = hiperparametros,
                    metric = "Accuracy",
                    trControl = control_train)
modelo_lda



# TRAINING : Error train

predicciones <- predict(object = modelo_lda,
                        newdata = datos_train_prep,
                        type = "prob")

AUC_LDA_train <- auc(response = datos_train_prep$pb, 
                     predictor = predicciones$pres) 
AUC_LDA_train




# TESTING : Error test

predicciones_raw <- predict(modelo_lda, newdata = datos_test_prep,
                            type = "raw")
confusionMatrix <- confusionMatrix(data = predicciones_raw, reference = datos_test_prep$pb,
                                   positive = "pres")

confusionMatrix



precision(data = predicciones_raw, reference = datos_test_prep$pb,
          positive = "pres")



# Error de test
error_test <- mean(predicciones_raw != datos_test_prep$pb)
paste("Model Error Test", round(error_test*100, 2), "%")



# ROC
# C?lculo de la curva
predicciones <- predict(object = modelo_lda,
                        newdata = datos_test_prep,
                        type = "prob")

curva_roc <- roc(response = datos_test_prep$pb, 
                 predictor = predicciones$pres) 

# tr <- threshold(curva_roc, 'spec_sens')
# Gr?fico de la curva
# plot(curva_roc)

d= cbind(curva_roc$specificities, curva_roc$sensitivities)
colnames(d) <- c("specificities","sensitivities")
d <- as.data.frame(d)
# pdf("ROC_v5.pdf")
p2 <- ggplot(d,aes(specificities, sensitivities), )+
  geom_line()+
  geom_point() + scale_x_reverse()
p2
# dev.off()

AUC_LDA_test <- auc(response = datos_test_prep$pb, 
                    predictor = predicciones$pres) 
AUC_LDA_test





lda_prec  <-  raster::predict(predicted,modelo_lda, progress="text",type = "prob", 
                              args=c("extrapolate=T", "doclamp=TRUE"))

lda_prec_corrected <- 1 - lda_prec




# write.csv(RFpm,"RFpm.csv", row.names = FALSE)

### graficando ggplot2
mapRF.ggplot2 <- lda_prec_corrected %>% as("SpatialPixelsDataFrame") %>% 
  as.data.frame()
e2 <- ggplot() + geom_tile(data = mapRF.ggplot2, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis(option="turbo") + xlab("Longitud") + ylab("Latitud")+ 
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_nautical) + 
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "gray45"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "gray45"),
        panel.background = element_rect(fill = "gray54",colour = "gray85",size = 0.5, linetype = "solid"))+ labs(fill = "Probabilidad") + 
  labs(subtitle= "Modelamiento de la distribución de"~~italic("Rubus sp.")~~"",title="LDA")
e2


Gradient Boosting


# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN

particiones  <- 5
repeticiones <- 10

# Hiperparámetros
hiperparametros <- expand.grid(interaction.depth = c(1, 2),
                               n.trees = c(500, 1000, 2000),
                               shrinkage = c(0.001, 0.01, 0.1),
                               n.minobsinnode = c(2, 5, 15))

set.seed(123)
seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)
for (i in 1:(particiones * repeticiones)) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO

control_train <- trainControl(method = "repeatedcv", number = particiones,
                              repeats = repeticiones, seeds = seeds,
                              returnResamp = "final", verboseIter = FALSE,
                              allowParallel = TRUE)

# AJUSTE DEL MODELO

set.seed(342)
modelo_boost <- train(pb ~ ., data = datos_train_prep,
                      method = "gbm",
                      tuneGrid = hiperparametros,
                      metric = "Accuracy",
                      trControl = control_train,
                      # Número de árboles ajustados
                      distribution = "adaboost",
                      verbose = FALSE)
modelo_boost



# REPRESENTACIÓN GRÁFICA

ggplot(modelo_boost, highlight = TRUE) +
  labs(title = "Evolución del accuracy del modelo Gradient Boosting") +
  guides(color = guide_legend(title = "shrinkage"),
         shape = guide_legend(title = "shrinkage")) +
  theme_bw() +
  theme(legend.position = "bottom")




# TRAINING : Error train

predicciones <- predict(object = modelo_boost,
                        newdata = datos_train_prep,
                        type = "prob")

AUC_GBM_train <- auc(response = datos_train_prep$pb, 
                     predictor = predicciones$pres) 
AUC_GBM_train




# TESTING : Error test

predicciones_raw <- predict(modelo_rf, newdata = datos_test_prep,
                            type = "raw")
confusionMatrix <- confusionMatrix(data = predicciones_raw, reference = datos_test_prep$pb,
                                   positive = "pres")

confusionMatrix



precision(data = predicciones_raw, reference = datos_test_prep$pb,
          positive = "pres")



# Error de test
error_test <- mean(predicciones_raw != datos_test_prep$pb)
paste("Model Error Test", round(error_test*100, 2), "%")



# ROC
# C?lculo de la curva
predicciones <- predict(object = modelo_boost,
                        newdata = datos_test_prep,
                        type = "prob")

curva_roc <- roc(response = datos_test_prep$pb, 
                 predictor = predicciones$pres) 

# tr <- threshold(curva_roc, 'spec_sens')
# Gr?fico de la curva
# plot(curva_roc)

d= cbind(curva_roc$specificities, curva_roc$sensitivities)
colnames(d) <- c("specificities","sensitivities")
d <- as.data.frame(d)
# pdf("ROC_v5.pdf")
p2 <- ggplot(d,aes(specificities, sensitivities), )+
  geom_line()+
  geom_point() + scale_x_reverse()
p2
# dev.off()

AUC_GBM_test <- auc(response = datos_test_prep$pb, 
                    predictor = predicciones$pres) 
AUC_GBM_test





boost_prec  <-  raster::predict(predicted,modelo_boost, progress="text",type = "prob", 
                                args=c("extrapolate=T", "doclamp=TRUE"))

boost_prec_corrected <- 1 - boost_prec




# write.csv(RFpm,"RFpm.csv", row.names = FALSE)

### graficando ggplot2
mapRF.ggplot2 <- boost_prec_corrected %>% as("SpatialPixelsDataFrame") %>% 
  as.data.frame()
f2 <- ggplot() + geom_tile(data = mapRF.ggplot2, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis(option="turbo") + xlab("Longitud") + ylab("Latitud")+ 
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_nautical) + 
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "gray45"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "gray45"),
        panel.background = element_rect(fill = "gray54",colour = "gray85",size = 0.5, linetype = "solid"))+ labs(fill = "Probabilidad") + 
  labs(subtitle= "Modelamiento de la distribución de"~~italic("Rubus sp.")~~"",title="Gradient Boosting")
f2



# Se trasforma el dataframe devuelto por resamples() para separar el nombre del
# modelo y las métricas en columnas distintas.
metricas_resamples <- resultados_resamples$values %>%
  gather(key = "modelo", value = "valor", -Resample) %>%
  separate(col = "modelo", into = c("modelo", "metrica"),
           sep = "~", remove = TRUE)
metricas_resamples %>% head()



metricas_resamples %>% 
  group_by(modelo, metrica) %>% 
  summarise(media = mean(valor)) %>%
  spread(key = metrica, value = media) %>%
  arrange(desc(Accuracy))



metricas_resamples %>%
  filter(metrica == "Accuracy") %>%
  group_by(modelo) %>% 
  summarise(media = mean(valor)) %>%
  ggplot(aes(x = reorder(modelo, media), y = media, label = round(media, 2))) +
  geom_segment(aes(x = reorder(modelo, media), y = 0,
                   xend = modelo, yend = media),
               color = "grey50") +
  geom_point(size = 7, color = "firebrick") +
  geom_text(color = "white", size = 2.5) +
  scale_y_continuous(limits = c(0, 1)) +
  # Accuracy basal
  geom_hline(yintercept = 0.62, linetype = "dashed") +
  annotate(geom = "text", y = 0.72, x = 8.5, label = "Accuracy basal") +
  labs(title = "Validación: Accuracy medio repeated-CV",
       subtitle = "Modelos ordenados por media",
       x = "modelo") +
  coord_flip() +
  theme_bw()



metricas_resamples %>% filter(metrica == "Accuracy") %>%
  group_by(modelo) %>% 
  mutate(media = mean(valor)) %>%
  ungroup() %>%
  ggplot(aes(x = reorder(modelo, media), y = valor, color = modelo)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  scale_y_continuous(limits = c(0, 1)) +
  # Accuracy basal
  geom_hline(yintercept = 0.62, linetype = "dashed") +
  annotate(geom = "text", y = 0.65, x = 8.5, label = "Accuracy basal") +
  theme_bw() +
  labs(title = "Validación: Accuracy medio repeated-CV",
       subtitle = "Modelos ordenados por media") +
  coord_flip() +
  theme(legend.position = "none")


# Test de Friedman para comparar el accuracy de los modelos


matriz_metricas <- metricas_resamples %>% filter(metrica == "Accuracy") %>%
  spread(key = modelo, value = valor) %>%
  select(-Resample, -metrica) %>% as.matrix()
friedman.test(y = matriz_metricas)



# Comparaciones múltiples con un test suma de rangos de Wilcoxon


metricas_accuracy <- metricas_resamples %>% filter(metrica == "Accuracy")
comparaciones  <- pairwise.wilcox.test(x = metricas_accuracy$valor, 
                                       g = metricas_accuracy$modelo,
                                       paired = TRUE,
                                       p.adjust.method = "holm")

# Se almacenan los p_values en forma de dataframe
comparaciones <- comparaciones$p.value %>%
  as.data.frame() %>%
  rownames_to_column(var = "modeloA") %>%
  gather(key = "modeloB", value = "p_value", -modeloA) %>%
  na.omit() %>%
  arrange(modeloA) 

comparaciones


# Error de test


predicciones <- extractPrediction(
  models = modelos,
  testX = datos_test_prep %>% select(-pb),
  testY = datos_test_prep$pb
)
predicciones %>% head()



metricas_predicciones <- predicciones %>%
  mutate(acierto = ifelse(obs == pred, TRUE, FALSE)) %>%
  group_by(object, dataType) %>%
  summarise(accuracy = mean(acierto))

metricas_predicciones %>%
  spread(key = dataType, value = accuracy) %>%
  arrange(desc(Test))



ggplot(data = metricas_predicciones,
       aes(x = reorder(object, accuracy), y = accuracy,
           color = dataType, label = round(accuracy, 2))) +
  geom_point(size = 8) +
  scale_color_manual(values = c("orangered2", "gray50")) +
  geom_text(color = "white", size = 3) +
  scale_y_continuous(limits = c(0, 1)) +
  # Accuracy basal
  geom_hline(yintercept = 0.62, linetype = "dashed") +
  annotate(geom = "text", y = 0.66, x = 8.5, label = "Accuracy basal") +
  coord_flip() +
  labs(title = "Accuracy de entrenamiento y test", 
       x = "modelo") +
  theme_bw() + 
  theme(legend.position = "bottom")

comparaciones %>% filter((modeloA == "rf") | (modeloA == "SVMradial" & modeloB == "rf"))



# AUC ---------------------------------------------------------------------



Value <- rbind(AUC_RF_test[1],AUC_RF_train[1],
               AUC_GBM_test[1],AUC_GBM_train[1],
               AUC_KNN_test[1],AUC_KNN_train[1],
               AUC_LDA_test[1],AUC_LDA_train[1],
               AUC_Log_test[1],AUC_Log_train[1],
               AUC_NB_test[1],AUC_NB_train[1],"0.8162545","0.9078527")
Model <- c("RF","RF","GBM","GBM","KNN","KNN","LDA","LDA","Log","Log","NB","NB","MaxEnt","MaxEnt")
Partition <- c("Testing","Training","Testing","Training","Testing","Training","Testing","Training","Testing","Training","Testing","Training","Testing","Training")
AUC_final <- cbind(Value,Model,Partition)
colnames(AUC_final) <- c("Value","Model","Partition")
AUC_final <- as.data.frame(AUC_final)
AUC_final$Value <- as.numeric(AUC_final$Value)
AUC_finalB <- AUC_final
# Guardar la tabla en un archivo CSV
write.csv(AUC_finalB, file = "AUC_finalLandrace.csv", row.names = FALSE)




ggplot(data = AUC_final,
       aes(x = reorder(Model, Value), y = Value,
           color = Partition, label = round(Value, 2))) +
  geom_point(size = 8) +
  scale_color_manual(values = c("#01665e", "#8c96c6")) +
  geom_text(color = "white", size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  # Accuracy basal
  # geom_hline(yintercept = 0.8, linetype = "dashed", colour="grey") +
  # annotate(geom = "text", y = 0.7, x = 6, label = "AUC basal") +
  coord_flip() +
  labs(title = "", 
       x = "Machine Learning Model", y="Area Under the Curve (AUC)") +
  theme_bw() + 
  labs(subtitle= ""~~italic("Polylepis sericea")~~"") + labs(title = "B")



# Post-modelamiento -------------------------------------------------------


Occur <- read_delim("Datos presencia/name_all_thin1.csv", delim = ',', col_names = T) |> as.data.frame() |> na.omit()
Occur <-  Occur[,-1]
ExtRast <- terra::extract(Prediction$lyr1, Occur) #Extract values from occs points
Occur1 <- sort(ExtRast$lyr1) #numeric order
PercentThreshold <- 0.3 #Excluding 80% of total presence data
RclVal  <-  Occur1[round(length(Occur1) * PercentThreshold) + 1]
RclVal


Current <- Prediction$lyr1 

M_b <- M >= RclVal
plot(M_b)
stack(M_b)
 raster::writeRaster(M_b, filename = "Projection_threshold/current_Total.tif", overwrite=T)


# Future mean threshold 

#Total area
GCM <- c('MI', 'MP', "MR")
SSP <- c('2', '3')
year <- c('50', '70')

for (GC in 1: length(GCM)) {
  for (Sp in 1: length(SSP)) {
    for (ye in 1: length(year)) {
      
      pat = paste0(SSP[Sp],'_',year[ye],'.asc') 
      
      L <- rast(list.files(path="Final_models/M_0.1_F_lq_Set_6_E/", pattern =pat, full.names = T)) 
      
      M <- mean(L)
      
      M_b <- M >= RclVal
      
      stack(M_b)
      
      raster::writeRaster(M_b, filename = paste0("Projection_threshold/",SSP[Sp],'_',year[ye],".tif"))
    }
  } 
}



# Suitability evaluation 

##Import shapefile South_Am
# shape <- 'Vectores/Lim_Nari.shp'
# Narino = vect(shape)
# crs(Narino) <- '+proj=longlat +datum=WGS84 +no_defs'


# Create current suitability shape ----------------------------------------

#Create niche current shape 
tif <- read_stars("Projection_threshold/current_Total.tif")
sf <- st_as_sf(tif, merge = T)
suit <- sf[sf$current_Total.tif == 1, ]

st_write(suit, 'Vectores/current_suitability_Total.shp')

shape <- 'Vectores/current_suitability_Total.shp' 
ca = vect(shape)

# Creating map suitability 

#Load shapefile
Narino <- st_read('Vectores/Lim_Nari.shp')


#Select shapefile

#Select total shapefile
#Total <- st_read('Final_maps/Vectors/current_suitability_Total.shp')
Bcocke <- st_read('Vectores/current_suitability_Total.shp')
sf <- st_as_sf(Bcocke, merge = T)
plot(sf)


# Niche map ---------------------------------------------------------------



# Proyecciones futuras ----------------------------------------------------

  wcF <- list.files("D:/AGROSAVIA_ 2020 _BUHO/PAPERS_2020/Batata ML/Analisis2023/2040ssp585", pattern = 'bio', full.names = TRUE)
envirem1f <- stack(wcFiles)
Predictionf <- predict(enviremf, BestModels, type = "cloglog")

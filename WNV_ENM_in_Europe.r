# Note on the environmental rasters:
	# - 4 models for the past and present projections (obsclim and counterclim): GSWP3-W5E5, 20CRv3, 20CRv3-ERA5, and 20CRv3-W5E5;
	#		and 6 periods: 1901-1919, 1920-1939, 1940-1959, 1960-1979, 1980-1999, and 2000-2019 for temperature, precipitation and
	#		relative humidity. Population and land cover data do not change between models
	# - 10 models for the future projections (bias-adjusted): MRI-ESM2-0, MPI-ESM1-2-HR, MIROC6, IPSL-CM6A-LR, GFDL-ESM4, EC-Earth3,
	#		CNRM-ESM2-1, CNRM-CM6-1, CanESM5, and UKESM1-0-LL; and 5 periods: 1995-2014, 2020-2039, 2040-2059, 2060-2079, and 2080-2099
	# - climatic variables: "hur" = near surface relative humidity, "pr" = precipitation, and "tas" = near surface air temperature;
	#		units: relative humidity in % (??, some values >100), temperature in Kelvins (°C+273.15), and precipitation in kg/m2/second

library(blockCV)
library(dismo)
library(exactextractr)
library(gbm)
library(geosphere)
library(lubridate)
library(maptools)
library(ncdf4)
library(ncf)
library(raster)
library(RColorBrewer)
library(rgdal)
library(rgeos)
library(sf)
library(vioplot)

savingPlots = FALSE
showingPlots = FALSE

# 1. Preparing the data frame with the WNV occurrence data + visualisation

nuts3 = crop(shapefile("MOOD_NUTs3_WillyW/VectornetDATAforMOODjan21.shp"), extent(-10,33,34.5,72)) # NUT3 polygons shapefile
nutsM = crop(shapefile("MOOD_NUTs3_WillyW/VectornetMAPforMOODjan21.shp"), extent(-10,33,34.5,72)) # modified NUT3 shapefile
nuts3 = subset(nuts3, !CountryISO%in%c("MA","DZ","TN","MT","TR","CI","MD","UA","RU","FO","IS","GE","BY","IS","GL","FO","CY","SJ"))
nutsM = subset(nutsM, !CountryISO%in%c("MA","DZ","TN","MT","TR","CI","MD","UA","RU","FO","IS","GE","BY","IS","GL","FO","CY","SJ"))
correspondences = shapefile("MOOD_NUTs3_WillyW/vnMOODdatamapcodejoinpoly.shp")@data[,c("DATLOCODE","MAPLOCODE")]
missingIDs1 = nuts3@data[which(!nuts3@data[,"LocationCo"]%in%correspondences[,"DATLOCODE"]),"LocationCo"]
missingIDs2 = nutsM@data[which(!nutsM@data[,"LocationCo"]%in%correspondences[,"MAPLOCODE"]),"LocationCo"]
data = read.csv("ECDC_WNV_cases.csv", head=T); dates = decimal_date(ymd(data[,"DateOfOnsetISOdate"]))
data = data[which(dates<2020),]; data = data[which(data[,"Imported"]!="Y"),]
data[which(data[,"PlaceOfInfectionEVD"]=="FR813"),] = "FRJ13"; data[which(data[,"PlaceOfNotification"]=="FR813"),] = "FRJ13"
data[which(data[,"PlaceOfInfectionEVD"]=="FR823"),] = "FRL03"; data[which(data[,"PlaceOfNotification"]=="FR823"),] = "FRL03"
data[which(data[,"PlaceOfInfectionEVD"]=="FR824"),] = "FRL04"; data[which(data[,"PlaceOfNotification"]=="FR824"),] = "FRL04"
data[which(data[,"PlaceOfInfectionEVD"]=="FR826"),] = "FRL06"; data[which(data[,"PlaceOfNotification"]=="FR826"),] = "FRL06"
data[which(data[,"PlaceOfInfectionEVD"]=="FR831"),] = "FRM01"; data[which(data[,"PlaceOfNotification"]=="FR831"),] = "FRM01"
data[which(data[,"PlaceOfInfectionEVD"]=="HU101"),] = "HU110"; data[which(data[,"PlaceOfNotification"]=="HU101"),] = "HU110"
data[which(data[,"PlaceOfInfectionEVD"]=="HU102"),] = "HU120"; data[which(data[,"PlaceOfNotification"]=="HU102"),] = "HU120"
	# "Classification" = "PROB": any person meeting the clinical criteria AND with at least one of the following two: 
					 #   (i) an epidemiological link, (ii) a laboratory test for a probable case
					 # = "CONF": any person meeting the laboratory criteria for case confirmation
					 # --> not assigning a "presence" to a polygon only associated with "PROB" cases		   
nuts3_data = as.data.frame(matrix(nrow=dim(nuts3@data)[1], ncol=4)); colnames(nuts3_data) = c("country","NUT3","observations","atLeastOneConfirmed")
nutsM_data = as.data.frame(matrix(nrow=dim(nutsM@data)[1], ncol=3)); colnames(nutsM_data) = c("country","NUTM","observations")
nuts3_data$country = nuts3@data[,"CountryNam"]; nuts3_data$NUT3 = nuts3@data[,"LocationCo"]; nuts3_data$observations = 0; nuts3_data$atLeastOneConfirmed = FALSE
nutsM_data$country = nutsM@data[,"CountryNam"]; nutsM_data$NUTM = nutsM@data[,"LocationCo"]; nutsM_data$observations = 0
missingIDs3 = unique(data[which((!data[,"PlaceOfInfectionEVD"]%in%nuts3_data[,"NUT3"])&(!data[,"PlaceOfNotification"]%in%nuts3_data[,"NUT3"])),"PlaceOfNotification"])
	# Missing IDs: "ES61", "ITF1", "ITH3", "ITG2", "EL30", "ITH4", "ITF5", and "AT13" --> do not correspond to NUT3 IDs
for (i in 1:dim(data)[1])
	{
		index = which(nuts3_data[,"NUT3"]==data[i,"PlaceOfInfectionEVD"])
		if (length(index) != 1) index = which(nuts3_data[,"NUT3"]==data[i,"PlaceOfNotification"])
		# if (length(index) != 1) print(c(i,data[i,"PlaceOfInfectionEVD"],data[i,"PlaceOfNotification"],data[i,"CountryName"]))
		if (length(index) == 1) nuts3_data[index,"observations"] = nuts3_data[index,"observations"] + 1
		if ((length(index) == 1)&(data[i,"Classification"]=="CONF")) nuts3_data[index,"atLeastOneConfirmed"] = TRUE
	}
nber_of_nut3_polygons_with_at_least_one_observation = sum(nuts3_data[,"observations"]>0) # 185
nber_of_nut3_polygons_with_at_least_one_confirmed_observation = sum((nuts3_data[,"observations"]>0)&(nuts3_data[,"atLeastOneConfirmed"]==TRUE)) # 169
nuts3_data[which(nuts3_data[,"atLeastOneConfirmed"]==FALSE),"observations"] = 0
for (i in 1:dim(nuts3_data)[1])
	{
		if (nuts3_data[i,"observations"] > 0)
			{
				if (nuts3_data[i,"NUT3"]%in%nutsM_data[,"NUTM"])
					{
						index = which(nutsM_data[,"NUTM"]==nuts3_data[i,"NUT3"])
						if (length(index) != 1)
							{
								print(i)
							}	else	{
								nutsM_data[index,"observations"] = nutsM_data[index,"observations"]+nuts3_data[i,"observations"]
							}
					}	else	{
						# print(nuts3_data[i,"country"])
						nutsM_ID = correspondences[which(correspondences[,"DATLOCODE"]==nuts3_data[i,"NUT3"]),"MAPLOCODE"]
						if ((length(nutsM_ID) > 0)&&(!is.na(nutsM_ID)))
							{
								index = which(nutsM_data[,"NUTM"]==nutsM_ID)
								nutsM_data[index,"observations"] = nutsM_data[index,"observations"]+nuts3_data[i,"observations"]
							}	else	{			
								print(nuts3_data[i,"NUT3"])
							}
					}
			}
	}
if (sum(nuts3_data[,"observations"]) != sum(nutsM_data[,"observations"])) print("Houston, we have a problem!")
nuts3_data[which(nuts3_data[,"country"]%in%c("Bosnia and Herzegovina","Switzerland")),"observations"] = NA
nutsM_data[which(nutsM_data[,"country"]%in%c("Bosnia and Herzegovina","Switzerland")),"observations"] = NA
contour = unionSpatialPolygons(nutsM, rep(1,length(nutsM)))
if (savingPlots)
	{
		pdf("WNV_cases_map1_NEW.pdf", width=8, height=5.7) # 1° version
		par(mfrow=c(1,2), oma=c(0.1,0.1,0.1,0.1), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		colourScale = colorRampPalette(brewer.pal(9,"YlOrBr"))(131)[21:121]
		cases_per_km2 = nuts3_data[,"observations"]/area(nuts3,unit="km")
		cases_per_km2[which(cases_per_km2>(10^-8))] = 10^-8
		cols = colourScale[((cases_per_km2/(10^-8))*100)+1]
		cols[which(cases_per_km2==0)] = "gray90"
		cols[which(is.na(cases_per_km2))] = "gray80"
		plot(contour, lwd=0.4, border="gray30", col=NA)
		plot(nuts3, col=cols, border="gray50", lwd=0.1, add=T)
		rast = raster(as.matrix(c(0,10^-8)))
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.080,0.095,0.70,0.96),
			 adj=3, axis.args=list(at=seq(2*(10^-9),10^-8,2*(10^-9)), labels=c("2 E-9","4 E-9","6 E-9","8 E-9","E-8"), cex.axis=0.6,
			 lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.7, col.axis="gray30", line=0, mgp=c(0,0.4,0)), alpha=1, side=3)
		cases_per_km2 = nutsM_data[,"observations"]/area(nutsM,unit="km")
		cases_per_km2[which(cases_per_km2>(10^-8))] = 10^-8
		cols = colourScale[((cases_per_km2/(10^-8))*100)+1]
		cols[which(cases_per_km2==0)] = "gray90"
		cols[which(is.na(cases_per_km2))] = "gray80"
		plot(contour, lwd=0.4, border="gray30", col=NA)
		plot(nutsM, col=cols, border="gray50", lwd=0.1, add=T)
		dev.off()

		population = brick(paste0("Environmental_rasters/ISIMIP3a/population_histsoc_0p5deg_annual_2000_2019_timmean.nc4"), varname="total-population")
		population_log_nuts3 = rep(NA, dim(nuts3@data)[1]); population_log_nutsM = rep(NA, dim(nutsM@data)[1])
		population_nuts3 = rep(NA, dim(nuts3@data)[1]); population_nutsM = rep(NA, dim(nutsM@data)[1])
		for (j in 1:length(nuts3))
			{
				maxArea = 0; polIndex = 0
				for (k in 1:length(nuts3@polygons[[j]]@Polygons))
					{
						if (maxArea < nuts3@polygons[[j]]@Polygons[[k]]@area)
							{
								maxArea = nuts3@polygons[[j]]@Polygons[[k]]@area; polIndex = k
							}
					}
				pol1 = nuts3@polygons[[j]]@Polygons[[polIndex]]; p = Polygon(pol1@coords)
				ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
				pol2 = sf::st_as_sfc(sps); st_crs(pol2) = crs(nutsM)
				if (j == 1) crs(population) = crs(pol2) 
				population_log_nuts3[j] = log10(exactextractr::exact_extract(population,pol2,fun="sum")+1)
				population_nuts3[j] = exactextractr::exact_extract(population,pol2,fun="sum")
			}
		for (j in 1:length(nutsM))
			{
				maxArea = 0; polIndex = 0
				for (k in 1:length(nutsM@polygons[[j]]@Polygons))
					{
						if (maxArea < nutsM@polygons[[j]]@Polygons[[k]]@area)
							{
								maxArea = nutsM@polygons[[j]]@Polygons[[k]]@area; polIndex = k
							}
					}
				pol1 = nutsM@polygons[[j]]@Polygons[[polIndex]]; p = Polygon(pol1@coords)
				ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
				pol2 = sf::st_as_sfc(sps); st_crs(pol2) = crs(nutsM)
				if (j == 1) crs(population) = crs(pol2) 
				population_log_nutsM[j] = log10(exactextractr::exact_extract(population,pol2,fun="sum")+1)
				population_nutsM[j] = exactextractr::exact_extract(population,pol2,fun="sum")
			}
		pdf("WNV_cases_map2_NEW.pdf", width=8, height=5.0) # 2° version
		par(mfrow=c(2,4), oma=c(0,0,0,0.6), mar=c(0,0,0,1.5), lwd=0.2, col="gray30")		
		colourScale = colorRampPalette(brewer.pal(9,"YlOrBr"))(131)[21:121]
		cases_per_km2 = nuts3_data[,"observations"]/area(nuts3,unit="km")
		cases_per_km2[which(cases_per_km2>(10^-8))] = 10^-8
		cols = colourScale[((cases_per_km2/(10^-8))*100)+1]
		cols[which(cases_per_km2==0)] = "gray90"
		cols[which(is.na(cases_per_km2))] = "gray70"
		plot(contour, lwd=0.4, border="gray30", col=NA)
		plot(nuts3, col=cols, border="gray50", lwd=0.1, add=T)
		mtext(expression(bold(A)), side=3, line=-1.6, at=-9, cex=0.70, col="gray30")
		mtext("WNV cases/km2", side=3, line=-2.7, at=0, cex=0.45, col="gray30")
		cases_per_km2 = nutsM_data[,"observations"]/area(nutsM,unit="km")
		cases_per_km2[which(cases_per_km2>(10^-8))] = 10^-8
		cols = colourScale[((cases_per_km2/(10^-8))*100)+1]
		cols[which(cases_per_km2==0)] = "gray90"
		cols[which(is.na(cases_per_km2))] = "gray70"
		plot(contour, lwd=0.4, border="gray30", col=NA)
		plot(nutsM, col=cols, border="gray50", lwd=0.1, add=T)
		mtext("WNV cases/km2", side=3, line=-2.7, at=0, cex=0.45, col="gray30")
		rast = raster(as.matrix(c(0,10^-8)))
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.840,0.855,0.05,0.50), adj=3,
			 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
		colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(131)[21:121]
		cases_per_10E5people = (nuts3_data[,"observations"]/population_nuts3)*(10^5)
		cases_per_10E5people[which(cases_per_10E5people>(100))] = 100
		cols = colourScale[((cases_per_10E5people/100)*100)+1]
		cols[which(cases_per_10E5people==0)] = "gray90"
		cols[which(is.na(cases_per_10E5people))] = "gray70"
		plot(contour, lwd=0.4, border="gray30", col=NA)
		plot(nuts3, col=cols, border="gray50", lwd=0.1, add=T)
		mtext(expression(bold(B)), side=3, line=-1.6, at=-9, cex=0.70, col="gray30")
		mtext("WNV cases/105 people", side=3, line=-2.7, at=0, cex=0.45, col="gray30")
		cases_per_10E5people = (nutsM_data[,"observations"]/population_nutsM)*(10^5)
		cases_per_10E5people[which(cases_per_10E5people>(100))] = 100
		cols = colourScale[((cases_per_10E5people/100)*100)+1]
		cols[which(cases_per_10E5people==0)] = "gray90"
		cols[which(is.na(cases_per_10E5people))] = "gray70"
		plot(contour, lwd=0.4, border="gray30", col=NA)
		plot(nutsM, col=cols, border="gray50", lwd=0.1, add=T)
		mtext("WNV cases/105 people", side=3, line=-2.7, at=0, cex=0.45, col="gray30")
		rast = raster(as.matrix(c(0,100)))
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.840,0.855,0.05,0.50), adj=3,
			 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
		dev.off()
	}

# 2. Training the BRT models for all the subsequent projections

envVariableNames = c("temperature_winter","temperature_spring","temperature_summer","temperature_inFall",
				   "precipitation_winter","tprecipitation_spring","precipitation_summer","precipitation_inFall",
				   "relative_humidity_winter","relative_humidity_spring","relative_humidity_summer","relative_humidity_inFall",
				   "primary_forest_areas","primary_non-forest_areas","secondary_forest_areas","secondary_non-forest_areas",
				   "croplands_all_categories","managed_pasture_and_rangeland","human_pop_density_log10")
nutsM_datas = list(); models_isimip3a = c("gswp3-w5e5","20crv3","20crv3-era5","20crv3-w5e5")
models_isimip3a_names = c("GSWP3-W5E5", "20CRv3", "20CRv3-ERA5", "20CRv3-W5E5")
if (!file.exists("Training_dFrames.rds"))
	{
		i = 1 # for the future projections, we only train one BRT model
		for (i in 1:length(models_isimip3a))
			{
				environmentalValues = matrix(nrow=dim(nutsM_data)[1], ncol=length(envVariableNames)); colnames(environmentalValues) = envVariableNames
				temperature = brick(paste0("Environmental_rasters/ISIMIP3a/tas_day_obsclim_historical_",models_isimip3a[i],"_2000_2019_ymonmean.nc"))
				temperature_winter = mean(temperature[[12]],temperature[[1]],temperature[[2]])-273.15 # conversion to Celcius degrees
				temperature_spring = mean(temperature[[3]],temperature[[4]],temperature[[5]])-273.15 # conversion to Celcius degrees
				temperature_summer = mean(temperature[[6]],temperature[[7]],temperature[[8]])-273.15 # conversion to Celcius degrees
				temperature_inFall = mean(temperature[[9]],temperature[[10]],temperature[[11]])-273.15 # conversion to Celcius degrees
				precipitation = brick(paste0("Environmental_rasters/ISIMIP3a/pr_day_obsclim_historical_",models_isimip3a[i],"_2000_2019_ymonmean.nc"))
				precipitation_winter = mean(precipitation[[12]],precipitation[[1]],precipitation[[2]])*60*60*24 # conversion to kg/m2/day
				precipitation_spring = mean(precipitation[[3]],precipitation[[4]],precipitation[[5]])*60*60*24 # conversion to kg/m2/day
				precipitation_summer = mean(precipitation[[6]],precipitation[[7]],precipitation[[8]])*60*60*24 # conversion to kg/m2/day
				precipitation_inFall = mean(precipitation[[9]],precipitation[[10]],precipitation[[11]])*60*60*24 # conversion to kg/m2/day
				relative_humidity = brick(paste0("Environmental_rasters/ISIMIP3a/hurs_day_obsclim_historical_",models_isimip3a[i],"_2000_2019_ymonmean.nc"))
				relative_humidity_winter = mean(relative_humidity[[12]],relative_humidity[[1]],relative_humidity[[2]])
				relative_humidity_spring = mean(relative_humidity[[3]],relative_humidity[[4]],relative_humidity[[5]])
				relative_humidity_summer = mean(relative_humidity[[6]],relative_humidity[[7]],relative_humidity[[8]])
				relative_humidity_inFall = mean(relative_humidity[[9]],relative_humidity[[10]],relative_humidity[[11]])
				land_cover = nc_open(paste0("Environmental_rasters/ISIMIP3a/landcover_annual_2000_2015_timmean.nc4"))
				population = brick(paste0("Environmental_rasters/ISIMIP3a/population_histsoc_0p5deg_annual_2000_2019_timmean.nc4"), varname="total-population")
				landCoverVariableIDs = names(land_cover$var); land_covers1 = list(); land_covers2 = list(); land_covers3 = list()
				landCoverVariableNames = as.character(read.csv("Environmental_rasters/Luse.csv")[1:12,2])
				for (j in 2:13)
					{
						land_covers1[[j-1]] = brick(paste0("Environmental_rasters/ISIMIP3a/landcover_annual_2000_2015_timmean.nc4"), varname=landCoverVariableIDs[j])
					}
				variable_codes = c("croplands","pastures","urbanAreas","primaryForest","primaryNonF","secondaryForest","secondaryNonF")
				variable_names = c("crops","pasture","urban land","forested primary land","non-forested primary land",
						  		   "potentially forested secondary land","potentially non-forested secondary land")
				for (j in 1:length(variable_names))
					{
						names = gsub("\\."," ",landCoverVariableNames); indices = which(landCoverVariableNames==variable_names[j])
						if (length(indices) == 0) indices = which(grepl(variable_names[j],names))
						if (variable_names[j] == "pasture") indices = c(indices, which(grepl("rangeland",names)))
						land_cover = land_covers1[[indices[1]]]; names(land_cover) = variable_codes[j]; # print(indices)
						if (length(indices) > 1)
							{
								for (k in 2:length(indices)) land_cover[] = land_cover[]+land_covers1[[indices[k]]][]
							}
						land_covers2[[j]] = land_cover[[1]]; land_covers3[[j]] = raster::aggregate(land_cover[[1]],2)
					}
				envVariables = list()			
				envVariables[[1]] = temperature_winter; envVariables[[2]] = temperature_spring
				envVariables[[3]] = temperature_summer; envVariables[[4]] = temperature_inFall
				envVariables[[5]] = precipitation_winter; envVariables[[6]] = precipitation_spring
				envVariables[[7]] = precipitation_summer; envVariables[[8]] = precipitation_inFall
				envVariables[[9]] = relative_humidity_winter; envVariables[[10]] = relative_humidity_spring
				envVariables[[11]] = relative_humidity_summer; envVariables[[12]] = relative_humidity_inFall
				envVariables[[13]] = land_covers2[[4]] # primary forest areas
				envVariables[[14]] = land_covers2[[5]] # primary non-forest areas
				envVariables[[15]] = land_covers2[[6]] # secondary forest areas
				envVariables[[16]] = land_covers2[[7]] # secondary non-forest areas
				envVariables[[17]] = land_covers2[[1]] # croplands (all catergories)
				envVariables[[18]] = land_covers2[[2]] # managed pasture + rangeland
				envVariables[[19]] = population # human population (not log-transformed)
				for (j in 1:length(envVariables)) envVariables[[j]] = crop(envVariables[[j]], nutsM, snap="out")
				for (j in 1:length(envVariables)) envVariables[[j]] = mask(envVariables[[j]], nutsM)
				areas_nuts_M_km = area(nutsM, unit="km")
				for (j in 1:length(nutsM))
					{
						maxArea = 0; polIndex = 0
						for (k in 1:length(nutsM@polygons[[j]]@Polygons))
							{
								if (maxArea < nutsM@polygons[[j]]@Polygons[[k]]@area)
									{
										maxArea = nutsM@polygons[[j]]@Polygons[[k]]@area; polIndex = k
									}
							}
						pol1 = nutsM@polygons[[j]]@Polygons[[polIndex]]; p = Polygon(pol1@coords)
						ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
						pol2 = sf::st_as_sfc(sps); st_crs(pol2) = crs(nutsM)
						for (k in 1:18)
							{
								if ((j == 1)&(k == 1)) crs(envVariables[[k]]) = crs(pol2) 
								environmentalValues[j,k] = exactextractr::exact_extract(envVariables[[k]], pol2, fun="mean")
							}
						if ((j == 1)&(k == 1)) crs(envVariables[[19]]) = crs(pol2) 
						environmentalValues[j,19] = (log10(exactextractr::exact_extract(envVariables[[19]],pol2,fun="sum")+1))/areas_nuts_M_km[j]
					}
				nutsM_datas[[i]] = cbind(nutsM_data, environmentalValues)
				if ((i == 1)&(savingPlots))
					{
						pdf("All_env_factors_1_NEW.pdf", width=8, height=4.8); par(mfrow=c(2,4), oma=c(0,0,0,0.6), mar=c(0,0,0,0.8), lwd=0.2, col="gray30") # 1° version
						colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]
						cols = colourScale[(((environmentalValues[,3]-min(environmentalValues[,3],na.rm=T))/(max(environmentalValues[,3],na.rm=T)-min(environmentalValues[,3],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,3],na.rm=T),max(environmentalValues[,3],na.rm=T))))
						mtext("Air temperature", side=3, line=-1.5, at=0, cex=0.50, col="gray30"); mtext("(summer, °C)", side=3, line=-2.2, at=0, cex=0.50, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.870,0.885,0.05,0.95), adj=3,
							 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(101)
						cols = colourScale[(((environmentalValues[,7]-min(environmentalValues[,7],na.rm=T))/(max(environmentalValues[,7],na.rm=T)-min(environmentalValues[,7],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,7],na.rm=T),max(environmentalValues[,7],na.rm=T))))
						mtext("Precipitation", side=3, line=-1.5, at=0, cex=0.50, col="gray30"); mtext("(summer, kg/m2/day)", side=3, line=-2.2, at=0, cex=0.50, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.870,0.885,0.05,0.95), adj=3,
							 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
						colourScale = c("#E5E5E5",colorRampPalette(c("gray97","chartreuse4"),bias=1)(121)[21:121])
						cols = colourScale[(((environmentalValues[,13]-min(environmentalValues[,13],na.rm=T))/(max(environmentalValues[,13],na.rm=T)-min(environmentalValues[,13],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,13],na.rm=T),max(environmentalValues[,13],na.rm=T))))
						mtext("Primary forested", side=3, line=-1.5, at=0, cex=0.50, col="gray30"); mtext("areas", side=3, line=-2.2, at=0, cex=0.50, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.870,0.885,0.05,0.95), adj=3,
							 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
						colourScale = c("#E5E5E5",colorRampPalette(c("gray97","olivedrab3"),bias=1)(121)[21:121])
						cols = colourScale[(((environmentalValues[,15]-min(environmentalValues[,15],na.rm=T))/(max(environmentalValues[,15],na.rm=T)-min(environmentalValues[,15],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,15],na.rm=T),max(environmentalValues[,15],na.rm=T))))
						mtext("Secondary forested", side=3, line=-1.5, at=0, cex=0.50, col="gray30"); mtext("areas", side=3, line=-2.2, at=0, cex=0.50, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.870,0.885,0.05,0.95), adj=3,
							 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"PuBu"))(151)[1:101]
						cols = colourScale[(((environmentalValues[,11]-min(environmentalValues[,11],na.rm=T))/(max(environmentalValues[,11],na.rm=T)-min(environmentalValues[,11],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,11],na.rm=T),max(environmentalValues[,11],na.rm=T))))
						mtext("Relative humidity", side=3, line=-1.5, at=0, cex=0.50, col="gray30"); mtext("(summer, %)", side=3, line=-2.2, at=0, cex=0.50, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.870,0.885,0.05,0.95), adj=3,
							 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
						colourScale = c("#E5E5E5",colorRampPalette(c("gray97","navajowhite4"),bias=1)(121)[21:121])
						cols = colourScale[(((environmentalValues[,17]-min(environmentalValues[,17],na.rm=T))/(max(environmentalValues[,17],na.rm=T)-min(environmentalValues[,17],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,17],na.rm=T),max(environmentalValues[,17],na.rm=T))))
						mtext("Croplands", side=3, line=-1.5, at=0, cex=0.50, col="gray30"); mtext("(all categories)", side=3, line=-2.2, at=0, cex=0.50, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.870,0.885,0.05,0.95), adj=3,
							 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
						colourScale = c("#E5E5E5",colorRampPalette(c("gray97","burlywood3"),bias=1)(121)[21:121])
						cols = colourScale[(((environmentalValues[,18]-min(environmentalValues[,18],na.rm=T))/(max(environmentalValues[,18],na.rm=T)-min(environmentalValues[,18],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,18],na.rm=T),max(environmentalValues[,18],na.rm=T))))
						mtext("Pastures and", side=3, line=-1.5, at=0, cex=0.50, col="gray30"); mtext("rangeland", side=3, line=-2.2, at=0, cex=0.50, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.870,0.885,0.05,0.95), adj=3,
							 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"BuPu"))(151)[21:121]; values = environmentalValues[,19]; values[values>(10^-8)] = (10^-8)
						cols = colourScale[(((values-min(values,na.rm=T))/(max(values,na.rm=T)-min(values,na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(values,na.rm=T),max(values,na.rm=T))))
						mtext("Human population", side=3, line=-1.5, at=0, cex=0.50, col="gray30"); mtext("density (log10/km2)", side=3, line=-2.2, at=0, cex=0.50, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.870,0.885,0.05,0.95), adj=3,
							 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
						dev.off()
						pdf("All_env_factors_2_NEW.pdf", width=8, height=5.8); par(mfrow=c(3,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30") # 2° version
						colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]
						cols = colourScale[(((environmentalValues[,1]-min(environmentalValues[,1:4],na.rm=T))/(max(environmentalValues[,1:4],na.rm=T)-min(environmentalValues[,1:4],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,1:4],na.rm=T),max(environmentalValues[,1:4],na.rm=T))))
						mtext("Air temperature", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(winter, °C)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]
						cols = colourScale[(((environmentalValues[,2]-min(environmentalValues[,1:4],na.rm=T))/(max(environmentalValues[,1:4],na.rm=T)-min(environmentalValues[,1:4],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,1:4],na.rm=T),max(environmentalValues[,1:4],na.rm=T))))
						mtext("Air temperature", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(spring, °C)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]
						cols = colourScale[(((environmentalValues[,3]-min(environmentalValues[,1:4],na.rm=T))/(max(environmentalValues[,1:4],na.rm=T)-min(environmentalValues[,1:4],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,1:4],na.rm=T),max(environmentalValues[,1:4],na.rm=T))))
						mtext("Air temperature", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(summer, °C)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]
						cols = colourScale[(((environmentalValues[,4]-min(environmentalValues[,1:4],na.rm=T))/(max(environmentalValues[,1:4],na.rm=T)-min(environmentalValues[,1:4],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,1:4],na.rm=T),max(environmentalValues[,1:4],na.rm=T))))
						mtext("Air temperature", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(fall, °C)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = c("#E5E5E5",colorRampPalette(c("gray97","chartreuse4"),bias=1)(121)[21:121])
						cols = colourScale[(((environmentalValues[,13]-min(environmentalValues[,13],na.rm=T))/(max(environmentalValues[,13],na.rm=T)-min(environmentalValues[,13],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,13],na.rm=T),max(environmentalValues[,13],na.rm=T))))
						mtext("Primary forested", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("areas", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = c("#E5E5E5",colorRampPalette(c("gray97","darkseagreen4"),bias=1)(121)[21:121])
						cols = colourScale[(((environmentalValues[,14]-min(environmentalValues[,14],na.rm=T))/(max(environmentalValues[,14],na.rm=T)-min(environmentalValues[,14],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,14],na.rm=T),max(environmentalValues[,14],na.rm=T))))
						mtext("Primary non-forested", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("areas", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(101)
						cols = colourScale[(((environmentalValues[,5]-min(environmentalValues[,5:8],na.rm=T))/(max(environmentalValues[,5:8],na.rm=T)-min(environmentalValues[,5:8],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,5:8],na.rm=T),max(environmentalValues[,5:8],na.rm=T))))
						mtext("Precipitation", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(winter, kg/m2/day)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(101)
						cols = colourScale[(((environmentalValues[,6]-min(environmentalValues[,5:8],na.rm=T))/(max(environmentalValues[,5:8],na.rm=T)-min(environmentalValues[,5:8],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,5:8],na.rm=T),max(environmentalValues[,5:8],na.rm=T))))
						mtext("Precipitation", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(spring, kg/m2/day)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(101)
						cols = colourScale[(((environmentalValues[,7]-min(environmentalValues[,5:8],na.rm=T))/(max(environmentalValues[,5:8],na.rm=T)-min(environmentalValues[,5:8],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,5:8],na.rm=T),max(environmentalValues[,5:8],na.rm=T))))
						mtext("Precipitation", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(summer, kg/m2/day)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(101)
						cols = colourScale[(((environmentalValues[,8]-min(environmentalValues[,5:8],na.rm=T))/(max(environmentalValues[,5:8],na.rm=T)-min(environmentalValues[,5:8],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,5:8],na.rm=T),max(environmentalValues[,5:8],na.rm=T))))
						mtext("Precipitation", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(fall, kg/m2/day)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = c("#E5E5E5",colorRampPalette(c("gray97","burlywood3"),bias=1)(121)[21:121])
						cols = colourScale[(((environmentalValues[,18]-min(environmentalValues[,18],na.rm=T))/(max(environmentalValues[,18],na.rm=T)-min(environmentalValues[,18],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,18],na.rm=T),max(environmentalValues[,18],na.rm=T))))
						mtext("Pastures and", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("rangeland", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = c("#E5E5E5",colorRampPalette(c("gray97","olivedrab3"),bias=1)(121)[21:121])
						cols = colourScale[(((environmentalValues[,15]-min(environmentalValues[,15],na.rm=T))/(max(environmentalValues[,15],na.rm=T)-min(environmentalValues[,15],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,15],na.rm=T),max(environmentalValues[,15],na.rm=T))))
						mtext("Secondary forested", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("areas", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"PuBu"))(151)[1:101]
						cols = colourScale[(((environmentalValues[,9]-min(environmentalValues[,9:12],na.rm=T))/(max(environmentalValues[,9:12],na.rm=T)-min(environmentalValues[,9:12],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,9:12],na.rm=T),max(environmentalValues[,9:12],na.rm=T))))
						mtext("Relative humidity", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(winter, %)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"PuBu"))(151)[1:101]
						cols = colourScale[(((environmentalValues[,10]-min(environmentalValues[,9:12],na.rm=T))/(max(environmentalValues[,9:12],na.rm=T)-min(environmentalValues[,9:12],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,9:12],na.rm=T),max(environmentalValues[,9:12],na.rm=T))))
						mtext("Relative humidity", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(spring, %)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"PuBu"))(151)[1:101]
						cols = colourScale[(((environmentalValues[,11]-min(environmentalValues[,9:12],na.rm=T))/(max(environmentalValues[,9:12],na.rm=T)-min(environmentalValues[,9:12],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,9:12],na.rm=T),max(environmentalValues[,9:12],na.rm=T))))
						mtext("Relative humidity", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(summer, %)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"PuBu"))(151)[1:101]
						cols = colourScale[(((environmentalValues[,12]-min(environmentalValues[,9:12],na.rm=T))/(max(environmentalValues[,9:12],na.rm=T)-min(environmentalValues[,9:12],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,9:12],na.rm=T),max(environmentalValues[,9:12],na.rm=T))))
						mtext("Relative humidity", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(fall, %)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = c("#E5E5E5",colorRampPalette(c("gray97","navajowhite4"),bias=1)(121)[21:121])
						cols = colourScale[(((environmentalValues[,17]-min(environmentalValues[,17],na.rm=T))/(max(environmentalValues[,17],na.rm=T)-min(environmentalValues[,17],na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(environmentalValues[,17],na.rm=T),max(environmentalValues[,17],na.rm=T))))
						mtext("Croplands", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(all categories)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						colourScale = colorRampPalette(brewer.pal(9,"BuPu"))(151)[21:121]; values = environmentalValues[,19]; values[values>(10^-8)] = (10^-8)
						cols = colourScale[(((values-min(values,na.rm=T))/(max(values,na.rm=T)-min(values,na.rm=T)))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
						rast = raster(as.matrix(c(min(values,na.rm=T),max(values,na.rm=T))))
						mtext("Human population", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("density (log10/km2)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
							 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
						dev.off()
					}
			}
		saveRDS(nutsM_datas, "Training_dFrames.rds")
	}	else	{
		nutsM_datas = readRDS("Training_dFrames.rds")
	}
samplingPtsMinDist = function(observations, minDist=500, nberOfPoints=5)
	{
		# function written by Jean Artois (source: Dhingra, Artois, et al. 2016, eLife)
		indices = rep(NA, nberOfPoints)
		selection_list = list(1:nrow(observations)) 
  		indices[1] = sample(1:dim(observations)[1], 1)
		dists = list(spDistsN1(as.matrix(observations), as.matrix(observations[indices[1],]), longlat=T))
		for (i in 2:nberOfPoints)
			{
    				selection = which(dists[[(i-1)]] > minDist)
    				if (length(selection) == 0)
    					{
    						stop("Restarts the function with a smaller minimum distance")
					}
    				selection_list[[i]] = selection
    				test = table(unlist(selection_list))
    				indices_minDist = as.numeric(names(which(test==i)))
    				indices[i] = sample(indices_minDist, 1)   
				dists[[i]] = spDistsN1(as.matrix(observations), as.matrix(observations[indices[i],]), longlat=T)
			}
		return(indices)
	}
foldSelection = function(observations, selectedPoints)
	{
		# function written by Jean Artois (source: Dhingra, Artois, et al. 2016, eLife)
		fold_selection = sapply(1:nrow(observations), function(i) which.min(spDistsN1(as.matrix(selectedPoints), as.matrix(observations[i,]), longlat=T)))
		return(fold_selection)
	}
training_new_models = FALSE
if (training_new_models)
	{
		for (i in 1:length(models_isimip3a))
			{
				data = data.frame(as.matrix(nutsM_datas[[i]][,3:dim(nutsM_datas[[i]])[2]]))
				centroids = coordinates(nutsM); indices = c()
				for (j in 1:dim(data)[1])
					{
						if (sum(!is.na(data[j,1:dim(data)[2]])) == dim(data)[2]) indices = c(indices, j)
					}
				data = data[indices,]; centroids = centroids[indices,]
				data[which(data[,"observations"]>0),"observations"] = 1
				plottingCorrelogram = FALSE
				if (plottingCorrelogram == TRUE)
					{
						correlogram = ncf::correlog(centroids[,1], centroids[,2], data[,"observations"], na.rm=T, increment=10, resamp=0, latlon=T)
						dev.new(width=4.5, height=3); par(mar=c(2.5,2.2,1.0,1.5))
						plot(correlogram$mean.of.class[-1], correlogram$correlation[-1], ann=F, axes=F, lwd=0.2, cex=0.5, col=NA, ylim=c(-0.4,1.0), xlim=c(0,3500))
						abline(h=0, lwd=0.5, col="red", lty=2)
						lines(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, col="gray30")
						points(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, cex=0.25, col="gray90", pch=16)
						points(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, cex=0.25, col="gray30", pch=1)
						axis(side=1, pos=-0.4, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,-0.05,0), at=seq(0,9000,1000))
						axis(side=2, pos=0, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.18,0), at=seq(-0.4,1,0.2))
						title(xlab="distance (km2)", cex.lab=0.7, mgp=c(0.3,0,0), col.lab="gray30")
						title(ylab="correlation", cex.lab=0.7, mgp=c(0.4,0,0), col.lab="gray30")
					}
				theRanges = c(500,500)*1000 # distance in meters
				nberOfReplicates = 10 # one replicate = one folds partition
				gbm.x = colnames(data)[2:dim(data)[2]]
				gbm.y = "observations"
				offset = NULL
				tree.complexity = 5 # "tc" = number of nodes in the trees
				learning.rate = 0.005 # "lr" = contribution of each tree to the growing model
				bag.fraction = 0.80 # proportion of data used to train a given tree
				site.weights = rep(1, dim(data)[1])
				var.monotone = rep(0, length(gbm.x))
				n.folds = 5
				prev.stratify = TRUE
				family = "bernoulli"
				n.trees = 100 # initial number of trees
				step.size = 10 # interval at which the predictive deviance is computed and logged
							   # (at each interval, the folds are successively used as test data set
							   # and the remaining folds as training data sets to compute the deviance)
				max.trees = 10000 # maximum number of trees that will be considered
				tolerance.method = "auto"
				tolerance = 0.001
				plot.main = TRUE
				plot.folds = FALSE
				verbose = TRUE
				silent = FALSE
				keep.fold.models = FALSE
				keep.fold.vector = FALSE
				keep.fold.fit = FALSE
				showingFoldsPlot = FALSE		
				brt_model_ccvs = list() # classic cross-validations (CCVs)
				brt_model_scvs = list() # spatial cross-validations (SCVs)
				AUCs = matrix(nrow=nberOfReplicates, ncol=2)
				colnames(AUCs) = c("CCV_AUC","SCV_AUC")
				for (j in 1:nberOfReplicates)
					{
						# BRT with classic (standard) cross-validation (CCV):
						n.trees = 100; fold.vector = NULL
						brt_model_ccvs[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
							var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
							verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
						dev.copy2pdf(file=paste0("All_the_BRT_models/",models_isimip3a[i],"_CCV_replicate_",j,".pdf")); dev.off()
						AUCs[j,"CCV_AUC"] = brt_model_ccvs[[j]]$cv.statistics$discrimination.mean # Mean test AUC (from the AUCs computed on each fold tested as test data in the CCV)
						object = brt_model_ccvs[[j]]; n.trees = brt_model_ccvs[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						projection = predict.gbm(object, data[,2:dim(data)[2]], n.trees, type, single.tree)
						# BRT with spatial (geographic) cross-validation (SCV) based on the blocks generation of Valavi et al. (2019, MEE):
						spdf = SpatialPointsDataFrame(centroids, data[,1:dim(data)[2]], proj4string=crs(nutsM))
						myblocks = spatialBlock(spdf, species="observations", rasterLayer=NULL, k=n.folds, theRange=theRanges[1], selection="random")
						fold.vector = myblocks$foldID; n.trees = 100
						brt_model_scvs[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
							var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
							verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
						dev.copy2pdf(file=paste0("All_the_BRT_models/",models_isimip3a[i],"_SCV_replicate_",j,".pdf")); dev.off()
						AUCs[j,"SCV_AUC"] = brt_model_scvs[[j]]$cv.statistics$discrimination.mean # Mean test AUC (from the AUCs computed on each fold tested as test data in the SCV)
						object = brt_model_scvs[[j]]; n.trees = brt_model_scvs[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						projection = predict.gbm(object, data[,2:dim(data)[2]], n.trees, type, single.tree)
					}
				saveRDS(brt_model_ccvs, paste0("All_the_BRT_models/",models_isimip3a[i],"_models_CCV.rds"))
				saveRDS(brt_model_scvs, paste0("All_the_BRT_models/",models_isimip3a[i],"_models_SCV.rds"))
				write.csv(AUCs, paste0("All_the_BRT_models/",models_isimip3a[i],"_CCV_SCV_AUCs.csv"), row.names=F, quote=F)
			}
	}
nutsM_projections = list(); maxV = 0
for (i in 1:length(models_isimip3a))
	{
		brt_model_scvs = readRDS(paste0("All_the_BRT_models/",models_isimip3a[i],"_models_SCV.rds"))
		df = as.data.frame(nutsM_datas[[i]]); colnames(df) = gsub("-","\\.",colnames(df))
		nutsM_projection1 = matrix(nrow=dim(df)[1], ncol=length(brt_model_scvs))
		for (j in 1:length(brt_model_scvs))
			{
				n.trees = brt_model_scvs[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
				nutsM_projection1[,j] = predict.gbm(brt_model_scvs[[j]], df, n.trees, type, single.tree)
			}
		nutsM_projection2 = matrix(nrow=dim(df)[1], ncol=1)
		for (j in 1:dim(nutsM_projection1)[1])
			{
				nutsM_projection2[j,1] = mean(nutsM_projection1[j,])
			}
		nutsM_projections[[i]] = nutsM_projection2
	}
for (i in 1:length(nutsM_projections))
	{
		if (maxV < max(nutsM_projections[[i]])) maxV = max(nutsM_projections[[i]])
	}
if (savingPlots)
	{
		pdf("Training_dFrames_NEW.pdf", width=8, height=5.8); par(mfrow=c(3,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")		
		colourScale = rev(colorRampPalette(brewer.pal(11,"RdBu"))(121)[11:111])
		for (i in 1:length(nutsM_projections))
			{
				cols = colourScale[(((nutsM_projections[[i]]-0)/(1-0))*100)+1]
				plot(nutsM, col=cols, border=NA, lwd=0.1); rast = raster(as.matrix(c(0,1)))
				mtext(models_isimip3a_names[i], side=3, line=-1.5, at=0, cex=0.50, col="gray30"); mtext("", side=3, line=-2.2, at=0, cex=0.50, col="gray30")
				if (i == 6)
					{
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.120,0.135,0.70,0.95), adj=3,
							 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
					}
			}
		dev.off()
	}

# 3. Compouting the relative influences and response curves

	# 3.1. Comparison of the relative influence of each environmental factor

computing_relative_influences = FALSE
if (computing_relative_influences)
	{
		relativeInfluences = matrix(0, nrow=length(envVariableNames), ncol=length(models_isimip3a))
		row.names(relativeInfluences) = envVariableNames; colnames(relativeInfluences) = models_isimip3a
		for (i in 1:length(models_isimip3a))
			{
				brt_model_scvs = readRDS(paste0("All_the_BRT_models/",models_isimip3a[i],"_models_SCV.rds"))
				for (j in 1:length(brt_model_scvs))
					{
						for (k in 1:length(envVariableNames))
							{
								relativeInfluences[k,i] = relativeInfluences[k,i] + summary(brt_model_scvs[[j]])[gsub("-","\\.",envVariableNames)[k],"rel.inf"]
							}
					}
				relativeInfluences[,i] = relativeInfluences[,i]/length(brt_model_scvs)
			}
		write.table(round(relativeInfluences,1), "Relative_influences.csv", quote=F, sep=",")
	}

	# 3.2. Comparison of the response curves for each environmental factor

envVariableValues_list = list()
for (i in 1:length(models_isimip3a))
	{
		data = data.frame(as.matrix(nutsM_datas[[i]][,3:dim(nutsM_datas[[i]])[2]]))
		centroids = coordinates(nutsM); indices = c()
		for (j in 1:dim(data)[1])
			{
				if (sum(!is.na(data[j,1:dim(data)[2]])) == dim(data)[2]) indices = c(indices, j)
			}
		data = data[indices,]; centroids = centroids[indices,]
		data[which(data[,"observations"]>0),"observations"] = 1
		data = data[which(data[,"observations"]==1),]
		envVariableValues = matrix(nrow=3, ncol=length(envVariableNames))
		row.names(envVariableValues) = c("median","minV","maxV")
		colnames(envVariableValues) = envVariableNames
		for (j in 1:length(envVariableNames))
			{
				minV = min(data[,gsub("-","\\.",envVariableNames)[j]], na.rm=T)
				maxV = max(data[,gsub("-","\\.",envVariableNames)[j]], na.rm=T)
				medianV = median(data[,gsub("-","\\.",envVariableNames)[j]], na.rm=T)
				envVariableValues[,j] = cbind(medianV, minV, maxV)
			}
		envVariableValues_list[[i]] = envVariableValues
	}
if (savingPlots)
	{
		pdf("Response_curves_NEW.pdf", width=8, height=5); par(mfrow=c(4,5), oma=c(1.3,1.5,1,0.5), mar=c(2.5,1,0.5,1), lwd=0.2, col="gray30")
		for (i in 1:length(envVariableNames))
			{
				projections_list = list(); dfs = list()
				for (j in 1:1)
					{
						valuesInterval = 0.1; valuesInterval = (envVariableValues_list[[j]]["maxV",i]-envVariableValues_list[[j]]["minV",i])/100
						df = data.frame(matrix(nrow=length(seq(envVariableValues_list[[j]]["minV",i],envVariableValues_list[[j]]["maxV",i],valuesInterval)),ncol=length(envVariableNames)))
						colnames(df) = gsub("-","\\.",envVariableNames)
						for (k in 1:length(envVariableNames))
							{
								valuesInterval = 0.1; valuesInterval = (envVariableValues_list[[j]]["maxV",k]-envVariableValues_list[[j]]["minV",k])/100
								if (i == k) df[,gsub("-","\\.",envVariableNames)[k]] = seq(envVariableValues_list[[j]]["minV",k],envVariableValues_list[[j]]["maxV",k],valuesInterval)
								if (i != k) df[,gsub("-","\\.",envVariableNames)[k]] = rep(envVariableValues_list[[j]]["median",k],dim(df)[1])
							}
						dfs[[j]] = df; projections = list()
						brt_model_scvs = readRDS(paste0("All_the_BRT_models/",models_isimip3a[j],"_models_SCV.rds"))
						for (k in 1:length(brt_model_scvs))
							{
								n.trees = brt_model_scvs[[k]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
								projection = predict.gbm(brt_model_scvs[[k]], newdata=df, n.trees, type, single.tree)
								if ((j == 1)&(k == 1))
									{
										minX = min(df[,gsub("-","\\.",envVariableNames)[i]]); maxX = max(df[,gsub("-","\\.",envVariableNames)[i]])
										minY = min(projection); maxY = max(projection)
									}	else	{
										if (minX > min(df[,gsub("-","\\.",envVariableNames)[i]])) minX = min(df[,gsub("-","\\.",envVariableNames)[i]])
										if (maxX < max(df[,gsub("-","\\.",envVariableNames)[i]])) maxX = max(df[,gsub("-","\\.",envVariableNames)[i]])
										if (minY > min(projection)) minY = min(projection)
										if (maxY < max(projection)) maxY = max(projection)
									}
								projections[[k]] = projection
							}
						projections_list[[j]] = projections
					}
				cols = c("red","chartreuse3")
				for (k in 1:length(brt_model_scvs))
					{
						for (j in 1:1)
							{
								if ((j == 1)&(k == 1))
									{
										plot(dfs[[j]][,gsub("-","\\.",envVariableNames)[i]],projections_list[[j]][[k]],col=cols[j],ann=F,axes=F,lwd=0.2,type="l",xlim=c(minX,maxX),ylim=c(minY,maxY))
									}	else	{
										lines(dfs[[j]][,gsub("-","\\.",envVariableNames)[i]],projections_list[[j]][[k]],col=cols[j],lwd=0.2)
									}
							}
					}
				axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.040, col.axis="gray30", mgp=c(0,0.09,0))
				axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.040, col.axis="gray30", mgp=c(0,0.25,0))
				title(ylab="predicted values", cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
				title(xlab=envVariableNames[i], cex.lab=0.9, mgp=c(1.1,0,0), col.lab="gray30")
				box(lwd=0.2, col="gray30")
			}
		dev.off()
	}

# 4. Performing the past and present ENM projections

scenarios = c("counterclim","obsclim"); savingFiles = FALSE
year_intervals = c("1901_1919","1920_1939","1940_1959","1960_1979","1980_1999","2000_2019")
if (!file.exists("ISIMIP3a_dFrames.rds"))
	{
		nutsM_datas1 = list()
		for (i in 1:length(models_isimip3a))
			{
				nutsM_datas2 = list()
				for (h in 1:length(scenarios))
					{
						nutsM_datas3 = list()
						for (g in 1:length(year_intervals))
							{
								environmentalValues = matrix(nrow=dim(nutsM_data)[1], ncol=length(envVariableNames)); colnames(environmentalValues) = envVariableNames
								temperature = brick(paste0("Environmental_rasters/ISIMIP3a/tas_day_",scenarios[h],"_historical_",models_isimip3a[i],"_",year_intervals[g],"_ymonmean.nc"))
								temperature_winter = mean(temperature[[12]],temperature[[1]],temperature[[2]])-273.15 # conversion to Celcius degrees
								temperature_spring = mean(temperature[[3]],temperature[[4]],temperature[[5]])-273.15 # conversion to Celcius degrees
								temperature_summer = mean(temperature[[6]],temperature[[7]],temperature[[8]])-273.15 # conversion to Celcius degrees
								temperature_inFall = mean(temperature[[9]],temperature[[10]],temperature[[11]])-273.15 # conversion to Celcius degrees
								precipitation = brick(paste0("Environmental_rasters/ISIMIP3a/pr_day_",scenarios[h],"_historical_",models_isimip3a[i],"_",year_intervals[g],"_ymonmean.nc"))
								precipitation_winter = mean(precipitation[[12]],precipitation[[1]],precipitation[[2]])*60*60*24 # conversion to kg/m2/day
								precipitation_spring = mean(precipitation[[3]],precipitation[[4]],precipitation[[5]])*60*60*24 # conversion to kg/m2/day
								precipitation_summer = mean(precipitation[[6]],precipitation[[7]],precipitation[[8]])*60*60*24 # conversion to kg/m2/day
								precipitation_inFall = mean(precipitation[[9]],precipitation[[10]],precipitation[[11]])*60*60*24 # conversion to kg/m2/day
								relative_humidity = brick(paste0("Environmental_rasters/ISIMIP3a/hurs_day_",scenarios[h],"_historical_",models_isimip3a[i],"_",year_intervals[g],"_ymonmean.nc"))
								relative_humidity_winter = mean(relative_humidity[[12]],relative_humidity[[1]],relative_humidity[[2]])
								relative_humidity_spring = mean(relative_humidity[[3]],relative_humidity[[4]],relative_humidity[[5]])
								relative_humidity_summer = mean(relative_humidity[[6]],relative_humidity[[7]],relative_humidity[[8]])
								relative_humidity_inFall = mean(relative_humidity[[9]],relative_humidity[[10]],relative_humidity[[11]])
								if (g < length(year_intervals)) land_cover = nc_open(paste0("Environmental_rasters/ISIMIP3a/landcover_annual_",year_intervals[g],"_timmean.nc4"))
								if (g == length(year_intervals)) land_cover = nc_open(paste0("Environmental_rasters/ISIMIP3a/landcover_annual_2000_2015_timmean.nc4"))
								population = brick(paste0("Environmental_rasters/ISIMIP3a/population_histsoc_0p5deg_annual_",year_intervals[g],"_timmean.nc4"), varname="total-population")
								landCoverVariableIDs = names(land_cover$var); land_covers1 = list(); land_covers2 = list(); land_covers3 = list()
								landCoverVariableNames = as.character(read.csv("Environmental_rasters/Luse.csv")[1:12,2])
								for (j in 2:13)
									{
										if (g < length(year_intervals)) fileName = paste0("Environmental_rasters/ISIMIP3a/landcover_annual_",year_intervals[g],"_timmean.nc4")
										if (g == length(year_intervals)) fileName = paste0("Environmental_rasters/ISIMIP3a/landcover_annual_2000_2015_timmean.nc4")
										land_covers1[[j-1]] = brick(fileName, varname=landCoverVariableIDs[j])
									}
								variable_codes = c("croplands","pastures","urbanAreas","primaryForest","primaryNonF","secondaryForest","secondaryNonF")
								variable_names = c("crops","pasture","urban land","forested primary land","non-forested primary land",
										  		   "potentially forested secondary land","potentially non-forested secondary land")
								for (j in 1:length(variable_names))
									{
										names = gsub("\\."," ",landCoverVariableNames); indices = which(landCoverVariableNames==variable_names[j])
										if (length(indices) == 0) indices = which(grepl(variable_names[j],names))
										if (variable_names[j] == "pasture") indices = c(indices, which(grepl("rangeland",names)))
										land_cover = land_covers1[[indices[1]]]; names(land_cover) = variable_codes[j]; # print(indices)
										if (length(indices) > 1)
											{
												for (k in 2:length(indices)) land_cover[] = land_cover[]+land_covers1[[indices[k]]][]
											}
										land_covers2[[j]] = land_cover[[1]]; land_covers3[[j]] = raster::aggregate(land_cover[[1]],2)
									}
								envVariables = list()			
								envVariables[[1]] = temperature_winter; envVariables[[2]] = temperature_spring
								envVariables[[3]] = temperature_summer; envVariables[[4]] = temperature_inFall
								envVariables[[5]] = precipitation_winter; envVariables[[6]] = precipitation_spring
								envVariables[[7]] = precipitation_summer; envVariables[[8]] = precipitation_inFall
								envVariables[[9]] = relative_humidity_winter; envVariables[[10]] = relative_humidity_spring
								envVariables[[11]] = relative_humidity_summer; envVariables[[12]] = relative_humidity_inFall
								envVariables[[13]] = land_covers2[[4]] # primary forest areas
								envVariables[[14]] = land_covers2[[5]] # primary non-forest areas
								envVariables[[15]] = land_covers2[[6]] # secondary forest areas
								envVariables[[16]] = land_covers2[[7]] # secondary non-forest areas
								envVariables[[17]] = land_covers2[[1]] # croplands (all catergories)
								envVariables[[18]] = land_covers2[[2]] # managed pasture + rangeland
								envVariables[[19]] = population # human population (not log-transformed)
								for (j in 1:length(envVariables)) envVariables[[j]] = crop(envVariables[[j]], nutsM, snap="out")
								for (j in 1:length(envVariables)) envVariables[[j]] = mask(envVariables[[j]], nutsM)
								areas_nuts_M_km = area(nutsM, unit="km")
								for (j in 1:length(nutsM))
									{
										maxArea = 0; polIndex = 0
										for (k in 1:length(nutsM@polygons[[j]]@Polygons))
											{
												if (maxArea < nutsM@polygons[[j]]@Polygons[[k]]@area)
													{
														maxArea = nutsM@polygons[[j]]@Polygons[[k]]@area; polIndex = k
													}
											}
										pol1 = nutsM@polygons[[j]]@Polygons[[polIndex]]; p = Polygon(pol1@coords)
										ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
										pol2 = sf::st_as_sfc(sps); st_crs(pol2) = crs(nutsM)
										for (k in 1:18)
											{
												if ((j == 1)&(k == 1)) crs(envVariables[[k]]) = crs(pol2) 
												environmentalValues[j,k] = exactextractr::exact_extract(envVariables[[k]], pol2, fun="mean")
											}
										if ((j == 1)&(k == 1)) crs(envVariables[[19]]) = crs(pol2) 
										environmentalValues[j,19] = (log10(exactextractr::exact_extract(envVariables[[19]],pol2,fun="sum")+1))/areas_nuts_M_km[j]
									}
								nutsM_datas3[[g]] = environmentalValues
							}
						nutsM_datas2[[h]] = nutsM_datas3
					}
				nutsM_datas1[[i]] = nutsM_datas2
			}
		saveRDS(nutsM_datas1, "ISIMIP3a_dFrames.rds")
	}	else	{
		nutsM_datas1 = readRDS("ISIMIP3a_dFrames.rds")
	}
nutsM_projections1 = list()
for (i in 1:length(models_isimip3a))
	{
		brt_model_scvs = readRDS(paste0("All_the_BRT_models/",models_isimip3a[i],"_models_SCV.rds"))
		nutsM_projections2 = list()
		for (h in 1:length(scenarios))
			{
				nutsM_projections3 = list()
				for (g in 1:length(year_intervals))
					{
						df = as.data.frame(nutsM_datas1[[i]][[h]][[g]]); colnames(df) = gsub("-","\\.",colnames(df))
						nutsM_projection4 = matrix(nrow=dim(df)[1], ncol=length(brt_model_scvs))
						for (j in 1:length(brt_model_scvs))
							{
								n.trees = brt_model_scvs[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
								nutsM_projection4[,j] = predict.gbm(brt_model_scvs[[j]], df, n.trees, type, single.tree)
							}
						nutsM_projections3[[g]] = nutsM_projection4
					}
				nutsM_projections2[[h]] = nutsM_projections3
			}
		nutsM_projections1[[i]] = nutsM_projections2
	}
if (savingPlots)
	{
		for (i in 1:length(nutsM_projections1))
			{
				pdf(paste0("ISIMIP3a_projections/ISIMIP3a_",models_isimip3a_names[i],".pdf"), width=8, height=5.8); par(mfrow=c(3,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")		
				colourScale = rev(colorRampPalette(brewer.pal(11,"RdBu"))(121)[11:111])
				for (h in 1:length(scenarios))
					{
						for (g in 1:length(year_intervals))
							{
								nutsM_projection = matrix(nrow=dim(nutsM_projections1[[i]][[h]][[g]])[1], ncol=1)
								for (j in 1:dim(nutsM_projection)[1])
									{
										nutsM_projection[j,1] = mean(nutsM_projections1[[i]][[h]][[g]][j,])
									}
								cols = colourScale[(((nutsM_projection[,1]-0)/(1-0))*100)+1]
								plot(contour, lwd=0.4, border="gray30", col=NA)
								plot(nutsM, col=cols, border=NA, lwd=0.1, add=T); rast = raster(as.matrix(c(0,1)))
								if ((h == 1)&(g == 1)) mtext(expression(bold(A)), side=3, line=-1.4, at=-8.5, cex=0.70, col="gray30")
								if ((h == 2)&(g == 1)) mtext(expression(bold(B)), side=3, line=-1.4, at=-8.5, cex=0.70, col="gray30")
								if (h == 1)
									{
										mtext(gsub("_","-",year_intervals[g]), side=3, line=-2.5, at=1, cex=0.50, col="gray30")
									}
								if ((h == 2)&(g == length(year_intervals)))
									{
										plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
											 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.2, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
									}								
							}
					}
				dev.off()
			}
		min_max_values_1 = matrix(nrow=12, ncol=2)
		for (i in 1:12)
			{
				minV = 9999; maxV = -9999
				for (j in 1:length(nutsM_datas1))
					{
						for (g in 1:length(year_intervals))
							{
								if (min(nutsM_datas1[[j]][[2]][[g]][,i],na.rm=T) < minV) minV = min(nutsM_datas1[[j]][[2]][[g]][,i],na.rm=T)
								if (max(nutsM_datas1[[j]][[2]][[g]][,i],na.rm=T) > maxV) maxV = max(nutsM_datas1[[j]][[2]][[g]][,i],na.rm=T)
							}
					}
				min_max_values_1[i,1] = minV; min_max_values_1[i,2] = maxV
			}
		for (i in 1:12)
			{
				pdf(paste0("ISIMIP3a_envVariables/Abs_values/ISIMIP3a_",colnames(nutsM_datas1[[1]][[2]][[1]])[i],"_values.pdf"), width=8, height=7.73)
				par(mfrow=c(4,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")	
				if (i %in% c(1:4)) colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]
				if (i %in% c(5:8)) colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(101)
				if (i %in% c(9:12)) colourScale = colorRampPalette(brewer.pal(9,"PuBu"))(151)[1:101]
				for (j in 1:length(nutsM_datas1))
					{
						for (g in 1:length(year_intervals))
							{
								cols = colourScale[(((nutsM_datas1[[j]][[2]][[g]][,i]-min_max_values_1[i,1])/(min_max_values_1[i,2]-min_max_values_1[i,1]))*100)+1]
								plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=cols, lwd=0.1, add=T)
								# if ((j == 1)&(g == 1)) mtext(expression(bold(A)), side=3, line=-1.4, at=-8.5, cex=0.70, col="gray30")
								# if ((j == 2)&(g == 1)) mtext(expression(bold(B)), side=3, line=-1.4, at=-8.5, cex=0.70, col="gray30")
								# if ((j == 3)&(g == 1)) mtext(expression(bold(C)), side=3, line=-1.4, at=-8.5, cex=0.70, col="gray30")
								# if ((j == 4)&(g == 1)) mtext(expression(bold(D)), side=3, line=-1.4, at=-8.5, cex=0.70, col="gray30")
								if (g == 1) mtext(models_isimip3a_names[j], side=3, line=-1.7, at=1, cex=0.50, col="gray30")
								if (j == 1) mtext(gsub("_","-",year_intervals[g]), side=3, line=-2.5, at=1, cex=0.50, col="gray30")
								if ((j == 4)&(g == length(year_intervals)))
									{
										rast = raster(as.matrix(cbind(min_max_values_1[i,1],min_max_values_1[i,2])))
										plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
											 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.2, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
									}								
							}						
					}
				dev.off()
			}
		min_max_values_2 = matrix(nrow=12, ncol=2)
		for (i in 1:12)
			{
				minV = 9999; maxV = -9999
				for (j in 1:length(nutsM_datas1))
					{
						for (g in 1:length(year_intervals))
							{
								if (min(nutsM_datas1[[j]][[2]][[g]][,i]-nutsM_datas1[[j]][[2]][[6]][,i],na.rm=T) < minV) minV = min(nutsM_datas1[[j]][[2]][[g]][,i]-nutsM_datas1[[j]][[2]][[6]][,i],na.rm=T)
								if (max(nutsM_datas1[[j]][[2]][[g]][,i]-nutsM_datas1[[j]][[2]][[6]][,i],na.rm=T) > maxV) maxV = max(nutsM_datas1[[j]][[2]][[g]][,i]-nutsM_datas1[[j]][[2]][[6]][,i],na.rm=T)
							}
					}
				if (abs(minV) < abs(maxV)) minV = -maxV
				if (abs(maxV) < abs(minV)) maxV = -minV
				min_max_values_2[i,1] = minV; min_max_values_2[i,2] = maxV
			}
		for (i in 1:12)
			{
				pdf(paste0("ISIMIP3a_envVariables/Differences/ISIMIP3a_",colnames(nutsM_datas1[[1]][[2]][[1]])[i],"_difference.pdf"), width=8, height=7.73)
				par(mfrow=c(4,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")	
				for (j in 1:length(nutsM_datas1))
					{
						for (g in 1:(length(year_intervals)-1))
							{
								colourScale = rev(colorRampPalette(brewer.pal(11,"PuOr"))(121)[11:111])
								differences = nutsM_datas1[[j]][[2]][[g]][,i]-nutsM_datas1[[j]][[2]][[6]][,i]
								cols = colourScale[(((differences-min_max_values_2[i,1])/(min_max_values_2[i,2]-min_max_values_2[i,1]))*100)+1]
								plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, add=T)
								if (j == 1) mtext(gsub("_","-",year_intervals[g]), side=3, line=-2.5, at=1, cex=0.50, col="gray30")
								if (g == 1) mtext(models_isimip3a_names[j], side=3, line=-1.7, at=1, cex=0.50, col="gray30")
								if ((j == length(nutsM_datas1))&(g == (length(year_intervals)-1)))
									{
										rast = raster(as.matrix(cbind(min_max_values_2[i,1],min_max_values_2[i,2])))
										plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
											 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.2, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
									}					
							}						
						if (i %in% c(1:4)) colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]
						if (i %in% c(5:8)) colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(101)
						if (i %in% c(9:12)) colourScale = colorRampPalette(brewer.pal(9,"PuBu"))(151)[1:101]
						cols = colourScale[(((nutsM_datas1[[j]][[2]][[1]][,i]-min_max_values_1[i,1])/(min_max_values_1[i,2]-min_max_values_1[i,1]))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=cols, lwd=0.1, add=T)			
						if (j == 1) mtext(gsub("_","-",year_intervals[length(year_intervals)]), side=3, line=-2.5, at=1, cex=0.50, col="gray30")
						if (j == length(nutsM_datas1))
							{
								rast = raster(as.matrix(cbind(min_max_values_1[i,1],min_max_values_1[i,2])))
								plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
									 axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.2, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
							}								
					}
				dev.off()
			}
	}
population_counts = matrix(nrow=dim(nutsM@data)[1], ncol=length(year_intervals))
colnames(population_counts) = year_intervals
for (i in 1:length(year_intervals))
	{
		population = brick(paste0("Environmental_rasters/ISIMIP3a/population_histsoc_0p5deg_annual_",year_intervals[i],"_timmean.nc4"), varname="total-population")
		for (j in 1:length(nutsM))
			{
				maxArea = 0; polIndex = 0
				for (k in 1:length(nutsM@polygons[[j]]@Polygons))
					{
						if (maxArea < nutsM@polygons[[j]]@Polygons[[k]]@area)
							{
								maxArea = nutsM@polygons[[j]]@Polygons[[k]]@area; polIndex = k
							}
					}
				pol1 = nutsM@polygons[[j]]@Polygons[[polIndex]]; p = Polygon(pol1@coords)
				ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
				pol2 = sf::st_as_sfc(sps); st_crs(pol2) = crs(nutsM)
				if (j == 1) crs(population) = crs(pol2) 
				population_counts[j,i] = exactextractr::exact_extract(population, pol2, fun="sum")
			}
	}
if (savingPlots)
	{
		colours1 = c(rgb(255,127,14,255,maxColorValue=255), rgb(31,119,180,255,maxColorValue=255))
		colours2 = c(rgb(255,127,14,125,maxColorValue=255), rgb(31,119,180,125,maxColorValue=255))
		cols1 = c(rep(colours1[1], length(year_intervals)), rep(colours1[2], length(year_intervals)))
		cols2 = c(rep(colours2[1], length(year_intervals)), rep(colours2[2], length(year_intervals)))
		thresholds = c(0.1, 0.5)
		for (i in 1:length(thresholds))
			{
				tab = matrix(nrow=10, ncol=12); colNames = c()
				for (j in 1:length(year_intervals))
					{
						colNames = c(colNames, paste0(year_intervals[j],"_counterclim"))
						for (k in 1:dim(nutsM_projections1[[1]][[1]][[j]])[2])
							{
								vS = population_counts[,j]
								vS[which(nutsM_projections1[[1]][[1]][[j]][,k]<thresholds[i])] = 0
								tab[k,j] = sum(vS)
							}
					}
				for (j in 1:length(year_intervals))
					{
						colNames = c(colNames, paste0(year_intervals[j],"_obsclim"))
						for (k in 1:dim(nutsM_projections1[[1]][[2]][[j]])[2])
							{
								vS = population_counts[,j]
								vS[which(nutsM_projections1[[1]][[2]][[j]][,k]<thresholds[i])] = 0
								tab[k,6+j] = sum(vS)
							}
					}
				colnames(tab) = colNames
				if (savingFiles) write.csv(tab, paste0("ISIMIP3a_model1_",gsub("\\.",",",as.character(thresholds[i])),".csv"), row.names=F, quote=F)
				pdf(paste0("ISIMIP3a_violin_",gsub("\\.",",",as.character(thresholds[i])),"_NEW.pdf"), width=4, height=2.2)
				par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1.5,1.5,0,0), mgp=c(1.2,0.75,0), lwd=0.2, col="gray30")
				if (i == 1) plot(0:1, 0:1, type="n", xlim=c(190.5,201.5), ylim=c(0.55*(10^8),1.9*(10^8)), axes=F, ann=F)
				if (i == 2) plot(0:1, 0:1, type="n", xlim=c(190.5,201.5), ylim=c(0.3*(10^7),5.7*(10^7)), axes=F, ann=F)
				vioplot(tab, ann=F, axes=F, at=rep(c(1910,1930,1950,1970,1990,2010)/10,2), border=cols1, col=cols2, use.cols=T, horizontal=F, lineCol=NA, rectCol=NA, colMed=NA, add=T)
				axis(1, at=c(1880,1910,1930,1950,1970,1990,2010,2030)/10, labels=c(1880,1910,1930,1950,1970,1990,2010,2030), mgp=c(0,0.11,0), lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.03, col="gray30", col.axis="gray30", col.lab="gray30")
				if (i == 1) axis(2, at=seq(0.3*(10^8),2.1*(10^8),3*(10^7)), labels=c(30,60,90,120,150,180,210), mgp=c(1,0.30,0), lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.03, col="gray30", col.axis="gray30", col.lab="gray30")
				if (i == 2) axis(2, at=seq(0,6*(10^7),10^7), labels=c(0,10,20,30,40,50,60), mgp=c(1,0.30,0), lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.03, col="gray30", col.axis="gray30", col.lab="gray30")
				dev.off()
			}
		pdf(paste0("ISIMIP3a_violin_SI_NEW.pdf"), width=8, height=3.6)
		par(mfrow=c(2,3), oma=c(0,1,0,0), mar=c(2,2.5,0,0), mgp=c(1.2,0.75,0), lwd=0.2, col="gray30")
		for (i in 1:length(thresholds))
			{
				for (h in 2:length(models_isimip3a))
					{
						tab = matrix(nrow=10, ncol=12); colNames = c()
						for (j in 1:length(year_intervals))
							{
								colNames = c(colNames, paste0(year_intervals[j],"_counterclim"))
								for (k in 1:dim(nutsM_projections1[[h]][[1]][[j]])[2])
									{
										vS = population_counts[,j]
										vS[which(nutsM_projections1[[h]][[1]][[j]][,k]<thresholds[i])] = 0
										tab[k,j] = sum(vS)
									}
								
							}
						for (j in 1:length(year_intervals))
							{
								colNames = c(colNames, paste0(year_intervals[j],"_obsclim"))
								for (k in 1:dim(nutsM_projections1[[h]][[2]][[j]])[2])
									{
										vS = population_counts[,j]
										vS[which(nutsM_projections1[[h]][[2]][[j]][,k]<thresholds[i])] = 0
										tab[k,6+j] = sum(vS)
									}
							}
						colnames(tab) = colNames
						if (savingFiles) write.csv(tab, paste0("ISIMIP3a_model",h,"_",gsub("\\.",",",as.character(thresholds[i])),".csv"), row.names=F, quote=F)
						if (i == 1) plot(0:1, 0:1, type="n", xlim=c(190.5,201.5), ylim=c(0.63*(10^8),1.83*(10^8)), axes=F, ann=F)
						if (i == 2) plot(0:1, 0:1, type="n", xlim=c(190.5,201.5), ylim=c(0.8*(10^7),6.2*(10^7)), axes=F, ann=F)
						vioplot(tab, ann=F, axes=F, at=rep(c(1910,1930,1950,1970,1990,2010)/10,2), border=cols1, col=cols2, use.cols=T, horizontal=F, lineCol=NA, rectCol=NA, colMed=NA, add=T)
						axis(1, at=c(1880,1910,1930,1950,1970,1990,2010,2030)/10, labels=c(1880,1910,1930,1950,1970,1990,2010,2030), mgp=c(0,0.06,0), lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.025, col="gray30", col.axis="gray30", col.lab="gray30")
						if (i == 1) axis(2, at=seq(0.5*(10^8),1.9*(10^8),2*(10^7)), labels=c(50,70,90,110,130,150,170,190), mgp=c(1,0.20,0), lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.02, col="gray30", col.axis="gray30", col.lab="gray30")
						if (i == 2) axis(2, at=seq(0,7*(10^7),10^7), labels=c(0,10,20,30,40,50,60,70), mgp=c(1,0.20,0), lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.02, col="gray30", col.axis="gray30", col.lab="gray30")
						if (h == 2) title(ylab="Population at risk (in million people)", cex.lab=0.80, mgp=c(1.5,0,0), col.lab="gray30")
						mtext(paste0(models_isimip3a_names[h]," - ecological suitability >",as.character(thresholds[i])), side=3, at=196, line=-0.8, cex=0.5, col="gray30")
					}
			}
		dev.off()
	}

# 5. Performing all the future ENM projections

models_isimip3b = c("canesm5","cnrm-cm6-1","cnrm-esm2-1","ec-earth3","gfdl-esm4","ipsl-cm6a-lr","miroc6","mpi-esm1-2-hr","mri-esm2-0","ukesm1-0-ll")
models_isimip3b_names = c("CanESM5","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL")
scenarios = c("ssp126","ssp370","ssp585"); year_intervals = c("2020_2039","2040_2059","2060_2079","2080_2099"); scenario_names = c("RCP 2.6 - SSP1","RCP 6.0 - SSP3","RCP 8.5 - SSP7")
if (!file.exists("ISIMIP3b_dFrames.rds"))
	{
		nutsM_datas1 = list()
		for (i in 1:length(models_isimip3b))
			{
				nutsM_datas2 = list()
				for (h in 1:length(scenarios))
					{
						nutsM_datas3 = list()
						for (g in 1:length(year_intervals))
							{
								environmentalValues = matrix(nrow=dim(nutsM_data)[1], ncol=length(envVariableNames)); colnames(environmentalValues) = envVariableNames
								temperature = brick(paste0("Environmental_rasters/ISIMIP3b/tas_day_bias-adjusted_",scenarios[h],"_",models_isimip3b[i],"_",year_intervals[g],"_ymonmean.nc"))
								temperature_winter = mean(temperature[[12]],temperature[[1]],temperature[[2]])-273.15 # conversion to Celcius degrees
								temperature_spring = mean(temperature[[3]],temperature[[4]],temperature[[5]])-273.15 # conversion to Celcius degrees
								temperature_summer = mean(temperature[[6]],temperature[[7]],temperature[[8]])-273.15 # conversion to Celcius degrees
								temperature_inFall = mean(temperature[[9]],temperature[[10]],temperature[[11]])-273.15 # conversion to Celcius degrees
								precipitation = brick(paste0("Environmental_rasters/ISIMIP3b/pr_day_bias-adjusted_",scenarios[h],"_",models_isimip3b[i],"_",year_intervals[g],"_ymonmean.nc"))
								precipitation_winter = mean(precipitation[[12]],precipitation[[1]],precipitation[[2]])*60*60*24 # conversion to kg/m2/day
								precipitation_spring = mean(precipitation[[3]],precipitation[[4]],precipitation[[5]])*60*60*24 # conversion to kg/m2/day
								precipitation_summer = mean(precipitation[[6]],precipitation[[7]],precipitation[[8]])*60*60*24 # conversion to kg/m2/day
								precipitation_inFall = mean(precipitation[[9]],precipitation[[10]],precipitation[[11]])*60*60*24 # conversion to kg/m2/day
								relative_humidity = brick(paste0("Environmental_rasters/ISIMIP3b/hurs_day_bias-adjusted_",scenarios[h],"_",models_isimip3b[i],"_",year_intervals[g],"_ymonmean.nc"))
								relative_humidity_winter = mean(relative_humidity[[12]],relative_humidity[[1]],relative_humidity[[2]])
								relative_humidity_spring = mean(relative_humidity[[3]],relative_humidity[[4]],relative_humidity[[5]])
								relative_humidity_summer = mean(relative_humidity[[6]],relative_humidity[[7]],relative_humidity[[8]])
								relative_humidity_inFall = mean(relative_humidity[[9]],relative_humidity[[10]],relative_humidity[[11]])
								land_cover = nc_open(paste0("Environmental_rasters/ISIMIP3b/landcover_",scenarios[h],"_annual_",year_intervals[g],"_timmean.nc4"))
								population = brick(paste0("Environmental_rasters/ISIMIP3b/population_bias-adjustedsoc_5min_annual_",year_intervals[g],"_timmean.nc4"), varname="number_of_people")
								landCoverVariableIDs = names(land_cover$var); land_covers1 = list(); land_covers2 = list(); land_covers3 = list()
								landCoverVariableNames = as.character(read.csv("Environmental_rasters/Luse.csv")[1:12,2])
								for (j in 2:13)
									{
										land_covers1[[j-1]] = brick(paste0("Environmental_rasters/ISIMIP3b/landcover_",scenarios[h],"_annual_",year_intervals[g],"_timmean.nc4"), varname=landCoverVariableIDs[j])
									}
								variable_codes = c("croplands","pastures","urbanAreas","primaryForest","primaryNonF","secondaryForest","secondaryNonF")
								variable_names = c("crops","pasture","urban land","forested primary land","non-forested primary land",
										  		   "potentially forested secondary land","potentially non-forested secondary land")
								for (j in 1:length(variable_names))
									{
										names = gsub("\\."," ",landCoverVariableNames); indices = which(landCoverVariableNames==variable_names[j])
										if (length(indices) == 0) indices = which(grepl(variable_names[j],names))
										if (variable_names[j] == "pasture") indices = c(indices, which(grepl("rangeland",names)))
										land_cover = land_covers1[[indices[1]]]; names(land_cover) = variable_codes[j]; # print(indices)
										if (length(indices) > 1)
											{
												for (k in 2:length(indices)) land_cover[] = land_cover[]+land_covers1[[indices[k]]][]
											}
										land_covers2[[j]] = land_cover[[1]]; land_covers3[[j]] = raster::aggregate(land_cover[[1]],2)
									}
								envVariables = list()			
								envVariables[[1]] = temperature_winter; envVariables[[2]] = temperature_spring
								envVariables[[3]] = temperature_summer; envVariables[[4]] = temperature_inFall
								envVariables[[5]] = precipitation_winter; envVariables[[6]] = precipitation_spring
								envVariables[[7]] = precipitation_summer; envVariables[[8]] = precipitation_inFall
								envVariables[[9]] = relative_humidity_winter; envVariables[[10]] = relative_humidity_spring
								envVariables[[11]] = relative_humidity_summer; envVariables[[12]] = relative_humidity_inFall
								envVariables[[13]] = land_covers2[[4]] # primary forest areas
								envVariables[[14]] = land_covers2[[5]] # primary non-forest areas
								envVariables[[15]] = land_covers2[[6]] # secondary forest areas
								envVariables[[16]] = land_covers2[[7]] # secondary non-forest areas
								envVariables[[17]] = land_covers2[[1]] # croplands (all catergories)
								envVariables[[18]] = land_covers2[[2]] # managed pasture + rangeland
								envVariables[[19]] = population # human population (not log-transformed)
								for (j in 1:length(envVariables)) envVariables[[j]] = crop(envVariables[[j]], nutsM, snap="out")
								for (j in 1:length(envVariables)) envVariables[[j]] = mask(envVariables[[j]], nutsM)
								areas_nuts_M_km = area(nutsM, unit="km")
								for (j in 1:length(nutsM))
									{
										maxArea = 0; polIndex = 0
										for (k in 1:length(nutsM@polygons[[j]]@Polygons))
											{
												if (maxArea < nutsM@polygons[[j]]@Polygons[[k]]@area)
													{
														maxArea = nutsM@polygons[[j]]@Polygons[[k]]@area; polIndex = k
													}
											}
										pol1 = nutsM@polygons[[j]]@Polygons[[polIndex]]; p = Polygon(pol1@coords)
										ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
										pol2 = sf::st_as_sfc(sps); st_crs(pol2) = crs(nutsM)
										for (k in 1:18)
											{
												if ((j == 1)&(k == 1)) crs(envVariables[[k]]) = crs(pol2) 
												environmentalValues[j,k] = exactextractr::exact_extract(envVariables[[k]], pol2, fun="mean")
											}
										if ((j == 1)&(k == 1)) crs(envVariables[[19]]) = crs(pol2) 
										environmentalValues[j,19] = (log10(exactextractr::exact_extract(envVariables[[19]],pol2,fun="sum")+1))/areas_nuts_M_km[j]
									}
								nutsM_datas3[[g]] = environmentalValues
							}
						nutsM_datas2[[h]] = nutsM_datas3
					}
				nutsM_datas1[[i]] = nutsM_datas2
			}
		saveRDS(nutsM_datas1, "ISIMIP3b_dFrames.rds")
	}	else	{
		nutsM_datas1 = readRDS("ISIMIP3b_dFrames.rds")
	}
nutsM_projections1 = list()
for (h in 1:length(scenarios))
	{
		nutsM_projections2 = list()
		for (g in 1:length(year_intervals))
			{
				for (i in 1:length(models_isimip3b))
					{
						df = as.data.frame(nutsM_datas1[[i]][[h]][[g]]); colnames(df) = gsub("-","\\.",colnames(df))
						if (i == 1) nutsM_projection3 = matrix(nrow=dim(df)[1], ncol=length(models_isimip3b))
						nutsM_projection4 = matrix(nrow=dim(df)[1], ncol=length(brt_model_scvs))
						for (j in 1:length(brt_model_scvs))
							{
								n.trees = brt_model_scvs[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
								nutsM_projection4[,j] = predict.gbm(brt_model_scvs[[j]], df, n.trees, type, single.tree)
							}
						nutsM_projection5 = matrix(nrow=dim(df)[1], ncol=1)
						for (j in 1:dim(nutsM_projection4)[1])
							{
								nutsM_projection5[j,1] = mean(nutsM_projection4[j,])
							}
						nutsM_projection3[,i] = nutsM_projection5
					}
				nutsM_projection6 = matrix(nrow=dim(nutsM_projection3)[1], ncol=1)
				for (i in 1:dim(nutsM_projection3)[1])
					{
						nutsM_projection6[i,1] = mean(nutsM_projection4[i,])
					}
				nutsM_projections2[[g]] = nutsM_projection6
			}
		nutsM_projections1[[h]] = nutsM_projections2
	}
if (savingPlots)
	{
		pdf(paste0("ISIMIP3b_average_NEW.pdf"), width=8, height=5.8); par(mfrow=c(3,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")		
		colourScale = rev(colorRampPalette(brewer.pal(11,"RdBu"))(121)[11:111])
		for (h in 1:length(scenarios))
			{
				plot.new(); plot.new()
				for (g in 1:length(year_intervals))
					{
						cols = colourScale[(((nutsM_projections1[[h]][[g]]-0)/(1-0))*100)+1]
						plot(contour, lwd=0.4, border="gray30", col=NA)
						plot(nutsM, col=cols, border=NA, lwd=0.1, add=T); rast = raster(as.matrix(c(0,1)))
						mtext(scenario_names[h], side=3, line=-1.5, at=1, cex=0.50, col="gray30")
						mtext(gsub("_","-",year_intervals[g]), side=3, line=-2.3, at=1, cex=0.50, col="gray30")
					}
			}
		dev.off()
	}
brt_model_scvs = readRDS(paste0("All_the_BRT_models/gswp3-w5e5_models_SCV.rds"))
nutsM_projections1 = list()
for (i in 1:length(models_isimip3b))
	{
		nutsM_projections2 = list()
		for (h in 1:length(scenarios))
			{
				nutsM_projections3 = list()
				for (g in 1:length(year_intervals))
					{
						df = as.data.frame(nutsM_datas1[[i]][[h]][[g]]); colnames(df) = gsub("-","\\.",colnames(df))
						nutsM_projection4 = matrix(nrow=dim(df)[1], ncol=length(brt_model_scvs))
						for (j in 1:length(brt_model_scvs))
							{
								n.trees = brt_model_scvs[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
								nutsM_projection4[,j] = predict.gbm(brt_model_scvs[[j]], df, n.trees, type, single.tree)
							}
						nutsM_projections3[[g]] = nutsM_projection4
					}
				nutsM_projections2[[h]] = nutsM_projections3
			}
		nutsM_projections1[[i]] = nutsM_projections2
	}
if (savingPlots)
	{
		for (i in 1:length(nutsM_projections1))
			{
				pdf(paste0("ISIMIP3b_projections/ISIMIP3b_",models_isimip3b_names[i],".pdf"), width=8, height=5.8); par(mfrow=c(3,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")		
				colourScale = rev(colorRampPalette(brewer.pal(11,"RdBu"))(121)[11:111])
				for (h in 1:length(scenarios))
					{
						plot.new(); plot.new()
						for (g in 1:length(year_intervals))
							{
								nutsM_projection = matrix(nrow=dim(nutsM_projections1[[i]][[h]][[g]])[1], ncol=1)
								for (j in 1:dim(nutsM_projection)[1])
									{
										nutsM_projection[j,1] = mean(nutsM_projections1[[i]][[h]][[g]][j,])
									}
								cols = colourScale[(((nutsM_projection[,1]-0)/(1-0))*100)+1]
								plot(contour, lwd=0.4, border="gray30", col=NA)
								plot(nutsM, col=cols, border=NA, lwd=0.1, add=T); rast = raster(as.matrix(c(0,1)))
								mtext(scenario_names[h], side=3, line=-1.5, at=1, cex=0.50, col="gray30")
								mtext(gsub("_","-",year_intervals[g]), side=3, line=-2.3, at=1, cex=0.50, col="gray30")
							}
					}
				dev.off()
			}
	}
population_counts = matrix(nrow=dim(nutsM@data)[1], ncol=length(year_intervals))
colnames(population_counts) = year_intervals
for (i in 1:length(year_intervals))
	{
		population = brick(paste0("Environmental_rasters/ISIMIP3b/population_bias-adjustedsoc_5min_annual_",year_intervals[i],"_timmean.nc4"), varname="number_of_people")
		for (j in 1:length(nutsM))
			{
				maxArea = 0; polIndex = 0
				for (k in 1:length(nutsM@polygons[[j]]@Polygons))
					{
						if (maxArea < nutsM@polygons[[j]]@Polygons[[k]]@area)
							{
								maxArea = nutsM@polygons[[j]]@Polygons[[k]]@area; polIndex = k
							}
					}
				pol1 = nutsM@polygons[[j]]@Polygons[[polIndex]]; p = Polygon(pol1@coords)
				ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
				pol2 = sf::st_as_sfc(sps); st_crs(pol2) = crs(nutsM)
				if (j == 1) crs(population) = crs(pol2) 
				population_counts[j,i] = exactextractr::exact_extract(population, pol2, fun="sum")
			}
	}
if (savingPlots)
	{
		colours1 = c(rgb(194,230,153,255,maxColorValue=255), rgb(120,198,121,255,maxColorValue=255), rgb(49,163,84,255,maxColorValue=255))
		colours2 = c(rgb(194,230,153,125,maxColorValue=255), rgb(120,198,121,125,maxColorValue=255), rgb(49,163,84,125,maxColorValue=255))
		cols1 = c(rep(colours1[1], length(year_intervals)), rep(colours1[2], length(year_intervals)), rep(colours1[3], length(year_intervals)))
		cols2 = c(rep(colours2[1], length(year_intervals)), rep(colours2[2], length(year_intervals)), rep(colours2[3], length(year_intervals)))
		for (i in 1:length(thresholds))
			{
				tab = matrix(nrow=10*10, ncol=(3*4)); colNames = c()
				for (j in 1:length(scenarios))
					{
						for (k in 1:length(year_intervals))
							{
								colNames = c(colNames, paste0(scenarios[j],"_",year_intervals[k]))
								for (l in 1:length(models_isimip3b))
									{
										for (m in 1:dim(nutsM_projections1[[l]][[j]][[k]])[2])
											{
												vS = population_counts[,k]
												vS[which(nutsM_projections1[[l]][[j]][[k]][,m]<thresholds[i])] = 0
												tab[((l-1)*dim(nutsM_projections1[[l]][[j]][[k]])[2])+m,((j-1)*length(year_intervals))+k] = sum(vS)
											}
									}
							}
					}
				colnames(tab) = colNames
				if (savingFiles) write.csv(tab, paste0("ISIMIP3b_allModels_",gsub("\\.",",",as.character(thresholds[i])),".csv"), row.names=F, quote=F)
				pdf(paste0("ISIMIP3b_violin_",gsub("\\.",",",as.character(thresholds[i])),"_NEW.pdf"), width=4, height=2.2)
				par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1.5,1.5,0,0), mgp=c(1.2,0.75,0), lwd=0.2, col="gray30")
				# if (i == 1) plot(0:1, 0:1, type="n", xlim=c(190.5,201.5), ylim=c(0.9*(10^8),1.9*(10^8)), axes=F, ann=F)
				if (i == 1) plot(0:1, 0:1, type="n", xlim=c(198.5,209.5), ylim=c(1.5*(10^8),4.2*(10^8)), axes=F, ann=F)
				if (i == 2) plot(0:1, 0:1, type="n", xlim=c(198.5,209.5), ylim=c(0,10.1*(10^7)), axes=F, ann=F)
				vioplot(tab, ann=F, axes=F, at=rep(c(2030,2050,2070,2090)/10,3), border=cols1, col=cols2, use.cols=T, horizontal=F, lineCol=NA, rectCol=NA, colMed=NA, add=T)
				axis(1, at=c(1970,1990,2010,2030,2050,2070,2090,2110,2130)/10, labels=c(1970,1990,2010,2030,2050,2070,2090,2110,2130), mgp=c(0,0.06,0), lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.025, col="gray30", col.axis="gray30", col.lab="gray30")
				if (i == 1) axis(2, at=seq(1.0*(10^8),4.5*(10^8),5*(10^7)), labels=c(100,150,200,250,300,350,400,450), mgp=c(1,0.20,0), lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.02, col="gray30", col.axis="gray30", col.lab="gray30")
				if (i == 2) axis(2, at=seq(-10^7,12*(10^7),2*(10^7)), labels=c(0,20,40,60,80,100,120), mgp=c(1,0.20,0), lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.02, col="gray30", col.axis="gray30", col.lab="gray30")
				if (h == 2) title(ylab="Population at risk (in million people)", cex.lab=0.80, mgp=c(1.5,0,0), col.lab="gray30")
				dev.off()
			}
	}

# 6. Compressing the environmental files for GitHub

directories = c("Environmental_rasters/ISIMIP3a","Environmental_rasters/ISIMIP3b"); wd = getwd()
for (i in 1:length(directories))
	{
		files = list.files(directories[i]); files = files[which(grepl("\\.nc",files))]
		setwd(paste0(wd,"/",directories[i]))
		for (j in 1:length(files))
			{
				system(paste0("tar -zcvf ",gsub("\\.nc","\\.tar.gz",gsub("\\.nc4","\\.tar.gz",files[j]))," ",files[j]))
			}
		setwd(wd)
	}
	

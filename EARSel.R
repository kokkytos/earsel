
setwd(dirname(rstudioapi::getSourceEditorContext()$path))#for Rscript inside RStudio


## Settings ------------------------------------------------------------------------
source('config/settings.R')




## Libraries ----------------------------------------------------------------------
library(raster)
library(MODIS)
library(beepr)
library(rgdal)

MODISoptions(MODISserverOrder = c("LPDAAC","LAADS")) #run lpdaacLogin(server = "LPDAAC") first, writes cred. in ~/.netrc, https://www.rdocumentation.org/packages/MODIS/versions/1.1.0/topics/lpdaacLogin

## Setup some VARS
NUTSSHP<-readOGR(NUTS_DIR,NUTS_SHP)

YEARSSEQ <-seq(as.Date(ndvistartdate), as.Date(ndvienddate), "year")
years<-format(as.Date(YEARSSEQ, format = "%Y-%m-%d"),"%Y")

years<-setdiff(years, c("1999", "2000","2002","2005")) #years to exclude,#tags:subset,years




## Custom functions --------------------------------------------------------------

# create mask file
createMask<-function(shp, NUTS_IDs,outproj){#better implementation with proj
  nuts<- subset(shp, NUTS_ID %in% NUTS_IDs)
  nuts$NUTS_ID2<-as.numeric(1.00)
  nuts <- sp::spTransform(nuts, CRS(outproj))
  return(nuts)
}


# Calculate annual mean
makeAnnualMean<-function(year){
  print(paste(length(year), " files to calculate mean"))
  s<-stack(year)
  annualNDVI <- mean(s, na.rm=TRUE) # for parallelization check https://stackoverflow.com/questions/41318164/r-faster-mean-calculation-in-raster-stack
  
  scalefactor<-0.0001 # see.Product page
  annualNDVI <- annualNDVI*scalefactor
}

# generate filename by year
genfilenameYear<-function(year, template, DIR){
  file<-sprintf(template, year)
  file<-file.path(DIR,file)
  file<-raster(file)
}

#function to calculate SANUI
CalcSanui<-function(ndvi,SeaWinds,ols  ){
  ndviPositive<-ndvi
  ndviPositive[ndvi<0]<-NA#set negative values to NA
  sanui <- (1-ndviPositive)*SeaWinds*(ols/satDN)#satDN=63
  
}


#todo:maybe useless
#'function:Get city extent from shapes
#'
#' @param dns 
#' @param NUTS_IDs, a vector with NUTS_IDs to subset shp
#' @return an extent object of the selected features
#' @export
#' @examples
getCityExtent<-function(shp, NUTS_IDs){
  nuts3subs<- subset(shp, NUTS_ID %in% NUTS_IDs)
  e<-extent(nuts3subs)
  ext<-list(extentobj = e, proj4string = shp@proj4string)
  return (ext)
}

#todo:maybe useless
#function:reprojects extent object
#'
#' @param ext The extent object to be reprojected
#' @param inproj The current CRS of extent obj, as proj4string 
#' @param outproj The new CRS of extent obj, as proj4string 
#' @return The reprojected extent object
#' @export
reprExt<-function(ext, inproj, outproj){
  ext <- as(ext, "SpatialPolygons")
  sp::proj4string(ext) <- inproj
  extout <- sp::spTransform(ext, CRS(outproj))
  extout <-extent(bbox(extout))
  return(extout)
}




fSol <- function(x,s) {
  r<-s[[x]]#the raster obj
  r[r[] ==  255] <- NA    # clouds (255) set to NoData (for orignal stable light) !!!!!!!!!!!!!!!!
  r[r[] <=    6] <- NA
  cellStats(r, stat='sum', na.rm=TRUE)
}


fAol<-function(r,mask) {
  r<-setValues(r,1)
  r<-mask(r, mask,  updatevalue=NA)#todo:useless , ols alread masked
  ncell(r[!is.na(r[])])# [todo:is it OK]
  }

#' Aggregates the population for a vector of NUTS3 units of a year
#'
#' @param df dataframe with population data 
#' @param nuts_ids the NUTS3 of a city (ids of GEO column of dataframe)
#' @param year the year for aggregation
#'
#' @return the population as integer
#' @export
#'
#' @examples getPopByYear(df,  c("EL301","EL302"), 2001))
SumVarByNutsidsYear<-function(df, nuts_ids, year){
  dfsub<-subset(df, GEO %in% nuts_ids & TIME== year)
  #check for NA or zero in data and write to log file
  if (any(is.na(dfsub$newValue)||dfsub$newValue==0)){
    nullv<-dfsub[dfsub$newValue==0 | is.na(dfsub$newValue),]
    write.csv(nullv, file = paste("na_zero_",year,"_",deparse(substitute(df)),"_",mycity,".csv"))
    print(paste("NA,zero values in Year:", year, ", Dataframe:",deparse(substitute(df)),", City:",mycity))
  }
  
  sum(dfsub$newValue)
}




#uncomment to use them with loop, check closing } at the end
finalfunction<-function(x){
  
  #for the city get: NUTS_ID, extent, rextent, name and plot
  city<-x #cities[[2]]#Reference by city index, [todo: replace cities[[1]] with x for apply loop]
  NUTS_IDs<-city$NUTS_ID#get a vector with nuts ids
  ext<-getCityExtent(NUTSSHP, NUTS_IDs)#city extent (all NUTS3)
  ext<-reprExt(ext$extentobj,ext$proj4string,moll) #reproject extent to moll
  
  mycity<-as.character(city[1,"city"]) #city name
  
  citymask<-createMask(NUTSSHP, NUTS_IDs,moll)# createMask returns SpatialPolygonsDataFrame
  
  plot(citymask,main=mycity)
  
  
  ## City Population ----------------------------------------------------------------------------
  
  popcsv <- read.csv(file=POP_CSV, header=TRUE, sep=",",stringsAsFactors = F)
  popcsv$newValue <- suppressWarnings(as.numeric(gsub(",", "", as.character(popcsv$Value))))#new column with population data as numeric 
  
  pop<-sapply(years, function(x) SumVarByNutsidsYear(df=popcsv,nuts_ids=NUTS_IDs, year=x)) #sum population for NUTS3 for each year
  
  
  
  
  ## City GDP ----------------------------------------------------------------------------
  
  gdpcsv <- read.csv(file=GDP_CSV, header=TRUE, sep=",",stringsAsFactors = F)
  gdpcsv$newValue <- suppressWarnings(as.numeric(gsub(",", "", as.character(gdpcsv$Value))))#new column with population data as numeric 
  
  gdp<-sapply(years, function(x) SumVarByNutsidsYear(df=gdpcsv,nuts_ids=NUTS_IDs, year=x))#sum GDP for NUTS3 for each year
  
  
  
  
  ## DMSP/OLS ----------------------------------------------------------------------------------
  
  #batch generate rasters
  ols<-lapply(years, FUN=genfilenameYear,template="FIL%s.tif", DIR=OLSDIR)
  
  #batch crop rasters
  ols<-lapply(ols, FUN = crop, y= ext, snap='near',datatype ="INT1U")#todo να φύγει από δω (?)
  
  #stack rasters
  ols<-brick(ols)
  ols[ols>satDN]<-satDN #set values>saturationDN to saturationDN(63)
  
  #mask
  ols<-mask(ols, citymask,  updatevalue=NA)
  #rts_ols<-rts(ols,YEARSSEQ)
  
  
  
  
  ## ---NDVI---------------------------------------------------------------------------------
  
  cll <- getCollection(product = ndviproduct, forceCheck = TRUE)
  ndvi<-MODIS::runGdal( # run MODIS:::checkTools('GDAL') to check for GDAL library
    job = mycity,
    product = ndviproduct,
    extent =  ols[[1]],#extent.wgs84,#todo:pass raster OLS
    SDSstring="1",
    collection = cll,
    begin = ndvistartdate,
    end = ndvienddate,
    outDirPath = "NDVI",
    overwrite= TRUE,
    checkIntegrity = TRUE,
    #quiet=TRUE,
    wait = 20
  )
  
  names(ndvi)<-c(mycity)#set city name
  names(ndvi[[mycity]])<-format(as.Date(as.character(names(ndvi[[mycity]])), format = "%Y-%m-%d"),"%Y") #convert %Y-%m-%d to %Y
  
  ndvi<-lapply(split(x = ndvi[[mycity]], f = names( ndvi[[mycity]])), unlist) #split and apply, for each year create a list of filenames
  
  ndvi<-ndvi[years]#subset, keep only working years. #todo:remove to keep all years,#tags:subset,years
  ndvibyear<-sapply(ndvi, makeAnnualMean)
  
  
  ndvibyear<-brick(ndvibyear)
  #mask ndvi
  ndvibyear<-mask(ndvibyear, citymask,  updatevalue=NA)
  
  
  
  
  ## SEAWINDS -------------------------------------------------------------------
  #batch generate rasters
  SW<-lapply(years, FUN=genfilenameYear,template="Global_quev_%s_JFM_PR.tif", DIR=SEAWINDSDIR)
  
  #batch reproject rasters
  SW<-lapply(SW, 
             FUN=projectRaster, 
             to = ols[[1]],
             method    = 'ngb',
             alignOnly = FALSE,
             over      = FALSE
  )
  
  #stack rasters
  SW<-brick(SW)
  SW<-SW/maxValue(SW)#normalize Seawinds
  
  #mask
  SW<-mask(SW, citymask,  updatevalue=NA)
  names(SW)<-years
  
  
  
  
  ## AoL - SoL - iSoL  --------------------------------------------------------------------------
  aol<-fAol(ols[[1]],citymask)
  sol<-sapply(1:nlayers(ols), function(x) fSol(x, s=ols))
  names(sol)<-years
  
  #as time series obj
  #sol<-ts(sol,
  #        start=c(min(as.numeric(years))),
  #        end=c(max(as.numeric(years))), 
  #        frequency=1) 
  
  isol<-sol/aol/maxValue(ols) #todo:replace 63 but first set ols>63 to 63
  #names(isol)<-years
  
  
  
  
  ## VANUI -----------------------------------------------------------------------
  vanui <- (1-ndvibyear)*(ols/satDN)#returns rastet brick, ndvibyear,ols are rbricks
  #vanui<-mask(vanui, citymask,  updatevalue=NA)#todo:useless , input bricks already masked
  vanui<-cellStats(vanui, stat='sum', na.rm=TRUE)
  #rts_vanui<-rts(vanui,YEARSSEQ)
  vanui<-vanui/aol #normalized
  
  
  
  
  ## SANUI ------------------------------------------------------------------------
  sanui<-CalcSanui(ndvibyear,SW,ols)
  sanui<-mask(sanui, citymask,  updatevalue=NA)#todo:useless , ols alread masked
  sanui<-cellStats(sanui, stat='sum', na.rm=TRUE)
  sanui<-sanui/aol #normalized
  #rts_sanui<-rts(sanui,YEARSSEQ)
  
  
  
  ## Correlations ---------------
  
  #dataframe to save results correlation
  results.df <- data.frame(matrix(ncol = 4, nrow = 0))
  names(results.df)<-c("City", "Var1", "Var2", "Pearson")
  
  cordf = data.frame( rep(mycity, times=6),
                      c(rep("isol", times=2),rep("vanui", times=2) , rep("sanui", times=2)) , 
                      rep(c("pop", "gdp"), times =3),
                      c(
                        cor(isol,pop),
                        cor(isol,gdp), 
                        cor(vanui,pop),
                        cor(vanui,gdp),
                        cor(sanui,pop),
                        cor(sanui,gdp)
                        
                      )
  )
  names(cordf)<-c("City", "Var1", "Var2", "Pearson")
  results.df<-rbind(results.df, cordf)
  # Write csv to df
  write.csv(results.df,file.path(OUTPUTDIR, paste(mycity,".cor.csv",sep = "")),row.names=TRUE)
  
  
  ## save data to csv -----------------------------------------------------------------
  datadf <- data.frame(matrix(ncol = length(years), nrow = 0))
  
  datadf<-rbind(datadf, pop)
  datadf<-rbind(datadf, gdp)
  datadf<-rbind(datadf, isol)
  datadf<-rbind(datadf, vanui)
  datadf<-rbind(datadf, sanui)
  
  names(datadf)<-years
  row.names(datadf) <-c("Popupation", "GDP", "iSol", "SumVanui", "SumSanui")
  write.csv(datadf,file.path(OUTPUTDIR, paste(mycity, ".data.csv",sep = "")),row.names=TRUE)
  
  
  ## Results ------------------------------------------------------------------------------------
  dimnames(mycity)<-NULL#todo
  res<-list(mycity,pop, gdp, ols, ndvibyear, SW, vanui,sanui,aol, sol, isol)
  names(res)<-c('city','Population', 'GDP', 'DMSP/OLS','NDVI', 'SeaWinds',
                'vanui','sanui', 'AoL', 'SoL', 'iSoL')
  return(res)#todo:check return vars 

} #also uncomment return(res)


## Apply to each city  -----------------------------------------------------------------
result<-lapply(cities,  finalfunction)









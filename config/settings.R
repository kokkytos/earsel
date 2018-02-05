setwd("config")

csv <- read.csv(file="cities.csv", header=TRUE, sep=",",stringsAsFactors = TRUE)
cities<-split(csv, csv$city)


setwd("../")
#set CRS definitions
wgs84<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
moll<-"+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


OUTPUTDIR<-"output"

OLSDIR<-'OLS'
OLS_FILE<-file.path(OLSDIR, 'FIL2013.tif')# DMSP/OLS file, source: dstath
satDN<-63 # saturation DN of DMSP/OLS


NUTS_DIR<-"NUTS/NUTS_2013_01M_SH/data"
NUTS_SHP<-"NUTS_RG_01M_2013"


POP_CSV<-"NUTS/pop/demo_r_pjanaggr3_1_Data.csv"
GDP_CSV<-"NUTS/gdp/nama_10r_3gdp_1_Data.csv"


SEAWINDSDIR<-'SEAWINDS'

#MODIS
ndvistartdate <-"2000-02-01" #MOD13A3 Temporal Coverage	February 18, 2000 -
ndvienddate <-"2009-12-31"


ndviproduct <- "MOD13A3" # Vegetation Indices Monthly L3 Global 1km: https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13a3
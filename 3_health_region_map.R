################################################################################
## Map regional variation in low-value ECG use 
################################################################################

## -----------------------------------------------------------------------------
## SECTION 1: SETUP
## -----------------------------------------------------------------------------
# Load required libraries
library(sf)  
library(rgdal) 
library(maptools)
library(goeveg)

# Import shapefiles
thailand.sf <- read_sf("~/Thailand Shapefiles/tha_admbnda_adm1_rtsd_20190221.shp")
thailand.spdf <- readOGR("~/tha_admbnda_adm1_rtsd_20190221.shp")

# Import regional prevalence data from script 1
lvc_prev <- read.csv("~/lvc_regional_prev_FY18-19.csv")

# Calculate coefficient of quartile deviation
q1 <- quantile(lvc_prev$perc, prob = 0.25)[[1]]
q3 <- quantile(lvc_prev$perc, prob = 0.75)[[1]]
cqd <- (q3-q1)/(q3+q1)

## -----------------------------------------------------------------------------
## SECTION 2: MAKE HEALTH REGION BASE MAP
## -----------------------------------------------------------------------------
# Convert provinces into NHSO-defined health regions
thailand.sf <- thailand.sf %>%
  mutate(nhso_region = case_when(ADM1_EN %in% c("Chiang Mai","Chiang Rai","Nan","Phayao","Phrae",
                                                "Mae Hong Son","Lampang","Lamphun") ~ 1,
                                 ADM1_EN %in% c("Phitsanulok","Phetchabun","Tak","Sukhothai","Uttaradit") ~ 2,
                                 ADM1_EN %in% c("Nakhon Sawan","Kamphaeng Phet","Chai Nat","Phichit","Uthai Thani") ~ 3,
                                 ADM1_EN %in% c("Nonthaburi","Saraburi","Lop Buri","Nakhon Nayok","Pathum Thani",
                                                "Phra Nakhon Si Ayutthaya","Ang Thong","Sing Buri") ~ 4,
                                 ADM1_EN %in% c("Ratchaburi","Kanchanaburi","Prachuap Khiri Khan","Phetchaburi",
                                                "Samut Songkhram","Nakhon Pathom","Suphan Buri","Samut Sakhon") ~ 5,
                                 ADM1_EN %in% c("Rayong","Chanthaburi","Chachoengsao","Chon Buri","Trat",
                                                "Prachin Buri","Samut Prakan","Sa Kaeo") ~ 6,
                                 ADM1_EN %in% c("Khon Kaen","Kalasin","Maha Sarakham", "Roi Et") ~ 7,
                                 ADM1_EN %in% c("Sakon Nakhon","Nakhon Phanom","Nong Khai",
                                                "Loei","Nong Bua Lam Phu", "Bueng Kan","Udon Thani") ~ 8,
                                 ADM1_EN %in% c("Nakhon Ratchasima", "Chaiyaphum","Buri Ram","Surin") ~ 9,
                                 ADM1_EN %in% c("Ubon Ratchathani","Mukdahan","Yasothon","Si Sa Ket",
                                                "Amnat Charoen") ~ 10,
                                 ADM1_EN %in% c("Surat Thani","Krabi","Chumphon","Phuket","Nakhon Si Thammarat",
                                                "Phangnga","Ranong") ~ 11,
                                 ADM1_EN %in% c("Songkhla","Trang","Narathiwat","Pattani","Phatthalung",
                                                "Yala","Satun") ~ 12,
                                 ADM1_EN == "Bangkok" ~ 13))


## -----------------------------------------------------------------------------
## SECTION 3: MAP LVC RATES
## -----------------------------------------------------------------------------

# Join prevalence data with shapefile
thailand.sf <- thailand.sf %>%
  left_join(lvc_prev, by = 'nhso_region')

# Plot map
map <- thailand.sf %>%
  group_by(nhso_region)  %>% 
  ggplot() +
  geom_sf(data = thailand.sf, aes(fill = perc, col = perc)) +
  theme_void() +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  labs(fill="Potential low-value care rate",
       col="Potential low-value care rate")

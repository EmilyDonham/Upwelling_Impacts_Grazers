#Written by Emily Donham 11/1/2020
#Plots processed pH data

######################################################################################################
######################################################################################################
library(plyr); library(dplyr);library(broom);library(ggplot2); library(lubridate);library(LoLinR); library(stringr)
library(lme4); library(lmerTest); library(multcomp); library(phytotools); library(googledrive); library(Rmisc)
library(tibble); library(ggpubr); library(wesanderson); library(tidyverse);library(vegan);
library(lsmeans); library(RLRsim); library(pracma)
######################################################################################################
######################################################################################################

rm(list = ls())

#Import files
L = data.frame(list.files(pattern = "SWC_*")) #Files that start with code for Stillwater Cove
names(L)[names(L) =="list.files.pattern....SWC_..."] <- "Dep"
as.character(L$Dep)

#Create empty dataframe
pHDatFinal = data.frame(Date_UTC = factor(),Time_UTC= factor(), Depth = integer(), 
                        Temperature = double(), pH_int = double(),
                        pH_extT = double(),  pH_extF = double(), OxyUM = double(), 
                        OxymgL = double(), Sal = double(), QC = integer(), DT = POSIXct(), 
                        i = integer())  #Initialize empty dataframe

#Loops through and concatenates data from separate deployments
for (i in 1:length(L$Dep)) {
  Dep = L$Dep[i]
  pHDat = read.csv(paste(L$Dep[i], sep=""),
                   skip = 4)
  pHDat$DT <- as.POSIXct(paste(pHDat$Date_UTC, pHDat$Time_UTC), format = "%m/%d/%Y %H:%M:%S")
  pHDat$i <- i
  pHDatFinal <- rbind(pHDatFinal, pHDat)
  }

pHDatGood <- subset(pHDatFinal, pHDatFinal$QC == 1)
pHDatGood <- subset(pHDatGood, pHDatGood$DT < "2020-10-15 00:00:00")
#Removes data that aren't good, but have passed minimum thresholds
pHDatGood$pH_int <- ifelse(pHDatGood$DT < '2016-03-11' | pHDatGood$DT > '2016-05-20', pHDatGood$pH_int, NaN)
pHDatGood$pH_int <- ifelse(pHDatGood$DT < '2017-03-15' | pHDatGood$DT > '2017-04-20', pHDatGood$pH_int, NaN)
pHDatGood$pH_int <- ifelse(pHDatGood$DT < '2017-12-09' | pHDatGood$DT > '2018-01-15', pHDatGood$pH_int, NaN)

#Create dataframe with average values of temp, pH, DO for plotting
Exp <- c("Experiment 1","Experiment 1","Experiment 1","Experiment 1","Experiment 1","Experiment 1",
         "Experiment 2","Experiment 2","Experiment 2","Experiment 2","Experiment 2","Experiment 2")
Temp <- c(12.17, 11.7, 11.04, 10.68, 10.4, 10.33, 13.8, 13.3, 12.6, 11.8, 11.7, 11.3)
pH <- c(7.79, 7.73, 7.65, 7.58, 7.52, 7.49, 8.01, 7.82, 7.75, 7.66, 7.62, 7.56)
DO <- c(7.64, 6.89, 5.86, 5.18, 4.78, 4.69, 8.82, 7.68, 6.62, 5.58, 5.43, 4.67)
Aves <- data.frame(Exp,Temp,pH,DO) #Combine into one dataframe

###################################################################################################
#Figure 1 - plot of time series data
a <- ggplot(pHDatGood, aes(x = DT, y = pH_int, group = i)) + 
    geom_line() +
    theme_classic() +
    labs(x="Date (MM/YY)", y = "pH", element_text(size = 20)) +
    theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
    scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
a

b <- ggplot(pHDatGood, aes(x = DT, y = Temperature, group = i)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
b

c <- ggplot(pHDatGood, aes(x = DT, y = OxymgL, group = i)) + 
  geom_line() +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x="Date (MM/YY)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0))
c

#Figure 1 - plot of time series data
l <- ggarrange(a, b, c, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, legend = "bottom", common.legend = TRUE)
l
ggsave(plot = l, file = "Stillwater_timeseries_Month.png", 
       type = "cairo-png",  bg = "white",
       width = 20, height = 30, units = "cm", dpi = 300)
###################################################################################################

#color palette
pal = wes_palette("Zissou1",6,type = "continuous")

#Create "date" column with common year to ease separating into seasons
pHDatGood <- pHDatGood %>%
  mutate(date=ymd_hm(format(pHDatGood$DT, "2016-%m-%d-%H:%M")))

#Get rid of NA data for scatterplots
pHDatGood_3 <- na.omit(pHDatGood)

#Filter into upwelling season
pHDatGood_UpwellingAll = pHDatGood_3 %>%
  filter(pHDatGood_3$date < '2016-10-01' & pHDatGood_3$date > '2016-04-01') 

#Filter into "non" upwelling season
pHDatGood_NonUpwellingAll = pHDatGood_3 %>%
  filter(pHDatGood_3$date < '2016-3-01' | pHDatGood_3$date > '2016-11-01') 

###################################################################################################
#Figure 2 - Create scatterplots of sensor data with average mesocosm experiment
e <- ggplot(pHDatGood_UpwellingAll, aes(x = pH_int, y = OxymgL)) + 
  geom_point(size = 0.7, position = "jitter") +
  geom_point(data = Aves, aes(x = pH, y = DO, color = Exp), size = 2) +
  scale_color_brewer(palette = "Dark2") +
  stat_smooth(method = "lm" , formula = y~x, se = FALSE, color = pal[4]) +
  labs(y="Dissolved Oxygen (mg/L)", x="pH") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = c(0.85, 0.15), legend.key.width = unit(0.5,'cm'))
e

f <- ggplot(pHDatGood_UpwellingAll, aes(x = pH_int, y = Temperature)) + 
  geom_point(size = 0.7, position = "jitter") +
  geom_point(data = Aves, aes(x = pH, y = Temp, color = Exp), size = 2) +
  scale_color_brewer(palette = "Dark2") +
  stat_smooth(method = "lm" , formula = y~x, se = FALSE, color = pal[4]) +
  labs(y="Temperature (\u00B0C)", x="pH") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = c(0.85, 0.15), legend.key.width = unit(0.5,'cm'))
f

g <- ggplot(pHDatGood_UpwellingAll, aes(x = OxymgL, y = Temperature)) + 
  geom_point(size = 0.7, position = "jitter") +
  geom_point(data = Aves, aes(x = DO, y = Temp, color = Exp), size = 2) +
  scale_color_brewer(palette = "Dark2") +
  stat_smooth(method = "lm" , formula = y~x, se = FALSE, color = pal[4]) +
  labs(y="Temperature (\u00B0C)", x="Dissolved Oxygen (mg/L)") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = c(0.85, 0.15), legend.key.width = unit(0.5,'cm'))
g

#Figure 2 - Create scatterplots of sensor data with average mesocosm experiment
p <- ggarrange(e, f, g, labels = c("(a)", "(b)", "(c)"), ncol = 3, nrow = 1, legend = "bottom", common.legend = TRUE)
p
ggsave(plot = p, file = "Stillwater_Reg_plots.png", 
       type = "cairo-png",  bg = "white",
       width = 30, height = 10, units = "cm", dpi = 300)
###################################################################################################

#Regressions of time series separated by "season"
pHTC <- lm(pHDatGood_UpwellingAll$Temperature ~ pHDatGood_UpwellingAll$pH_int)
summary(pHTC)
pHDO <- lm(pHDatGood_UpwellingAll$OxymgL ~ pHDatGood_UpwellingAll$pH_int)
summary(pHDO)
DOTC <- lm(pHDatGood_UpwellingAll$Temperature ~ pHDatGood_UpwellingAll$OxymgL)
summary(DOTC)

pHTC <- lm(pHDatGood_NonUpwellingAll$Temperature ~ pHDatGood_NonUpwellingAll$pH_int)
summary(pHTC)
pHDO <- lm(pHDatGood_NonUpwellingAll$OxymgL ~ pHDatGood_NonUpwellingAll$pH_int)
summary(pHDO)
DOTC <- lm(pHDatGood_NonUpwellingAll$Temperature ~ pHDatGood_NonUpwellingAll$OxymgL)
summary(DOTC)


#Calculate offset of DO at given pH
DOexp <- 11.36266*Aves$pH -83.3112
DOoff <- Aves$DO - DOexp
DOoff

Aves <- data.frame(Aves, DOexp, DOoff) #Combine into one dataframe
#Calculate average offset in DO for each experiment
AvesSum <- summarySE(data=Aves, measurevar = "DOoff",
                      groupvars = c("Exp"), na.rm = TRUE)


###################################################################################################
#Figure S1 - Create scatterplots of sensor data by season
d <- ggplot(pHDatGood_UpwellingAll, aes(x = pH_int, y = OxymgL, color = Temperature)) + 
  geom_point(size = 1.5, position = "jitter") +
  scale_color_gradientn(colours = pal, name = "Temp (\u00B0C)", limits = c(8, 18)) +
  labs(y="Dissolved Oxygen (mg/L)", x="pH") +
  theme_classic() +
  #ggtitle("Upwelling") +
  ylim(2, 10) +
  xlim(7.6, 8.1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
d <- d + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5))+ scale_y_continuous(limits = c(2,10), breaks = c(3,6,9))
d

e <- ggplot(pHDatGood_NonUpwellingAll, aes(x = pH_int, y = OxymgL, color = Temperature)) + 
  geom_point(size = 1.5, position = "jitter") +
  scale_color_gradientn(colours = pal, name = "Temp (\u00B0C)", limits = c(8, 18)) +
  labs(y="", x="pH") +
  theme_classic() +
  ylim(2, 10) +
  xlim(7.6, 8.1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
e <- e + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5))  + scale_y_continuous(limits = c(2,10), breaks = c(3,6,9))
e

#Figure S1 - Create scatterplots of sensor data by season
q <- ggarrange(d, e, labels = c("(a)", "(b)"), ncol = 2, nrow = 1, legend = "bottom", common.legend = TRUE)
q

ggsave(plot = q, file = "Up_NonUP_scatter.png", 
       type = "cairo-png",  bg = "white",
       width = 20, height = 10, units = "cm", dpi = 300)
###################################################################################################


#Respirometry code Written by Umi Hoshijima ~ 09/2018
#Edited by Emily Donham with additonal analyses ~ 02/2020

# LoLinR is not on CRAN, so install with: 
#install_github('colin-olito/LoLinR')
######################################################################################################
######################################################################################################
library(plyr); library(dplyr);library(broom);library(ggplot2); library(lubridate);library(LoLinR); library(stringr)
library(lme4); library(lmerTest); library(multcomp); library(phytotools); library(googledrive); library(Rmisc)
library(tibble); library(ggpubr); library(wesanderson); library(tidyverse);library(vegan);
library(lsmeans); library(RLRsim); library(gridExtra); library(ggfortify); library(effects); library(remef)

######################################################################################################
######################################################################################################

rm(list = ls())

# Volume of vials used; S = small, L = large
Svial_vol_ml = 99.957
Lvial_vol_ml = 229.7059

# This is going to read in one of the respirometry sheets found in biological data: 
datUrchin = read.csv('RespAllSF_QC.csv')

#hm changes it into "time object"
# seconds converts into second since midnight
# as.integer just makes it an integer
TIME8 = as.integer(seconds(hm(datUrchin$TIME1)))

# The respirometry data was entered just how it was taken - with multiple columns of oxygen. I want this in 
#"long form" - that is, with 1 column for time, another for oxygen, another for temp. 
dat_long = data.frame(Vial = NA,Critter.ID=NA, Time = NA, Oxygen = NA, Temp = NA, Run = NA, TRT = NA, CS = NA)  #Initialize empty dataframe

for (b in seq(from=7, to=25, by=3)){ #rbind the proper columns, so column 4 (ID) and the others. 
  # This takes the time from the column, and changes it to "seconds since capped"
  datUrchin[,b] = as.integer(seconds(hm(datUrchin[,b])))-TIME8
  temp = datUrchin[,c(2,4,((b):(b+2)),1, 40, 41)] # This is how to do 4,5,6,7,8 -- 4,9,10,11,12
  # This adds a "vial" column - which you can ignore if you already have one. 
  colnames(temp) = colnames(dat_long) #column names for rbind
  dat_long = rbind(dat_long, temp)
}

# gets rid of any data without a corresponding oxygen value
dat_long = dat_long[complete.cases(dat_long$Oxygen),]

# Now going onto LoLinR to use local linear regression to calc respiration rate: 
dat_long$combined <- paste(dat_long$Critter.ID,dat_long$Run)
u <-data.frame(unique(dat_long$combined))
f = data.frame(Slope = NA,Sig=NA, Temp = NA, SAMPLE = NA, TIMESTEP = NA)  #Initialize empty dataframe

for (i in 1:nrow(u)) { 
  x <- dplyr::select(filter(dat_long, combined == u[i,]),c(Vial, Critter.ID, Time, Oxygen, Temp, Run))
  mean_temp = mean(x$Temp, na.rm = TRUE)
  regression = rankLocReg(xall = x$Time, yall = x$Oxygen, alpha = .7)
  summary(regression)
  a = regression$allRegs
  f[i,] <- c(a[1,]$b1 * -1,a[1,]$b1LoCI > 0 | a[1,]$b1UpCI<0, mean_temp, x$Critter.ID[1], x$Run[1])
}

#Add some identifying columns to the results of LoLin
dat_sum1 <- merge(datUrchin[,1:6],f, by=c("SAMPLE","TIMESTEP"))

#Create columns for calculating out rates
dat_sum1$resp_mg_sec <- NA
dat_sum1$resp_g_sec <- NA
dat_sum1$resp_mol_sec <- NA
dat_sum1$resp_umol_min <- NA

#Make slope numeric, chr for some reason...
dat_sum1$Slope = as.numeric(dat_sum1$Slope)

#Filter for controls only
Controls1 = dat_sum1 %>%
  filter(CS == "C") 

#Calc mean of controls
Con_Sum1 <- summarySE(data=Controls1, measurevar = "Slope",
                      groupvars = c("TIMESTEP","RUN"), na.rm = TRUE)
Con_Sum2 <- summarySE(data=Controls1, measurevar = "Slope",
                      groupvars = c("TIMESTEP"), na.rm = TRUE)
#Subset only samples
Samp1 = dat_sum1 %>%
  filter(CS == "S")

#Adds run specific ave control to each sample
Samp1 <- merge(Con_Sum1, Samp1, by.x = c("TIMESTEP", "RUN"), by.y = c("TIMESTEP", "RUN"))

#Correct slope of samples by subracting the ave slope of controls in a given run
Samp1$Slope.y = Samp1$Slope.y-Samp1$Slope.x

# Slopes are in mg/L per second. convert to something else: 
# firstly multiply by volume of vial in liters to convert to mg per second: 
for (i in 1:nrow(Samp1)) {
  if(Samp1$VTYPE[i]=='S') {
    Samp1$resp_mg_sec[i] = Samp1$Slope.y[i] * (Svial_vol_ml/1000)
    Samp1$resp_g_sec[i] = Samp1$resp_mg_sec[i] / 1000
    Samp1$resp_mol_sec[i] = Samp1$resp_g_sec[i] / 31.998
    Samp1$resp_umol_min[i] = Samp1$resp_mol_sec[i] * 1000000 * 60
  }else{
    Samp1$resp_mg_sec[i] = Samp1$Slope.y[i] * (Lvial_vol_ml/1000)
    Samp1$resp_g_sec[i] = Samp1$resp_mg_sec[i] / 1000
    Samp1$resp_mol_sec[i] = Samp1$resp_g_sec[i] / 31.998
    Samp1$resp_umol_min[i] = Samp1$resp_mol_sec[i] * 1000000 * 60
  }
}   

#Filter out all samples
Samp2 = datUrchin %>%
  filter(CS == "S") 
Samp2 = Samp2[complete.cases(Samp2[,7]),]
Samp2 = Samp2[,c(1,4,26:40)]
dat_sum1 <- merge(Samp1, Samp2, by=c("TIMESTEP", "SAMPLE"))

# Calc mass specific metabolism
dat_sum1$resp_corr_2 = dat_sum1$resp_umol_min/dat_sum1$WW..g.

# Plots wt by mass specific metabolic rate
ggplot(dat_sum1, aes(x = WW..g., y = resp_corr_2)) + geom_point() + theme_bw() +
  labs(y = 'mass specific respiration rate (umol/hr/g)', x = 'wet weight (g)') 

# Log transform, remember will be weird if neg
dat_sum1$log_resp = log(dat_sum1$resp_umol_min)
dat_sum1$log_massWW = log(dat_sum1$WW..g.)

# Log-log plot, should be increasing
ggplot(dat_sum1, aes(x = log_massWW, y = log_resp))+
  geom_point()+
  geom_smooth(method = 'lm') + 
  theme_bw()+
  labs(x = 'log(mass)', y = 'log(respiration rate)')

# Finds the slope of the log-log plot which will be used to linearize data as per metabolic theory
lmer_resp_corr_2 = lm(log_resp ~ log_massWW*TRT, data = dat_sum1)
summary(lmer_resp_corr_2)
coef(lmer_resp_corr_2)
coef1 = 0.841916 

mean_mass = mean(dat_sum1$WW..g., na.rm=TRUE)
  
# Transforming data to remove the effect of size on metabolic rate
dat_sum1$resp_corr_2 = (dat_sum1$resp_umol_min/dat_sum1$WW..g.) * 
             (dat_sum1$WW..g./mean_mass)^(1-coef1)
  
#Linearized plot of weight versus respiration rate
  ggplot(dat_sum1, aes(x = WW..g., y = resp_corr_2))+geom_point()+geom_smooth(method = 'lm') +
    theme_bw()+
    labs(x = 'Wet weight (g)', y = 'respiration rate (corr), umol/hr')

######################################################################################################
######################################################################################################
##Now moving on to Tegula!
######################################################################################################
######################################################################################################
  
#Remove all the data that are not necessary
rm(list = setdiff(ls(), c("dat_sum1", "datUrchin", "coef1","Svial_vol_ml")))

# Volume of vials used, all small for tegula
Svial_vol_ml = 99.957

#Import snail data now
datSnail = read.csv('RespAllTP_QC.csv')

#hm changes it into "time object"
# seconds convers into second since midnight
# as.integer just makes it an integer
TIME8 = as.integer(seconds(hm(datSnail$TIME1)))

#Need to alter data structure, so first create blank dataframe
dat_long = data.frame(Vial = NA,Critter.ID=NA, Time = NA, Oxygen = NA, Temp = NA, Run = NA, TRT = NA, CS = NA)  #Initialize empty dataframe

for (b in seq(from=7, to=25, by=3)){ #rbind the proper columns, so column 4 (ID) and the others. 
  # This takes the time from the column, and changes it to "seconds since capped"
  datSnail[,b] = as.integer(seconds(hm(datSnail[,b])))-TIME8
  temp = datSnail[,c(2,4,((b):(b+2)),1, 40, 41)] # This is how to do 4,5,6,7,8 -- 4,9,10,11,12
  # This adds a "vial" column - which you can ignore if you already have one. 
  colnames(temp) = colnames(dat_long) #column names for rbind
  dat_long = rbind(dat_long, temp)
}

# gets rid of any data without a corresponding oxygen value
dat_long = dat_long[complete.cases(dat_long$Oxygen),]

# Now going onto LoLinR: 
dat_long$combined <- paste(dat_long$Critter.ID,dat_long$Run)
dat_long$Oxygen <- as.numeric(dat_long$Oxygen)
dat_long$Temp <- as.numeric(dat_long$Temp)
u <-data.frame(unique(dat_long$combined))
f = data.frame(Slope = NA,Sig=NA, Temp = NA, SAMPLE = NA, TIMESTEP = NA)  #Initialize empty dataframe

for (i in 1:nrow(u)) { 
  x <- dplyr::select(filter(dat_long, combined == u[i,]),c(Vial, Critter.ID, Time, Oxygen, Temp, Run))
  mean_temp = mean(x$Temp, na.rm = TRUE)
  regression = rankLocReg(xall = x$Time, yall = x$Oxygen, alpha = .7)
  summary(regression)
  a = regression$allRegs
  f[i,] <- c(a[1,]$b1 * -1,a[1,]$b1LoCI > 0 | a[1,]$b1UpCI<0, mean_temp, x$Critter.ID[1], x$Run[1])
}

dat_sum2 <- merge(datSnail[,1:40],f, by=c("SAMPLE","TIMESTEP"))
dat_sum2 <- dat_sum2[-c(0,9:39)]

dat_sum2$resp_mg_sec <- NA
dat_sum2$resp_g_sec <- NA
dat_sum2$resp_mol_sec <- NA
dat_sum2$resp_umol_min <- NA

#Make slope numeric, chr for some reason...
dat_sum2$Slope = as.numeric(dat_sum2$Slope)

#Filter for controls only
Controls2 = dat_sum2 %>%
  filter(CS == "C") 

#Calc mean of controls
Con_Sum2 <- summarySE(data=Controls2, measurevar = "Slope",
                     groupvars = c("TIMESTEP","RUN"), na.rm = TRUE)
Con_Sum3 <- summarySE(data=Controls2, measurevar = "Slope",
                      groupvars = c("TIMESTEP"), na.rm = TRUE)

Samp = dat_sum2 %>%
  filter(CS == "S")

Samp <- merge(Con_Sum2, Samp, by.x = c("TIMESTEP", "RUN"), by.y = c("TIMESTEP", "RUN"))

#Correct slope of samples by subracting the ave slope of controls in a run
Samp$Slope.y = Samp$Slope.y-Samp$Slope.x

# Slopes are in mg/L per second. convert to something else: ####
# firstly multiply by volume of vial in liters to convert to mg per second: 
for (i in 1:nrow(Samp)) {
  if(Samp$VTYPE[i]=='S') {
    Samp$resp_mg_sec[i] = Samp$Slope.y[i] * (Svial_vol_ml/1000)
    Samp$resp_g_sec[i] = Samp$resp_mg_sec[i] / 1000
    Samp$resp_mol_sec[i] = Samp$resp_g_sec[i] / 31.998
    Samp$resp_umol_min[i] = Samp$resp_mol_sec[i] * 1000000 * 60
  }else{
    Samp$resp_mg_sec[i] = Samp$Slope.y[i] * (Lvial_vol_ml/1000)
    Samp$resp_g_sec[i] = Samp$resp_mg_sec[i] / 1000
    Samp$resp_mol_sec[i] = Samp$resp_g_sec[i] / 31.998
    Samp$resp_umol_min[i] = Samp$resp_mol_sec[i] * 1000000 * 60
  }
}   

Samp2 = datSnail %>%
  filter(CS == "S") 
Samp2 = Samp2[complete.cases(Samp2[,7]),]
Samp2 = Samp2[,c(1,4,26:40)]
dat_sum2 <- merge(Samp, Samp2, by=c("TIMESTEP", "SAMPLE"))

names(dat_sum2)[names(dat_sum2) == "TRT.x"] <- "TRT"
dat_sum2I <- subset(dat_sum2, dat_sum2$TIMESTEP ==1)

#Correct respirometry data by weight measurements
dat_sum2$resp_corr_2 = dat_sum2$resp_umol_min/dat_sum2$WW..g.

ggplot(dat_sum2, aes(x = WW..g., y = resp_corr_2)) + geom_point() + theme_bw() +
  labs(y = 'mass specific respiration rate (umol/hr/g)', x = 'wet weight (g)') 

#Remove resp rates <0
dat_sum2_0 <- subset(dat_sum2, dat_sum2$resp_corr_2 > 0)

#Log transforming resp rate and mass
dat_sum2_0$log_resp = log(dat_sum2_0$resp_umol_min)
dat_sum2_0$log_massWW = log(dat_sum2_0$WW..g.) 

ggplot(dat_sum2_0, aes(x = log_massWW, y = log_resp))+
  geom_point()+
  geom_smooth(method = 'lm') + 
  theme_bw()+
  labs(x = 'log(mass)', y = 'log(respiration rate)')

#Calc slope of regression of log mass and log resp to get coefficient to correct subsequent
#respiration rate data to remove mass
lmer_resp_corr_2 = lm(log_resp ~ log_massWW*TRT, data = dat_sum2_0)
summary(lmer_resp_corr_2)
coef(lmer_resp_corr_2)
coef2 = 0.64133 # using wet weight

mean_mass = mean(dat_sum2_0$WW..g.)

#Doing mass correction
dat_sum2_0$resp_corr_2 = (dat_sum2_0$resp_umol_min/dat_sum2_0$WW..g.) * 
          (dat_sum2_0$WW..g./mean_mass)^(1-coef2)

rm(list = setdiff(ls(), c("dat_sum1", "dat_sum2","dat_sum2_0","datUrchin","datSnail","coef1","coef2")))

######################################################################################################
######################################################################################################
### On to water chemistry
######################################################################################################
######################################################################################################

chem = read.csv('DiscreteSamples_All_Chronic.csv') # Read in discrete sample data file, already passed discrete samples through CO2SYS
chem$DT <- mdy_hm(paste(chem$DATE, chem$TIME)) # Creating R date/time stamp
chem$DT <- force_tz(as.POSIXct(chem$DT, origin = as.POSIXct("1970-01-01", TZ = "America/Los_Angeles"), 
        TZ = "America/Los_Angeles"), tzone = "America/Los_Angeles") #Super weird, but couldn't get it to not be in UTC, so had to do this
chem$Tco2 <- chem$HCO3.out..mmol.kgSW. + chem$CO3.out..mmol.kgSW. + chem$CO2.out..mmol.kgSW. # Calculating total CO2

# Summary stats, changed the measurevar to get all values for TABLE 1
pHave <- summarySE(data=chem, measurevar = "Tco2",
                    groupvars = c("TRT","TIMEPOINT"))
pH2ave <- summarySE(data=pHave, measurevar = "Tco2",
                    groupvars = c("TRT"))

DO = read.csv('DO_processed_Chronic.csv') # Read in Vernier probes dissolved oxygen data (FROM HEADER BUCKETS)
# Create R date/time
DO$DT <- as.POSIXct((DO$SDN - 719529)*86400, origin = "1970-01-01", TZ = "America/Los_Angeles")
DO$DT <- format(DO$DT, format = '%Y-%m-%d %H:%M')
pH = read.csv('pH_processed_Chronic.csv') # Read in Durafet probe pH data (FROM HEADER BUCKETS)
pH$DT <- force_tz(as.POSIXct((pH$SDN - 719529)*86400, origin = "1970-01-01", TZ = "America/Los_Angeles"), tzone = "America/Los_Angeles")
pH$DT <- format(pH$DT, format = '%Y-%m-%d %H:%M')

# Taking pH data and converting to long format
pH_long <- data.frame(rbind(cbind(rep('H1',nrow(pH)), pH$DT, pH$pH_1, pH$TC_1), cbind(rep('H2',nrow(pH)), pH$DT, pH$pH_2, pH$TC_2),
                 cbind(rep('H3',nrow(pH)), pH$DT, pH$pH_3, pH$TC_3), cbind(rep('H4',nrow(pH)), pH$DT, pH$pH_4, pH$TC_4),
                 cbind(rep('H5',nrow(pH)), pH$DT, pH$pH_5, pH$TC_5), cbind(rep('H6',nrow(pH)), pH$DT, pH$pH_6, pH$TC_6)))
pH_long <- rename(pH_long, c("H" = "X1", "Date" = "X2", "pH" = "X3", "Temp" = "X4")) #Renmae columns
pH_long$Date <- as.POSIXct(pH_long$Date, origin = "1970-01-01", TZ = "America/Los_Angeles")
options(digits = 4)
pH_long$pH <- as.numeric(pH_long$pH) # Converting to numeric from factor is odd...
pH_long$Temp <- as.numeric(pH_long$Temp)

# Taking DO data and converting to long format
DO_long <- data.frame(rbind(cbind(rep('H1',nrow(DO)), DO$SDN, DO$mgL_2, DO$TC_2), cbind(rep('H2',nrow(DO)), DO$SDN, DO$mgL_3, DO$TC_3),
                            cbind(rep('H3',nrow(DO)), DO$SDN, DO$mgL_4, DO$TC_4), cbind(rep('H4',nrow(DO)), DO$SDN, DO$mgL_5, DO$TC_5),
                            cbind(rep('H5',nrow(DO)), DO$SDN, DO$mgL_6, DO$TC_6), cbind(rep('H6',nrow(DO)), DO$SDN, DO$mgL_7, DO$TC_7)))
               
DO_long <- rename(DO_long, c("H" = "X1", "Date" = "X2", "DO" = "X3", "Temp" = "X4")) #Renmae columns
DO_long$Date <- as.POSIXct((DO$SDN - 719529)*86400, origin = "1970-01-01", TZ = "America/Los_Angeles")
DO_long$DO <- as.numeric(DO_long$DO) # Converting to numeric from factor is odd...
DO_long$Temp <- as.numeric(DO_long$Temp)

# Now we'll calibrate the durafets by calculating an offset from discrete samples throughout our experiment
Headers = chem %>%    # Subset only header samples 
    filter(SH == "H") 
Headers$ID <- factor(Headers$ID)    #Need to reset the levels of the factors, for looping
u <- data.frame(unique(Headers$TRT))    #This will be used to loop through all headers
u <- rename(u, c("ID" = "unique.Headers.TRT."))
u$ID <- factor(u$ID)
offset <- data.frame(H1 = NA, H2 = NA, H3 = NA, H4 = NA, H5 = NA, H6 = NA)
for (i in 1:nrow(u)) {
  temp = Headers %>%
    filter(ID == u$ID[i])
    for (j in 1:nrow(temp)) {
      temp2 = pH_long %>%
        filter(H == u$ID[i])
      ind <- which.min(abs(temp$DT[j]-temp2$Date))
      offset[j, i] <- temp2$pH[ind] - temp$pH.out[j]
}
}   
offset <- cbind(offset, unique(Headers$DATE)) # combine dates to time series of discrete samples
offset <-rename(offset, c("Date" = "unique(Headers$DATE)"))

# Need to restructure to plot
offsetP <- data.frame(rbind(cbind(rep('H1',nrow(offset)), offset$H1, offset$Date), cbind(rep('H2',nrow(offset)), offset$H2, offset$Date),
                            cbind(rep('H3',nrow(offset)), offset$H3, offset$Date), cbind(rep('H4',nrow(offset)), offset$H4, offset$Date),
                            cbind(rep('H5',nrow(offset)), offset$H5, offset$Date), cbind(rep('H6',nrow(offset)), offset$H6, offset$Date)))
offsetP$X2 <- as.numeric(offsetP$X2)

# Plot offsets over time for all headers
i <- ggplot(offsetP, aes(x = factor(X3), y = X2, color = X1)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~offsetP$X1) +
  scale_color_manual(values = wes_palette("Zissou1", 6, "continuous")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="Durafet Calibration",y="Offset", x = "Time Point")
i

# Calc ave offset for each durafet probe
offsetM <- summarySE(data=offsetP, measurevar = "X2",
                    groupvars = "X1")
for(i in 1:nrow(pH_long)) {
  if(pH_long$H[i] =='H1') {
    pH_long$pH[i] = pH_long$pH[i] - offsetM$X2[1]
  }else if (pH_long$H[i] =='H2'){
    pH_long$pH[i] = pH_long$pH[i] - offsetM$X2[2]
  }else if (pH_long$H[i] =='H3'){
    pH_long$pH[i] = pH_long$pH[i] - offsetM$X2[3]
  }else if (pH_long$H[i] =='H4'){
    pH_long$pH[i] = pH_long$pH[i] - offsetM$X2[4]
  }else if (pH_long$H[i] =='H5'){
    pH_long$pH[i] = pH_long$pH[i] - offsetM$X2[5]
  }else if (pH_long$H[i] =='H6'){
    pH_long$pH[i] = pH_long$pH[i] - offsetM$X2[6]
  }else {
    pH_long$pH[i] = pH_long$pH[i]
  }
}   

# Now add the source water durafets, we didn't take discrete samples of these since they aren't in the model, 
# but helpful to get an idea of what the mixing was like
pH_mix <- data.frame(rbind(cbind(rep('H7',nrow(pH)), pH$DT, pH$pH_7, pH$TC_7), 
                           cbind(rep('H8',nrow(pH)), pH$DT, pH$pH_8, pH$TC_8)))
pH_mix <- rename(pH_mix, c("H" = "X1", "Date" = "X2", "pH" = "X3", "Temp" = "X4"))
pH_mix$Date <- as.POSIXct(pH_mix$Date, origin = "1970-01-01", TZ = "America/Los_Angeles")
pH_mix$pH <- as.numeric(pH_mix$pH) # Converting to numeric from factor is odd...
pH_mix$Temp <- as.numeric(pH_mix$Temp)

pHcal <- rbind(pH_long, pH_mix)
pHcal <- subset(pHcal, pHcal[ , 3] > 7.2) #Remove values under 7.2 since these are not good data
pHcal <- subset(pHcal, !is.na(pHcal[ , 3])) #Remove NANs

#Remove time period where system was being backflushed and sensor were briefly out of the water
pHcal <- subset(pHcal, pHcal$Date >= as.POSIXct("2019-11-10 12:00:00", origin = "1970-01-01", tz = "America/Los_Angeles") 
              | pHcal$Date <= as.POSIXct("2019-11-10 10:30:00", origin = "1970-01-01", tz = "America/Los_Angeles"))
pHcal$H <- factor(pHcal$H)

# Bring in YSI measurements to plot points over durafet time series
YSI = read.csv('YSI_Fall2019.csv')
# Convert time
YSI$DT <- force_tz(as.POSIXct(with(YSI, mdy(Date) + hm(Time)), origin = as.POSIXct("1970-01-01", 
           TZ = "America/Los_Angeles")), tzone = "America/Los_Angeles")
pHcalsub <- pHcal[pHcal$H == "H1"|pHcal$H == "H2"|pHcal$H == "H3"
                  |pHcal$H == "H4"|pHcal$H == "H5"|pHcal$H == "H6",]

DO_longsub <- DO_long[DO_long$H == "H1"|DO_long$H == "H2"|DO_long$H == "H3"
                  |DO_long$H == "H4"|DO_long$H == "H5"|DO_long$H == "H6",]

######################################################################################################
######################################################################################################
## Let's do some PCA regressions!
######################################################################################################
######################################################################################################

## Now calc PCA on average variables
pHcalsub <- pHcal[pHcal$H %in% c("H1","H2","H3","H4","H5","H6"),]
#Calculate summary stats for pH, temp and DO to do PCA on
pHave <- summarySE(data=pHcalsub, measurevar = "pH", 
                    groupvars = "H")
pHTCave <- summarySE(data=pHcalsub, measurevar = "Temp",
                    groupvars = "H")
DO_long <- na.omit(DO_long)
DOave <- summarySE(data=DO_long, measurevar = "DO",
                   groupvars = "H")
#Combine variable for PCA regression
PCAvar <- data.frame(cbind(pHave$H,pHave$pH,pHave$sd,pHTCave$Temp,pHTCave$sd,DOave$DO,DOave$sd))
names(PCAvar)[1:7] <- c("Headers","pH","pH SD","Temperature","Temperature SD","DO","DO SD") # Change column names
PCAvar$Headers <- c("H1", "H2", "H3", "H4", "H5", "H6")

TS.PCA <- prcomp(PCAvar[2:7], scale=TRUE) #Run PCA of pH, DO, tempDO, temppH
summary(TS.PCA) #Look at the amount of variance explained
bi1 <- biplot(TS.PCA, col = wes_palette("Zissou1", 2))
bi1

Load <- data.frame(TS.PCA$x) #Turn loadings into dataframe
TrtPCA <- data.frame(c(PCAvar, Load)) #Combine PCA loadings with OG data

#Making preliminary plots of data across treaments based on "PCs"
PC1Reg <- summarySE(data=TrtPCA, measurevar = "PC1",
                    groupvars = "Headers")
PC2Reg <- summarySE(data=TrtPCA, measurevar = "PC2",
                    groupvars = "Headers")
PC3Reg <- summarySE(data=TrtPCA, measurevar = "PC3",
                    groupvars = "Headers")
PC4Reg <- summarySE(data=TrtPCA, measurevar = "PC4",
                    groupvars = "Headers")

# Adding average PC1 score to urchin and snail data for plotting
for(i in 1:nrow(pHcalsub)) {
  if(pHcalsub$H[i] =="H1") {
    pHcalsub$PC1[i] = PC1Reg$PC1[1]
  }else if (pHcalsub$H[i] =="H2"){
    pHcalsub$PC1[i] = PC1Reg$PC1[2]
  }else if (pHcalsub$H[i] =="H3"){
    pHcalsub$PC1[i] = PC1Reg$PC1[3]
  }else if (pHcalsub$H[i] =="H4"){
    pHcalsub$PC1[i] = PC1Reg$PC1[4]
  }else if (pHcalsub$H[i] =="H5"){
    pHcalsub$PC1[i] = PC1Reg$PC1[5]
  }else if (pHcalsub$H[i] =="H6"){
    pHcalsub$PC1[i] = PC1Reg$PC1[6]
  }else {
    pHcalsub$PC1[i] = NA
  }
}   

for(i in 1:nrow(YSI)) {
  if(YSI$H[i] =="H1") {
    YSI$PC1[i] = PC1Reg$PC1[1]
  }else if (YSI$H[i] =="H2"){
    YSI$PC1[i] = PC1Reg$PC1[2]
  }else if (YSI$H[i] =="H3"){
    YSI$PC1[i] = PC1Reg$PC1[3]
  }else if (YSI$H[i] =="H4"){
    YSI$PC1[i] = PC1Reg$PC1[4]
  }else if (YSI$H[i] =="H5"){
    YSI$PC1[i] = PC1Reg$PC1[5]
  }else if (YSI$H[i] =="H6"){
    YSI$PC1[i] = PC1Reg$PC1[6]
  }else {
    YSI$PC1[i] = NA
  }
}   

pal = wes_palette("Zissou1",6,type = "continuous")

# Some plots of time series data
i <- ggplot(pHcalsub, aes(x = Date, y = pH), group = factor(H)) +
  geom_line(aes(x = Date, y = pH, color = factor(H)), alpha = 1/3) +
  theme_classic() +
  scale_color_manual(values = rev(wes_palette("Zissou1", 6, "continuous"))) +
  geom_point(data = YSI, mapping = 
               aes(x = DT, y = pH, color = factor(H)), size = 1) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "2 week", expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "right", legend.key.width = unit(0.5,'cm')) +
  ylim(7.4, 7.9) +
  labs(y="pH") + 
  labs(color="Header")
i <- i + guides(color = guide_legend(nrow=1)) + scale_y_continuous(breaks = c(7.4,7.6,7.8))
i

j <- ggplot(pHcalsub, aes(x = Date, y = Temp, color = factor(H)), group = factor(H)) +
  geom_line(alpha = 1/2.5) +
  theme_classic() +
  scale_color_manual(values = rev(wes_palette("Zissou1", 6, "continuous"))) +
  geom_point(data = YSI, mapping = 
               aes(x = DT, y = Temp, color = factor(H)), size = 1) +
  theme(legend.title = element_blank(), legend.position = "right", legend.key.width = unit(0.5,'cm')) +
  scale_x_datetime(date_breaks = "2 week", expand = c(0, 0)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(9, 15) +
  labs(color="Header") +
  labs(y="Temperature (\u00B0C)")
j <- j + guides(color = guide_legend(nrow=1))
j

k <- ggplot(DO_longsub, aes(x = Date, y = DO, color = factor(H)), group = factor(H)) +
  geom_line(alpha = 1/2.5) +
  theme_classic() +
  scale_color_manual(values = rev(wes_palette("Zissou1", 6, "continuous"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="Dissolved Oxygen (mg/L)", x="Date") +
  scale_x_datetime(date_breaks = "2 week", date_labels = "%m/%d",
                   expand = c(0, 0)) +
  labs(color="Header") +
  geom_point(data = YSI, mapping = 
               aes(x = DT, y = DO1, color = factor(H)), size = 1) +
  theme(legend.title = element_blank(), legend.position = "right", legend.key.width = unit(0.5,'cm')) +
  ylim(3.5, 10) 
k <- k + guides(color = guide_legend(nrow=1))  + scale_y_continuous(breaks = c(3,6,9))

k

l <- ggarrange(i, j, k, labels = c("(a)", "(b)", "(c)"), common.legend = TRUE, legend = "bottom", ncol = 1, nrow = 3)
l

#Create FIGURE S2
ggsave(plot = l, file = "FigS2.png", 
       type = "cairo-png",  bg = "white",
       width = 15, height = 25, units = "cm", dpi = 300)



# Adding average PC1 score to urchin and snail data for plotting
for(i in 1:nrow(dat_sum1)) {
  if(dat_sum1$TRT[i] ==1) {
    dat_sum1$PC1[i] = PC1Reg$PC1[1]
  }else if (dat_sum1$TRT[i] ==2){
    dat_sum1$PC1[i] = PC1Reg$PC1[2]
  }else if (dat_sum1$TRT[i] ==3){
    dat_sum1$PC1[i] = PC1Reg$PC1[3]
  }else if (dat_sum1$TRT[i] ==4){
    dat_sum1$PC1[i] = PC1Reg$PC1[4]
  }else if (dat_sum1$TRT[i] ==5){
    dat_sum1$PC1[i] = PC1Reg$PC1[5]
  }else if (dat_sum1$TRT[i] ==6){
    dat_sum1$PC1[i] = PC1Reg$PC1[6]
  }else {
    dat_sum1$PC1[i] = NA
  }
}   

for(i in 1:nrow(dat_sum2_0)) {
  if(dat_sum2_0$TRT[i] ==1) {
    dat_sum2_0$PC1[i] = PC1Reg$PC1[1]
  }else if (dat_sum2_0$TRT[i] ==2){
    dat_sum2_0$PC1[i] = PC1Reg$PC1[2]
  }else if (dat_sum2_0$TRT[i] ==3){
    dat_sum2_0$PC1[i] = PC1Reg$PC1[3]
  }else if (dat_sum2_0$TRT[i] ==4){
    dat_sum2_0$PC1[i] = PC1Reg$PC1[4]
  }else if (dat_sum2_0$TRT[i] ==5){
    dat_sum2_0$PC1[i] = PC1Reg$PC1[5]
  }else if (dat_sum2_0$TRT[i] ==6){
    dat_sum2_0$PC1[i] = PC1Reg$PC1[6]
  }else {
    dat_sum2_0$PC1[i] = NA
  }
}  

for(i in 1:nrow(dat_sum1)) {
  if(dat_sum1$TRT[i] ==1) {
    dat_sum1$PC2[i] = PC2Reg$PC2[1]
  }else if (dat_sum1$TRT[i] ==2){
    dat_sum1$PC2[i] = PC2Reg$PC2[2]
  }else if (dat_sum1$TRT[i] ==3){
    dat_sum1$PC2[i] = PC2Reg$PC2[3]
  }else if (dat_sum1$TRT[i] ==4){
    dat_sum1$PC2[i] = PC2Reg$PC2[4]
  }else if (dat_sum1$TRT[i] ==5){
    dat_sum1$PC2[i] = PC2Reg$PC2[5]
  }else if (dat_sum1$TRT[i] ==6){
    dat_sum1$PC2[i] = PC2Reg$PC2[6]
  }else {
    dat_sum1$PC2[i] = NA
  }
}   

for(i in 1:nrow(dat_sum2_0)) {
  if(dat_sum2_0$TRT[i] ==1) {
    dat_sum2_0$PC2[i] = PC2Reg$PC2[1]
  }else if (dat_sum2_0$TRT[i] ==2){
    dat_sum2_0$PC2[i] = PC2Reg$PC2[2]
  }else if (dat_sum2_0$TRT[i] ==3){
    dat_sum2_0$PC2[i] = PC2Reg$PC2[3]
  }else if (dat_sum2_0$TRT[i] ==4){
    dat_sum2_0$PC2[i] = PC2Reg$PC2[4]
  }else if (dat_sum2_0$TRT[i] ==5){
    dat_sum2_0$PC2[i] = PC2Reg$PC2[5]
  }else if (dat_sum2_0$TRT[i] ==6){
    dat_sum2_0$PC2[i] = PC2Reg$PC2[6]
  }else {
    dat_sum2_0$PC2[i] = NA
  }
} 


for(i in 1:nrow(dat_sum1)) {
  if(dat_sum1$TRT[i] ==1) {
    dat_sum1$PC3[i] = PC3Reg$PC3[1]
  }else if (dat_sum1$TRT[i] ==2){
    dat_sum1$PC3[i] = PC3Reg$PC3[2]
  }else if (dat_sum1$TRT[i] ==3){
    dat_sum1$PC3[i] = PC3Reg$PC3[3]
  }else if (dat_sum1$TRT[i] ==4){
    dat_sum1$PC3[i] = PC3Reg$PC3[4]
  }else if (dat_sum1$TRT[i] ==5){
    dat_sum1$PC3[i] = PC3Reg$PC3[5]
  }else if (dat_sum1$TRT[i] ==6){
    dat_sum1$PC3[i] = PC3Reg$PC3[6]
  }else {
    dat_sum1$PC3[i] = NA
  }
}   

for(i in 1:nrow(dat_sum2_0)) {
  if(dat_sum2_0$TRT[i] ==1) {
    dat_sum2_0$PC3[i] = PC3Reg$PC3[1]
  }else if (dat_sum2_0$TRT[i] ==2){
    dat_sum2_0$PC3[i] = PC3Reg$PC3[2]
  }else if (dat_sum2_0$TRT[i] ==3){
    dat_sum2_0$PC3[i] = PC3Reg$PC3[3]
  }else if (dat_sum2_0$TRT[i] ==4){
    dat_sum2_0$PC3[i] = PC3Reg$PC3[4]
  }else if (dat_sum2_0$TRT[i] ==5){
    dat_sum2_0$PC3[i] = PC3Reg$PC3[5]
  }else if (dat_sum2_0$TRT[i] ==6){
    dat_sum2_0$PC3[i] = PC3Reg$PC3[6]
  }else {
    dat_sum2_0$PC3[i] = NA
  }
} 



for(i in 1:nrow(dat_sum1)) {
  if(dat_sum1$TRT[i] ==1) {
    dat_sum1$PC4[i] = PC4Reg$PC4[1]
  }else if (dat_sum1$TRT[i] ==2){
    dat_sum1$PC4[i] = PC4Reg$PC4[2]
  }else if (dat_sum1$TRT[i] ==3){
    dat_sum1$PC4[i] = PC4Reg$PC4[3]
  }else if (dat_sum1$TRT[i] ==4){
    dat_sum1$PC4[i] = PC4Reg$PC4[4]
  }else if (dat_sum1$TRT[i] ==5){
    dat_sum1$PC4[i] = PC4Reg$PC4[5]
  }else if (dat_sum1$TRT[i] ==6){
    dat_sum1$PC4[i] = PC4Reg$PC4[6]
  }else {
    dat_sum1$PC4[i] = NA
  }
}   

for(i in 1:nrow(dat_sum2_0)) {
  if(dat_sum2_0$TRT[i] ==1) {
    dat_sum2_0$PC4[i] = PC4Reg$PC4[1]
  }else if (dat_sum2_0$TRT[i] ==2){
    dat_sum2_0$PC4[i] = PC4Reg$PC4[2]
  }else if (dat_sum2_0$TRT[i] ==3){
    dat_sum2_0$PC4[i] = PC4Reg$PC4[3]
  }else if (dat_sum2_0$TRT[i] ==4){
    dat_sum2_0$PC4[i] = PC4Reg$PC4[4]
  }else if (dat_sum2_0$TRT[i] ==5){
    dat_sum2_0$PC4[i] = PC4Reg$PC4[5]
  }else if (dat_sum2_0$TRT[i] ==6){
    dat_sum2_0$PC4[i] = PC4Reg$PC4[6]
  }else {
    dat_sum2_0$PC4[i] = NA
  }
} 


#Mass correct grazing rate data
dat_sum1$bin <- rep(NA, nrow(dat_sum1))
for(x in c('A', 'B', 'C')) dat_sum1$bin[grep(x, dat_sum1$SAMPLE)] <- x
dat_sum1$graze_corr <- dat_sum1$GrazingCorr_gWWperday/dat_sum1$WW..g.

dat_sum2$bin <- rep(NA, nrow(dat_sum2))
for(x in c('A', 'B', 'C')) dat_sum2$bin[grep(x, dat_sum2$SAMPLE)] <- x
dat_sum2$graze_corr <- dat_sum2$GrazingCorr_gWWperday/dat_sum2$WW..g.

dat_sum2_0$bin <- rep(NA, nrow(dat_sum2_0))
for(x in c('A', 'B', 'C')) dat_sum2_0$bin[grep(x, dat_sum2_0$SAMPLE)] <- x
dat_sum2_0$graze_corr <- dat_sum2_0$GrazingCorr_gWWperday/dat_sum2_0$WW..g.


######################################################################################################
######################################################################################################
#Now on to short experiment
######################################################################################################
######################################################################################################

Svial_vol_ml = 99.957
Lvial_vol_ml = 229.7059

# This is going to read in one of the respirometry sheets: 
datUrchinS = read.csv('RespAllSF_Short.csv')

#hm changes it into "time object"
# seconds convers into second since midnight
# as.integer just makes it an integer
TIME8 = as.integer(seconds(hm(datUrchinS$TIME1)))

# The respirometry data was entered just how it was taken - with multiple columns of oxygen. I want this in 
#"long form" - that is, with 1 column for time, another for oxygen, another for temp. 
# if this is confusing, compare "dat_long" I make with the following code and hopefully this will make sense! 
# you could just as easily do this in excel but I wanted to do it in code as I had quite a bit of data! 
dat_long = data.frame(Vial = NA,Critter.ID=NA, Time = NA, Oxygen = NA, Temp = NA, Run = NA, TRT = NA, CS = NA)  #Initialize empty dataframe

for (b in seq(from=7, to=25, by=3)){ #rbind the proper columns, so column 4 (ID) and the others. 
  # This takes the time from the column, and changes it to "seconds since capped"
  datUrchinS[,b] = as.integer(seconds(hm(datUrchinS[,b])))-TIME8
  temp = datUrchinS[,c(2,4,((b):(b+2)),1, 40, 41)] # This is how to do 4,5,6,7,8 -- 4,9,10,11,12
  # This adds a "vial" column - which you can ignore if you already have one. 
  colnames(temp) = colnames(dat_long) #column names for rbind
  dat_long = rbind(dat_long, temp)
}

# gets rid of any data without a corresponding oxygen value
dat_long = dat_long[complete.cases(dat_long$Oxygen),]

# Now going onto LoLinR: 
dat_long$combined <- paste(dat_long$Critter.ID,dat_long$Run)
u <-data.frame(unique(dat_long$combined))
f = data.frame(Slope = NA,Sig=NA, Temp = NA, SAMPLE = NA, TIMESTEP = NA)  #Initialize empty dataframe

for (i in 1:nrow(u)) { 
  x <- dplyr::select(filter(dat_long, combined == u[i,]),c(Vial, Critter.ID, Time, Oxygen, Temp, Run))
  mean_temp = mean(x$Temp, na.rm = TRUE)
  regression = rankLocReg(xall = x$Time, yall = x$Oxygen, alpha = .7)
  summary(regression)
  a = regression$allRegs
  f[i,] <- c(a[1,]$b1 * -1,a[1,]$b1LoCI > 0 | a[1,]$b1UpCI<0, mean_temp, x$Critter.ID[1], x$Run[1])
}


#Add some identifying columns to the results of LoLin
dat_sum3 <- merge(datUrchinS[,1:6],f, by=c("SAMPLE","TIMESTEP"))

#Create columns for calculating out rates
dat_sum3$resp_mg_sec <- NA
dat_sum3$resp_g_sec <- NA
dat_sum3$resp_mol_sec <- NA
dat_sum3$resp_umol_min <- NA

#Make slope numeric, chr for some reason...
dat_sum3$Slope = as.numeric(dat_sum3$Slope)

#Filter for controls only
Controls1 = dat_sum3 %>%
  filter(CS == "C") 

#Calc mean of controls
Con_Sum1 <- summarySE(data=Controls1, measurevar = "Slope",
                      groupvars = c("TIMESTEP","RUN"), na.rm = TRUE)

#Subset only samples
Samp1 = dat_sum3 %>%
  filter(CS == "S")

#Adds run specific ave control to each sample
Samp1 <- merge(Con_Sum1, Samp1, by.x = c("TIMESTEP", "RUN"), by.y = c("TIMESTEP", "RUN"))

#Correct slope of samples by subracting the ave slope of controls in a run
Samp1$Slope.y = Samp1$Slope.y-Samp1$Slope.x

# Slopes are in mg/L per second. convert to something else: ####
# firstly multiply by volume of vial in liters to convert to mg per second: 
for (i in 1:nrow(Samp1)) {
  if(Samp1$VTYPE[i]=='S') {
    Samp1$resp_mg_sec[i] = Samp1$Slope.y[i] * (Svial_vol_ml/1000)
    Samp1$resp_g_sec[i] = Samp1$resp_mg_sec[i] / 1000
    Samp1$resp_mol_sec[i] = Samp1$resp_g_sec[i] / 31.998
    Samp1$resp_umol_min[i] = Samp1$resp_mol_sec[i] * 1000000 * 60
  }else{
    Samp1$resp_mg_sec[i] = Samp1$Slope.y[i] * (Lvial_vol_ml/1000)
    Samp1$resp_g_sec[i] = Samp1$resp_mg_sec[i] / 1000
    Samp1$resp_mol_sec[i] = Samp1$resp_g_sec[i] / 31.998
    Samp1$resp_umol_min[i] = Samp1$resp_mol_sec[i] * 1000000 * 60
  }
}   

#Filter out all samples
Samp2 = datUrchinS %>%
  filter(CS == "S") 
Samp2 = Samp2[complete.cases(Samp2[,7]),]
Samp2 = Samp2[,c(1,4,26:40)]
dat_sum3 <- merge(Samp1, Samp2, by=c("TIMESTEP", "SAMPLE"))

# Calc mass specific metabolism
dat_sum3$resp_corr = dat_sum3$resp_umol_min/dat_sum3$WW..g.

# Plots wt by mass specific metabolic rate
ggplot(dat_sum3, aes(x = WW..g., y = resp_corr)) + geom_point() + theme_bw() +
  labs(y = 'mass specific respiration rate (umol/hr/g)', x = 'wet weight (g)') 

# Log transform, remember will be weird if neg
dat_sum3$log_resp = log(dat_sum3$resp_umol_min)
dat_sum3$log_massWW = log(dat_sum3$WW..g.)

# Log-log plot, should be increasing
ggplot(dat_sum3, aes(x = log_massWW, y = log_resp))+
  geom_point()+
  geom_smooth(method = 'lm') + 
  theme_bw()+
  labs(x = 'log(mass)', y = 'log(respiration rate)')

# Finds the slope of the log-log plot which will be used to linearize data as per metabolic theory
lmer_resp_corr_2 = lm(log_resp ~ log_massWW*TRT, data = dat_sum3)
summary(lmer_resp_corr_2)
coef(lmer_resp_corr_2)
coef1 = 0.781437
#Using wet weight
#coef = 0.74524 #Umi's AK coef

mean_mass = mean(dat_sum3$WW..g.)

# Transforming data to remove the effect of size on metabolic rate
dat_sum3$resp_corr = (dat_sum3$resp_umol_min/dat_sum3$WW..g.) * 
  (dat_sum3$WW..g./mean_mass)^(1-coef1)

######################################################################################################
######################################################################################################
# On to snails...
######################################################################################################
######################################################################################################

Svial_vol_ml = 99.957
Lvial_vol_ml = 229.7059

# This is going to read in one of the respirometry sheets: 
datSnailS = read.csv('RespAllTP_Short_edit.csv')

#hm changes it into "time object"
# seconds convers into second since midnight
# as.integer just makes it an integer
TIME8 = as.integer(seconds(hm(datSnailS$TIME1)))

# The respirometry data was entered just how it was taken - with multiple columns of oxygen. I want this in 
#"long form" - that is, with 1 column for time, another for oxygen, another for temp. 
# if this is confusing, compare "dat_long" I make with the following code and hopefully this will make sense! 
# you could just as easily do this in excel but I wanted to do it in code as I had quite a bit of data! 
dat_long = data.frame(Vial = NA,Critter.ID=NA, Time = NA, Oxygen = NA, Temp = NA, Run = NA, TRT = NA, CS = NA)  #Initialize empty dataframe

for (b in seq(from=7, to=25, by=3)){ #rbind the proper columns, so column 4 (ID) and the others. 
  # This takes the time from the column, and changes it to "seconds since capped"
  datSnailS[,b] = as.integer(seconds(hm(datSnailS[,b])))-TIME8
  temp = datSnailS[,c(2,4,((b):(b+2)),1, 40, 41)] # This is how to do 4,5,6,7,8 -- 4,9,10,11,12
  # This adds a "vial" column - which you can ignore if you already have one. 
  colnames(temp) = colnames(dat_long) #column names for rbind
  dat_long = rbind(dat_long, temp)
}

# gets rid of any data without a corresponding oxygen value
dat_long = dat_long[complete.cases(dat_long$Oxygen),]

# Now going onto LoLinR: 
dat_long$combined <- paste(dat_long$Critter.ID,dat_long$Run)
u <-data.frame(unique(dat_long$combined))
f = data.frame(Slope = NA,Sig=NA, Temp = NA, SAMPLE = NA, TIMESTEP = NA)  #Initialize empty dataframe

for (i in 1:nrow(u)) { 
  x <- dplyr::select(filter(dat_long, combined == u[i,]),c(Vial, Critter.ID, Time, Oxygen, Temp, Run))
  mean_temp = mean(x$Temp, na.rm = TRUE)
  regression = rankLocReg(xall = x$Time, yall = x$Oxygen, alpha = .7)
  summary(regression)
  a = regression$allRegs
  f[i,] <- c(a[1,]$b1 * -1,a[1,]$b1LoCI > 0 | a[1,]$b1UpCI<0, mean_temp, x$Critter.ID[1], x$Run[1])
}


#Add some identifying columns to the results of LoLin
dat_sum4 <- merge(datSnailS[,1:6],f, by=c("SAMPLE","TIMESTEP"))

#Create columns for calculating out rates
dat_sum4$resp_mg_sec <- NA
dat_sum4$resp_g_sec <- NA
dat_sum4$resp_mol_sec <- NA
dat_sum4$resp_umol_min <- NA

#Make slope numeric, chr for some reason...
dat_sum4$Slope = as.numeric(dat_sum4$Slope)

#Filter for controls only
Controls1 = dat_sum4 %>%
  filter(CS == "C") 

#Calc mean of controls
Con_Sum1 <- summarySE(data=Controls1, measurevar = "Slope",
                      groupvars = c("TIMESTEP","RUN"), na.rm = TRUE)

#Subset only samples
Samp1 = dat_sum4 %>%
  filter(CS == "S")

#Adds run specific ave control to each sample
Samp1 <- merge(Con_Sum1, Samp1, by.x = c("TIMESTEP", "RUN"), by.y = c("TIMESTEP", "RUN"))

#Correct slope of samples by subracting the ave slope of controls in a run
Samp1$Slope.y = Samp1$Slope.y-Samp1$Slope.x

# Slopes are in mg/L per second. convert to something else: ####
# firstly multiply by volume of vial in liters to convert to mg per second: 
for (i in 1:nrow(Samp1)) {
  if(Samp1$VTYPE[i]=='S') {
    Samp1$resp_mg_sec[i] = Samp1$Slope.y[i] * (Svial_vol_ml/1000)
    Samp1$resp_g_sec[i] = Samp1$resp_mg_sec[i] / 1000
    Samp1$resp_mol_sec[i] = Samp1$resp_g_sec[i] / 31.998
    Samp1$resp_umol_min[i] = Samp1$resp_mol_sec[i] * 1000000 * 60
  }else{
    Samp1$resp_mg_sec[i] = Samp1$Slope.y[i] * (Lvial_vol_ml/1000)
    Samp1$resp_g_sec[i] = Samp1$resp_mg_sec[i] / 1000
    Samp1$resp_mol_sec[i] = Samp1$resp_g_sec[i] / 31.998
    Samp1$resp_umol_min[i] = Samp1$resp_mol_sec[i] * 1000000 * 60
  }
}   

#Filter out all samples
Samp2 = datSnailS %>%
  filter(CS == "S") 
Samp2 = Samp2[complete.cases(Samp2[,7]),]
Samp2 = Samp2[,c(1,4,26:40)]
dat_sum4 <- merge(Samp1, Samp2, by=c("TIMESTEP", "SAMPLE"))

# Calc mass specific metabolism
dat_sum4$resp_corr = dat_sum4$resp_umol_min/dat_sum4$WW..g.

# Plots wt by mass specific metabolic rate
ggplot(dat_sum4, aes(x = WW..g., y = resp_corr)) + geom_point() + theme_bw() +
  labs(y = 'mass specific respiration rate (umol/hr/g)', x = 'wet weight (g)') 

# Log transform, remember will be weird if neg

dat_sum4$log_resp = log(dat_sum4$resp_umol_min)
dat_sum4$log_massWW = log(dat_sum4$WW..g.)
dat_sum4 <- subset(dat_sum4, dat_sum4$resp_corr > 0)

# Log-log plot, should be increasing
ggplot(dat_sum4, aes(x = log_massWW, y = log_resp))+
  geom_point()+
  geom_smooth(method = 'lm') + 
  theme_bw()+
  labs(x = 'log(mass)', y = 'log(respiration rate)')

# Finds the slope of the log-log plot which will be used to linearize data as per metabolic theory
lmer_resp_corr_2 = lm(log_resp ~ log_massWW*TRT, data = dat_sum4)
summary(lmer_resp_corr_2)
coef(lmer_resp_corr_2)
coef1 = 0.65843 #Using wet weight

mean_mass = mean(dat_sum4$WW..g.)

# Transforming data to remove the effect of size on metabolic rate
dat_sum4$resp_corr = (dat_sum4$resp_umol_min/dat_sum4$WW..g.) * 
  (dat_sum4$WW..g./mean_mass)^(1-coef1)

ggplot(dat_sum4, aes(x = TRT, y = resp_corr, group = TRT))+geom_boxplot() + theme_bw() + 
  labs(title = 'Respiration Rate, corrected with log-log fit')






######################################################################################################
######################################################################################################
### On to water chemistry, should probably switch to processing DS data in R 
######################################################################################################
######################################################################################################

chem2 = read.csv('DiscreteSamples_All_Acute.csv')
chem2$DT <- mdy_hm(paste(chem2$DATE, chem2$TIME)) # Creating R date/time stamp
chem2$DT <- force_tz(as.POSIXct(chem2$DT, origin = as.POSIXct("1970-01-01", TZ = "America/Los_Angeles"), 
                               TZ = "America/Los_Angeles"), tzone = "America/Los_Angeles") #Super weird, but couldn't get it to not be in UTC, so had to do this
chem2$Tco2 <- chem2$HCO3.out..mmol.kgSW. + chem2$CO3.out..mmol.kgSW. + chem2$CO2.out..mmol.kgSW. # Calculating total CO2

# Summary stats TABLE 1; change measurevar for each carbon species
pHave <- summarySE(data=chem2, measurevar = "WAr.out",
                   groupvars = c("TRT","TIMEPOINT"))
pH2ave <- summarySE(data=pHave, measurevar = "WAr.out",
                    groupvars = c("TRT"))

DO2 = read.csv('DO_processed_Acute.csv')
# Create R date/time
DO2$DT <- as.POSIXct((DO2$SDN - 719529)*86400, origin = "1970-01-01", TZ = "America/Los_Angeles")
DO2$DT <- format(DO2$DT, format = '%Y-%m-%d %H:%M')
pH2 = read.csv('pH_processed_Acute.csv')
pH2$DT <- force_tz(as.POSIXct((pH2$SDN - 719529)*86400, origin = "1970-01-01", TZ = "America/Los_Angeles"), tzone = "America/Los_Angeles")
pH2$DT <- format(pH2$DT, format = '%Y-%m-%d %H:%M')

# Taking pH data and converting to long format
pH_long2 <- data.frame(rbind(cbind(rep('H1',nrow(pH2)), pH2$DT, pH2$pH_1, pH2$TC_1), cbind(rep('H2',nrow(pH2)), pH2$DT, pH2$pH_2, pH2$TC_2),
                            cbind(rep('H3',nrow(pH2)), pH2$DT, pH2$pH_3, pH2$TC_3), cbind(rep('H4',nrow(pH2)), pH2$DT, pH2$pH_4, pH2$TC_4),
                            cbind(rep('H5',nrow(pH2)), pH2$DT, pH2$pH_5, pH2$TC_5), cbind(rep('H6',nrow(pH2)), pH2$DT, pH2$pH_6, pH2$TC_6)))
pH_long2 <- pH_long2 %>%  #Rename columns
  rename("H" = "X1", "Date" = "X2", "pH" = "X3", "Temp" = "X4")
pH_long2$Date <- as.POSIXct(pH_long2$Date, origin = "1970-01-01", TZ = "America/Los_Angeles")
pH_long2$pH <- as.numeric(pH_long2$pH) # Converting to numeric from factor is odd...
pH_long2$Temp <- as.numeric(pH_long2$Temp)

DO_long2 <- data.frame(rbind(cbind(rep('H1',nrow(DO2)), DO2$SDN, DO2$mgL_2, DO2$TC_2), cbind(rep('H2',nrow(DO2)), DO2$SDN, DO2$mgL_3, DO2$TC_3),
                            cbind(rep('H3',nrow(DO2)), DO2$SDN, DO2$mgL_4, DO2$TC_4), cbind(rep('H4',nrow(DO2)), DO2$SDN, DO2$mgL_5, DO2$TC_5),
                            cbind(rep('H5',nrow(DO2)), DO2$SDN, DO2$mgL_6, DO2$TC_6), cbind(rep('H6',nrow(DO2)), DO2$SDN, DO2$mgL_7, DO2$TC_7)))

DO_long2 <- DO_long2 %>%  #Rename columns
  rename("H" = "X1", "Date" = "X2", "DO" = "X3", "Temp" = "X4")
DO_long2$Date <- as.POSIXct((DO2$SDN - 719529)*86400, origin = "1970-01-01", TZ = "America/Los_Angeles")
DO_long2$DO <- as.numeric(DO_long2$DO)# Converting to numeric from factor is odd...
DO_long2$Temp <- as.numeric(DO_long2$Temp)


# Now we'll calibrate the durafets by calculating an offset from discrete samples
# throughout our experiment
Headers2 = chem2 %>%    # Subset only header samples 
  filter(SH == "H") 
Headers2$ID <- factor(Headers2$ID)    #Need to reset the levels of the factors, for looping
u2 <- data.frame(unique(Headers2$TRT))    #This will be used to loop through all headers
u2 <- rename(u2, "ID" = "unique.Headers2.TRT.")
u2$ID <- factor(u2$ID)
offset2 <- data.frame(H1 = NA, H2 = NA, H3 = NA, H4 = NA, H5 = NA, H6 = NA)
for (i in 1:nrow(u2)) {
  temp2 = Headers2 %>%
    filter(ID == u2$ID[i])
  for (j in 1:nrow(temp2)) {
    temp3 = pH_long2 %>%
      filter(H == u2$ID[i])
    ind <- which.min(abs(temp2$DT[j]-temp3$Date))
    offset2[j, i] <- temp3$pH[ind] - temp2$pH.out[j]
  }
}   
offset2 <- cbind(offset2, unique(Headers2$DATE)) # combine dates to time series of discrete samples
offset2 <-rename(offset2, "Date" = "unique(Headers2$DATE)")
# Need to restructure to plot
offsetP2 <- data.frame(rbind(cbind(rep('H1',nrow(offset2)), offset2$H1, offset2$Date), cbind(rep('H2',nrow(offset2)), offset2$H2, offset2$Date),
                            cbind(rep('H3',nrow(offset2)), offset2$H3, offset2$Date), cbind(rep('H4',nrow(offset2)), offset2$H4, offset2$Date),
                            cbind(rep('H5',nrow(offset2)), offset2$H5, offset2$Date), cbind(rep('H6',nrow(offset2)), offset2$H6, offset2$Date)))
offsetP2$X2 <- as.numeric(offsetP2$X2)
# Plot offsets over time for all headers
i <- ggplot(offsetP2, aes(x = factor(X3), y = X2, color = X1)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~X1) +
  scale_color_manual(values = wes_palette("Zissou1", 6, "continuous")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="Durafet Calibration",y="Offset", x = "Time Point")
i


# Calc ave offset, for now, may want to remove outliers eventually...
offsetM2 <- summarySE(data=offsetP2, measurevar = "X2",
                     groupvars = "X1")
for(i in 1:nrow(pH_long2)) {
  if(pH_long2$H[i] =='H1') {
    pH_long2$pH[i] = pH_long2$pH[i] - offsetM2$X2[1]
  }else if (pH_long2$H[i] =='H2'){
    pH_long2$pH[i] = pH_long2$pH[i] - offsetM2$X2[2]
  }else if (pH_long2$H[i] =='H3'){
    pH_long2$pH[i] = pH_long2$pH[i] - offsetM2$X2[3]
  }else if (pH_long2$H[i] =='H4'){
    pH_long2$pH[i] = pH_long2$pH[i] - offsetM2$X2[4]
  }else if (pH_long2$H[i] =='H5'){
    pH_long2$pH[i] = pH_long2$pH[i] - offsetM2$X2[5]
  }else if (pH_long2$H[i] =='H6'){
    pH_long2$pH[i] = pH_long2$pH[i] - offsetM2$X2[6]
  }else {
    pH_long2$pH[i] = pH_long2$pH[i]
  }
}   
# Now add the source water durafets, we didn't take discrete samples of these
pH_mix2 <- data.frame(rbind(cbind(rep('H7',nrow(pH2)), pH2$DT, pH2$pH_7, pH2$TC_7), 
                           cbind(rep('H8',nrow(pH2)), pH2$DT, pH2$pH_8, pH2$TC_8)))
pH_mix2 <- pH_mix2 %>%  #Rename columns
  rename("H" = "X1", "Date" = "X2", "pH" = "X3", "Temp" = "X4")
pH_mix2$Date <- as.POSIXct(pH_mix2$Date, origin = "1970-01-01", TZ = "America/Los_Angeles")
pH_mix2$pH <- as.numeric(pH_mix2$pH) # Converting to numeric from factor is odd...
pH_mix2$Temp <- as.numeric(pH_mix2$Temp)

pHcal2 <- rbind(pH_long2, pH_mix2)
pHcal2 <- subset(pHcal2, pHcal2[ , 3] > 7.2)
pHcal2 <- subset(pHcal2, !is.na(pHcal2[ , 3]))
pHcal2 <- subset(pHcal2, pHcal2$Date >= as.POSIXct("2019-11-10 12:00:00", origin = "1970-01-01", tz = "America/Los_Angeles") 
                | pHcal2$Date <= as.POSIXct("2019-11-10 10:30:00", origin = "1970-01-01", tz = "America/Los_Angeles"))
pHcal2$H <- factor(pHcal2$H)

# Bring in YSI measurements to plot points over durafet time series
YSI2 = read.csv('YSI_Summer2020.csv')
# Convert time
YSI2$DT <- force_tz(as.POSIXct(with(YSI2, mdy(DATE) + hm(Time)), origin = as.POSIXct("1970-01-01", 
                                                                                   TZ = "America/Los_Angeles")), tzone = "America/Los_Angeles")
YSIcorr <- merge(chem2, YSI2, by = c('ID','DATE'))
YSIcorr$Diff <- YSIcorr$pH.out-YSIcorr$pH
YSIoff <- summarySE(data=YSIcorr, measurevar = "Diff",
                   groupvars = c("TRT"))

for(i in 1:nrow(YSI2)) {
  if(YSI2$ID[i] =='1A') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[1]
  }else if (YSI2$ID[i] =='1B') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[1]
  }else if (YSI2$ID[i] =='1C') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[1]
  }else if (YSI2$ID[i] =='2A') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[2]
  }else if (YSI2$ID[i] =='2B') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[2]
  }else if (YSI2$ID[i] =='2C') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[2]
  }else if (YSI2$ID[i] =='3A') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[3]
  }else if (YSI2$ID[i] =='3B') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[3]
  }else if (YSI2$ID[i] =='3C') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[3]
  }else if (YSI2$ID[i] =='4A') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[4]
  }else if (YSI2$ID[i] =='4B') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[4]
  }else if (YSI2$ID[i] =='4C') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[4]
  }else if (YSI2$ID[i] =='5A') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[5]
  }else if (YSI2$ID[i] =='5B') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[5]
  }else if (YSI2$ID[i] =='5C') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[5]
  }else if (YSI2$ID[i] =='6A') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[6]
  }else if (YSI2$ID[i] =='6B') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[6]
  }else if (YSI2$ID[i] =='6C') {
    YSI2$pH[i] = YSI2$pH[i] + YSIoff$Diff[6]
  }else {
    YSI2$pH[i] = YSI2$pH[i]
  }
}   

######################################################################################################
######################################################################################################
## Let's do some PCA regressions!
######################################################################################################
######################################################################################################

## Now calc PCA on average variables
pHcalsub2 <- pHcal2[pHcal2$H %in% c("H1","H2","H3","H4","H5","H6"),]
pHave2 <- summarySE(data=pHcalsub2, measurevar = "pH",
                   groupvars = "H")
pHTCave2 <- summarySE(data=pHcalsub2, measurevar = "Temp",
                     groupvars = "H")
DO_long2 <- na.omit(DO_long2)
DOave2 <- summarySE(data=DO_long2, measurevar = "DO",
                   groupvars = "H")

PCAvar2 <- data.frame(cbind(pHave2$H,pHave2$pH,pHave2$sd,pHTCave2$Temp,pHTCave2$sd,DOave2$DO,DOave2$sd))
names(PCAvar2)[1:7] <- c("Headers","pH","pH SD","Temperature","Temperature SD","DO","DO SD")
PCAvar2$Headers <- c("H1", "H2", "H3", "H4", "H5", "H6")

TS.PCA2 <- prcomp(PCAvar2[2:7], scale=TRUE) #Run PCA of pH, DO, tempDO, temppH
summary(TS.PCA2) #Look at the amount of variance explained

bi1 <- autoplot(TS.PCA, label = TRUE, label.vjust = 2, colour = rev(pal), size = 4, loadings = TRUE, 
                loadings.colour = "black", loadings.label = TRUE, loadings.label.vjust = 0.5, 
                loadings.label.colour = "black", loadings.label.hjust = 1)
bi1 <- bi1 + theme(panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_rect(colour = "black", fill = NA))
bi1 <- bi1 +  scale_y_continuous(breaks = c(-0.4, 0, 0.4))

bi2 <- autoplot(TS.PCA2, label = TRUE, label.vjust = 2, colour = rev(pal), size = 4, loadings = TRUE, 
                loadings.colour = "black", loadings.label = TRUE, loadings.label.vjust = 0.5, 
                loadings.label.colour = "black", loadings.label.hjust = 1)
bi2 <- bi2 + theme(panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_rect(colour = "black", fill = NA))
bi2

#FIGURE 3 plot
l <- ggarrange(bi1,bi2, labels = c("(a)", "(b)"), common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 1)
ggsave(plot = l, file = "Fig3.pdf", bg = "white",
       width = 20, height = 10, units = "cm", dpi = 300)


Load2 <- data.frame(TS.PCA2$x) #Turn loadings into dataframe
TrtPCA2 <- data.frame(c(PCAvar2, Load)) #Combine PCA loadings with OG data


#Making preliminary plots of data across treaments based on "PCs"
PC1Reg2 <- summarySE(data=TrtPCA2, measurevar = "PC1",
                    groupvars = "Headers")
PC2Reg2 <- summarySE(data=TrtPCA2, measurevar = "PC2",
                    groupvars = "Headers")
PC3Reg2 <- summarySE(data=TrtPCA2, measurevar = "PC3",
                    groupvars = "Headers")
PC4Reg2 <- summarySE(data=TrtPCA2, measurevar = "PC4",
                    groupvars = "Headers")

# Adding PC1 & PC2 scores to urchin and snail data for plotting & stats
for(i in 1:nrow(dat_sum3)) {
  if(dat_sum3$TRT[i] ==1) {
    dat_sum3$PC1[i] = PC1Reg2$PC1[1]
  }else if (dat_sum3$TRT[i] ==2){
    dat_sum3$PC1[i] = PC1Reg2$PC1[2]
  }else if (dat_sum3$TRT[i] ==3){
    dat_sum3$PC1[i] = PC1Reg2$PC1[3]
  }else if (dat_sum3$TRT[i] ==4){
    dat_sum3$PC1[i] = PC1Reg2$PC1[4]
  }else if (dat_sum3$TRT[i] ==5){
    dat_sum3$PC1[i] = PC1Reg2$PC1[5]
  }else if (dat_sum3$TRT[i] ==6){
    dat_sum3$PC1[i] = PC1Reg2$PC1[6]
  }else {
    dat_sum3$PC1[i] = NA
  }
}   

for(i in 1:nrow(dat_sum4)) {
  if(dat_sum4$TRT[i] ==1) {
    dat_sum4$PC1[i] = PC1Reg2$PC1[1]
  }else if (dat_sum4$TRT[i] ==2){
    dat_sum4$PC1[i] = PC1Reg2$PC1[2]
  }else if (dat_sum4$TRT[i] ==3){
    dat_sum4$PC1[i] = PC1Reg2$PC1[3]
  }else if (dat_sum4$TRT[i] ==4){
    dat_sum4$PC1[i] = PC1Reg2$PC1[4]
  }else if (dat_sum4$TRT[i] ==5){
    dat_sum4$PC1[i] = PC1Reg2$PC1[5]
  }else if (dat_sum4$TRT[i] ==6){
    dat_sum4$PC1[i] = PC1Reg2$PC1[6]
  }else {
    dat_sum4$PC1[i] = NA
  }
}  

for(i in 1:nrow(dat_sum3)) {
  if(dat_sum3$TRT[i] ==1) {
    dat_sum3$PC2[i] = PC2Reg2$PC2[1]
  }else if (dat_sum3$TRT[i] ==2){
    dat_sum3$PC2[i] = PC2Reg2$PC2[2]
  }else if (dat_sum3$TRT[i] ==3){
    dat_sum3$PC2[i] = PC2Reg2$PC2[3]
  }else if (dat_sum3$TRT[i] ==4){
    dat_sum3$PC2[i] = PC2Reg2$PC2[4]
  }else if (dat_sum3$TRT[i] ==5){
    dat_sum3$PC2[i] = PC2Reg2$PC2[5]
  }else if (dat_sum3$TRT[i] ==6){
    dat_sum3$PC2[i] = PC2Reg2$PC2[6]
  }else {
    dat_sum3$PC2[i] = NA
  }
}   

for(i in 1:nrow(dat_sum4)) {
  if(dat_sum4$TRT[i] ==1) {
    dat_sum4$PC2[i] = PC2Reg2$PC2[1]
  }else if (dat_sum4$TRT[i] ==2){
    dat_sum4$PC2[i] = PC2Reg2$PC2[2]
  }else if (dat_sum4$TRT[i] ==3){
    dat_sum4$PC2[i] = PC2Reg2$PC2[3]
  }else if (dat_sum4$TRT[i] ==4){
    dat_sum4$PC2[i] = PC2Reg2$PC2[4]
  }else if (dat_sum4$TRT[i] ==5){
    dat_sum4$PC2[i] = PC2Reg2$PC2[5]
  }else if (dat_sum4$TRT[i] ==6){
    dat_sum4$PC2[i] = PC2Reg2$PC2[6]
  }else {
    dat_sum4$PC2[i] = NA
  }
} 

# Adding average PC1 score to urchin and snail data for plotting
for(i in 1:nrow(pHcalsub2)) {
  if(pHcalsub2$H[i] =="H1") {
    pHcalsub2$PC1[i] = PC1Reg2$PC1[1]
  }else if (pHcalsub2$H[i] =="H2"){
    pHcalsub2$PC1[i] = PC1Reg2$PC1[2]
  }else if (pHcalsub2$H[i] =="H3"){
    pHcalsub2$PC1[i] = PC1Reg2$PC1[3]
  }else if (pHcalsub2$H[i] =="H4"){
    pHcalsub2$PC1[i] = PC1Reg2$PC1[4]
  }else if (pHcalsub2$H[i] =="H5"){
    pHcalsub2$PC1[i] = PC1Reg2$PC1[5]
  }else if (pHcalsub2$H[i] =="H6"){
    pHcalsub2$PC1[i] = PC1Reg2$PC1[6]
  }else {
    pHcalsub2$PC1[i] = NA
  }
}   

for(i in 1:nrow(YSI2)) {
  if(YSI2$H[i] =="H1") {
    YSI2$PC1[i] = PC1Reg2$PC1[1]
  }else if (YSI2$H[i] =="H2"){
    YSI2$PC1[i] = PC1Reg2$PC1[2]
  }else if (YSI2$H[i] =="H3"){
    YSI2$PC1[i] = PC1Reg2$PC1[3]
  }else if (YSI2$H[i] =="H4"){
    YSI2$PC1[i] = PC1Reg2$PC1[4]
  }else if (YSI2$H[i] =="H5"){
    YSI2$PC1[i] = PC1Reg2$PC1[5]
  }else if (YSI2$H[i] =="H6"){
    YSI2$PC1[i] = PC1Reg2$PC1[6]
  }else {
    YSI2$PC1[i] = NA
  }
}   


# Some plots of time series data
m <- ggplot(pHcalsub2, aes(x = Date, y = pH)) +
  geom_line(aes(x = Date, y = pH, color = factor(H)), alpha = 1/3) +
  theme_classic() +
  scale_color_manual(values = rev(wes_palette("Zissou1", 6, "continuous"))) +
  geom_point(data = YSI2, mapping = 
               aes(x = DT, y = pH, color = factor(H)), size = 1) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "2 week", expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  #theme(legend.position = "none") +
  #scale_fill_gradientn(colours = pal, name = "Upwelling Intensity (PC1)") +
  theme(legend.title = element_blank(), legend.position = "bottom", legend.key.width = unit(0.5,'cm')) +
  ylim(7.5, 8.1) +
  labs(color = "Header") +
  labs(y="pH")
m 
m <- m + guides(color = guide_legend(nrow=1))

j <- ggplot(pHcalsub2, aes(x = Date, y = Temp, color = factor(H))) +
  geom_line(aes(x = Date, y = Temp, color = factor(H)), alpha = 1/3) +
  theme_classic() +
  scale_color_manual(values = rev(wes_palette("Zissou1", 6, "continuous"))) +
  geom_point(data = YSI2, aes(x = DT, y = Temp, color = factor(H)), size = 1) +
  #theme(legend.position = "none") +
  #scale_fill_gradientn(colours = pal, name = "Upwelling Intensity (PC1)") +
  theme(legend.title = element_blank(), legend.position = "bottom", legend.key.width = unit(0.5,'cm')) +
  scale_x_datetime(date_breaks = "2 week", expand = c(0, 0)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(10, 15) +
  labs(color = "Header") +
  labs(y="Temperature (\u00B0C)")
j 
j <- j + guides(color = guide_legend(nrow=1))+ scale_y_continuous(breaks = c(10,12,14,16))

#Resample to 15 min to ease plotting and remove lines for data gaps
H1 <- DO_long2[which(DO_long2$H=="H1"),]
H1$Date <- round_date(H1$Date, "second")
H1_full <- data.frame(Date = seq(H1$Date[1], H1$Date[nrow(H1)], by=900))
H1_full$Date <- ceiling_date(H1_full$Date, "second")
H1 <- merge(H1, H1_full, all=TRUE)
H1$H <- "H1"

H2 <- DO_long2[which(DO_long2$H=="H2"),]
H2$Date <- round_date(H2$Date, "second")
H2_full <- data.frame(Date = seq(H2$Date[1], H2$Date[nrow(H2)], by=900))
H2_full$Date <- ceiling_date(H2_full$Date, "second")
H2 <- merge(H2, H2_full, all=TRUE)
H2$H <- "H2"

H3 <- DO_long2[which(DO_long2$H=="H3"),]
H3$Date <- round_date(H3$Date, "second")
H3_full <- data.frame(Date = seq(H3$Date[1], H3$Date[nrow(H3)], by=900))
H3_full$Date <- ceiling_date(H3_full$Date, "second")
H3 <- merge(H3, H3_full, all=TRUE)
H3$H <- "H3"

H4 <- DO_long2[which(DO_long2$H=="H4"),]
H4$Date <- round_date(H4$Date, "second")
H4_full <- data.frame(Date = seq(H4$Date[1], H4$Date[nrow(H4)], by=900))
H4_full$Date <- ceiling_date(H4_full$Date, "second")
H4 <- merge(H4, H4_full, all=TRUE)
H4$H <- "H4"

H5 <- DO_long2[which(DO_long2$H=="H5"),]
H5$Date <- round_date(H5$Date, "second")
H5_full <- data.frame(Date = seq(H5$Date[1], H5$Date[nrow(H5)], by=900))
H5_full$Date <- ceiling_date(H5_full$Date, "second")
H5 <- merge(H5, H5_full, all=TRUE)
H5$H <- "H5"

H6 <- DO_long2[which(DO_long2$H=="H6"),]
H6$Date <- round_date(H6$Date, "second")
H6_full <- data.frame(Date = seq(H6$Date[1], H6$Date[nrow(H6)], by=900))
H6_full$Date <- ceiling_date(H6_full$Date, "second")
H6 <- merge(H6, H6_full, all=TRUE)
H6$H <- "H6"

H1$Delete <- NA
H2$Delete <- NA
H3$Delete <- NA
H4$Delete <- NA
H5$Delete <- NA
H6$Delete <- NA

for(i in 2:(nrow(H1)-1)) {
  if (is.na(H1$DO[i]) == TRUE & is.na(H1$DO[i+1]) == FALSE
      & is.na(H1$DO[i-1]) == FALSE) {
    H1$Delete[i] <- 'D'
  } else {
    H1$Delete[i] <- 'F'
  }
} 
H1 <- H1[which(H1$Delete=="F"),]
 
for(i in 2:(nrow(H2)-1)) {
  if (is.na(H2$DO[i]) == TRUE & is.na(H2$DO[i+1]) == FALSE
      & is.na(H2$DO[i-1]) == FALSE) {
    H2$Delete[i] <- 'D'
  } else {
    H2$Delete[i] <- 'F'
  }
} 
H2 <- H2[which(H2$Delete=="F"),]

for(i in 2:(nrow(H3)-1)) {
  if (is.na(H3$DO[i]) == TRUE & is.na(H3$DO[i+1]) == FALSE
      & is.na(H3$DO[i-1]) == FALSE) {
    H3$Delete[i] <- 'D'
  } else {
    H3$Delete[i] <- 'F'
  }
} 
H3 <- H3[which(H3$Delete=="F"),]

for(i in 2:(nrow(H4)-1)) {
  if (is.na(H4$DO[i]) == TRUE & is.na(H4$DO[i+1]) == FALSE
      & is.na(H4$DO[i-1]) == FALSE) {
    H4$Delete[i] <- 'D'
  } else {
    H4$Delete[i] <- 'F'
  }
} 
H4 <- H4[which(H4$Delete=="F"),]

for(i in 2:(nrow(H5)-1)) {
  if (is.na(H5$DO[i]) == TRUE & is.na(H5$DO[i+1]) == FALSE
      & is.na(H5$DO[i-1]) == FALSE) {
    H5$Delete[i] <- 'D'
  } else {
    H5$Delete[i] <- 'F'
  }
} 
H5 <- H5[which(H5$Delete=="F"),]

for(i in 2:(nrow(H6)-1)) {
  if (is.na(H6$DO[i]) == TRUE & is.na(H6$DO[i+1]) == FALSE
      & is.na(H6$DO[i-1]) == FALSE) {
    H6$Delete[i] <- 'D'
  } else {
    H6$Delete[i] <- 'F'
  }
} 
H6 <- H6[which(H6$Delete=="F"),]

DOtest <- rbind(H1,H2,H3,H4,H5,H6)


k <- ggplot(DOtest, aes(x = Date, y = DO, color = factor(H)), group = factor(H)) +
  geom_line(alpha = 1/3) +
  theme_classic() +
  #theme(legend.position = "none") +
  scale_color_manual(values = rev(wes_palette("Zissou1", 6, "continuous"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="Dissolved Oxygen (mg/L)", x="Date") +
  scale_x_datetime(date_breaks = "1 week", date_labels = "%m/%d",
                   expand = c(0, 0)) +
  geom_point(data = YSI2, mapping = 
               aes(x = DT, y = DO1, color = factor(H)), size = 1) +
  #scale_fill_gradientn(colours = pal, name = "Upwelling Intensity") +
  theme(legend.title = element_blank(), legend.position = "bottom", legend.key.width = unit(0.5,'cm')) +
  labs(color = "Header") +
  ylim(3.75, 10) 
k 
k <- k + guides(color = guide_legend(nrow=1)) + scale_y_continuous(breaks = c(5,7,9))

#Figure S3
l <- ggarrange(m, j, k, labels = c("(a)", "(b)", "(c)"), common.legend = TRUE, legend = "bottom", ncol = 1, nrow = 3)
l
ggsave(plot = l, file = "FigS3.png", 
       type = "cairo-png",  bg = "white",
       width = 15, height = 25, units = "cm", dpi = 300)

dat_sum3$bin <- rep(NA, nrow(dat_sum3))
for(x in c('A', 'B', 'C')) dat_sum3$bin[grep(x, dat_sum3$SAMPLE)] <- x
dat_sum3$graze_corr <- dat_sum3$GrazingCorr_gWWperday/dat_sum3$WW..g.

dat_sum4$bin <- rep(NA, nrow(dat_sum4))
for(x in c('A', 'B', 'C')) dat_sum4$bin[grep(x, dat_sum4$SAMPLE)] <- x
dat_sum4$graze_corr <- dat_sum4$GrazingCorr_gWWperday/dat_sum4$WW..g.

##################################################################################

rm(list = setdiff(ls(), c("dat_sum1","dat_sum2_0","dat_sum3","dat_sum4","datUrchin","datSnail","datUrchinS","daySnailS")))

#Add in days for when measurement was made
for(i in 1:(nrow(dat_sum1))) {
if(dat_sum1$TIMESTEP[i] == 1){
  dat_sum1$Days[i] = 0
} else if (dat_sum1$TIMESTEP[i] == 2) {
  dat_sum1$Days[i] = 29
} else if (dat_sum1$TIMESTEP[i] == 3) {
  dat_sum1$Days[i] = 60
} else 
  dat_sum1$Days[i] = 85
}

for(i in 1:(nrow(dat_sum2_0))) {
  if(dat_sum2_0$TIMESTEP[i] == 1){
    dat_sum2_0$Days[i] = 0
  } else if (dat_sum2_0$TIMESTEP[i] == 2) {
    dat_sum2_0$Days[i] = 29
  } else if (dat_sum2_0$TIMESTEP[i] == 3) {
    dat_sum2_0$Days[i] = 57
  } else 
    dat_sum2_0$Days[i] = 82
}

for(i in 1:(nrow(dat_sum3))) {
  if(dat_sum3$TIMESTEP[i] == 1){
    dat_sum3$Days[i] = 0
  } else 
    dat_sum3$Days[i] = 3
}

for(i in 1:(nrow(dat_sum4))) {
  if(dat_sum4$TIMESTEP[i] == 1){
    dat_sum4$Days[i] = 0
  } else 
    dat_sum4$Days[i] = 3
}

#Select just data for timesteps not = 1 (initial)
dat_sum1_234 <- subset(dat_sum1, dat_sum1$TIMESTEP == 2 | dat_sum1$TIMESTEP == 3 | dat_sum1$TIMESTEP == 4)
dat_sum2_234 <- subset(dat_sum2_0, dat_sum2_0$TIMESTEP == 2 | dat_sum2_0$TIMESTEP == 3 | dat_sum2_0$TIMESTEP == 4)


#Models for resp for S. fran chronic exp
SFrespAll <- lmer(resp_corr_2 ~ PC1 + Days + (1|bin), data = dat_sum1_234)
SFrespAll2 <- lmer(resp_corr_2 ~ PC1 + PC2 + Days + (1|bin), data = dat_sum1_234)
summary(SFrespAll)
anova(SFrespAll)
AIC(SFrespAll, SFrespAll2)
emmeans(SFrespAll, list(pairwise ~ Days), adjust = "tukey")

#Convert negative grazing rate to zero
#Models for grazing rate of S. fran chronic exp.
dat_sum1_234$graze_corr <- ifelse(dat_sum1_234$graze_corr < 0, 0, dat_sum1_234$graze_corr)
SFgrazeAll <- lmer(graze_corr ~ PC1 + Days + (1|bin), data = dat_sum1_234)
SFgrazeAll2 <- lmer(graze_corr ~ PC1 + PC2 + Days + (1|bin), data = dat_sum1_234)
summary(SFgrazeAll)
anova(SFgrazeAll)
AIC(SFgrazeAll, SFgrazeAll2)

t <- ggplot(dat_sum1_234, aes(x=PC1, y=resp_corr_2, group = as.factor(PC1), color = as.factor(TIMESTEP))) +
  geom_point() +
  facet_wrap(~TIMESTEP) +
  theme_classic()
t

#Convert negative grazing rate to zero
dat_sum2_234$graze_corr <- ifelse(dat_sum2_234$graze_corr < 0, 0, dat_sum2_234$graze_corr)
#Models for T. pulligo resp chronic exp
TPrespAll <- lmer(resp_corr_2 ~ PC1 + Days + (1|bin), data = dat_sum2_234)
TPrespAll2 <- lmer(resp_corr_2 ~ PC1 + PC2 + Days + (1|bin), data = dat_sum2_234)
summary(TPrespAll)
anova(TPrespAll)
AIC(TPrespAll, TPrespAll2)
emmeans(TPrespAll, list(pairwise ~ Days), adjust = "tukey")

#Models for grazing rate of T. pulligo chronic exp.
TPgrazeAll <- lmer(graze_corr ~ PC1 + Days + (1|bin), data = dat_sum2_234)
TPgrazeAll2 <- lmer(graze_corr ~ PC1 + PC2 + Days + (1|bin), data = dat_sum2_234)
summary(TPgrazeAll)
anova(TPgrazeAll)
AIC(TPgrazeAll, TPgrazeAll2)

t <- ggplot(dat_sum2_234, aes(x=PC1, y=graze_corr, group = as.factor(PC1), color = as.factor(TIMESTEP))) +
  geom_point() +
  facet_wrap(~TIMESTEP) +
  theme_classic()
t

#Models for resp and grazing responses S. fran acute exp
dat_sum3$graze_corr <- ifelse(dat_sum3$graze_corr < 0, 0, dat_sum3$graze_corr)
SFrespS <- lmer(resp_corr ~ PC1 * Days + (1|bin), data = dat_sum3)
SFrespS2 <- lmer(resp_corr ~ PC1 + PC2 + Days + (1|bin), data = dat_sum3)
summary(SFrespS)
anova(SFrespS)
AIC(SFrespS, SFrespS2)
emmeans(SFrespS, list(pairwise ~ Days), adjust = "tukey")

SFgrazeSAll <- lmer(graze_corr ~ PC1 + Days + (1|bin), data = dat_sum3)
SFgrazeSAll2 <- lmer(graze_corr ~ PC1 + PC2 + Days + (1|bin), data = dat_sum3)
summary(SFgrazeSAll)
anova(SFgrazeSAll)
AIC(SFgrazeSAll, SFgrazeSAll2)
emmeans(SFgrazeSAll, list(pairwise ~ Days), adjust = "tukey")

#Models for resp and grazing responses T. pulligo acute exp
dat_sum4$graze_corr <- ifelse(dat_sum4$graze_corr < 0, 0, dat_sum4$graze_corr)
TPrespS <- lmer(resp_corr ~ PC1 * Days + (1|bin), data = dat_sum4)
TPrespS2 <- lmer(resp_corr ~ PC1 + PC2 + Days + (1|bin), data = dat_sum4)
summary(TPrespS)
anova(TPrespS)
AIC(TPrespS, TPrespS2)

TPgrazeSAll <- lmer(graze_corr ~ PC1 * Days + (1|bin), data = dat_sum4)
TPgrazeSAll2 <- lmer(graze_corr ~ PC1 + PC2 + Days + (1|bin), data = dat_sum4)
summary(TPgrazeSAll)
anova(TPgrazeSAll)
AIC(TPgrazeSAll, TPgrazeSAll2)

#Subset initial and final data points for both S. fran and T. pulligo chronic exp
dat_sum1_1 <- subset(dat_sum1, dat_sum1$TIMESTEP == 1)
dat_sum2_1 <- subset(dat_sum2_0, dat_sum2_0$TIMESTEP == 1)

dat_sum1_4 <- subset(dat_sum1, dat_sum1$TIMESTEP == 4)
dat_sum2_4 <- subset(dat_sum2_0, dat_sum2_0$TIMESTEP == 4)

#merge initial and final data by id so that we can calculate change in wet weight and buoyant wt
dat_sum1_WWBW <- merge(dat_sum1_1, dat_sum1_4, by = "SAMPLE")
dat_sum2_WWBW <- merge(dat_sum2_1, dat_sum2_4, by = "SAMPLE")

#Calculate out the % change in buoyant wt and wet weight
dat_sum1_WWBW$rGRW <- (dat_sum1_WWBW$WW..g..y - dat_sum1_WWBW$WW..g..x)/dat_sum1_WWBW$WW..g..x
dat_sum2_WWBW$rGRW <- (dat_sum2_WWBW$WW..g..y - dat_sum2_WWBW$WW..g..x)/dat_sum2_WWBW$WW..g..x

dat_sum1_WWBW$rGRB <- (dat_sum1_WWBW$BW..g..y - dat_sum1_WWBW$BW..g..x)/dat_sum1_WWBW$BW..g..x
dat_sum2_WWBW$rGRB <- (dat_sum2_WWBW$BW..g..y - dat_sum2_WWBW$BW..g..x)/dat_sum2_WWBW$BW..g..x
#Models for growth and calcification responses for Experiment 1
#Start with S. fran
SFgrowthWW <- lmer(rGRW ~ PC1.x + (1|bin.x), data = dat_sum1_WWBW)
SFgrowthWW2 <- lmer(rGRW ~ PC1.x + PC2.x + (1|bin.x), data = dat_sum1_WWBW)
summary(SFgrowthWW)
anova(SFgrowthWW)
AIC(SFgrowthWW, SFgrowthWW2)

SFgrowthBW <- lmer(rGRB ~ PC1.x + (1|bin.x), data = dat_sum1_WWBW)
SFgrowthBW2 <- lmer(rGRB ~ PC1.x + PC2.x + (1|bin.x), data = dat_sum1_WWBW)
summary(SFgrowthBW)
anova(SFgrowthBW)
AIC(SFgrowthBW, SFgrowthBW2)

#Now on to P. pulligo
TPgrowthWW <- lmer(rGRW ~ PC1.x + (1|bin.x), data = dat_sum2_WWBW)
TPgrowthWW2 <- lmer(rGRW ~ PC1.x + PC2.x + (1|bin.x), data = dat_sum2_WWBW)
summary(TPgrowthWW)
anova(TPgrowthWW)
AIC(TPgrowthWW, TPgrowthWW2)

TPgrowthBW <- lmer(rGRB ~ PC1.x + (1|bin.x), data = dat_sum2_WWBW)
TPgrowthBW2 <- lmer(rGRB ~ PC1.x + PC2.x + (1|bin.x), data = dat_sum2_WWBW)
summary(TPgrowthBW)
anova(TPgrowthBW)
AIC(TPgrowthBW, TPgrowthBW2)

####################################################################################
#Creating partial residual plots of responses for figures!#############################

#Create partial residual plot of S. fran respiration rate (chronic) removing effect of timestep & random effects
y_PC1Resp <- remef(SFrespAll, fix = "Days", ran = "all")
SFLPC1Resp <- dat_sum1_234[!is.na(dat_sum1_234$resp_corr_2), ]
SFLPC1Resp$resp_PR <- y_PC1Resp
SF_C2 <- SFLPC1Resp[which(SFLPC1Resp$TIMESTEP == 2),]
LM1 <- lm(resp_PR ~ PC1, data = SF_C2)
summary(LM1)
SF_C3 <- SFLPC1Resp[which(SFLPC1Resp$TIMESTEP == 3),]
LM2 <- lm(resp_PR ~ PC1, data = SF_C3)
summary(LM2)
SF_C4 <- SFLPC1Resp[which(SFLPC1Resp$TIMESTEP == 4),]
LM3 <- lm(resp_PR ~ PC1, data = SF_C4)
summary(LM3)

library(extrafont)
font_import()
windowsFonts()
t <- ggplot(SFLPC1Resp, aes(x=PC1, y=resp_PR, group = as.factor(TIMESTEP), 
                           color = as.factor(TIMESTEP), fill = as.factor(TIMESTEP))) +
  geom_point(alpha = 0.25) + 
  stat_smooth(method = "lm" , formula = y~x, se = TRUE) +
  scale_color_brewer(name = "Month", labels = c("1","2","3"), palette = "Dark2") +
  scale_fill_brewer(name = "Month", labels = c("1","2","3"), palette = "Dark2") +
  ylab(bquote('Partial Residual Respiration ('*mu~'mol' ~hr^-1 ~ind^-1*')'))+ 
  xlab(element_blank()) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85)) 
  #theme(text=element_text(size=12, family="Times New Roman"))
t

#Create partial residual plot of T. pulligo respiration rate (chronic) removing effect of timestep & random effects
y_PC1Resp <- remef(TPrespAll, fix = "Days", ran = "all")
TPLPC1Resp <- dat_sum2_234[!is.na(dat_sum2_234$resp_corr_2), ]
TPLPC1Resp$resp_PR <- y_PC1Resp
TP_C2 <- TPLPC1Resp[which(TPLPC1Resp$TIMESTEP == 2),]
LM1 <- lm(resp_PR ~ PC1, data = TP_C2)
summary(LM1)
TP_C3 <- TPLPC1Resp[which(TPLPC1Resp$TIMESTEP == 3),]
LM2 <- lm(resp_PR ~ PC1, data = TP_C3)
summary(LM2)
TP_C4 <- TPLPC1Resp[which(TPLPC1Resp$TIMESTEP == 4),]
LM3 <- lm(resp_PR ~ PC1, data = TP_C4)
summary(LM3)

a <- ggplot(TPLPC1Resp, aes(x=PC1, y=resp_PR, group = as.factor(TIMESTEP), 
                            color = as.factor(TIMESTEP), fill = as.factor(TIMESTEP))) +
  stat_smooth(method = "lm" , formula = y~x, se = TRUE) +
  geom_point(alpha = 0.25) + 
  scale_color_brewer(name = "Month", labels = c("1","2","3"), palette = "Dark2") +
  scale_fill_brewer(name = "Month", labels = c("1","2","3"), palette = "Dark2") +
  ylab(bquote('Partial Residual Respiration ('*mu~'mol' ~hr^-1 ~ind^-1*')'))+ 
  xlab(element_blank()) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85)) 
 # theme(text=element_text(size=12, family="Times New Roman"))
a 

#Create partial residual plot of S. fran grazing rate (chronic) removing effect of timestep & random effects
y_PC1graze <- remef(SFgrazeAll, fix = "Days", ran = "all")
SFLPC1graze <- dat_sum1_234[!is.na(dat_sum1_234$graze_corr), ]
SFLPC1graze$graze_PR <- y_PC1graze
SF_C2 <- SFLPC1graze[which(SFLPC1graze$TIMESTEP == 2),]
LM1 <- lm(graze_PR ~ PC1, data = SF_C2)
summary(LM1)
SF_C3 <- SFLPC1graze[which(SFLPC1graze$TIMESTEP == 3),]
LM2 <- lm(graze_PR ~ PC1, data = SF_C3)
summary(LM2)
SF_C4 <- SFLPC1graze[which(SFLPC1graze$TIMESTEP == 4),]
LM2 <- lm(graze_PR ~ PC1, data = SF_C4)
summary(LM2)

r <- ggplot(SFLPC1graze, aes(x=PC1, y=graze_PR, group = as.factor(TIMESTEP), 
                            color = as.factor(TIMESTEP), fill = as.factor(TIMESTEP))) +
  geom_point(alpha = 0.25) + 
  stat_smooth(method = "lm" , formula = y~x, se = TRUE) +
  scale_color_brewer(name = "Month", labels = c("1","2","3"), palette = "Dark2") +
  scale_fill_brewer(name = "Month", labels = c("1","2","3"), palette = "Dark2") +
  ylab(bquote('Partial Residual Grazing (g kelp wet wt'~d^-1*')'))+ 
  xlab("PC1") +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85)) 
  #theme(text=element_text(size=12, family="Times New Roman"))
r <- r + scale_y_continuous(breaks = c(0.050,0.125,0.20))
r

#Create partial residual plot of T. pulligo grazing rate (chronic) removing effect of timestep & random effects
y_PC1graze2 <- remef(TPgrazeAll, fix = "Days", ran = "all")
TPLPC1graze <- dat_sum2_234[!is.na(dat_sum2_234$graze_corr), ]
TPLPC1graze$graze_PR <- y_PC1graze2
TP_C2 <- TPLPC1graze[which(TPLPC1graze$TIMESTEP == 2),]
LM5 <- lm(graze_PR ~ PC1, data = TP_C2)
summary(LM5)
TP_C3 <- TPLPC1graze[which(TPLPC1graze$TIMESTEP == 3),]
LM2 <- lm(graze_PR ~ PC1, data = TP_C3)
summary(LM2)
TP_C4 <- TPLPC1graze[which(TPLPC1graze$TIMESTEP == 4),]
LM2 <- lm(graze_PR ~ PC1, data = TP_C4)
summary(LM2)

b <- ggplot(TPLPC1graze, aes(x=PC1, y=graze_PR, group = as.factor(TIMESTEP), 
                             color = as.factor(TIMESTEP), fill = as.factor(TIMESTEP))) +
  geom_point(alpha = 0.25) + 
  stat_smooth(method = "lm" , formula = y~x, se = TRUE) +
  scale_color_brewer(name = "Month", labels = c("1","2","3"), palette = "Dark2") +
  scale_fill_brewer(name = "Month", labels = c("1","2","3"), palette = "Dark2") +
  ylab(bquote('Partial Residual Grazing (g kelp wet wt'~d^-1*')'))+ 
  xlab("PC1") +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85)) 
  #theme(text=element_text(size=12, family="Times New Roman"))

b <- b + scale_y_continuous(breaks = c(0.00,0.06,0.12))
b

#Create partial residual plot of S. fran repiration rate (acute) removing effect of timestep & random effects
y_PC1Resp2 <- remef(SFrespS, fix = "Days", ran = "all")
SFSC1Resp <- dat_sum3[!is.na(dat_sum3$resp_corr), ]
SFSC1Resp$resp_PR <- y_PC1Resp2
SF_0 <- SFSC1Resp[which(SFSC1Resp$TIMESTEP == 1),]
LM1 <- lm(resp_PR ~ PC1, data = SF_0)
summary(LM1)
SF_72 <- SFSC1Resp[which(SFSC1Resp$TIMESTEP == 2),]
LM2 <- lm(resp_PR ~ PC1, data = SF_72)
summary(LM2)
s <- ggplot(SFSC1Resp, aes(x=PC1, y=resp_PR, group = as.factor(TIMESTEP), 
                            color = as.factor(TIMESTEP), fill = as.factor(TIMESTEP))) +
  geom_point(alpha = 0.25) + 
  stat_smooth(method = "lm" , formula = y~x, se = TRUE) +
  scale_color_brewer(name = "Days", labels = c("0","3"), palette = "Dark2") +
  scale_fill_brewer(name = "Days", labels = c("0","3"), palette = "Dark2") +
  ylab(element_blank())+ 
  xlab(element_blank()) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85))
  #theme(text=element_text(size=12, family="Times New Roman"))

s

#Create partial residual plot of T. pulligo respiration rate (acute) removing effect of timestep & random effects
y_PC1Resp2 <- remef(TPrespS, fix = "Days", ran = "all")
TPSC1Resp <- dat_sum4[!is.na(dat_sum4$resp_corr), ]
TPSC1Resp$resp_PR <- y_PC1Resp2

e <- ggplot(TPSC1Resp, aes(x=PC1, y=resp_PR, group = as.factor(TIMESTEP), 
                           color = as.factor(TIMESTEP), fill = as.factor(TIMESTEP))) +
  geom_point() +
  scale_color_brewer(name = "Days", labels = c("0","3"), palette = "Dark2") +
  scale_fill_brewer(name = "Days", labels = c("0","3"), palette = "Dark2") +
  ylab(bquote('Partial Residual Respiration ('*mu~'mol' ~hr^-1 ~ind^-1*')'))+ 
  xlab(element_blank()) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85)) 
  #theme(text=element_text(size=12, family="Times New Roman"))

e

#Create partial residual plot of S. fran grazing rate (acute) removing effect of timestep & random effects
y_PC1graze2 <- remef(SFgrazeSAll, fix = "Days", ran = "all")
SFSC1grazeS <- dat_sum3[!is.na(dat_sum3$graze_corr), ]
SFSC1grazeS$graze_PR <- y_PC1graze2
SF_0 <- SFSC1grazeS[which(SFSC1grazeS$TIMESTEP == 1),]
LM1 <- lm(graze_PR ~ PC1, data = SF_0)
summary(LM1)
SF_72 <- SFSC1grazeS[which(SFSC1grazeS$TIMESTEP == 2),]
LM2 <- lm(graze_PR ~ PC1, data = SF_72)
summary(LM2)
max <- max(SFSC1grazeS$graze_PR)
min <- min(SFSC1grazeS$graze_PR)

q <- ggplot(SFSC1grazeS, aes(x=PC1, y=graze_PR, group = as.factor(TIMESTEP), 
                           color = as.factor(TIMESTEP), fill = as.factor(TIMESTEP))) +
  geom_point(alpha = 0.25) + 
  stat_smooth(method = "lm" , formula = y~x, se = TRUE) +
  scale_color_brewer(name = "Days", labels = c("0","3"), palette = "Dark2") +
  scale_fill_brewer(name = "Days", labels = c("0","3"), palette = "Dark2") +
  ylab(element_blank())+ 
  xlab("PC1") +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85)) 
  #theme(text=element_text(size=12, family="Times New Roman"))

q <- q + scale_y_continuous(breaks = c(0.000,0.075,0.150))
q

#Create partial residual plot of T. pulligo grazing rate (acute) removing effect of timestep & random effects
y_PC1graze2 <- remef(TPgrazeSAll, fix = "Days", ran = "all")
TPSC1grazeS <- dat_sum4[!is.na(dat_sum4$graze_corr), ]
TPSC1grazeS$graze_PR <- y_PC1graze2

f <- ggplot(TPSC1grazeS, aes(x=PC1, y=graze_PR, group = as.factor(TIMESTEP), 
                             color = as.factor(TIMESTEP), fill = as.factor(TIMESTEP))) +
  geom_point() +
  scale_color_brewer(name = "Days", labels = c("0","3"), palette = "Dark2") +
  scale_fill_brewer(name = "Days", labels = c("0","3"), palette = "Dark2") +
  ylab(bquote('Partial Residual Grazing (g kelp wet wt'~d^-1*')'))+ 
  xlab("PC1") +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85)) 
  #theme(text=element_text(size=12, family="Times New Roman"))

f 

#Create partial residual plot of S. fran growth rate removing random effects
y_PC1WW2 <- remef(SFgrowthWW, ran = "all")
SFSC1WW<- dat_sum1_WWBW[!is.na(dat_sum1_WWBW$rGRW), ]
SFSC1WW$gW <- y_PC1WW2
LM1 <- lm(gW ~ PC1.x, data = SFSC1WW)
summary(LM1)


m <- ggplot(SFSC1WW, aes(x=PC1.x, y=gW, color = CS.x, fill = CS.x)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm" , formula = y~x, se = TRUE) +
  scale_color_manual(values = "#7570B3") +
  scale_fill_manual(values = "#7570B3") +
  ylab(bquote('Partial Residual Growth (%'~Delta*' in wet wt)'))+ 
  xlab(element_blank()) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("PC1") 
  #theme(text=element_text(size=12, family="Times New Roman")) 

m

#Create partial residual plot of T. pulligo growth rate removing random effects
y_PC1WW2 <- remef(TPgrowthWW, ran = "all")
TPSC1WW<- dat_sum2_WWBW[!is.na(dat_sum2_WWBW$rGRW), ]
TPSC1WW$gW <- y_PC1WW2
LM1 <- lm(gW ~ PC1.x, data = TPSC1WW)
summary(LM1)


d <- ggplot(TPSC1WW, aes(x=PC1.x, y=gW,color = CS.x, fill = CS.x)) +
  geom_point() +
  scale_color_manual(values = "#7570B3") +
  scale_fill_manual(values = "#7570B3") +
  ylab(bquote('Partial Residual Growth (%'~Delta*' in wet wt)'))+ 
  xlab(element_blank()) +
  theme_classic() +
  theme(legend.position = "none") 
  #theme(text=element_text(size=12, family="Times New Roman"))

d <- d + scale_y_continuous(breaks = c(-0.15,-0.05,0.05,0.15,0.25))
d

#Create partial residual plot of S. fran calcification rate removing random effects
y_PC1BW2 <- remef(SFgrowthBW, ran = "all")
SFSC1BW<- dat_sum1_WWBW[!is.na(dat_sum1_WWBW$rGRB), ]
SFSC1BW$gB <- y_PC1BW2
LM1 <- lm(gB ~ PC1.x, data = SFSC1BW)
summary(LM1)

n <- ggplot(SFSC1BW, aes(x=PC1.x, y=gB, color = CS.x, fill = CS.x)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm" , formula = y~x, se = TRUE) +
  scale_color_manual(values = "#7570B3") +
  scale_fill_manual(values = "#7570B3") +
  ylab(bquote('Partial Residual Calcification (%'~Delta*' in buoyant wt)'))+ 
  xlab("PC1") +
  theme_classic() +
  theme(legend.position = "none") 
  #theme(text=element_text(size=12, family="Times New Roman"))

n <- n + scale_y_continuous(breaks = c(0.00,0.6,1.2,1.8))

#Create partial residual plot of T. pulligo calcification rate removing random effects
y_PC1BW2 <- remef(TPgrowthBW, ran = "all")
TPSC1BW<- dat_sum2_WWBW[!is.na(dat_sum2_WWBW$rGRB), ]
TPSC1BW$gB <- y_PC1BW2
LM1 <- lm(gB ~ PC1.x, data = TPSC1BW)
summary(LM1)

c <- ggplot(TPSC1BW, aes(x=PC1.x, y=gB,color = CS.x, fill = CS.x)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm" , formula = y~x, se = TRUE) +
  scale_color_manual(values = "#7570B3") +
  scale_fill_manual(values = "#7570B3") +
  ylab(bquote('Partial Residual Calcification (%'~Delta*' in buoyant wt)'))+ 
  xlab("PC1") +
  theme_classic() +
  theme(legend.position = "none") 
 # theme(text=element_text(size=12, family="Times New Roman"))

c <- c + scale_y_continuous(breaks = c(0.00,0.12,0.24))
c

#Figure 4
l <- ggarrange(t,s,r,q, labels = c("(a)", "(b)", "(c)", "(d)"), ncol = 2, nrow = 2)
l
#Figure 5
TP <- ggarrange(a,b,c, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3)
TP
#Figure 6
SFgrowth <- ggarrange(m, n, labels = c("(a)","(b)"), ncol = 1, nrow = 2)
SFgrowth
#Figure S4
Sup <- ggarrange(e,f,d, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3)
Sup

ggsave(plot = l, file = "Fig4.png", 
       type = "cairo-png",  bg = "transparent",
       width = 30, height = 20, units = "cm", dpi = 300)

ggsave(plot = TP, file = "Fig5.png", 
       type = "cairo-png",  bg = "transparent",
       width = 15, height = 35, units = "cm", dpi = 300)

ggsave(plot = SFgrowth, file = "Fig6.png", 
       type = "cairo-png",  bg = "transparent",
       width = 15, height = 23.5, units = "cm", dpi = 300)

ggsave(plot = Sup, file = "FigS4.png", 
       type = "cairo-png",  bg = "transparent",
       width = 15, height = 35, units = "cm", dpi = 300)










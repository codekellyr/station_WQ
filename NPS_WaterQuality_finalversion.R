################################################################################
#   Compare TPO4 at N01 and SRS1C to NE1, independent of seasonal variability  #
################################################################################

# K. McCaffrey, last edited 30 Aug 2023

library(readxl)
library(lubridate)
library(ggplot2)
library(jtools)
library(zoo)
library(gridExtra)
library(imputeTS)
library(forecast)
library(stlplus)
library(tseries)
library(FSA)

#### Load and format data, General Exploration ####
data<-read_excel("n01_ne1_srs1c_TP_lat_long.xlsx", sheet="n01_ne1_srs1c_TP_lat_long")
names(data)<-c("Project", "Station", "Date", "Time", "SampleCode", "Sample",
               "Coll", "UD", "DIS", "Depth.m", "Depth.unit", "DCS", "Latitude",
               "Longitude", "TPO4")
# TPO4 - total phosphate as P (mg/L)
# DIS - discharge code (1/0). 1 is discharging at structure on L29
# DSC - depth to consolidate substrate in meters
# Depth - water depth in meters

points<-read_excel("n01_ne1_srs1c_TP_lat_long.xlsx", sheet="stationCoordinatesLatLongUTM")

# format data columns
data$Project<-as.factor(data$Project)
data$Station<-as.factor(data$Station)
data$Time<-format(data$Time, format="%H:%M")
data$Date<-as.Date(as.character(data$Date), "%Y%m%d")

str(data)

# plot(data$Date, data$DCS) 
# plot(data$Date, data$Depth.m)
# plot(data$Date, data$TPO4)

# plot total P by station, through time
plot<-ggplot(data, aes(Date, TPO4))+
  geom_line()+
  ggtitle("")+
  xlab("Date")+
  ylab("Total Phosphate (mg/L)") +
  scale_x_date(breaks=seq(as.Date("2015-01-01"), as.Date("2023-12-31"), by="2 years"), 
               date_labels="%Y", 
               minor_breaks="6 months")+
  theme_nice(base_size=12)+
  theme(axis.text.x=element_text(angle=35))
plot2<-plot+facet_grid(.~Station)
raw_P_time_plot<-plot2+geom_vline(xintercept = data$Date[data$DIS==1], color="red")
# png("raw_P_time_plot.png", width =6.5, height =5, units="in", res=300)
raw_P_time_plot
# dev.off()

# plot water depth by station, through time
plot<-ggplot(data, aes(Date, Depth.m))+
  geom_line()+
  ggtitle("")+
  xlab("Date")+
  ylab("Water Depth (m)") +
  scale_x_date(breaks=seq(as.Date("2015-01-01"), as.Date("2023-12-31"), by="2 years"), 
               date_labels="%Y", 
               minor_breaks="6 months")+
  theme_nice(base_size=12)+
  theme(axis.text.x=element_text(angle=35))
plot2<-plot+facet_grid(.~Station)
raw_depth_time_plot<-plot2+geom_vline(xintercept = data$Date[data$DIS==1], color="red")
# png("raw_depth_time_plot.png", width =6.5, height =5, units="in", res=300)
raw_depth_time_plot
# dev.off()

# plot DCS by station, through time
plot<-ggplot(data, aes(Date, DCS))+
  geom_line()+
  ggtitle("")+
  xlab("Date")+
  ylab("Depth to Consolidated Sediment (m)") +
  scale_x_date(breaks=seq(as.Date("2015-01-01"), as.Date("2023-12-31"), by="2 years"), 
               date_labels="%Y", 
               minor_breaks="6 months")+
  theme_nice(base_size=12)+
  theme(axis.text.x=element_text(angle=35))
plot2<-plot+facet_grid(.~Station)
raw_DCS_time_plot<-plot2+geom_vline(xintercept = data$Date[data$DIS==1], color="red")
# png("raw_DCS_time_plot.png", width =6.5, height =5, units="in", res=300)
raw_DCS_time_plot
# dev.off()

# split data by station
NE1<-subset(data, data$Station=="NE1")
N01<-subset(data, data$Station=="N01")
SRS1C<-subset(data, data$Station=="SRS1C")

# basic summary stats
summary(NE1$TPO4); sd(NE1$TPO4)
summary(SRS1C$TPO4); sd(SRS1C$TPO4)
summary(N01$TPO4); sd(N01$TPO4)

# test for normality, examine distribution
hist(N01$TPO4); shapiro.test(N01$TPO4) # p < 0.001
hist(NE1$TPO4); shapiro.test(NE1$TPO4) # p < 0.001
hist(SRS1C$TPO4); shapiro.test(SRS1C$TPO4) # p < 0.001

# also test depths for normality
shapiro.test(N01$Depth.m)
shapiro.test(NE1$Depth.m)
shapiro.test(SRS1C$Depth.m)

shapiro.test(N01$DCS)
shapiro.test(NE1$DCS)
shapiro.test(SRS1C$DCS)

# non-parametric - is Depth correlated with DCS?
cor.test(N01$Depth.m,N01$DCS, method="spearman")
cor.test(NE1$Depth.m,NE1$DCS, method="spearman")
cor.test(SRS1C$Depth.m,SRS1C$DCS, method="spearman")

# is DCS correlated with TPO4?
cor.test(N01$DCS,N01$TPO4, method="spearman")
cor.test(NE1$DCS,NE1$TPO4, method="spearman")
cor.test(SRS1C$DCS,SRS1C$TPO4, method="spearman")

# is Depth correlated with TPO4?
cor.test(N01$Depth.m, N01$TPO4, method="spearman")
cor.test(NE1$Depth.m, NE1$TPO4, method="spearman")
cor.test(SRS1C$Depth.m, SRS1C$TPO4, method="spearman")

#### Create timeseries, examine missing data and outliers ####
range(N01$Date); range(NE1$Date); range(SRS1C$Date) # 8 years of data, roughly monthly

# convert TPO4 to a ts object (missing values filled as NA)
N01_ts<-data.frame(date=N01$Date, TPO4=N01$TPO4)
N01_ts<-read.zoo(N01_ts, FUN=as.yearmon)
N01_ts<-as.ts(N01_ts)

NE1_ts<-data.frame(date=NE1$Date, TPO4=NE1$TPO4)
NE1_ts<-read.zoo(NE1_ts, FUN=as.yearmon)
NE1_ts<-as.ts(NE1_ts)

SRS1C_ts<-data.frame(date=SRS1C$Date, TPO4=SRS1C$TPO4)
SRS1C_ts<-read.zoo(SRS1C_ts, FUN=as.yearmon)
SRS1C_ts<-as.ts(SRS1C_ts)

tsp(N01_ts) # 2015, 2023, 12
tsp(NE1_ts)
tsp(SRS1C_ts)

# what proportion of data is missing at each station?
# a large proportion, or multiple consecutive missing values may distort seasonality
statsNA(N01_ts) # 17.5% missing values in time series, 10 gaps, max size 5, avg 1.7
statsNA(NE1_ts) # 21.6% missing values in time series, 12 gaps, max size 4, avg 1.75
statsNA(SRS1C_ts) # 24.7% missing values in time series, 10 gaps, max size 8, avg 2.4

ggplot_na_distribution(N01_ts)
ggplot_na_distribution(NE1_ts)
ggplot_na_distribution(SRS1C_ts) 

# are there outliers in the data?
boxplot(N01$TPO4)$out # 4 outliers
boxplot(NE1$TPO4)$out # 10 outliers
boxplot(SRS1C$TPO4)$out # 4 outliers
# maybe - look for outliers after accounting for seasonal variation
tsoutliers(N01_ts) # 4 outliers
tsoutliers(NE1_ts) # 2 outliers
tsoutliers(SRS1C_ts) # 1 outlier

#### Create dataset versions with and without outlier replacement and impute NA ####

## test methods of imputation
eg1<-ggplot_na_imputations(x_with_na=N01_ts, x_with_imputations=na_kalman(N01_ts, model="StructTS"),
                           title="Kalman StructTS",
                           ylab = "Total Phosphate(mg/L)")
eg2<-ggplot_na_imputations(x_with_na=N01_ts, x_with_imputations=na_kalman(N01_ts, model="auto.arima"),
                           title="Kalman auto.arima",
                           ylab = "Total Phosphate(mg/L)") # NO
eg3<-ggplot_na_imputations(x_with_na=N01_ts, x_with_imputations=na_interpolation(N01_ts, option="spline", method="natural"),
                           title="Interpolation natural splines",
                           ylab = "Total Phosphate(mg/L)")

eg1
eg2
eg3
# png("impute_plot_N01.png", width =11, height =6, units="in", res=300)
grid.arrange(eg1, eg2, eg3, nrow=1)
# dev.off()

# I think kalman is most realistic (kalman can be good for ts with seasonality and trend)
ggplot_na_imputations(x_with_na=NE1_ts, x_with_imputations=na_kalman(NE1_ts, model="StructTS"))
ggplot_na_imputations(x_with_na=NE1_ts, x_with_imputations=na_kalman(NE1_ts, model="auto.arima")) # NO
ggplot_na_imputations(x_with_na=NE1_ts, x_with_imputations=na_interpolation(NE1_ts, option="spline", method="natural"))

ggplot_na_imputations(x_with_na=SRS1C_ts, x_with_imputations=na_kalman(SRS1C_ts, model="StructTS"))
ggplot_na_imputations(x_with_na=SRS1C_ts, x_with_imputations=na_kalman(SRS1C_ts, model="auto.arima")) # NO
ggplot_na_imputations(x_with_na=SRS1C_ts, x_with_imputations=na_interpolation(SRS1C_ts, option="spline", method="natural"))

## complete imputation
N01_ts_int<-na_kalman(N01_ts, model="StructTS")
NE1_ts_int<-na_kalman(NE1_ts, model="StructTS")
SRS1C_ts_int<-na_kalman(SRS1C_ts, model="StructTS")

## for dataset with modified outliers, replace outliers then impute
N01_ts_clean<-tsclean(N01_ts, replace.missing = F) # modify outliers with linear interpolation, do not replace missing here
N01_ts_clean<-na_kalman(N01_ts_clean, model="StructTS")

NE1_ts_clean<-tsclean(NE1_ts, replace.missing=F)
NE1_ts_clean<-na_kalman(NE1_ts_clean, model="StructTS")

SRS1C_ts_clean<-tsclean(SRS1C_ts, replace.missing=F)
SRS1C_ts_clean<-na_kalman(SRS1C_ts_clean, model="StructTS")



#### seasonal timeseries decomposition - with outliers ####

## do acf plot
ggAcf(N01_ts_int, lag.max=97) # data appears to have seasonality
ggAcf(NE1_ts_int, lag.max=97) # maybe just an outlier in Nov?
ggAcf(SRS1C_ts_int, lag.max=97) # maybe just outliers in Feb and Dec? appears to have some seasonality, not passing CI
# with white noise, expect 95% spikes to be inside

## fit STL window size
# for seasonal - assuming seasonal trend is changing on shorter term basis as water management changes
# also reasonable because depth is changing thorugh itme
plot_seasonal(stlplus(N01_ts_int, s.window=7))
plot_seasonal(stlplus(N01_ts_int, s.window=13)) # I think 13 looks good
plot_seasonal(stlplus(N01_ts_int, s.window=15))
plot_seasonal(stlplus(N01_ts_int, s.window=17))

plot_seasonal(stlplus(NE1_ts_int, s.window=7))
plot_seasonal(stlplus(NE1_ts_int, s.window=13)) # I think 13 looks good
plot_seasonal(stlplus(NE1_ts_int, s.window=15))
plot_seasonal(stlplus(NE1_ts_int, s.window=17))

plot_seasonal(stlplus(SRS1C_ts_int, s.window=7))
plot_seasonal(stlplus(SRS1C_ts_int, s.window=13)) # I think 13 looks good
plot_seasonal(stlplus(SRS1C_ts_int, s.window=15))
plot_seasonal(stlplus(SRS1C_ts_int, s.window=17))

# set window to 13, otherwise use default settings
N01_stl<-stlplus(N01_ts_int, s.window=13)
NE1_stl<-stlplus(NE1_ts_int, s.window=13)
SRS1C_stl<-stlplus(SRS1C_ts_int, s.window=13)

# png("N01_stl_plot_with_outliers.png", width =6, height =6, units="in", res=300)
plot(N01_stl, main="N01", xlab="Date", ylab="Total Phosphate (mg/L)") # remaining spikes in remainder correspond to more extreme values in raw data
# dev.off()

# png("NE1_stl_plot_with_outliers.png", width =6, height =6, units="in", res=300)
plot(NE1_stl, main="NE1", xlab="Date", ylab="Total Phosphate (mg/L)") # remaining spikes in remainder correspond to more extreme values in raw data
# dev.off()

# png("SRS1C_stl_plot_with_outliers.png", width =6, height =6, units="in", res=300)
plot(SRS1C_stl, main="SRS1C", xlab="Date", ylab="Total Phosphate (mg/L)") # remaining spikes in remainder correspond to more extreme values in raw data
# dev.off()

plot_trend(N01_stl) # plot trend and remainder with trend overlay
plot_trend(NE1_stl) 
plot_trend(SRS1C_stl)

plot_rembycycle(N01_stl) # plot remainder component by cycle
plot_rembycycle(NE1_stl)
plot_rembycycle(SRS1C_stl)

# is remainder stationary?
adf.test(N01_stl$data$remainder) # p < 0.01, stationary
adf.test(NE1_stl$data$remainder) # p < 0.01, stationary
adf.test(SRS1C_stl$data$remainder) # p < 0.01, stationary

# acf plots on remainder
# only going to use original data to find deseas and remain here, so not using imputed values
N01_tab<-data.frame(raw=c(N01_ts), impute=N01_stl$data$raw, 
                    deseas=(c(N01_ts)-N01_stl$data$seasonal), 
                    remain=(c(N01_ts)-(N01_stl$data$seasonal+N01_stl$data$trend)), 
                    date=as.yearmon(N01_stl$time))
NE1_tab<-data.frame(raw=c(NE1_ts), impute=NE1_stl$data$raw, 
                    deseas=(c(NE1_ts)-NE1_stl$data$seasonal), 
                    remain=(c(NE1_ts)-(NE1_stl$data$seasonal+NE1_stl$data$trend)), 
                    date=as.yearmon(NE1_stl$time))
SRS1C_tab<-data.frame(raw=c(SRS1C_ts), impute=SRS1C_stl$data$raw, 
                      deseas=(c(SRS1C_ts)-SRS1C_stl$data$seasonal), 
                      remain=(c(SRS1C_ts)-(SRS1C_stl$data$seasonal+SRS1C_stl$data$trend)), 
                      date=as.yearmon(SRS1C_stl$time))

plot(N01_tab$date, N01_tab$deseas, type="l")
plot(NE1_tab$date, NE1_tab$deseas, type="l")
plot(SRS1C_tab$date, SRS1C_tab$deseas, type="l")

test<-N01_tab[,c("date", "remain")]
test$date<-as.Date(test$date)
test<-read.zoo(test, FUN=as.yearmon)
N01_ts_rem<-as.ts(test)

test<-NE1_tab[,c("date", "remain")]
test$date<-as.Date(test$date)
test<-read.zoo(test, FUN=as.yearmon)
NE1_ts_rem<-as.ts(test)

test<-SRS1C_tab[,c("date", "remain")]
test$date<-as.Date(test$date)
test<-read.zoo(test, FUN=as.yearmon)
SRS1C_ts_rem<-as.ts(test)

#acf plot
ggAcf(N01_ts_rem, na.action=na.pass, lag.max=97) # just over 5% cross CI, no obvious pattern, likely extreme values
ggAcf(NE1_ts_rem, na.action=na.pass, lag.max=97) # just over or under 5% cross no obvious pattern, likely extreme values
ggAcf(SRS1C_ts_rem, na.action=na.pass, lag.max=97) # maybe some pattern left, but not regular

# combine for plotting
NE1_tab$month<-as.factor(month(NE1_tab$date))
N01_tab$month<-as.factor(month(N01_tab$date))
SRS1C_tab$month<-as.factor(month(SRS1C_tab$date))

NE1_tab$station<-rep("NE1")
N01_tab$station<-rep("N01")
SRS1C_tab$station<-rep("SRS1C")

deseas_tab<-rbind(NE1_tab, N01_tab, SRS1C_tab)

# ggplot(data=NE1_tab, aes(x=month, y=deseas))+
#   stat_boxplot(geom="errorbar")+
#   geom_boxplot()+
#   xlab("Month")+
#   ylab("TPO4 (mg/L)")
# ggplot(data=N01_tab, aes(x=month, y=deseas))+
#   stat_boxplot(geom="errorbar")+
#   geom_boxplot()+
#   xlab("Month")+
#   ylab("TPO4 (mg/L)")
# ggplot(data=SRS1C_tab, aes(x=month, y=deseas))+
#   stat_boxplot(geom="errorbar")+
#   geom_boxplot()+
#   xlab("Month")+
#   ylab("TPO4 (mg/L)")

ggplot(data=deseas_tab, aes(x=month, y=deseas, fill=station))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  xlab("Month")+
  ylab("TPO4 (mg/L)")

# compare to raw
ggplot(data=deseas_tab, aes(x=month, y=raw, fill=station))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  xlab("Month")+
  ylab("TPO4 (mg/L)")

# kruskal-wallis test, compare months per station
kruskal.test(NE1_tab$deseas, NE1_tab$month) # no sig diff
kruskal.test(N01_tab$deseas, N01_tab$month) # no sig diff
kruskal.test(SRS1C_tab$deseas, SRS1C_tab$month) # no sig diff
# no sig diff by month

# compare deases between stations
kruskal.test(deseas_tab$deseas, deseas_tab$month) # no sig diff
kruskal.test(deseas_tab$deseas, deseas_tab$station) # significant
dunnTest(deseas~station, data=deseas_tab, method="bonferroni") 
# Comparison         Z      P.unadj        P.adj
# 1   N01 - NE1 10.730608 7.311342e-27 2.193403e-26
# 2 N01 - SRS1C  3.611423 3.045213e-04 9.135639e-04
# 3 NE1 - SRS1C -6.921502 4.468806e-12 1.340642e-11

# all stations are different from each other

deseas_tab$station2<-as.factor(deseas_tab$station)
deseas_tab$station2<-factor(deseas_tab$station2, levels=c("NE1", "SRS1C", "N01"))
ggplot(data=deseas_tab, aes(x=station2, y=deseas, fill=station))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(notch=T)+
  xlab("Station")+
  ylab("TPO4 (mg/L)")+
  theme_nice()

#### seasonal timeseries decomposition - replaced outliers ####
# see if acf plots of remainder after removal of seasonality and trend improves

## do acf plot
ggAcf(N01_ts_clean, lag.max=97)+theme_nice() # 
ggAcf(NE1_ts_clean, lag.max=97) 
ggAcf(SRS1C_ts_clean, lag.max=97) # not clearing CI
# with white noise, expect 95% spikes to be inside

## choose window
plot_seasonal(stlplus(N01_ts_clean, s.window=7))
plot_seasonal(stlplus(N01_ts_clean, s.window=13)) # I think 13 looks good
plot_seasonal(stlplus(N01_ts_clean, s.window=15))
plot_seasonal(stlplus(N01_ts_clean, s.window=17))

plot_seasonal(stlplus(NE1_ts_clean, s.window=7))
plot_seasonal(stlplus(NE1_ts_clean, s.window=13)) # I think 13 looks good
plot_seasonal(stlplus(NE1_ts_clean, s.window=15))
plot_seasonal(stlplus(NE1_ts_clean, s.window=17))

plot_seasonal(stlplus(SRS1C_ts_clean, s.window=7))
plot_seasonal(stlplus(SRS1C_ts_clean, s.window=13)) # I think 13 looks good
plot_seasonal(stlplus(SRS1C_ts_clean, s.window=15))
plot_seasonal(stlplus(SRS1C_ts_clean, s.window=17))

# set window to 13
N01_stl_clean<-stlplus(N01_ts_clean, s.window=13)
NE1_stl_clean<-stlplus(NE1_ts_clean, s.window=13)
SRS1C_stl_clean<-stlplus(SRS1C_ts_clean, s.window=13)

# png("N01_stl_plot_without_outliers.png", width =6, height =6, units="in", res=300)
plot(N01_stl_clean, main="N01", xlab="Date", ylab="Total Phosphate (mg/L)") # remaining spikes in remainder correspond to more extreme values in raw data
# dev.off()
# png("NE1_stl_plot_without_outliers.png", width =6, height =6, units="in", res=300)
plot(NE1_stl_clean, main="NE1", xlab="Date", ylab="Total Phosphate (mg/L)") # remaining spikes in remainder correspond to more extreme values in raw data
# dev.off()
# png("SRS1C_stl_plot_without_outliers.png", width =6, height =6, units="in", res=300)
plot(SRS1C_stl_clean, main="SRS1C", xlab="Date", ylab="Total Phosphate (mg/L)") # remaining spikes in remainder correspond to more extreme values in raw data
# dev.off()

# remaining spikes in remainder correspond to more extreme values in raw data

plot_trend(N01_stl_clean) # plot trend and remainder with trend overlay
plot_trend(NE1_stl_clean) 
plot_trend(SRS1C_stl_clean)

plot_rembycycle(N01_stl_clean) # plot remainder component by cycle
plot_rembycycle(NE1_stl_clean)
plot_rembycycle(SRS1C_stl_clean)

# is remainder stationary?
adf.test(N01_stl_clean$data$remainder) # p < 0.01, stationary
adf.test(NE1_stl_clean$data$remainder) # p < 0.01, stationary
adf.test(SRS1C_stl_clean$data$remainder) # p < 0.01, stationary

# acf plots on remainder
# only going to use original data to find deseas and remain here, so not using imputed values
N01_tab_clean<-data.frame(raw=c(N01_ts), impute=N01_stl_clean$data$raw, 
                          deseas=(c(N01_ts)-N01_stl_clean$data$seasonal), 
                          remain=(c(N01_ts)-(N01_stl_clean$data$seasonal+N01_stl_clean$data$trend)), 
                          date=as.yearmon(N01_stl_clean$time))
NE1_tab_clean<-data.frame(raw=c(NE1_ts), impute=NE1_stl_clean$data$raw, 
                          deseas=(c(NE1_ts)-NE1_stl_clean$data$seasonal), 
                          remain=(c(NE1_ts)-(NE1_stl_clean$data$seasonal+NE1_stl_clean$data$trend)), 
                          date=as.yearmon(NE1_stl_clean$time))
SRS1C_tab_clean<-data.frame(raw=c(SRS1C_ts), impute=SRS1C_stl_clean$data$raw, 
                            deseas=(c(SRS1C_ts)-SRS1C_stl_clean$data$seasonal), 
                            remain=(c(SRS1C_ts)-(SRS1C_stl_clean$data$seasonal+SRS1C_stl_clean$data$trend)), 
                            date=as.yearmon(SRS1C_stl_clean$time))

plot(N01_tab_clean$date, N01_tab_clean$deseas, type="l")
plot(NE1_tab_clean$date, NE1_tab_clean$deseas, type="l")
plot(SRS1C_tab_clean$date, SRS1C_tab_clean$deseas, type="l")

test<-N01_tab_clean[,c("date", "remain")]
test$date<-as.Date(test$date)
test<-read.zoo(test, FUN=as.yearmon)
N01_ts_rem_clean<-as.ts(test)

test<-NE1_tab_clean[,c("date", "remain")]
test$date<-as.Date(test$date)
test<-read.zoo(test, FUN=as.yearmon)
NE1_ts_rem_clean<-as.ts(test)

test<-SRS1C_tab_clean[,c("date", "remain")]
test$date<-as.Date(test$date)
test<-read.zoo(test, FUN=as.yearmon)
SRS1C_ts_rem_clean<-as.ts(test)

#acf plot
ggAcf(N01_ts_rem_clean, na.action=na.pass, lag.max=97)+theme_nice() # just over 5% cross CI, no obvious pattern, likely extreme values
ggAcf(NE1_ts_rem_clean, na.action=na.pass, lag.max=97)+theme_nice() # just over or under 5% cross no obvious pattern, likely extreme values
ggAcf(SRS1C_ts_rem_clean, na.action=na.pass, lag.max=97)+theme_nice() # maybe some pattern left, but not regular
# better! SRS1c is only remaining one that may have some issues

# combine for plotting
NE1_tab_clean$month<-as.factor(month(NE1_tab_clean$date))
N01_tab_clean$month<-as.factor(month(N01_tab_clean$date))
SRS1C_tab_clean$month<-as.factor(month(SRS1C_tab_clean$date))

NE1_tab_clean$station<-rep("NE1")
N01_tab_clean$station<-rep("N01")
SRS1C_tab_clean$station<-rep("SRS1C")

deseas_tab_clean<-rbind(NE1_tab_clean, N01_tab_clean, SRS1C_tab_clean)

# ggplot(data=NE1_tab_clean, aes(x=month, y=deseas))+
#   stat_boxplot(geom="errorbar")+
#   geom_boxplot()+
#   xlab("Month")+
#   ylab("TPO4 (mg/L)")
# ggplot(data=N01_tab_clean, aes(x=month, y=deseas))+
#   stat_boxplot(geom="errorbar")+
#   geom_boxplot()+
#   xlab("Month")+
#   ylab("TPO4 (mg/L)")
# ggplot(data=SRS1C_tab_clean, aes(x=month, y=deseas))+
#   stat_boxplot(geom="errorbar")+
#   geom_boxplot()+
#   xlab("Month")+
#   ylab("TPO4 (mg/L)")

ggplot(data=deseas_tab_clean, aes(x=month, y=deseas, fill=station))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  xlab("Month")+
  ylab("TPO4 (mg/L)")

# compare to raw
ggplot(data=deseas_tab_clean, aes(x=month, y=raw, fill=station))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot()+
  xlab("Month")+
  ylab("TPO4 (mg/L)")

shapiro.test(NE1_tab_clean$deseas)
shapiro.test(N01_tab_clean$deseas)
shapiro.test(SRS1C_tab_clean$deseas)

# kruskal-wallis test, compare months per station
kruskal.test(NE1_tab_clean$deseas, NE1_tab_clean$month) # no sig diff
kruskal.test(N01_tab_clean$deseas, N01_tab_clean$month) # no sig diff
kruskal.test(SRS1C_tab_clean$deseas, SRS1C_tab_clean$month) # no sig diff
# no sig diff by month

# compare deases between stations
kruskal.test(deseas_tab_clean$deseas, deseas_tab_clean$month) # no sig diff
kruskal.test(deseas_tab_clean$deseas, deseas_tab_clean$station) # significant
dunnTest(deseas~station, data=deseas_tab_clean, method="bonferroni") 
# Comparison         Z      P.unadj        P.adj
# 1   N01 - NE1 10.509969 7.771922e-26 2.331577e-25
# 2 N01 - SRS1C  3.323149 8.900746e-04 2.670224e-03
# 3 NE1 - SRS1C -6.990564 2.737833e-12 8.213499e-12

# all stations are different from each other
deseas_tab_clean$station2<-as.factor(deseas_tab_clean$station)
deseas_tab_clean$station2<-factor(deseas_tab_clean$station2, levels=c("NE1", "SRS1C", "N01"))
# png("station_comparison_plot.png", width =6, height =6, units="in", res=300)
ggplot(data=deseas_tab_clean, aes(x=station2, y=deseas))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(notch=T, fill="grey70")+
  xlab("Station")+
  ylab("Total Phosphate (mg/L)")+
  theme(base_size=12)+
  theme_nice()
# dev.off()

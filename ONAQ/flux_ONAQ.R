setwd("C:/Users/u1297807/Box/NEON transpiration/Keeling/scripts")
#explore level 4 flux data for ONAQ
library(neonUtilities)
library(BiocManager)
library(rhdf5)
library(fractal)
library(ggplot2)
library(dygraphs)
library(xts)
source("../source/multiplot.R")

options(stringsAsFactors=F)

zipsByProduct(dpID="DP1.00024.001", package="basic",
              site=c("ONAQ"),
              savepath="C:/Users/u1297807/Box/NEON transpiration/Keeling/data",
              check.size=F,
              avg="30")

flux<-stackEddy(filepath = "C:/Users/u1297807/Box/NEON transpiration/Keeling/data/filesToStack00200", level="dp04")
names(flux)
head(flux$ONAQ)

stackByTable(filepath = "C:/Users/u1297807/Box/NEON transpiration/Keeling/data/filesToStack00024",
                  folder=T)
par<-read.csv(file = "C:/Users/u1297807/Box/NEON transpiration/Keeling/data/filesToStack00024/stackedFiles/PARPAR_30min.csv")

#look up terms
term<-unlist(strsplit(names(flux$ONAQ), split=".", fixed=T))
flux$objDesc[which(flux$objDesc$Object %in% term),]
flux$variables

#times are in GMT
timeB <- substring(flux$ONAQ$timeBgn, 1, nchar(flux$ONAQ$timeBgn)-4)
timeB <- strptime(timeB, format="%Y-%m-%dT%H:%M:%S", tz="GMT")
timeB <- as.POSIXct(timeB)
flux$ONAQ<-cbind(timeB, flux$ONAQ)

par$dt<-as.POSIXct(par$startDateTime, format="%Y-%m-%dT%H:%M:%SZ", tz="GMT")

###employ filtering criteria 
#NEE: this is -50/50 umol/m2/s
#LH: -200/500 W/m2 
#S: -500/1400 W/m2

flux$ONAQ$NEE2<-ifelse(flux$ONAQ$data.fluxCo2.nsae.flux>-50&flux$ONAQ$data.fluxCo2.nsae.flux<50, flux$ONAQ$data.fluxCo2.nsae.flux, NA)
flux$ONAQ$LH2<-ifelse(flux$ONAQ$data.fluxH2o.nsae.flux>-200&flux$ONAQ$data.fluxH2o.nsae.flux<500, flux$ONAQ$data.fluxH2o.nsae.flux, NA)
flux$ONAQ$S2<-ifelse(flux$ONAQ$data.fluxTemp.nsae.flux>-500&flux$ONAQ$data.fluxTemp.nsae.flux<1400, flux$ONAQ$data.fluxTemp.nsae.flux, NA)

#apply a moving window 7-point median filter
flux$ONAQ$NEE<-medianFilter(flux$ONAQ$NEE2, order=7)
flux$ONAQ$LH<-medianFilter(flux$ONAQ$LH2, order=7)
flux$ONAQ$S<-medianFilter(flux$ONAQ$S2, order=7)

#when is first and last non-NaN values
ind<-which(!is.nan(flux$ONAQ$data.fluxH2o.nsae.flux))
flux2<-flux$ONAQ[ind[1]:tail(ind,1),]

#add top PAR sensor
for(i in 1:nrow(flux2)){
  flux2$par[i]<-par$PARMean[which(flux2$timeB[i]==par$dt&par$verticalPosition==40)]
}

fig1a<-ggplot(flux2)+
  geom_line(aes(x=timeB, y=NEE), col="black")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y")+
  scale_y_continuous(expression(paste("NEE (", mu, "mol", " C", O[2], " ", m^-2, " ", s^-1, ")")))+
  theme_bw()+
  theme(axis.title.x = element_blank())

fig1ab<-ggplot(flux2)+
  geom_line(aes(x=timeB, y=par), col="black")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y")+
  scale_y_continuous(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1, ")")))+
  theme_bw()+
  theme(axis.title.x = element_blank())

fig1b<-ggplot(flux2)+
  geom_line(aes(x=timeB, y=LH), col="black")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y")+
  scale_y_continuous(expression(paste("LH (W ", m^-2, ")")))+
  theme_bw()+
  theme(axis.title.x = element_blank())
fig1c<-ggplot(flux2)+
  geom_line(aes(x=timeB, y=S), col="black")+
  scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y")+
  scale_y_continuous(expression(paste("SH (W ", m^-2, ")")))+
  theme_bw()+
  theme(axis.title.x = element_blank())

jpeg(filename = "../plots/NEE_par_ts.jpeg", height=5, width=10, units="in", res=600)
multiplot(fig1a, fig1ab, cols = 1)
dev.off()
#dygraph
flux3<-xts(flux2[,-1*c(2:36,37:38)], order.by=flux2$timeB)
dygraph(flux3)%>%dyRangeSelector()

###plot light response curve
fig2<-ggplot(flux2)+
  geom_point(aes(x=par, y=NEE))+
  scale_y_continuous(expression(paste("NEE (", mu, "mol", " C", O[2], " ", m^-2, " ", s^-1, ")")))+
  scale_x_continuous(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1, ")")))+
  theme_bw()

jpeg(filename = "../plots/lightresponse.jpeg", height=4, width=6, units="in", res=600)
print(fig2)
dev.off()
#bring in soil moisture + VPD
load(file="../../env/NEON/RH_high.Rdata")
#match RH_high to flux2 data
flux2$VPD<-c()
for(i in 1:nrow(flux2)){
  flux2$VPD[i]<-RH_high$RHMean[which(flux2$timeB[i]==RH_high$dt)]
}


###further cleaning
#separate  out day and night NEE
fluxn<-subset(flux2, flux2$par<=1)
fluxd<-subset(flux2, flux2$par>1)

ggplot(fluxd)+
  geom_point(aes(x=par, y=NEE), col="black")+
  scale_y_continuous(expression(paste("NEE (", mu, "mol", " C", O[2], " ", m^-2, " ", s^-1, ")")))+
  scale_x_continuous(expression(paste("PAR (", mu, "mol ", m^-2, " ", s^-1, ")")))+
  theme_bw()+
  theme()

ggplot(fluxd)+
  geom_point(aes(x=VPD, y=NEE), col="black")+
  scale_y_continuous(expression(paste("NEE (", mu, "mol", " C", O[2], " ", m^-2, " ", s^-1, ")")))+
  scale_x_continuous(expression(paste("VPD (kPa)")))+
  theme_bw()+
  theme()

##hourly average of par

par$hour <- format(par$dt, '%H')
par.hour<-aggregate(par$PARMean, by=list(par$hour), FUN=mean, na.rm=T)
aggregate(par$PARExpUncert, by=list(par$hour==5:10), FUN=mean, na.rm=T)


# NEON_eddy_PULL.R
# Author: A. Kagawa-Viviani
# Date: 10/16/2019 (NEON workshop)
# Description: NEON NEE data pull
#    1) pull nee, par, etc
#    2) stack nee (eddy) and par
#    3) write to .csv 
# Notes: need to identify wd, stations, periods of interest first and provide
#   as arguments to new functions pull2csv.NEE() and pull2csv.PAR()
#   ** I did not attempt to pull multiple sites so not sure how this will
#    perform as a bulk operation. Probably need to specify savepath as pointing 
#    to different folders by site (and make folders?)

library(neonUtilities)
library(rhdf5)  # also library(h5)

setwd("C:/Users/Aurora/OneDrive/NEON/")

## Sparkle will have identified desired sites, periods (start and end),
##    data products (default: eddy bundle, PAR, others TBD)

# eddy: NEON.DP4.00200
pull2csv.NEE<-function(dp="DP4.00200.001",
                       site,      # 4-letter code, i.e. "WREF"
                       startdate, # "YYYY-MM"
                       enddate,   # "YYYY-MM"
                       savepath){ # "Workshop_Day1/WREF/"
  zipsByProduct(dpID=dp, package="basic",
                site=site, 
                startdate=startdate, 
                enddate=enddate,
                savepath=savepath, 
                check.size=T)
  # Eddy flux bundle
  flux<-stackEddy(filepath=paste0(savepath,"filesToStack00200"), level="dp04")
  #return(flux)
  write.csv(flux[[site]], file=paste0(savepath, site, "_",
                                     startdate, "_", enddate, 
                                     ".csv"))
}

# PAR:  NEON.DP1.00024 (1, 300 min avg); NEON.DP2.00005 (gap-filled)
pull2csv.PAR<-function(dp="DP1.00024.001",
                       site,      # 4-letter code, i.e. "WREF"
                       startdate, # "YYYY-MM"
                       enddate,   # "YYYY-MM"
                       savepath){ # "Workshop_Day1/WREF/"
  zipsByProduct(dpID=dp, package="basic",
                site=site, 
                startdate=startdate, 
                enddate=enddate,
                savepath=savepath, 
                check.size=T)
  # DP1
  stackByTable(paste0(savepath,"filesToStack00024"), 
               folder=T, saveUnzippedFiles=T)
}


flux<-pull2csv.NEE(site="GUAN", 
                  startdate="2014-12", 
                  enddate="2019-12",
                  savepath="Workshop_Day1/")
par<-pull2csv.PAR(site="GUAN", 
                  startdate="2014-12", 
                  enddate="2019-12",
                  savepath="Workshop_Day1/")


### Messing around and checking data
par<-read.csv(paste0(savepath,"/filesToStack00024/stackedFiles/PARPAR_30min.csv"),
              stringsAsFactors = F)
timeB <- substring(par$endDateTime, 1, nchar(par$endDateTime)-1) # remove
timeB <- strptime(timeB, format="%Y-%m-%dT%H:%M:%S", tz="GMT")
timeB <- as.POSIXct(timeB)

plot(par$outPARMean~timeB)




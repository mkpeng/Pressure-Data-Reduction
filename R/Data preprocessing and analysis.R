

library(tidyverse)
library(lubridate)


# Process the raw pressure data sampled at 5s and saved as Rdata --------

######## Sampled at 5 seconds as the raw input for the processing. The sampled at 
setwd("F:/CRIO A-Excluded Patients for Mingkai/CRIO A-Excluded Patients for Mingkai/New folder")

filesname_1 <- list.files()
index <- unlist(lapply(filesname_1,nchar))
index <- which(index > 20)
filesname_1 <- filesname_1[index]
n <- length(filesname_1)
filesname_1 <- c(filesname_1[n],filesname_1[1:(n-1)])

id_218_pp <- NULL
for(i in 1:length(filesname_1)){
        new_name <- paste("id_218_",i,".csv",sep="")
        cleanFiles(filesname_1[i],newfile = new_name)
        #### Read the pressure data
}

########## Read the pressure data, combine them and save it for the next stage of use
pressure_data_1 <- read_csv(file = "id_218_1.csv",skip = 8,n_max=5664,col_names = FALSE)
pressure_data_2 <- read_csv(file = "id_218_2.csv",skip = 8,n_max=5664,col_names = FALSE)

id_218_pp <- cbind(pressure_data_1,pressure_data_2)
dim(id_218_pp)
save(id_218_pp,file = "Pressure data/id_218_pressure_data.RData")

############ calculate the correlation coefficient at the different sampling frequency
setwd("F:/CRIO A-Excluded Patients for Mingkai/Paper draft on pressure data/Data/Pressure data")
load("id_218_pressure_data.RData")

id_218_corr <- NULL
sample_freq <- c(5,30,60,120,240)
for(i in sample_freq){
        max_col <- dim(id_218_pp)[2]
        index <- seq(1,max_col,i/5)
        data_sub <- id_218_pp[,index]
        id_218 <- image_corr(batch = data_sub,threshold = 0.95,cut_pp = 5)
        id_218$freq <- i
        id_218_corr <- rbind(id_218_corr,id_218)
}

save(id_218_corr,file = "id_218_corr.RData")

####### The same process is repeated for all the patients

# The time not on bed and active in a selected 24 hours period ------------------------------------------

####################  ID: 217
load("timestamp/id_270.RData")
id_270_timestap$seq <- 1:dim(id_270_timestap)[1]
### load the correlation coefficent
load("id_270.RData")
table(id_270_corr$freq)


id_270_corr <- inner_join(id_270_corr,id_270_timestap)
id_270_corr <- filter(id_270_corr,freq %in% c(30,60,120))

######### Set the 24 hour period 
date1 <- as.POSIXct("2016-02-01 16:00:00",tz="UTC")
date2 <- as.POSIXct("2016-02-02 16:00:00",tz ="UTC")
int <- interval(date1, date2)
index <- id_270_corr$time_final %within% int
id_270_corr <- id_270_corr[index,]

output_270 <- frame_summary(data_corr=id_270_corr,static_period=120,threshold=c(seq(0.75,0.98,0.01)),id="270")

time_onbed <- output_270$time_onbed
time_not_onbed <- output_270$time_not_onbed

ggplot(time_onbed,aes(x=corr_threshold,y=time/3600,group=freq,color=as.factor(freq))) +geom_line()+geom_point()
ggplot(time_onbed,aes(x=corr_threshold,y=num_position,group=freq,color=as.factor(freq))) +geom_line()+geom_point()



###### ID: 218
###### ID: 218
load("timestamp/id_218.RData")
dim(id_218_timestap)
summary(id_218_timestap)

load("Pressure data/id_218_corr.Rdata")
dim(id_218_corr)
table(id_218_corr$freq)
id_218_corr$seq <- (id_218_corr$seq-1)*id_218_corr$freq+1

id_218_corr <- inner_join(id_218_corr,id_218_timestap)
id_218_corr <- filter(id_218_corr,freq %in% c(30,60,120))

######### Set the 24 hour period - ID: 218
date1 <- as.POSIXct("2015-11-09 15:00:00",tz="UTC")
date2 <- as.POSIXct("2015-11-10 15:00:00",tz ="UTC")
int <- interval(date1, date2)
index <- id_218_corr$time_final %within% int
id_218_corr <- id_218_corr[index,]


output_218 <- frame_summary(data_corr=id_218_corr,static_period=120,threshold=c(seq(0.75,0.98,0.01)),id="218")

time_onbed <- output_218$time_onbed
time_not_onbed <- output_218$time_not_onbed

ggplot(time_onbed,aes(x=corr_threshold,y=time/3600,group=freq,color=as.factor(freq))) +geom_line()+geom_point()
ggplot(time_onbed,aes(x=corr_threshold,y=num_position,group=freq,color=as.factor(freq))) +geom_line()+geom_point()


###### ID: 146
load("timestamp/id_146.RData")
dim(id_146_timestap)

load("Pressure data/id_146_corr.Rdata")
dim(id_146_corr)
id_146_corr$seq <- (id_146_corr$seq-1)*id_146_corr$freq+1

id_146_corr <- inner_join(id_146_corr,id_146_timestap)
id_146_corr <- filter(id_146_corr,freq %in% c(30,60,120))

######### Set the 24 hour period - ID: 146
date1 <- as.POSIXct("2015-08-01 15:00:00",tz="UTC")
date2 <- as.POSIXct("2015-08-02 15:00:00",tz ="UTC")
int <- interval(date1, date2)
index <- id_146_corr$time_final %within% int
id_146_corr <- id_146_corr[index,]

output_146 <- frame_summary(data_corr=id_146_corr,static_period=120,threshold=c(seq(0.75,0.98,0.01)),id="146")

time_onbed <- output_146$time_onbed
time_not_onbed <- output_146$time_not_onbed

ggplot(time_onbed,aes(x=corr_threshold,y=time/3600,group=freq,color=as.factor(freq))) +geom_line()+geom_point()
ggplot(time_onbed,aes(x=corr_threshold,y=num_position,group=freq,color=as.factor(freq))) +geom_line()+geom_point()




###### ID: 420
###### ID: 420
load("timestamp/id_420.RData")
dim(id_420_timestap)

load("Pressure data/id_420_corr.Rdata")
dim(id_420_corr)
table(id_420_corr$freq)
id_420_corr$seq <- (id_420_corr$seq-1)*id_420_corr$freq+1

id_420_corr <- inner_join(id_420_corr,id_420_timestap)
id_420_corr <- filter(id_420_corr,freq %in% c(30,60,120))

######### Set the 24 hour period - ID: 420
date1 <- as.POSIXct("2016-11-08 15:00:00",tz="UTC")
date2 <- as.POSIXct("2016-11-09 15:00:00",tz ="UTC")
int <- interval(date1, date2)
index <- id_420_corr$time_final %within% int
id_420_corr <- id_420_corr[index,]


output_420 <- frame_summary(data_corr=id_420_corr,static_period=120,threshold=c(seq(0.75,0.98,0.01)),id="420")

time_onbed <- output_420$time_onbed
time_not_onbed <- output_420$time_not_onbed

ggplot(time_onbed,aes(x=corr_threshold,y=time/3600,group=freq,color=as.factor(freq))) +geom_line()+geom_point()
ggplot(time_onbed,aes(x=corr_threshold,y=num_position,group=freq,color=as.factor(freq))) +geom_line()+geom_point()

########## Final output to save all the summary data on time on bed and 
save(output_270,output_146,output_218,output_420,file = "4p_summary.RData")

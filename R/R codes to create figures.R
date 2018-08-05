
###### Clean the working space nad 
rm(list=ls(all=TRUE))
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(lubridate)

# 2. Changes image: where the correlation coefficient shows the different positions --------

### Use the data from ID: 270 for create the plot
setwd("F:/CRIO A-Excluded Patients for Mingkai/Paper draft on pressure data/Data")

### Read the coordinate for the data points
coordinate <- read.csv("cooridnate.csv")

###### Load the data
### Two files: raw data: pp_data, and corr file: id_270_pp
##### Procssed data for patient id: 270
### pp_data: pressure data for each frames
### id_270_pp: Output from the correlation evaluation between frames
load("id_270_pressure_data.RData")
### Frame used for illustration
id_plot <- c(45000,45500,4500)/5
id_plot

section2 <- pp_data[,id_plot]
colnames(section2) <- paste("frame",1:3,sep="")

### Print the correlation coefficient between frames
for(i in 2:3){
  temp <- section2[,c(1,i)]
  temp <- temp[temp>5,]
  print(cor.test(temp[,1],temp[,2]))
}


### from wide to long
section2 <- bind_cols(section2,coordinate)
section2_long <- gather(section2, frame, pressure, frame1:frame3)
section2_long <- filter(section2_long,pressure>5)

frame_name <- c("frame1" ="Frame 1","frame2" = "Frame 2","frame3" = "Frame 3")
#### Heat map plot to show the pressure value distribution
tiff(file = "Figure 2 a .tiff", width = 6, height = 6.5, units = "in", res = 300,compression = "lzw")
ggplot(section2_long,aes(x=axis_x,y=axis_y,col=pressure))+geom_point()+
  scale_y_continuous(breaks=c(0,12,24,36,48),limits = c(0,48),
                     labels = c(0,12,24,36,48))+
  scale_x_continuous(breaks=c(-120,-100,-80,-60,-40,-20,0),limits = c(-120,0))+
  facet_wrap(~frame,ncol=1,labeller = as_labeller(frame_name))+theme_minimal()+
  theme(panel.grid.major = element_line(colour = "gray80"))
dev.off()


############ Scatter plot between frame 1 and frame 2
temp <- section2[,]
temp <- temp[temp>5,]
cor.test(temp[,1],temp[,2])
temp <- filter(temp,!is.na(temp$frame1))
tiff(file = "Figure 2 b .tiff", width = 2.5, height = 2.5, units = "in", res = 300,compression = "lzw")
ggplot(temp,aes(x=temp[,2],y=temp[,1]))+geom_point(size=0.1)+
  geom_abline(slope = 1,intercept = 0)+
  xlim(0,60)+ylim(0,60)+theme_bw()+xlab("Frame 2")+ylab("Frame 1")+
  annotate("text", x = 20, y = 60, label = "paste(rho, \"=0.977\")", parse = TRUE)
dev.off()


##########  Correlation plot between Frame 2 and 3
temp <- section2
temp <- temp[temp[,2]>5&temp[,3],]
cor.test(temp[,3],temp[,2])
temp <- filter(temp,!is.na(temp$frame2))

tiff(file = "Figure 2 c .tiff", width = 2.5, height = 2.5, units = "in", res = 300,compression = "lzw")
ggplot(temp,aes(x=temp[,2],y=temp[,3]))+geom_point(size=0.1)+
  geom_abline(slope = 1,intercept = 0)+
  xlim(0,60)+ylim(0,60)+theme_bw()+xlab("Frame 2")+ylab("Frame 3")+
  annotate("text", x = 20, y = 60, label = "paste(rho, \"=0.332\")", parse = TRUE)
dev.off()

# 3. Image comparison at different correlation threshold ------------------
# load("id_270_pressure_data.RData")

tiff(file = "Correlation plot-0.99.tiff", width = 10, height = 2.7,
     units = "in", res = 300,compression = "lzw")
ggplot(frame_similar,aes(x=axis_x,y=axis_y,col=log(m1)))+geom_point()+
  scale_y_continuous(breaks=c(0,12,24,36,48),limits = c(0,48),
                     labels = c(0,12,24,36,48))+
  scale_x_continuous(breaks=c(-120,-100,-80,-60,-40,-20,0),limits = c(-120,0))+
  facet_wrap(~type)+theme_minimal()+
  theme(panel.grid.major = element_line(colour = "gray80"))
dev.off()






# 4. Boxplot: show the changes of cell difference and pressure difference  ----------------

############## relationship of correlation coefficient with other statistics to describe frame similarity
load("id_270.RData")
########### Results from the different sampling frequency
table(id_270_corr$freq)

### Sampling frequency at 5s
id_270_corr_freq5 <- filter(id_270_corr,freq==5)
id_270_corr_freq5$corr_group <- floor(id_270_corr_freq5$corr_coef*100)

id_270_corr_freq5$corr_group_cluster <- ifelse(id_270_corr_freq5$corr_group< 85,84,id_270_corr_freq5$corr_group)
id_270_corr_freq5 <- filter(id_270_corr_freq5,total_cell > 500)

tiff(file = "Figure 4(a) - Box plot on cell difference.tiff", width = 7, height = 3, units = "in", res = 500,compression = "lzw")
ggplot(id_270_corr_freq5,aes(x=as.factor(corr_group_cluster/100),y=diff_cell))+
  geom_boxplot(outlier.size = 0.2)+
  xlab("Correlation coefficient")+
  ylab("Sum of cell difference \n between frames")+theme_bw()
dev.off()

tiff(file = "Figure 4(b) - Box plot on pressure difference.tiff", width = 7, height = 3, units = "in", res = 500,compression = "lzw")
ggplot(id_270_corr_freq5,aes(x=as.factor(corr_group_cluster/100),y=diff_sum))+geom_boxplot(outlier.size = 0.2)+
  xlab("Correlation coefficient")+
  ylab("Sum of absolute pressure difference \n between frames")+theme_bw()
dev.off()


# 6. Upside down: ECG to show the sampling frequency  ---------------------

id_270_corr <- NULL
sample_freq <- c(5,10,30,60,120,240)
sample_freq <- 60

for(i in sample_freq){
  max_col <- dim(pp_data)[2]
  index <- seq(1,max_col,i/5)
  data_sub <- pp_data[,index]
  id_270 <- image_corr(batch = data_sub,threshold = 0.95,cut_pp = 5)
  id_270$freq <- i
  id_270_corr <- rbind(id_270_corr,id_270)
}

save(id_270_corr,file = "id_270.RData")

id_270_corr$seq <- (id_270_corr$seq-1)*id_270_corr$freq+1
id_270_corr_at60s <- id_270_corr

save(id_270_corr_at60s,file = "id_270_at60s.RData")

#################### load the timestamp at every second
load("timestamp/id_270.RData")
id_270_timestap$seq <- 1:dim(id_270_timestap)[1]

id_270_corr_1 <- inner_join(id_270_corr,id_270_timestap)

######### Set the 24 hour period
library(lubridate)

date1 <- as.POSIXct("2016-02-01 15:00:00",tz="UTC")
date2 <- as.POSIXct("2016-02-02 15:00:00",tz ="UTC")
int <- interval(date1, date2)
index <- id_270_corr_1$time_final %within% int

id_270_corr_sub <- id_270_corr_1[index,]

library(scales)

id_270_corr_sub <- filter(id_270_corr_sub,freq!=10)
group_name <-  c("5"="Sampled at every 5s", "30"="Sampled at every 30s","60"="Sampled at every 60s", "120"="Sampled at every 120s",
                                    "240"="Sampled at every 240s")

tiff(file = "Figure 6 - ECG plot_test.tiff", width = 24, height = 12, units = "in", res = 500,compression = "lzw")
ggplot(id_270_corr_sub,aes(x=time_final,y=corr_coef))+geom_line(size = 0.5)+
  facet_wrap(~freq,ncol=1,labeller = as_labeller(group_name))+theme_minimal()+
  scale_x_datetime(breaks = seq(date1,date2, "6 hours"),labels = date_format("%m/%d %H:%M"))+
  scale_y_continuous(breaks=seq(0,1,0.25))+xlab("Time")+ylab("Correlation coefficient")+
  theme(panel.grid.minor=element_blank(),text = element_text(size=30)) + scale_size_manual(values=0.1, guide=FALSE)
dev.off()




# Figure on comparing correlation coefficient with other statistic --------


setwd("/Volumes/MINGKAI 1/Paper draft on pressure data/Data")

### load the timestamp
load("timestamp/id_270.RData")
id_270_corr_at60s <- filter(id_270_corr,freq==60)
id_270_timestap$seq <- 1:dim(id_270_timestap)[1]

id_270_corr_at60s <- inner_join(id_270_corr_at60s,id_270_timestap)

######### Set the 24 hour period
date1 <- as.POSIXct("2016-02-01 16:00:00",tz="UTC")
date2 <- as.POSIXct("2016-02-02 16:00:00",tz ="UTC")
int <- interval(date1, date2)
index <- id_270_corr_at60s$time_final %within% int

id_270_corr_sub <- id_270_corr_at60s[index,]
id_270_corr_sub <- id_270_corr_sub %>% mutate(mean = total_pressure/total_cell) %>%
        select(seq,time_final,mean,total_cell,number_cell,corr_coef)
### from wide to long

id_270_corr_sub <- gather(id_270_corr_sub,type,value,mean:corr_coef)


ggplot(id_270_corr_sub,aes(x=time_final,y=value))+geom_line()+facet_wrap(~type,scales="free_y",ncol = 1)


group_name <-  c("mean"="Mean pressure value",
                 "total_cell"="Number of activated sensors",
                 "number_cell"="Number of sensors with pressure value > 40mmHg",
                 "120"="Sampled at every 120s",
                 "corr_coef"="Correlation coefficient")

library(scales)
tiff(file = "Figure 8 - ECG plot_test.tiff", width = 24, height = 12, units = "in", res = 500,compression = "lzw")
ggplot(id_270_corr_sub,aes(x=time_final,y=value))+geom_line(size = 0.5)+
        facet_wrap(~type,ncol=1,labeller = as_labeller(group_name),scales = "free_y")+theme_minimal()+
        scale_x_datetime(breaks = seq(date1,date2, "6 hours"),labels = date_format("%m/%d %H:%M"))+
        xlab("Time")+ylab("")
        theme(panel.grid.minor=element_blank(),text = element_text(size=30)) + scale_size_manual(values=0.1, guide=FALSE)
dev.off()


jpeg(file = "Figure 8 - ECG plot_test.jpeg", width = 20, height = 10, units = "in", res = 500)
ggplot(id_270_corr_sub,aes(x=time_final,y=value))+geom_line(size = 0.5)+
        facet_wrap(~type,ncol=1,labeller = as_labeller(group_name),scales = "free_y")+theme_minimal()+
        scale_x_datetime(breaks = seq(date1,date2, "6 hours"),labels = date_format("%m/%d %H:%M"))+
        xlab("Time")+ylab("")
        theme(panel.grid.minor=element_blank(),text = element_text(size=30)) + scale_size_manual(values=0.1, guide=FALSE)
dev.off()


# time period not on the bed and active -------------------------------------------------------



#### Time not on bed

load("4p_summary.RData")
time_not_onbed <- bind_rows(output_270$time_not_onbed,output_146$time_not_onbed,output_218$time_not_onbed,output_420$time_not_onbed)
time_not_onbed$id <- paste("Patient ID:",time_not_onbed$id)
time_not_onbed$freq <- paste(time_not_onbed$freq,"s",sep = "")
time_not_onbed$freq <- factor(time_not_onbed$freq,levels=c("30s","60s","120s"))

tiff(file = "Figure 7 a - time not on bed.tiff", width = 6, height = 2.5, units = "in", res = 500,compression = "lzw")
ggplot(data=time_not_onbed, aes(x=id, y=time_out/60, fill=as.factor(freq))) +
  geom_bar(stat="identity", position=position_dodge())+xlab("Patients") +ylab("Time not on bed (minutes)")+
  theme_bw()+theme(legend.title=element_blank())
dev.off()


###### Time on bed

time_onbed <- bind_rows(output_270$time_onbed,output_146$time_onbed,output_218$time_onbed,output_420$time_onbed)
time_onbed$id <- paste("Patient ID:",time_onbed$id)
time_onbed$freq <- paste(time_onbed$freq,"s",sep = "")
time_onbed$freq <- factor(time_onbed$freq,levels=c("30s","60s","120s"))

tiff(file = "Figure 7 b - time being active - 24h.tiff", width = 7, height = 3, units = "in", res = 500,compression = "lzw")
ggplot(time_onbed,aes(x=corr_threshold,y=time/3600,group=freq,color=freq))+geom_line()+
  facet_wrap(~id,ncol = 2)+
  theme_bw()+xlab("Threshold of correlation coefficient to detect position change")+
  ylab("Total time with position change (minutes)")+theme(legend.title=element_blank())
dev.off()


tiff(file = "Figure 7 c - number of position change - 24hs.tiff", width = 7, height = 3, units = "in", res = 500,compression = "lzw")
ggplot(time_onbed,aes(x=corr_threshold,y=num_position,group=freq,color=freq,linetype=freq))+geom_line()+facet_wrap(~id,ncol = 2)+
  theme_bw()+xlab("Threshold of correlation coefficient to detect position change")+
  ylab("Number of position changes")+theme(legend.title=element_blank())
dev.off()




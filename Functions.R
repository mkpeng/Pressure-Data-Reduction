
# Functions to be loaded before the processing ----------------------------


# Pre-Process the data each frame ---------------------------------------------
cleanFiles<-function(file,newfile){
        writeLines(iconv(readLines(file,skipNul = TRUE)),newfile)
}


# Pressure Correlation Function to compare adjcent pressure data ----------

##### Correlation between adjacent frames

image_corr <- function(batch=patch,threshold=0.9,cut_pp=5){
        number_images <- dim(batch)[2] - 1
        mat_cor <- matrix(0,number_images,7)
        for(i in 1:number_images){
                temp <- cbind(batch[,i],batch[,(i+1)])
                temp <- data.frame(temp)
                colnames(temp) <- c("m1","m2")
                temp_1 <- filter(temp,m1 > cut_pp |m2 > cut_pp)
                ### correlation coefficient 
                mat_cor[i,2] <- cor(temp$m1,temp$m2)
                ### Difference of pressure between two frames
                mat_cor[i,3]<- sum(abs(temp_1$m1-temp_1$m2))
                temp <- cbind(temp,coordinate)
                #### Total pressure on the mattress
                #### Number of number being contacted
                mat_cor[i,6] <- sum(temp$m1 >=cut_pp)
                #### Difference of value between two frames
                temp_2 <- filter(temp,m1 >= cut_pp & m2 >= cut_pp)
                temp_1 <- filter(temp,m1 > cut_pp |m2 > cut_pp)
                mat_cor[i,5] <- sum(temp_1$m1)
                
                #### the number of cell difference with a threshold of 5
                mat_cor[i,4]<- dim(temp_1)[1] - dim(temp_2)[1]
                mat_cor[i,7] <- dim(temp_1)[1]
        }
        
        mat_cor[,1]<- 1:number_images
        yes_no_1f <- ifelse(mat_cor[,2] > threshold,1,0)
        mat_cor <- cbind(mat_cor,yes_no_1f)
        mat_cor <- data.frame(mat_cor)
        colnames(mat_cor) <- c("seq","corr_coef","diff_sum","diff_cell","total_pressure","total_cell","number_cell","corr_cut")
        return(mat_cor)
}



# frame_summary function --------------------------------------------------

################ This function summarize the time not on bed and time being active with continuous movment
########## Input:
#### data_corr: correlation coefficient data at the different sampling frequencies
#### static_period: if no movement within 120, it is defined as statitic period
#### threshold: threshold of correlaiton coefficient used to determine which there is movement between frames
#### id: patient identifier


frame_summary <- function(data_corr=id_270_corr,static_period=120,threshold=0.95,id="id_270"){
        data_corr <- arrange(data_corr,freq,seq)
        ### period not on bed
        data_corr_notonbed <- filter(data_corr,total_cell <= 500)
        data_corr_notonbed <- data_corr_notonbed %>% group_by(freq) %>% mutate(time_lag = lag(time_final))
        index <- which(is.na(data_corr_notonbed$time_lag))
        data_corr_notonbed$time_lag[index] <- data_corr_notonbed$time_final[index]
        data_corr_notonbed$time_diff <- as.numeric(difftime(data_corr_notonbed$time_final,data_corr_notonbed$time_lag,units="secs"))
        
        data_corr_notonbed <- mutate(data_corr_notonbed,time_diff_new =ifelse(time_diff > freq,0,time_diff),
                                     new_position= ifelse(time_diff_new==0,1,0))
        data_corr_notonbed <- data_corr_notonbed %>% group_by(freq) %>% mutate(new_position = cumsum(new_position))
        
        ### time period 1 
        time_not_onbed_1 <- data_corr_notonbed %>% group_by(freq) %>% summarise(notonbed = sum(time_diff_new))
        ### add the boundary time
        time_not_onbed_2 <- data_corr_notonbed %>% group_by(freq,new_position) %>% summarise(times = n())
        time_not_onbed_2 <- mutate(time_not_onbed_2,time_2 = ifelse(times > 1,freq,freq/2))
        time_not_onbed_2 <- time_not_onbed_2 %>% group_by(freq) %>% summarise(time_2 = sum(time_2))
        
        ## Total time
        time_not_onbed <- inner_join(time_not_onbed_1,time_not_onbed_2)
        time_not_onbed <- mutate(time_not_onbed,time_out = notonbed+ time_2)
        
        ### time being active, moving 
        time_onbed <- NULL
        for (i in threshold){
                for(j in static_period){
                        ####  select the position below threshold and positon on the bed
                        data_corr_onbed <- filter(data_corr,total_cell > 500,corr_coef < i)
                        
                        data_corr_onbed <- data_corr_onbed %>% group_by(freq) %>% mutate(time_lag = lag(time_final))
                        index <- which(is.na(data_corr_onbed$time_lag))
                        data_corr_onbed$time_lag[index] <- data_corr_onbed$time_final[index]
                        data_corr_onbed$time_diff <- as.numeric(difftime(data_corr_onbed$time_final,data_corr_onbed$time_lag,units="secs"))
                        #  data_corr_onbed <- filter(data_corr_onbed,time_diff <= freq+1)
                        #### static period might be smaller than sampling frequency
                        index_2 <- (data_corr_onbed$freq+1) ==data_corr_onbed$time_diff
                        data_corr_onbed$time_diff[index_2] <- data_corr_onbed$freq[index_2]
                        
                        data_corr_onbed <- mutate(data_corr_onbed,position_change = ifelse(time_diff > j & time_diff > freq,1,0))
                        
                        data_corr_onbed <- mutate(data_corr_onbed,time_diff=ifelse(time_diff > j & time_diff > freq,0,time_diff))
                        data_corr_onbed$position_change[index] <- 1
                        data_corr_onbed <- mutate(data_corr_onbed,position_change = cumsum(position_change))
                        
                        data_corr_onbed_summary <- data_corr_onbed %>% group_by(freq,position_change) %>% 
                                summarise(time_1 = sum(time_diff),times=n())
                        data_corr_onbed_summary <- mutate(data_corr_onbed_summary,time_1 = ifelse(times > 1, time_1+freq,time_1+freq/2))
                        data_corr_onbed_final <- data_corr_onbed_summary %>% group_by(freq) %>% summarise(num_position=max(position_change),time = sum(time_1))
                        data_corr_onbed_final$corr_threshold <- i
                        data_corr_onbed_final$static_period <- j
                        time_onbed <- rbind(time_onbed,data_corr_onbed_final)
                }
        }
        time_onbed$id <- id
        #time_not_onbed$id <- id
        output <- list(time_onbed = time_onbed,time_not_onbed = time_not_onbed)
        return(output)
}



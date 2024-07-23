get_needed_information <- function(){
  area <<- readline(prompt = "What is the area you want to explore? Belgian part of the North Sea, 
                  Northern Adriatic Sea or another place? \n possible answers: bpns, nas, other ")
  
  if (area == "bpns"){
    # Load LifeWatch (LW) data with LW package
    station_region <<- readline("What is your region of interest? nearshore, midshore or offshore?  ")
    station_code <<- readline("What is de code of your station? For example 130 ")
    station_list <- c("LW01", "LW02", "435", "W07bis", "W08", "W09", "W10", "421", "130", "700",
                        "780", "330", "230", "710", "ZG02", "120", "215")
    station_coord <- data.frame("lat" = c(2.256, 2.556, 2.7903, 3.0125, 2.35, 2.7, 2.4167, 
                                      2.45, 2.9054, 3.221, 3.0573, 2.8091, 2.8504, 3.1383, 2.5007, 2.7025, 2.6108),
                            "long" = c(51.5687 , 51.8, 51.5807, 51.588, 51.4583, 51.75,
                                       51.6833, 51.4805, 51.2706, 51.377, 51.4714, 51.4341, 51.3087, 51.4412, 51.3352, 51.1861, 51.2749))
                            
    if (station_code %in% station_list){
               station_lat <<- station_coord$lat[which(station_list %in% station_code)]       # coordinates station: latitude 
               station_lon <<- station_coord$long[which(station_list %in% station_code)]       # coordinates station: longitude
     }
    startdate <<- as.Date(readline("What is first date of your data? format YYYY-MM-DD  "))
    stopdate <<- as.Date(readline("What is last date of your data? format YYYY-MM-DD  "))  
  }
  if (area == "nas"){
    # Load LifeWatch (LW) data with LW package
    station_region <<- readline("What is your region of interest? nearshore, midshore or offshore?  ")
    station_code <<- readline("What is de code of your station? For example 130 ")
    startdate <<- as.Date(readline("What is first date of your data? format YYYY-MM-DD  "))
    stopdate <<- as.Date(readline("What is last date of your data? format YYYY-MM-DD  "))
    station_lat <<- as.numeric(readline("Please provide coordinates for your station: Latitude:  ")) # coordinates station: latitude 
    station_lon <<- as.numeric(readline("Please provide coordinates for your station: Longitude:  ")) # coordinates station: longitude
  }
  if (area == "other"){
    # Load LifeWatch (LW) data with LW package
    station_region <<- readline("What is your region of interest? nearshore, midshore or offshore?  ")
    station_code <<- readline("What is de code of your station? For example 130 ")
    startdate <<- as.Date(readline("What is first date of your data? format YYYY-MM-DD  "))
    stopdate <<- as.Date(readline("What is last date of your data? format YYYY-MM-DD  "))
    station_lat <<- as.numeric(readline("Please provide coordinates for your station: Latitude:  ")) # coordinates station: latitude 
    station_lon <<- as.numeric(readline("Please provide coordinates for your station: Longitude:  ")) # coordinates station: longitude
  }
    obs_years <<- seq(lubridate::year(startdate), lubridate::year(stopdate), 1)
    startdate_obs_cal <<- as.Date(paste0(obs_years[1], "-01-01"))
    stopdate_obs_cal <<- as.Date(paste0(obs_years[round((length(obs_years)/3)*2)-1], "-12-31"))
    startdate_obs_val <<- as.Date(paste0(obs_years[round((length(obs_years)/3)*2)], "-01-01"))
    stopdate_obs_val <<- as.Date(paste0(obs_years[length(obs_years)], "-12-31"))
}

number_of_iterations <- function(){
  numSimulations <<- as.numeric(readline(prompt = "How many iterations do you want to run?"))
}

get_inputdata <- function(station_region, station_selection, startdate, stopdate){
  ## GAM models
  # Load packages for creating GAM models
  packages <- list("mgcv",
          "lubridate",
          "ggplot2",
          "stats",
          "xts",
          "lwdataexplorer",
          "plyr",
          "dplyr")
  lapply(packages, library, character.only=TRUE)
  
  # Load functions for creating GAM models
  source(paste0(wd,"Rscripts/GAM/functions_GAM_v1_SP.R"))
  source(paste0(wd,"Rscripts/GAM/LW_DataAccess_functions_for_GAM.R"))


  
  
  days_run = as.numeric(difftime(as.Date(stopdate) , as.Date(startdate))) + 1 # how many days are in the run = diff (start - stop) + 1
 
 
  
        
  station_data <- get_abiotic_pigment(start = startdate, stop = stopdate, station_selection = station_selection, station_region = station_region)

  #This line is for visualization purposes of the performance of the models
  #par(mfrow=c(2,2))
  
  #GAM PO4
  #The selected model has to be defined as mod_p_select
  ####################################################################
  
  po4 <- station_data[!is.na(station_data$PO4),]
  day_station<-yday(po4$Date)
  year_station<-lubridate::year(po4$Date)
  
  models <- list()
  models[[1]] <- test_model(po4$PO4, day_station, year_station, 3, 3)
  models[[2]] <- test_model(po4$PO4, day_station, year_station, 4, 4)
  models[[3]] <- test_model(po4$PO4, day_station, year_station, 5, 5)
  models[[4]] <- test_model(po4$PO4, day_station, year_station, 6, 6)
  AICs <- c(models[[1]]$aic,
            models[[2]]$aic,
            models[[3]]$aic,
            models[[4]]$aic)
  
  mod_p_select <- models[[which.min(AICs)]]
  
  #Comparison observations vs. predictions (by points)
  #PO4
#  plot_po4_1 <- comparisonPoint(day_station, year_station, po4$PO4, mod_p_select)
#  plot_po4_1 + ylim(0,max(po4$PO4)) + xlim(0,max(po4$PO4)) + labs(x = "PO4 observed", y = "PO4 predicted") 
 
  #GAM NH4
  #The selected model has to be defined as mod_nh4_select
  ####################################################################
  nh4 <- station_data[!is.na(station_data$NH4),]
  day_station<-yday(nh4$Date)
  year_station<-lubridate::year(nh4$Date)
  
  models <- list()
  models[[1]] <- test_model(nh4$NH4, day_station, year_station, 3, 3) 
  models[[2]] <- test_model(nh4$NH4, day_station, year_station, 4, 4)
  models[[3]] <- test_model(nh4$NH4, day_station, year_station, 5, 5)
  models[[4]] <- test_model(nh4$NH4, day_station, year_station, 6, 6)
  AICs <- c(models[[1]]$aic,
            models[[2]]$aic,
            models[[3]]$aic,
            models[[4]]$aic)
  
  mod_nh4_select <- models[[which.min(AICs)]]
  
  #Comparison observations vs. predictions
#  plot_nh4_1 <- comparisonPoint(day_station, year_station, nh4$NH4, mod_nh4_select)
#  plot_nh4_1 + ylim(0,6) + xlim(0,6) + labs(x = "NH4 observed", y = "NH4 predicted") 

  #GAM NO2
  #The selected model has to be defined as mod_no2_select
  ####################################################################
  no2 <- station_data[!is.na(station_data$NO2),]
  day_station<-yday(no2$Date)
  year_station<-lubridate::year(no2$Date)
  
  models <- list()
  models[[1]] <- test_model(no2$NO2, day_station, year_station, 3, 3)
  models[[2]] <- test_model(no2$NO2, day_station, year_station, 4, 4) 
  models[[3]] <- test_model(no2$NO2, day_station, year_station, 5, 5) 
  models[[4]] <- test_model(no2$NO2, day_station, year_station, 6, 6)
  AICs <- c(models[[1]]$aic,
            models[[2]]$aic,
            models[[3]]$aic,
            models[[4]]$aic)
  
  mod_no2_select <- models[[which.min(AICs)]]
  
  #Comparison observations vs. predictions
#  plot_no2_1 <- comparisonPoint(day_station, year_station, no2$NO2, mod_no2_select)
#  plot_no2_1 + ylim(0,1) + xlim(0,1) + labs(x = "NO2 observed", y = "NO2 predicted") 

  #GAM NO3
  #The selected model has to be defined as mod_no3_select
  ####################################################################
  no3 <- station_data[!is.na(station_data$NO3),]
  day_station<-yday(no3$Date)
  year_station<-lubridate::year(no3$Date)
  
  models <- list()
  models[[1]] <- test_model(no3$NO3, day_station, year_station, 3, 3)
  models[[2]] <- test_model(no3$NO3, day_station, year_station, 4, 4)
  models[[3]] <- test_model(no3$NO3, day_station, year_station, 5, 5) 
  models[[4]] <- test_model(no3$NO3, day_station, year_station, 6, 6)
  AICs <- c(models[[1]]$aic,
            models[[2]]$aic,
            models[[3]]$aic,
            models[[4]]$aic)
  
  mod_no3_select <- models[[which.min(AICs)]]
  
  #Comparison observations vs. predictions
#  plot_no3_1 <- comparisonPoint(day_station, year_station, no3$NO3, mod_no3_select)
#  plot_no3_1 + ylim(0,20) + xlim(0,20) + labs(x = "NO3 observed", y = "NO3 predicted") 

  #GAM for Silica
  #The selected model has to be defined as mod_si_select
  #####################################################################################
  sil <- station_data[!is.na(station_data$SiO4),]
  day_station<-yday(sil$Date)
  year_station<-lubridate::year(sil$Date)
  
  models <- list()
  models[[1]] <- test_model(sil$SiO4, day_station, year_station, 3, 3)
  models[[2]] <- test_model(sil$SiO4, day_station, year_station, 4, 4)
  models[[3]] <- test_model(sil$SiO4, day_station, year_station, 5, 5)
  models[[4]] <- test_model(sil$SiO4, day_station, year_station, 6, 6) 
  AICs <- c(models[[1]]$aic,
            models[[2]]$aic,
            models[[3]]$aic,
            models[[4]]$aic)
  
  mod_si_select <- models[[which.min(AICs)]]
  
  #Comparison per point
  #SiO4
#  plot_si_1 <- comparisonPoint(day_station, year_station, sil$SiO4, mod_si_select)
#  plot_si_1 + ylim(0,6) + xlim(0,6) + labs(x = "SiO4 observed", y = "SiO4 predicted") 

  #GAM for Temperature
  #The selected model has to be defined as mod_temp_select
  #####################################################################################
  if (station_region[1] != "offshore"){
    temperature <- station_data[!is.na(station_data$Temperature),]
    day_station<-yday(temperature$Date)
    year_station<-lubridate::year(temperature$Date)
    
    models <- list()
    models[[1]] <- test_model(temperature$Temperature, day_station, year_station, 3, 3)
    models[[2]] <- test_model(temperature$Temperature, day_station, year_station, 4, 4) 
    models[[3]] <- test_model(temperature$Temperature, day_station, year_station, 5, 5)
    models[[4]] <- test_model(temperature$Temperature, day_station, year_station, 6, 6) 
    models[[5]] <- test_model(temperature$Temperature, day_station, year_station, 7, 7)
    AICs <- c(models[[1]]$aic,
              models[[2]]$aic,
              models[[3]]$aic,
              models[[4]]$aic,
              models[[5]]$aic)
    
    mod_temp_select <- models[[which.min(AICs)]]
    
    #Comparison per point
    #temperature
    #plot_temp <- comparisonPoint(day_station, year_station, temperature$Temperature, mod_temp_select)
    #plot_temp + ylim(0,25) + xlim(0,25) + labs(x = "Temperature observed", y = "Temperature predicted") 
  }

  
  ## Create daily time series
## Create prediction for NPZD model
##dummy years from 2011 to 2013 will be a copy of the data from 2014 or 2012 depending on availability of data
######################################################################################
predictions <- data.frame(matrix(nrow = days_run, ncol = 8))

names(predictions)[1] <- "day"
names(predictions)[2] <- "year"
names(predictions)[3] <- "DIN"
names(predictions)[4] <- "Temp"
names(predictions)[5] <- "po4"
names(predictions)[6] <- "sio4"
names(predictions)[7] <- "Date"
names(predictions)[8] <- "month"



#Temperature data from MVB



if (station_region[1] == "offshore"){
  temp_data  <- get_MVB_SST(start = startdate, stop = stopdate)
  colnames(temp_data) <- c("Date", "temperature")
  
  #Copy dates from MVB data
  predictions$Date <- seq(ymd(startdate), ymd(stopdate), by = 'day')
  predictions$day <- yday(predictions$Date)
  predictions$year <- lubridate::year(predictions$Date)
  predictions$month <- lubridate::month(predictions$Date)
  
  #Add temperature data from MVB directly, here there are not GAM applied
  predictions$Temp <- temp_data$temperature

  #Add values from GAM models
  ##Dummy years are a copy of 2014, new reference dates data frame is created to make the predictions
  #DIN data are only available form 2014 onward
  if (startdate < as.Date("2014-01-01")){
    temp_predictions <- predictions[,c(1,2,7)]
    temp_predictions$year[c(1:1096)] <- rep(2014,times =1096)
  }

  #DIN
  predictions[3]<-predict.gam(mod_nh4_select, newdata = temp_predictions[,c(1,2)], type = "response") +
  predict.gam(mod_no2_select, newdata = temp_predictions[,c(1,2)], type = "response") +
  predict.gam(mod_no3_select, newdata = temp_predictions[,c(1,2)], type = "response")

  #PO4
  #Observations of PO4 and SiO4 started in 2012 for offshore stations
  #Therefore, the predictions repeat the year 2012 to fill in the year 2011
  temp_predictions <- predictions[,c(1,2,7)]
  temp_predictions$year[c(1:365)] <- rep(2012,times = 365)
  predictions[5]<-predict.gam(mod_p_select, newdata = temp_predictions[,c(1,2)], type = "response") 

  #SiO4
  predictions[6]<-predict.gam(mod_si_select, newdata = temp_predictions[,c(1,2)], type = "response") 

}else{
  predictions$Date <- seq(ymd(startdate), ymd(stopdate), by = 'day')
  predictions$day <- yday(predictions$Date)
  predictions$year <- lubridate::year(predictions$Date)
  predictions$month <- lubridate::month(predictions$Date)
  
  #Add values from GAM models
  
  #DIN
  #There is no data available from 2011 to 2013 for DIN (see station data)
  #Therefore, 2014 is used as reference year for the dummy years
  #It is repeated 1096 times = 3 years * 365 + 1 day of leap year (2012)
  
  #Create a temporal data frame to calculate the predictions based on the best GAM 
  #select day, year and exact date
  if (startdate < as.Date("2014-01-01")){
    temp_predictions <- predictions[,c(1,2,7)]
    temp_predictions$year[c(1:1096)] <- rep(2014,times =1096)
    #Add NH4, NO3, NO2 values based on the exact dates
    predictions[3]<-predict.gam(mod_nh4_select, newdata = temp_predictions[,c(1,2)], type = "response") +
      predict.gam(mod_no2_select, newdata = temp_predictions[,c(1,2)], type = "response") +
      predict.gam(mod_no3_select, newdata = temp_predictions[,c(1,2)], type = "response")
  }

  #Add NH4, NO3, NO2 values based on the exact dates
  predictions[3]<-predict.gam(mod_nh4_select, newdata = predictions[,c(1,2)], type = "response") +
  predict.gam(mod_no2_select, newdata = predictions[,c(1,2)], type = "response") +
  predict.gam(mod_no3_select, newdata = predictions[,c(1,2)], type = "response")


  #Temperature
  #There is temperature data on 2011, therefore it is not necessary to do make a temporal predictions data frame
  predictions[4]<-predict.gam(mod_temp_select, newdata = predictions[,c(1,2)], type = "response") 

  #PO4

  predictions[5]<-predict.gam(mod_p_select, newdata = predictions[,c(1,2)], type = "response") 

  #SiO4

  predictions[6]<-predict.gam(mod_si_select, newdata = predictions[,c(1,2)], type = "response") 
}


##Monthly data visualization
########################################################
#predictions$month <- month(predictions$Date)
#station_data$month <- month(station_data$Date)

#comparison_graphs <- compareGraphs(station_data, predictions, days_run)
#DIN
#comparison_graphs[[1]]

#PO4
#comparison_graphs[[2]]

#SiO4
#comparison_graphs[[3]]

#Let's look at our time series.

#compare predictions as time series
#comparison_time <- comparisonTimeSeries(predictions, station_data)
#DIN
#comparison_time[[1]]

#PO4
#comparison_time[[2]]

#SiO4
#comparison_time[[3]]





# adjust negative DIN values to small concentration.
if (length(which(predictions$DIN < 0 )) != 0) {
neg_values_DIN <- which(predictions$DIN < 0 )

min_values <- station_data$DIN[order(station_data$DIN, decreasing=F)]
min_values <- na.omit(min_values)[1:20]
min_DIN <- min(min_values)
max_DIN <- max(min_values)

predictions$DIN[neg_values_DIN] <- runif(length(neg_values_DIN), min_DIN, max_DIN)
}

# adjust negative values to small concentration.
if (length(which(predictions$po4 < 0 )) != 0) {
  neg_values_po4 <- which(predictions$po4 < 0 )
  
  min_values <- station_data$PO4[order(station_data$PO4, decreasing=F)]
  min_values <- na.omit(min_values)[1:20]
  min_po4 <- min(min_values)
  max_po4 <- max(min_values)
  
  predictions$po4[neg_values_po4] <- runif(length(neg_values_po4), min_po4, max_po4)
}

# adjust negative values to small concentration.
if (length(which(predictions$sio4 < 0 )) != 0) {
  neg_values_sio4 <- which(predictions$sio4 < 0 )
  
  min_values <- station_data$SiO4[order(station_data$SiO4, decreasing=F)]
  min_values <- na.omit(min_values)[1:20]
  min_sio4 <- min(min_values)
  max_sio4 <- max(min_values)
  
  predictions$sio4[neg_values_sio4] <- runif(length(neg_values_sio4), min_sio4, max_sio4)
}


#Save predictions as csv
write.csv(predictions, file = paste0(wd,"Input data/inputData_npzd_station", station_region[1],".csv"),row.names = F)

  return(predictions)
}

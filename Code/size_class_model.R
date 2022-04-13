library(devtools)
library(ssPopModel)
library(plotly)
library(tidyverse)
library("lattice")
library('gridExtra')

#======================================
# Settings
#======================================

# To change with your path:
path = "C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper"

res.folder = file.path('Results', 'growth_rates')
fig.folder = file.path('Figures', 'growth_rates')

setwd(path)

# Phytoplankton group: (to choose in ('ORGNANO', 'ORGPICOPRO', 'REDNANO', 'REDPICOEUK',\
#'REDPICOPRO')
pfg = 'ORGNANO'

colors = hcl.colors(30, "Zissou 1")

# Need for the 24 hours of the day and the first hour
# of the next day to estimate 24 hourly growth rates
nb.hours = 24
E.points.per.hour = 6

## Biomass conversion coefficients
if (pfg %in% c('ORGPICOPRO', 'REDPICOPRO', 'REDPICOEUK')) {
  a = 0.26
  b = 0.86
} else {
  a = 0.433
  b = 0.863
}

#=======================================
# Load PAR and events data
#=======================================

# Load PAR data separately (not the same sampling frequency)
E <- read.csv(file.path('Data','par','par_10min.csv')) 
colnames(E) = c('time','par')
E.days = format(strptime(E$time,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d')

#=======================================
# Load PAR and events data
#=======================================

events = read.csv(file.path('Results','event_T_anomalies_suited.csv'))
events = events[,c('window_start', 'window_end')]
n.events = dim(events)[1]

#=======================================
# Iterate through the events
#=======================================

# Create the results storage
all.dates.params = c()
NPP.size.ts = c()
mu.ts = c()
loss.ts = c()

start <- Sys.time()
for (event.idx in 1:n.events){

  # Event metadata
  event = events[event.idx,]
  days = seq(as.Date(event$window_start), as.Date(event$window_end), by="days")
  n.days = length(days)
  
  # Import phytoplankton data
  pfg.path = file.path('Data', 'INSITU', 'biovolume_classes', pfg,
                       paste0(as.character(days[1]), '.csv'))
  
  if (file.exists(pfg.path)){
    df <- read_csv(pfg.path) 
  }else{ # If no corresponding phytoplankton data, skip the event
    print(paste0('event', as.character(days[1]),' has no phytoplankton data'))
    next
  }
  
  # Format the phytoplankton data
  dist <- add_column(df, pop = pfg, .after = "time")
  dist <- add_column(dist, par = NA, .after = 'pop')
  dfh <- transform_PSD(dist, time.step="1 hour", interval.to.geomean=T)

  # Get all the dates of the series
  data.days = format(strptime(dfh$time,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d')
  E.days = format(strptime(E$time,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d')
  
  # Estimate per day growth rate
  for (day.idx in 1:n.days){
    print(day.idx)
    
    # Phyto data
    distribution.day <- dfh[data.days == days[day.idx],-c(2, 3)]
    distribution.next.day <- dfh[data.days == days[day.idx + 1],-c(2, 3)]
    distribution <- bind_rows(distribution.day, distribution.next.day[c(1),])
    
    if (length(distribution$time) != 25){
      print(length(distribution$time))
      print('Lacks data for that day')
      next
    }
    
    # Par data 
    E.data = E[E.days == days[day.idx],]
    
    # TO DO: Change the break in continue once has shown to be effective
    if (sum(is.na(E.data$par)) > 0){
      print('NaNs in E')
      next
    }
    
    if (sum(is.na(distribution)) > 0){
      print('Missing distrib values')
      next
    }
    
    if (nrow(E.data) == 0){
      print('No PAR DATA')
      next
    }
    #=================================
    # Smooth the distributions
    #=================================
    
    # Smooth the size distribution
    dist <- round(t(apply(distribution[,-c(1)],1, function(x) smooth.spline(x, all.knots=TRUE)$y)))
    dist[which(dist < 1)] <- 0                              
    distribution[,-c(1)] <- dist
    
    ma <- round(apply(distribution[,-c(1, 2)], 1, function(x) sum(x,na.rm=T)))                   
    Vdist <- as.matrix(distribution[,-c(1, 2)]/ma)   
    
    plot_ly(z= Vdist) %>% 
      add_surface() %>%
      layout(scene = list(xaxis = list(autorange = "reversed")))
    
    # No need for PAR smoothing
    
    #=================================
    # Model estimation
    #=================================
    
    volbins <- as.numeric(colnames(distribution)[-c(1)])
    resol <- 10
    dt <- resol/60
    output <- ssPopModel::determine_opt_para(distribution,E.data$par,resol)
    
    #=================================
    # Format the output
    #=================================
    
    #*************************
    # Growth rate
    #*************************
    
    params <- output$parameters
    PSD <- output$PSD
    mu_N <- diff(log(rowSums(PSD[,-c(1)], na.rm=T))) / as.numeric(diff(PSD$time))
    d.mu_N <- sum(mu_N, na.rm = T) 
    
    mu_N.obs = diff(log(rowSums(distribution[,-c(1)], na.rm=T))) / as.numeric(diff(distribution$time))
    d.mu_N.obs <- sum(mu_N.obs, na.rm=T)
    
    loss_N = mu_N - mu_N.obs
    
    #*************************
    # Frequency distribution
    #*************************
    
    s <- round(apply(distribution[,-c(1)], 1, function(x) sum(x,na.rm=T)))                   
    Nproj <- as.matrix(PSD[,-c(1)])
    s2 <- round(apply(Nproj, 1, function(x) sum(x,na.rm=T)))
    Vproj <- as.matrix(Nproj/s2)
    
    #=================================
    # Save the results of the day
    #=================================
    
    #****************************
    # Observed vs predicted plot
    #****************************
    
    hours = format(strptime(distribution$time,'%Y-%m-%d %H:%M:%S'),'%H')
    hours = as.numeric(hours)
    day = format(strptime(distribution$time,'%Y-%m-%d %H:%M:%S'),'%Y-%m-%d')
    
    d.data = distribution[1:nb.hours,-c(1)]/s
    m = dim(d.data)[2]
    
    d.data$hour = hours[1:nb.hours]
    d.data = d.data %>% gather(SizeClass, freq,
                               colnames(d.data)[1]:colnames(d.data)[m])
    d.data$SizeClass = as.numeric(d.data$SizeClass)
    
    plot1 <- levelplot(freq ~ hour*SizeClass, data=d.data ,
                       main="Observed", col.regions = colors)
    
    d.pred = as.data.frame(Vproj[1:nb.hours,]) 
    d.pred$hour = hours[1:nb.hours]
    
    d.pred = d.pred %>% gather(SizeClass, freq,
                               colnames(d.pred)[1]:colnames(d.pred)[m])
    d.pred$SizeClass = as.numeric(d.pred$SizeClass)
    
    plot2 <- levelplot(freq ~ hour*SizeClass, data=d.pred,
                       main="Predicted", col.regions = colors)
    
    jpeg(file = file.path(fig.folder, 'obs_vs_pred', pfg, paste0(day[1],".jpg")))
    plot = grid.arrange(plot1, plot2, ncol=2, top = unique(day)[1])
    dev.off()
    
    #****************************
    # Parameters
    #****************************
    
    params = cbind(day[1], params)
    all.dates.params = rbind(all.dates.params, params)
    
    # Store them at each iteration (if a run fails)
    write.table(all.dates.params, 
                col.names = T,
                file.path(res.folder, 'best_params',
                          paste0(pfg,".csv")),
                sep = ',',
                row.names = F)
    
    #****************************
    # Division series
    #****************************
    
    # Hack here on the nb of values of mu_N
    mu.day = data.frame(distribution$time[1:nb.hours], mu_N)# Legacy: [1:nb.hours]
    colnames(mu.day) = c('date', 'mu')
    mu.ts = rbind(mu.ts, mu.day)
    
    # Store them at each iteration (if a run fails)
    write.table(mu.ts, 
                col.names = T,
                file.path(res.folder, 'mu', paste0(pfg, ".csv")),
                sep = ',',
                row.names = F)
    
    #****************************
    # Division series
    #****************************
    
    # Hack here on the nb of values of mu_N
    loss.day = data.frame(distribution$time[1:nb.hours], loss_N)# [1:nb.hours]
    colnames(loss.day) = c('date', 'loss')
    loss.ts = rbind(loss.ts, loss.day)
    
    # Store them at each iteration (if a run fails)
    write.table(loss.ts, 
                col.names = T,
                file.path(res.folder, 'loss',
                          paste0(pfg, ".csv")),
                sep = ',',
                row.names = F)
    
    
    #****************************
    # Division plot 
    #****************************
    
    jpeg(file = file.path(fig.folder, 'cell_cycles', pfg,
                          paste0(day[1],".jpg")))
    plot(mu_N, frame = FALSE, main = day[1], ylab="Div rate (h-1)", xlab='time (hours)')
    dev.off()
    
  }
}

end <- Sys.time()
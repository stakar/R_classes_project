#Preparing packages
library(JADE)
library(dplyr)
library(signal)
library(matrixStats)
library(MASS)


#narrow-band bank of filters

filtering <- function(data,low,high,freq){
  b.low <- low/(freq*0.5)
  b.high <- high/(freq*0.5)
  but <- signal::butter(n=4,W=c(b.low,b.high),type='pass')
  filtered <- signal::filtfilt(but$b,but$a,data)
  return(filtered)
}

filtmat <- function(electrodes,low,high,freq){
  
  el <- matrix(nrow = nrow(electrodes),ncol = 
                 ncol(electrodes)) 
  for (n in 1:ncol(electrodes)){
    tmp <- filtering(electrodes[,n],low,high,freq)
    el[,n] <- tmp
  }
  return(el)
}

whole_filtering <- function(electrodes){
  
  el.8<- filtmat(electrodes,7.9,8.1,256)
  el.14 <- filtmat(electrodes,13.9,14.1,256)
  el.28 <-filtmat(electrodes,27.9,28.1,256)
  
  return(list(el.8,el.14,el.28))
}

#variance analyzer

savgol <- function(va.matrix){
  sgol <- list()
  for (z in 1:length(va)) {
    tmp <- matrix(nrow = nrow(va[[z]]),
                  ncol = ncol(va[[z]]))
    
    for (n in 1:ncol(va[[z]])){
      tmp[,n] <- sgolayfilt(va[[z]][,n],
                            p=2,
                            n=length(va[[z]][,n])-1)
    }
    sgol[[z]] <- tmp
  }
  return(sgol)
}

#channel integrator

integrator <- function(sg.list){
  integrated <- matrix(nrow = nrow(sg.list[[1]]),
                       ncol = ncol(sg.list[[1]]))
  for (n in 1:length(sg.list)) {
    integrated[,n] <- rowMeans(sg.list[[n]])
  }
  return(integrated)
}

#signal normalizer

normalizer <- function(int.list){
  norm.list <- matrix(nrow = nrow(int.list),
                      ncol = ncol(int.list))
  for (n in 1:ncol(int.list)) {
    norm.list[,n] <- int.list[,n]/rowSums(int.list)
  }
  return(norm.list)
}

#Bakardjian system

sysBak <- function(path,list.electrodes,chunk){
  
  data <- read.csv(path,sep = ',',header = FALSE)
  data1 <- apply(data, FUN = as.numeric,MARGIN = c(1,2))
  data<-data1[(5*256):(20*256),]
  data <- data[((chunk-1)*256):(chunk*256),]
  bss <- AMUSE(data)
  electrodes <- bss$S[,list.electrodes]
  filtered <- whole_filtering(electrodes = electrodes)
  va <- list()
  for (n in 1:length(filtered)){
    va[[n]] <- abs(filtered[[n]])
  }
  sgfiltered <- savgol(va)
  integ <- integrator(sgfiltered)
  norm <- normalizer(integ)
  
  return(norm)
}
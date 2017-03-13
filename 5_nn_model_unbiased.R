rm(list=ls()) #Clear the R environment
library(jpeg) #Package to read a jpeg image 
library("EBImage") # for resize
directory = 'G:\\TomTOm\\exp2\\'  #My local directory #'C:\\Users\\jjaaae\\Documents\\Cyient Insights\\TomTom'
img_dir <- paste(directory,"\\imageset",sep="")
label_dir <- paste(directory,"\\prob_maps",sep="")
images = list.files(img_dir)
labels = list.files(label_dir)

nonlinear <- function(x,b){
  if (b)
    return(x*(1-x))
  else
    return(1/(1+exp(-(x+1))))
} 

k = 11 #what is k?  It seems to be width of a border around label images, or something similar based on label code below
l_rate = 0.01  #learning rate
Hn = 10
img_scale = 0.5

#k*(k+1)*10 uniformly distributed numbers between 0 and 0.5 in 10 column matrix - Rounded to 3 decimal places
coff_1 <- matrix(round(runif((k*k+1)*10, 0, 0.05),3),nrow = k*k+1)

#(Hn+1) uniformly distributed numbers between 0 and 0.5 in 1 column matrix - Rounded to 3 decimal places
coff_2 <- matrix(round(runif((Hn+1)*1, 0, 0.05),3),nrow = Hn+1)

#f=1
for (itr in 1:300) {
  lrate = l_rate/sqrt(itr)
  for (f in sample(c(1:20))) {
    #LOAD THE IMAGE DATA
    img = readJPEG(paste(img_dir,"\\",images[2*f-1],sep = "")) #read jpeg file (reading file every iteration, why?)
    # This was here originally: img = 0.2126*img[,,1] + 0.7152*img[,,2] + 0.0722*img[,,3] # Convert to grayscale
    #img <- channel(img, "gray") #Why not use something like this?
    img = resize(img, dim(img)[1]*img_scale, dim(img)[2]*img_scale)
    
    #LOAD THE LABEL MASK DATA
    label = readJPEG(paste(label_dir,"\\",labels[f],sep = ""))#read label image
    label = resize(label,dim(label)[1]*img_scale,dim(label)[2]*img_scale)     #rescale to match the image sizes
    label = ifelse(label>0.1,1,0)                             #apply threshold to mask labels
    label[-(((k+1)/2):(nrow(label)-(k-1)/2)),] = -1  #setting a bunch of middle rows to be = to -1, why?
    label[,-(((k+1)/2):(ncol(label)-(k-1)/2))] = -1  #setting a bunch of middle columns to be = to -1, why?
    
    A = which(label==1,TRUE)  #Get row/col of locations where label = 1
    B = which(label==0,TRUE)  #Get row/col of locations where label = 0
    B = B[sample(1:nrow(B),nrow(A)),]  #Get random sample from B that matches dimensions of A
    A = rbind(A,B)                     #Add B to the bottom of A
    A = A[sample(1:nrow(A)),]          #Randomly shuffle the row order of A/B combined matrix
    cost = 0
    
    #WHAT IS THIS FOR LOOP TRYING TO ACCOMPLISH?
    for (j in 1:(dim(A)[1]%/%100)) {  #j cycles from 1 to hundreds unit of dimension of A
      fst = (100*(j-1)+1)             #multiplier to get j in multiples of 100, but starting at 1 (1, 101, 201,...)
      
      #makes 1 column matrix with 1 label value repeated in 100 rows
      truelabel = matrix(label[A[fst,1],A[fst,2]], ncol = 1, nrow = 100)
      
      #
      L1 = matrix(img[(A[fst,1]-(k-1)/2):(A[fst,1]+(k-1)/2),(A[fst,2]-(k-1)/2):(A[fst,2]+(k-1)/2)],ncol = k*k)
      for (i in (100*(j-1)+2):(100*j)) {
        L1 = rbind(L1,matrix(img[(A[i,1]-(k-1)/2):(A[i,1]+(k-1)/2),(A[i,2]-(k-1)/2):(A[i,2]+(k-1)/2)],ncol = k*k))
        truelabel[(i-1)%%100+1] = label[A[i,1],A[i,2]]
      }
      L1 = cbind(matrix(1,nrow = nrow(L1),ncol = 1),L1)  #Add column of 1's as first column of L1
      L2 = nonlinear(L1 %*% coff_1,0)
      L2 = cbind(matrix(1,nrow = nrow(L2),ncol = 1),L2)  #Add column of 1's as first column of L2
      L3 = nonlinear(L2 %*% coff_2,0)
      err = (L3 - truelabel)
      cost = cost + sum(err^2)
      delta3 = err
      delta2 = (delta3%*%t(coff_2))*nonlinear(L2,1)
      delta2 = delta2[,-1]
      
      coff_2 = coff_2 - lrate*t(L2)%*%delta3
      coff_1 = coff_1 - lrate*t(L1)%*%delta2
      # print(c(itr,f,j,round(cost,2)))
    }
    print(c(itr,f,round(cost,2)))
  }
}

write.csv(coff_2, file=paste(directory,"coff_2.csv",sep = ""))
write.csv(coff_1, file=paste(directory,"coff_1.csv",sep = ""))

# rows = ((k-1)/2+1):(nrow(img)-(k-1)/2)
# bats = split(rows, ceiling(1:length(rows)/5))


rm(list=ls())
gc()
library(data.table)
library(jpeg)
library(EBImage)
directory = 'G:\\TomTOm\\exp2\\'
images = list.files(paste(directory,"imageset\\",sep=""))
labels = list.files(paste(directory,"prob_maps\\",sep=""))


nonlinear <- function(x,b){
  if (b)
    return(x*(1-x))
  else
    return(1/(1+exp(-(x+1))))
}

k=11
coff_2 = fread(paste(directory,"coff_2.csv",sep = ""))
coff_2[,1]<-NULL
coff_2 = as.matrix(coff_2)
coff_1 = fread(paste(directory,"coff_1.csv",sep = ""))
coff_1[,1]<-NULL
coff_1 = as.matrix(coff_1)

f=21
for (f in 21:22) {
  img = readJPEG(paste(directory,"imageset\\",images[2*f-1],sep = ""))#read jpeg file
  img = 0.2126*img[,,1] + 0.7152*img[,,2] + 0.0722*img[,,3] # Convert to grayscale
  img = resize(img,dim(img)[1]/2,dim(img)[2]/2)
  
  pred = matrix(0,nrow = dim(img)[1],ncol = dim(img)[2])
  for (i in ((k-1)/2+1):(nrow(img)-(k-1)/2) ) {
    for (j in ((k-1)/2+1):(ncol(img)-(k-1)/2) ) {
      dp = matrix(img[(i-(k-1)/2):(i+(k-1)/2),(j-(k-1)/2):(j+(k-1)/2)],ncol = k*k)
      dp = cbind(matrix(1,nrow = nrow(dp),ncol = 1),dp)
      Hlay = nonlinear(dp %*% coff_1,0)
      Hlay = cbind(matrix(1,nrow = nrow(Hlay),ncol = 1),Hlay)
      pred[i,j] = nonlinear(Hlay %*% coff_2,0) 
    }
    print(c(f,i))
  }
  # pred = normalize(pred)
  par(mar=c(1,1,1,1))
  plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
  rasterImage(pred,0,0,1,1)
  
  # label = readJPEG(paste(directory,"prob_maps\\",labels[f],sep = ""))#read label image
  # label = resize(label,dim(label)[1]/2,dim(label)[2]/2)
}


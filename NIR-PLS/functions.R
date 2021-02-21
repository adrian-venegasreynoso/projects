library(pracma)
library(Hmisc)
library(squash)
my_savgol <- function(spectra, points, forder, dorder)
{
        filtered_data <- matrix(0, nrow(spectra), ncol(spectra))

        for(i in 1:nrow(spectra))
                filtered_data[i,] <- savgol(spectra[i,], points, forder, dorder)

return(filtered_data)
}

snv <- function(spectra)
{
        spectra <- scale(t(spectra), center=TRUE, scale=TRUE)
        return(t(spectra))
}

spectra_plot_y <- function(mydata, wavelengths, ordered_y)
{
         #Create color map:
  xlim <- c(wavelengths[1],wavelengths[length(wavelengths)])
  ylim <- c(ordered_y[1], ordered_y[length(ordered_y)])
  map <- makecmap(ordered_y)                                                        
  mycol <- cmap(ordered_y, map = map)
  par(font=2,las=1,mar = c(5,4,4,10) + 0.1)
  
  #Plot spectra:
  matplot(wavelengths,t(mydata),font.axis=2,main="",          
          col=mycol,lty=1, xlab="",ylab="",type="l", xlim=xlim)
  minor.tick(nx=2, ny=2,tick.ratio=0.75)                                          
  par(mar = c(5,4,4,6) + 0.1)
  title(xlab="Wavelength (nm)",ylab="Absorbance (a.u.)",font.lab=2)
  
  #Plot color map:
  vkey(map, title = "SSC (ºBrix)", stretch=1.25, side=2, skip=2, x=1100, y=0.0075)

}
spectra_plot <- function(mydata, wavelengths)
{
         #Create color map:
  xlim <- c(wavelengths[1],wavelengths[length(wavelengths)])
  par(font=2,las=1,mar = c(5,4,4,10) + 0.1)
  
  #Plot spectra:
  matplot(wavelengths,t(mydata),font.axis=2,main="",          
          lty=1, xlab="",ylab="",type="l", xlim=xlim)
  minor.tick(nx=2, ny=2,tick.ratio=0.75)                                          
  par(mar = c(5,4,4,6) + 0.1)
  title(xlab="Wavelength (nm)",ylab="Absorbance (log[1/R])",font.lab=2)
}


chop_spectra <- function(spectra, wavelength, wl1, wl2)
{
        i_min <- which(wavelength >= wl1)[1]
        i_max <- which(wavelength >= wl2)[1]

        spectra <- spectra[,i_min:i_max]

        return(spectra)
#        wavelength <<- wavelength[i_min,i_max]
}

my_savgol=function(x,width,order,deriv)
{
  ##insert spdiags function
  spdiags <-function (arg1,arg2,arg3,arg4){
    B <- arg1 
    if (is.matrix(arg2))
      d <- matrix(arg2,dim(arg2)[1]*dim(arg2)[2],1)
    else
      d <- arg2
    p <- length(d) 
    A <- sparseMatrix(i = 1:arg3, j = 1:arg3, x = 0, dims= c(arg3,arg4))
    m <- dim(A)[1] 
    n <- dim(A)[2] 
    len<-matrix(0,p+1,1)
    for (k in 1:p)
      len[k+1] <- len[k]+length(max(1,1-d[k]):min(m,n-d[k])) 
    a <- matrix(0, len[p+1],3) 
    for (k in 1:p)
    {
      i <- t(max(1,1-d[k]):min(m,n-d[k])) 
      a[(len[k]+1):len[k+1],] <- c(i, i+d[k], B[(i+(m>=n)*d[k]),k]) 
    }
    res1 <- sparseMatrix(i = a[,1], j = a[,2], x = a[,3], dims = c(m,n))
    res1 <- apply(res1, 2, as.numeric)
    return (res1)
  }
  ##end spdiags function

  m=nrow(x)
  n=ncol(x)
  w=max(3,1+2*round((width-1)/2) )
  o=min(c(max(0,round(order)),5,w-1))
  d=min(max(0,round(deriv)),o)
  p=(w-1)/2
  xc=((-p:p)%*%matrix(1,1,1+o))^(t(matrix(1,1,w))%*%(0:o))
  we=qr.solve(xc,diag(w))
    options(warn=-1)
  b=apply((matrix(1,d,1)%*%matrix(1:(o+1-d),1,(o+1-d))+t(matrix(0:(d-1),1,(d)))%*%matrix(1,1,o+1-d)),2,prod)
  gg=matrix(1,n,1)%*%we[(d+1),]*b[1]
    library(Matrix)
  di=spdiags(gg,p:(-p),n,n)
    options(warn=0)
  w1=diag(b,nrow=length(b))%*%we[(d+1):(o+1),]
  di[1:w,1:(p+1)]=t(xc[1:(p+1),1:(1+o-d)]%*%w1)
  di[(n-w+1):n,(n-p):n]=t(xc[(p+1):w,1:(1+o-d)]%*%w1)
  result=x%*%di
}

#Piecewise Direct Standardization (PDS) algorithm:

#INPUT:   masterSpectra = Spectra acquired with the master instrument (matrix).
#         slaveSpectra = Spectra acquired with the slave instrument (matrix).
#         MWsize = Half size of the moving window (integer).
#         Ncomp = Number of latent variables used in the PLS model (integer).
#         wavelength = wavelength (numeric vector).

#OUTPUT:  P = the PDS transfer matrix.

PDS<-function(masterSpectra, slaveSpectra, MWsize, Ncomp, wavelength){
  
require(pls)

#Loop Initialization:
i<-MWsize
k<-i-1
#Creation of an empty P matrix:
P<-matrix(0,nrow=ncol(masterSpectra),ncol=ncol(masterSpectra)-(2*i)+2)
InterceptReg<-c()

while(i<=(ncol(masterSpectra)-k)){
  
  #PLS regression:
  fit<- plsr(masterSpectra[,i] ~ as.matrix(slaveSpectra[,(i-k):(i+k)]),
             ncomp=Ncomp, scale=F, method="oscorespls")
  
  #Extraction of the regression coefficients:
  coefReg<-as.numeric(coef(fit, ncomp=Ncomp, intercept = TRUE))
  InterceptReg<-c(InterceptReg,coefReg[1])
  coefReg<-coefReg[2:length(coefReg)]
  
  #Add coefficients to the transfer matrix:
  P[(i-k):(i+k),i-k]<-t(coefReg)
  rm(coefReg,fit)
  i<-i+1
  
  #Diplay progression:
  cat("\r",paste(round(i/ncol(masterSpectra)*100)," %",sep=""))}

P<-data.frame(matrix(0,nrow=ncol(masterSpectra),ncol=k), P,
              matrix(0,nrow=ncol(masterSpectra),ncol=k))
InterceptReg<-c(rep(0,k),InterceptReg,rep(0,k)) 

Output<-list(P = P , Intercept = InterceptReg)
return(Output)}

transfer_cal <- function(master, slave, wavelength, val) #Val is added
{
#Plot the master and slave spectra:
par(font=2,las=1,mar = c(5,5,4,2) + 0.1)                   
matplot(wavelength,t(master),font.axis=2,main="",        
        col="blue",lty=1, xlab="",ylab="",type="l",lwd=1)
par(new=TRUE)                                             
matplot(wavelength,t(slave),font.axis=2,main="",          
        col="forestgreen",lty=1, xlab="",ylab="",type="l",lwd=1)
minor.tick(nx=2, ny=2,tick.ratio=0.75)
par(mar = c(5,4,4,2) + 0.1)
title( xlab="Wavelength (nm)",ylab="Absorbance (Log(1/R))",cex.lab=1,font.lab=2)
legend("topleft",bty="n", legend=c("Master data","Slave data"),
       col=c("blue","forestgreen"),lty=c(1),cex=1,lwd=2)
Ncomp<-2  
MWsize<-2

#Compute the transfer matrix P:
Pmat<-PDS(master, slave, MWsize, Ncomp, wavelength)
print(length(Pmat$Intercept))
#Standardization of the slave data (using P and the intercept):
SlaveCor<-as.matrix(val)%*%as.matrix(Pmat$P)
SlaveCor<-sweep(SlaveCor, 2, as.numeric(t(Pmat$Intercept)), "+")


#Plot the master and slave spectra:
par(font=2,las=1,mar = c(5,5,4,2) + 0.1)                   
matplot(wavelength,t(master),font.axis=2,main="",        
        col="blue",lty=1, xlab="",ylab="",type="l",lwd=1)
par(new=TRUE)                                           
matplot(wavelength,t(val),font.axis=2,main="",          
        col="forestgreen",lty=1, xlab="",ylab="",type="l",lwd=1)
par(new=TRUE)                                                
matplot(wavelength,t(SlaveCor),font.axis=2,main="",     
        col="red",lty=1, xlab="",ylab="",type="l",lwd=1)
minor.tick(nx=2, ny=2,tick.ratio=0.75)
par(mar = c(5,4,4,2) + 0.1)
title( xlab="Wavelength (nm)",ylab="Absorbance (Log(1/R))",cex.lab=1,font.lab=2)    
legend("topleft",bty="n", legend=c("Master data","Slave data", "Stand. slave data"),
       col=c("blue","forestgreen","red"),lty=c(1),cex=1,lwd=3)

return(SlaveCor)
}


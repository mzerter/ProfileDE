CollocInferPlots = function(coefs,pars,lik,proc,times=NULL,data=NULL,
      cols=NULL,datacols=NULL,datanames=NULL,ObsPlot=TRUE,DerivPlot=TRUE)
{
  if(is.null(cols)){ cols = 1:ncol(coefs) }

 timevec = proc$more$qpts

 traj = proc$bvals$bvals%*%coefs
 colnames(traj) = proc$more$names
 
 dtraj = proc$bvals$dbvals%*%coefs

 ftraj = proc$more$fn(timevec,traj,pars,proc$more$more)
 
 otraj = lik$more$fn(timevec,traj,pars,lik$more$more)
 
 
 
 if(ObsPlot){
  if(is.null(times)){ times = 1:nrow(data) }
  if(is.null(datacols)){  
    if(ncol(otraj) == ncol(traj)){ datacols = cols}
    else{ datacols = 1:ncol(otraj) }
  }
  
  X11()
  matplot(timevec,otraj,type='l',lwd=2,xlab='time',ylab='Observations',
        cex.lab=1.5,cex.axis=1.5,col=datacols)
  if(!is.null(data)){
    if(is.null(datanames)){ datanames = substring(colnames(data),1,1) }
    if(is.null(datanames)){ datanames = 1:ncol(data) }
    if(is.null(times)){ times = 1:nrow(data) }
    matplot(times,data,col=datacols,add=TRUE,pch=datanames)
  }
 }
  
 if(DerivPlot){
  X11()
  matplot(timevec,ftraj,type='l',lwd=2,xlab='time',ylab='f(x): -, dx: --',
        cex.lab=1.5,cex.axis=1.5,col=cols,lty=1) 
  matplot(timevec,dtraj,type='l',lwd=2,lty=2,add=TRUE,col=cols)
  abline(h = 0)
  legend(x='topright',legend=proc$more$names,lwd=2,lty=1,col=cols)
 
  X11()
  par(mfrow=c(2,1),mar=c(4,4,1,1))
  matplot(timevec,dtraj-ftraj,type='l',lwd=2,xlab='time',ylab='dx-f(x)',
        cex.lab=1.5,cex.axis=1.5,col=cols,lty=1)
  abline(h=0)     
  matplot(timevec,traj,type='l',lwd=2,xlab='time',ylab='x',
        cex.lab=1.5,cex.axis=1.5,col=cols,lty=1)    
  legend(x='topright',legend=proc$more$names,lwd=2,lty=1,col=cols)
 }

 return(list(timevec=timevec,traj=traj,dtraj=dtraj,ftraj=ftraj,otraj=otraj))

}
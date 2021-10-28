R=function(d) {
  result=matrix(0,ncol=N,nrow=N);
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      result[i,j]=prod(d^((as.vector(abs(D[ind[i],]-D[ind[j],])))^1.9));
      
    }
  }
  result=result+t(result)+(1+10^(-5))*diag(N);
  return(result);
}



N=length(y1);
ind=1:N;

m3=dim(D3)[1]; m4=dim(D4)[1]; m5=dim(D5)[1];
y.exp=c(D3[,6],D4[,6],D5[,6]);
m=length(y.exp);
rel=function(x) {
  result=(x-means[2])/sds[2];
  return(as.numeric(result));
}


m=m3+m4+m5;

R.true=d.400[4]*R(d.400[1:3])+d.400[5]*diag(N); E=eigen(R.true);

w=as.vector(E$vectors%*%diag(1/E$values)%*%t(E$vectors)%*%as.matrix(y1[ind]));

R.pred=function(d,X) {
  M=dim(X)[1];
  
  result=matrix(0,ncol=N,nrow=M);
  for (i in 1:N) {
    for (j in 1:M) {
      result[j,i]=prod(d^((as.vector(abs(D[ind[i],]-X[j,])))^2));
      
    }
  }
  return(result);
}            

m.x=function(X) {
  R.new=R.pred(d.400[1:3],X);
  y1.pred=d.400[4]*as.vector(R.new%*%as.matrix(w));
  return(y1.pred)
}

#glass
x1=seq(D3[1,5],D3[m3,5]+.05,.01);
x2=seq(D4[1,5],D4[m4,5]+.05,.05);
x3=seq(D5[1,5],D5[m5,5]+.05,.05);


X1=cbind(rep(Mode(MC[,1]),length(x1)),cbind(rep(rel((Mode(MC[,4]+2*MC[,7+m+2])/2)),length(x1)),x1));
X2=cbind(rep(Mode(MC[,2]),length(x2)),cbind(rep(rel((Mode(MC[,4]+15*MC[,7+m+2])/15)),length(x2)),x2));
X3=cbind(rep(Mode(MC[,3]),length(x3)),cbind(rep(rel((Mode(MC[,4]+50*MC[,7+m+2])/50)),length(x3)),x3));
X4=cbind(rep(min(-1.28),length(x3)),cbind(rep(rel(0.005),length(x3)),x3)); #g.ps=.7, g.pc=.8, g.san=-1, g.pp=-1.4 (.01)

par(mfrow=c(1,1));
x1.t=x1*sds[3]+means[3];
m.x1=m.x(X1)*sds[4]+means[4];
plot(x1.t,m.x1,type='l',lwd=1,ylim=c(min(y.exp*sds[4]+means[4])-.3,max(y.exp*sds[4]+means[4])+.3),col='blue',xlab='Filler Volume Fraction',ylab='Elastic Modulus (GPa)',main='Glass PP (Holdout)'); 
points(D3[,5]*sds[3]+means[3],sds[4]*y.exp[1:m3]+means[4],pch='+',cex=1.5,col='blue')
arrows(x0=D3[,5]*sds[3]+means[3],y0=(sds[4]*y.exp[1:m3]+means[4])*.95,x1=D3[,5]*sds[3]+means[3],y1=(sds[4]*y.exp[1:m3]+means[4])*1.05,code=3,angle=90,length=.05,col='blue')
#points(D.cal[1:5,7]*sds[3]+means[3],m.cal[1:5],col='blue',pch='x',cex=1)

x2.t=x2*sds[3]+means[3];
m.x2=m.x(X2)*sds[4]+means[4];
lines(x2.t,m.x2,type='l',lwd=1,col='green',ylim=c(3,6)); points(D4[,5]*sds[3]+means[3],y.exp[(m3+1):(m3+m4)]*sds[4]+means[4],pch='*',cex=2,col='green')
arrows(x0=D4[,5]*sds[3]+means[3],y0=(sds[4]*y.exp[(m3+1):(m3+m4)]+means[4])*.95,x1=D4[,5]*sds[3]+means[3],y1=(sds[4]*y.exp[(m3+1):(m3+m4)]+means[4])*1.05,code=3,angle=90,length=.05,col='green')
#points(D.cal[6:9,7],y.cal[6:9],col='green',pch='+',cex=2)


x3.t=x3*sds[3]+means[3];
m.x3=m.x(X3)*sds[4]+means[4];

lines(x3.t,m.x3,type='l',lwd=1,col='orange',ylim=c(3,6)); points(D5[,5]*sds[3]+means[3],y.exp[(m3+m4+1):m]*sds[4]+means[4],pch=15,cex=1,col='orange')
arrows(x0=D5[,5]*sds[3]+means[3],y0=(sds[4]*y.exp[(m3+m4+1):m]+means[4])*.95,x1=D5[,5]*sds[3]+means[3],y1=(sds[4]*y.exp[(m3+m4+1):m]+means[4])*1.05,code=3,angle=90,length=.05,col='orange')
#points(D.cal[10:12,7],y.cal[10:12],col='orange',pch='x',cex=2)


m.x4=m.x(X4)*sds[4]+means[4];
lines(x3.t,m.x4,type='l',lwd=1,col='black')

legend('topleft',legend=c("4 micron","30 micron","100 micron","No Interphase"),col=c('blue','green','orange','black'),lty=rep(1,3),cex=1)



#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################



D3=D1; D4=D2;
m3=dim(D3)[1]; m4=dim(D4)[1]; m5=dim(D5)[1];
y.exp=c(D3[,6],D4[,6],D5[,6]);
m=length(y.exp);
rel=function(x) {
  result=(x-means[2])/sds[2];
  return(as.numeric(result));
}
N=length(y1);
ind=1:N;




m=m3+m4+m5;

R.true=d.400[4]*R(d.400[1:3])+d.400[5]*diag(N); E=eigen(R.true);

w=as.vector(E$vectors%*%diag(1/E$values)%*%t(E$vectors)%*%as.matrix(y1[ind]));

R.pred=function(d,X) {
  M=dim(X)[1];
  
  result=matrix(0,ncol=N,nrow=M);
  for (i in 1:N) {
    for (j in 1:M) {
      result[j,i]=prod(d^((as.vector(abs(D[ind[i],]-X[j,])))^2));
      
    }
  }
  return(result);
}            

m.x=function(X) {
  R.new=R.pred(d.400[1:3],X);
  y1.pred=d.400[4]*as.vector(R.new%*%as.matrix(w));
  return(y1.pred)
}


x1=seq(D3[1,5],D3[m3,5]+.05,.05);
x2=seq(D4[1,5],D4[m4,5]+.05,.05);


X1=cbind(rep(Mode(MC[,1]),length(x1)),cbind(rep(rel((Mode(MC[,4]+.0175*MC[,7+m+2])/.0175)),length(x1)),x1));
X2=cbind(rep(Mode(MC[,2]),length(x2)),cbind(rep(rel((Mode(MC[,4]+.2*MC[,7+m+2])/.2)),length(x2)),x2));

par(mfrow=c(1,1));
x1.t=x1*sds[3]+means[3];
m.x1=m.x(X1)*sds[4]+means[4];
plot(x1.t,m.x1,type='l',lwd=2,ylim=c(3,7.5),col='blue',xlab='Filler Volume Fraction',ylab='Elastic Modulus (GPa)',main='Alumina SAN'); points(D3[,5]*sds[3]+means[3],sds[4]*y.exp[1:m3]+means[4],pch='+',cex=1.5,col='blue')
arrows(x0=D3[,5]*sds[3]+means[3],y0=(sds[4]*y.exp[1:m3]+means[4])*.95,x1=D3[,5]*sds[3]+means[3],y1=(sds[4]*y.exp[(1:m3)]+means[4])*1.05,code=3,angle=90,length=.05,col='blue')



x2.t=x2*sds[3]+means[3];
m.x2=m.x(X2)*sds[4]+means[4];
lines(x2.t,m.x2,type='l',lwd=2,col='green'); points(D4[,5]*sds[3]+means[3],y.exp[(m3+1):(m3+m4)]*sds[4]+means[4],pch='*',cex=2,col='green')
arrows(x0=D4[,5]*sds[3]+means[3],y0=(sds[4]*y.exp[(m3+1):(m3+m4)]+means[4])*.95,x1=D4[,5]*sds[3]+means[3],y1=(sds[4]*y.exp[(m3+1):(m3+m4)]+means[4])*1.05,code=3,angle=90,length=.05,col='green')


X4=cbind(rep(min(-.7),length(x2)),cbind(rep(rel(0.0),length(x2)),x2));
m.x4=m.x(X4)*sds[4]+means[4];
lines(x2.t,m.x4,type='l',lwd=1,col='black')

legend('topleft',legend=c(".035 micron",".4 micron","No Interphase"),col=c('blue','green','black'),lty=rep(1,3),cex=1)


cred["g.pp",]=HPDinterval(as.mcmc(MC.7[,dim(MC.7)[2]]))

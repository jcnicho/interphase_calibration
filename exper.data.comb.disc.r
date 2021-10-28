##exper.data.combined
N=length(y1);
ind=1:N;


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
      result[j,i]=prod(d^((as.vector(abs(D[ind[i],]-X[j,c(1,2,7)])))^2));
      
    }
  }
  return(result);
}            

m.x=function(X) {
  R.new=R.pred(d.400[1:3],X);
  y1.pred=d.400[4]*as.vector(R.new%*%as.matrix(w));
  return(y1.pred)
}

R.disc=function(d,x) {
  m=length(x);
  result=matrix(0,ncol=m,nrow=m);
  for (i in 1:m) {
    for (j in 1:m) {
      result[i,j]=d^(abs(x[i]-x[j])^1.9);
    }
  }
  return(result);
}

###z[1:3] Interphase Modulus, z[4] Abs Thickness, z[5] tau1, z[6] tau2, z[7] sig2
g1=4; g2=4; g3=.01; lam=.001; alpha=5; beta=5;

logf=function(k,z,mx,E) {
  z[5:7]=exp(z[5:7]); z[7+m+1]=log.inv(z[7+m+1])
  if (k==1) {
    result=sum(dnorm(mx[1:m3],mean=y.exp[1:m3]-z[(7+1):(7+m3)],sd=sqrt(z[7]),log=TRUE));
    result=result+log(dtruncnorm(z[1],a=b.1[1],b=b.1[2],mean=c(-1),sd=sqrt(z[5])));
  } else if (k==1.1) {
    result=sum(dnorm(mx[(m3+1):(m3+m4)],mean=y.exp[(m3+1):(m3+m4)]-z[(7+m3+1):(7+m3+m4)],sd=sqrt(z[7]),log=TRUE));
    result=result+log(dtruncnorm(z[2],a=b.1[1],b=b.1[2],mean=-1,sd=sqrt(z[5])))
  } else if (k==1.2) {
    result=sum(dnorm(mx[(m3+m4+1):m],mean=y.exp[(m3+m4+1):m]-z[(7+m3+m4+1):(7+m)],sd=sqrt(z[7]),log=TRUE))
    result=result+log(dtruncnorm(z[3],a=b.1[1],b=b.1[2],mean=-1,sd=sqrt(z[5])))
  } else if (k==1.3) {
    result=sum(dnorm(mx,mean=y.exp-z[(7+1):(7+m)],sd=sqrt(z[7]),log=TRUE));
    result=result+log(dnorm(z[4],mean=0,sd=sqrt(z[6])))+log(dnorm(z[7+m+2],mean=0,sd=sqrt(z[6])))
  } else if (k==2) {
    result=sum(log(dtruncnorm(z[1:3],a=b.1[1],b=b.1[2],mean=c(-.5,-1,-1),sd=sqrt(z[5]))))+dinvgamma(z[5],shape=g1,scale=g1,log=TRUE);
  } else if (k==3) {
    result=log(dnorm(z[4],mean=0,sd=sqrt(z[6])))+log(dnorm(z[7+m+2],mean=0,sd=sqrt(z[6])))+dinvgamma(z[6],shape=g2,scale=g2,log=TRUE);
  } else if (k==4) {
    result=sum(dnorm(mx,mean=y.exp-z[(7+1):(7+m)],sd=sqrt(z[7])),log=TRUE)+dinvgamma(z[7],shape=g3,scale=g3,log=TRUE);
  } else if (k==4.1) {
    result=sum(dnorm(mx[1:m3],mean=y.exp[1:m3]-z[(7+1):(7+m3)],sd=sqrt(z[7]),log=TRUE));
    b=diag(1/sqrt(E[[1]]$values))%*%t(E[[1]]$vectors)%*%z[(7+1):(7+m3)];
    result=result-m3/2*log(2*pi*lam)-1/2*sum(log(abs(E[[1]]$values)))-1/(2*lam)*sum(b^2);
  } else if (k==4.2) {
    result=sum(dnorm(mx[(m3+1):(m3+m4)],mean=y.exp[(m3+1):(m3+m4)]-z[(7+m3+1):(7+m3+m4)],sd=sqrt(z[7]),log=TRUE));
    b=diag(1/sqrt(E[[2]]$values))%*%t(E[[2]]$vectors)%*%z[(7+1+m3):(7+m3+m4)];
    result=result-m4/2*log(2*pi*lam)-1/2*sum(log(abs(E[[2]]$values)))-1/(2*lam)*sum(b^2);
  } else if (k==4.3) {
    result=sum(dnorm(mx[(m3+m4+1):m],mean=y.exp[(m3+m4+1):m]-z[(7+m3+m4+1):(7+m)],sd=sqrt(z[7]),log=TRUE));
    b=diag(1/sqrt(E[[3]]$values))%*%t(E[[3]]$vectors)%*%z[(7+1+m3+m4):(7+m)];
    result=result-m5/2*log(2*pi*lam)-1/2*sum(log(abs(E[[3]]$values)))-1/(2*lam)*sum(b^2);
  } else {
    result=dbeta(z[7+m+1],alpha,beta,log=TRUE);
    b=diag(1/sqrt(E[[1]]$values))%*%t(E[[1]]$vectors)%*%z[(7+1):(7+m3)];
    result=result-m3/2*log(2*pi*lam)-1/2*sum(log(abs(E[[1]]$values)))-1/(2*lam)*sum(b^2);
    b=diag(1/sqrt(E[[2]]$values))%*%t(E[[2]]$vectors)%*%z[(7+1+m3):(7+m3+m4)];
    result=result-m4/2*log(2*pi*lam)-1/2*sum(log(abs(E[[2]]$values)))-1/(2*lam)*sum(b^2);
    b=diag(1/sqrt(E[[3]]$values))%*%t(E[[3]]$vectors)%*%z[(7+1+m3+m4):(7+m)];
    result=result-m5/2*log(2*pi*lam)-1/2*sum(log(abs(E[[3]]$values)))-1/(2*lam)*sum(b^2);
    
  }
  return(result);
}


c0=c(rep(.1,3),.0001,rep(1,3),rep(.0001,3),.1);
acc=rep(0,11);
z=rnorm(7+m+2,mean=c(-1,-1,-1,.1,0,0,-1,rep(0,m),0,0),sd=0);
M=15000; burn=10000;
m3=dim(D3)[1]; m4=dim(D4)[1]; m5=dim(D5)[1];
X1=cbind(as.matrix(rep(z[1],m3)),cbind(as.matrix(rep(rel(2*(z[4]+z[7+m+2]*2)/4),m3)),D3[,1:5]));
X2=cbind(as.matrix(rep(z[2],m4)),cbind(as.matrix(rep(rel(2*(z[4]+z[7+m+2]*15)/30),m4)),D4[,1:5]));
X3=cbind(as.matrix(rep(z[3],m5)),cbind(as.matrix(rep(rel(2*(z[4]+z[7+m+2]*50)/100),m5)),D5[,1:5]));
start=Sys.time()
mx.old=c(m.x(X1),m.x(X2),m.x(X3));
stop=Sys.time()
stop-start;
A=R.disc(log.inv(z[7+m+1]),D3[,5]); E1=eigen(A,symmetric=TRUE);
A=R.disc(log.inv(z[7+m+1]),D4[,5]); E2=eigen(A,symmetric=TRUE);
A=R.disc(log.inv(z[7+m+1]),D5[,5]); E3=eigen(A,symmetric=TRUE);
E=list(E1,E2,E3);


start=Sys.time();
for (i in 1:burn) {
  
  zn=rnorm(7+m+2,mean=0,sd=sqrt(c(c0[1:7],rep(c0[8],m3),rep(c0[9],m4),rep(c0[10],m5),c0[11],c0[4])));
  t1=z[1]+zn[1]; 
  t2=z[2]+zn[2]; 
  t3=z[3]+zn[3]; 
  t4=z[4]+zn[4];
  t5=z[7+m+2]+zn[7+m+2];
  
  ###do relative thickness and interphase modulus separately
  #X1=cbind(as.matrix(rep(temp[1],m3)),cbind(as.matrix(rep(temp[4],m3)),D3[,1:5]));
  #X2=cbind(as.matrix(rep(temp[2],m4)),cbind(as.matrix(rep(temp[4]/30,m4)),D4[,1:5]));
  #X3=cbind(as.matrix(rep(temp[3],m5)),cbind(as.matrix(rep(temp[4]/100,m5)),D5[,1:5]));
  #mx.new=c(m.x(X1),m.x(X2),m.x(X3));
  
  #a=logf(1.5,temp,mx.new)-logf(1.5,z,mx.old);
  #if (a>=log(runif(1,0,1))) {
  #  z=temp;
  #  acc[5]=acc[5]+1;
  #  mx.old=mx.new;
  #}
  
  
  
  
  X1=cbind(as.matrix(rep(t1,m3)),cbind(as.matrix(rep(rel(2*(z[4]+z[7+m+2]*2)/4),m3)),D3[,1:5]));
  X2=cbind(as.matrix(rep(t2,m4)),cbind(as.matrix(rep(rel(2*(z[4]+z[7+m+2]*15)/30),m4)),D4[,1:5]));
  X3=cbind(as.matrix(rep(t3,m5)),cbind(as.matrix(rep(rel(2*(z[4]+z[7+m+2]*50)/100),m5)),D5[,1:5]));
  if (t1<b.1[1] || t1>b.1[2]) {
    sorry=1;
  } else {
    mx.new=c(m.x(X1),mx.old[(m3+1):m]);
    temp=c(t1,z[2:(7+m+1)],z[7+m+2]);
    a=logf(1,temp,mx.new,E)-logf(1,z,mx.old,E);
    if (a>=log(runif(1,0,1))) {
      z[1]=temp[1];
      acc[1]=acc[1]+1;
      mx.old[1:m3]=mx.new[1:m3];
    }
  }
  
  if (t2<b.1[1] || t2>b.1[2]) {
    sorry=1;
  } else {
    mx.new=c(mx.old[1:m3],m.x(X2),mx.old[(m3+m4+1):(m3+m4+m5)]);
    temp=c(z[1],t2,z[3:(1+m+7)],z[7+m+2]);
    a=logf(1.1,temp,mx.new,E)-logf(1.1,z,mx.old,E);
    
    if (a>=log(runif(1,0,1))) {
      z[2]=temp[2];
      acc[2]=acc[2]+1;
      mx.old[(m3+1):(m3+m4)]=mx.new[(m3+1):(m3+m4)];
    }
  }
  
  if (t3<b.1[1] || t3>b.1[2]) {
    sorry=1;
  } else {
    mx.new=c(mx.old[1:(m3+m4)],m.x(X3));
    temp=c(z[1:2],t3,z[4:(1+m+7)],z[7+m+2]);
    a=logf(1.2,temp,mx.new,E)-logf(1.2,z,mx.old,E);
    if (a>=log(runif(1,0,1))) {
      z[3]=temp[3];
      acc[3]=acc[3]+1;
      mx.old[(m3+m4+1):m]=mx.new[(m3+m4+1):m];
    }
  }
  
  if (rel((t4+2*t5)/2)<=b.2[1] || rel((t4+2*t5)/2)>=b.2[2] || rel((t4+15*t5)/15)<=b.2[1] || rel((t4+15*t5)/15)>=b.2[2] || rel((t4+50*t5)/50)<=b.2[1] || rel((t4+50*t5)/50)>=b.2[2]) {
    sorry=1;
  } else {
    
    X1=cbind(as.matrix(rep(z[1],m3)),cbind(as.matrix(rep(rel(2*(t4+2*t5)/4),m3)),D3[,1:5]));
    X2=cbind(as.matrix(rep(z[2],m4)),cbind(as.matrix(rep(rel(2*(t4+15*t5)/30),m4)),D4[,1:5]));
    X3=cbind(as.matrix(rep(z[3],m5)),cbind(as.matrix(rep(rel(2*(t4+50*t5)/100),m5)),D5[,1:5]));
    mx.new=c(m.x(X1),m.x(X2),m.x(X3));
    
    
    temp=c(z[1:3],t4,z[5:(7+m+1)],t5);
    
    a=logf(1.3,temp,mx.new,E)-logf(1.3,z,mx.old,E);
    if (a>=log(runif(1,0,1))) {
      z=temp;
      acc[4]=acc[4]+1;
      mx.old=mx.new;
    }
  }
  
  
  temp=c(z[1:4],zn[5]+z[5],z[6:(1+m+7)],z[7+m+2]);
  a=logf(2,temp,mx.old,E)-logf(2,z,mx.old,E)-temp[5]+z[5];
  if (a >=log(runif(1,0,1))) {
    z=temp;
    acc[5]=acc[5]+1;
  }
  
  temp=c(z[1:5],zn[6]+z[6],z[7:(1+m+7)],z[7+m+2]);
  a=logf(3,temp,mx.old,E)-logf(3,z,mx.old,E)-temp[6]+z[6];
  if (a >=log(runif(1,0,1))) {
    z=temp;
    acc[6]=acc[6]+1;
  }
  
  temp=c(z[1:6],z[7]+zn[7],z[8:(7+m+1)],z[7+m+2]);
  a=logf(4,temp,mx.new,E)-logf(4,z,mx.old,E)-temp[7]+z[7];
  if (a >=log(runif(1,0,1))) {
    z=temp;
    acc[7]=acc[7]+1;
  }
  
  temp[7]=z[7];
  temp[(7+1):(7+m3)]=as.vector(E[[1]]$vectors%*%diag(sqrt(abs(E[[1]]$values)))%*%as.matrix(zn[(7+1):(7+m3)]))+z[(7+1):(7+m3)];
  a=logf(4.1,temp,mx.old,E)-logf(4.1,z,mx.old,E);
  if (a>log(runif(1,0,1))) {
    z=temp;
    acc[8]=acc[8]+1;
  }
  
  temp[(7+1):(7+m3)]=z[(7+1):(7+m3)]
  temp[(7+1+m3):(7+m3+m4)]=as.vector(E[[2]]$vectors%*%diag(sqrt(abs(E[[2]]$values)))%*%as.matrix(zn[(7+1+m3):(7+m3+m4)]))+z[(7+1+m3):(7+m3+m4)];
  a=logf(4.2,temp,mx.old,E)-logf(4.2,z,mx.old,E);
  if (a>log(runif(1,0,1))) {
    z=temp;
    acc[9]=acc[9]+1;
  }
  
  temp[(7+1):(7+m3+m4)]=z[(7+1):(7+m3+m4)]
  temp[(7+1+m3+m4):(7+m)]=as.vector(E[[3]]$vectors%*%diag(sqrt(abs(E[[3]]$values)))%*%as.matrix(zn[(7+1+m3+m4):(7+m)]))+z[(7+1+m3+m4):(7+m)];
  a=logf(4.3,temp,mx.old,E)-logf(4.3,z,mx.old,E);
  if (a>log(runif(1,0,1))) {
    z=temp;
    acc[10]=acc[10]+1;
  }
  
  temp[(7+1+m3+m4):(7+m)]=z[(7+1+m3+m4):(7+m)];
  temp[7+m+1]=z[7+m+1]+zn[7+m+1];
  A=R.disc(log.inv(temp[7+m+1]),D3[,5]); E1=eigen(A,symmetric=TRUE);
  A=R.disc(log.inv(temp[7+m+1]),D4[,5]); E2=eigen(A,symmetric=TRUE);
  A=R.disc(log.inv(temp[7+m+1]),D5[,5]); E3=eigen(A,symmetric=TRUE);
  E.new=list(E1,E2,E3);
  a=logf(5,temp,mx.old,E.new)-logf(5,z,mx.old,E)+log(exp(-z[7+m+1])/(1+exp(-z[7+m+1]))^2)-log(exp(-temp[7+m+1])/(1+exp(-temp[7+m+1]))^2);
  if (a>log(runif(1,0,1))) {
    z[7+m+1]=temp[7+m+1];
    acc[11]=acc[11]+1;
    E=E.new;
  }
  
  
  
  if(i %% 1000==0) {
    c0=c0+(acc/1000>.2)*.5*c0-(acc/1000<.15)*.5*c0;
    print(acc/1000);
    acc=rep(0,11);
    stop=Sys.time();
    print(stop-start);
    start=stop;
  }
}
stop=Sys.time()
stop-start;


MC=matrix(0,ncol=(7+m+2),nrow=(M-burn));

for (i in 1:(M-burn)) {
  zn=rnorm(7+m+2,mean=0,sd=sqrt(c(c0[1:7],rep(c0[8],m3),rep(c0[9],m4),rep(c0[10],m5),c0[11],c0[4])));
  t1=z[1]+zn[1]; 
  t2=z[2]+zn[2]; 
  t3=z[3]+zn[3]; 
  t4=z[4]+zn[4];
  t5=z[7+m+2]+zn[7+m+2]; 
  
  ###do relative thickness and interphase modulus separately
  #X1=cbind(as.matrix(rep(temp[1],m3)),cbind(as.matrix(rep(temp[4],m3)),D3[,1:5]));
  #X2=cbind(as.matrix(rep(temp[2],m4)),cbind(as.matrix(rep(temp[4]/30,m4)),D4[,1:5]));
  #X3=cbind(as.matrix(rep(temp[3],m5)),cbind(as.matrix(rep(temp[4]/100,m5)),D5[,1:5]));
  #mx.new=c(m.x(X1),m.x(X2),m.x(X3));
  
  #a=logf(1.5,temp,mx.new)-logf(1.5,z,mx.old);
  #if (a>=log(runif(1,0,1))) {
  #  z=temp;
  #  acc[5]=acc[5]+1;
  #  mx.old=mx.new;
  #}
  
  
  
  
  X1=cbind(as.matrix(rep(t1,m3)),cbind(as.matrix(rep(rel(2*(z[4]+z[7+m+2]*2)/4),m3)),D3[,1:5]));
  X2=cbind(as.matrix(rep(t2,m4)),cbind(as.matrix(rep(rel(2*(z[4]+z[7+m+2]*15)/30),m4)),D4[,1:5]));
  X3=cbind(as.matrix(rep(t3,m5)),cbind(as.matrix(rep(rel(2*(z[4]+z[7+m+2]*50)/100),m5)),D5[,1:5]));
  if (t1<b.1[1] || t1>b.1[2]) {
    sorry=1;
  } else {
    mx.new=c(m.x(X1),mx.old[(m3+1):m]);
    temp=c(t1,z[2:(7+m+1)],z[7+m+2]);
    a=logf(1,temp,mx.new,E)-logf(1,z,mx.old,E);
    if (a>=log(runif(1,0,1))) {
      z[1]=temp[1];
      acc[1]=acc[1]+1;
      mx.old[1:m3]=mx.new[1:m3];
    }
  }
  
  if (t2<b.1[1] || t2>b.1[2]) {
    sorry=1;
  } else {
    mx.new=c(mx.old[1:m3],m.x(X2),mx.old[(m3+m4+1):(m3+m4+m5)]);
    temp=c(z[1],t2,z[3:(1+m+7)],z[7+m+2]);
    a=logf(1.1,temp,mx.new,E)-logf(1.1,z,mx.old,E);
    
    if (a>=log(runif(1,0,1))) {
      z[2]=temp[2];
      acc[2]=acc[2]+1;
      mx.old[(m3+1):(m3+m4)]=mx.new[(m3+1):(m3+m4)];
    }
  }
  
  if (t3<b.1[1] || t3>b.1[2]) {
    sorry=1;
  } else {
    mx.new=c(mx.old[1:(m3+m4)],m.x(X3));
    temp=c(z[1:2],t3,z[4:(1+m+7)],z[7+m+2]);
    a=logf(1.2,temp,mx.new,E)-logf(1.2,z,mx.old,E);
    if (a>=log(runif(1,0,1))) {
      z[3]=temp[3];
      acc[3]=acc[3]+1;
      mx.old[(m3+m4+1):m]=mx.new[(m3+m4+1):m];
    }
  }
  
  if (rel((t4+2*t5)/2)<=b.2[1] || rel((t4+2*t5)/2)>=b.2[2] || rel((t4+15*t5)/15)<=b.2[1] || rel((t4+15*t5)/15)>=b.2[2] || rel((t4+50*t5)/50)<=b.2[1] || rel((t4+50*t5)/50)>=b.2[2]) {
    sorry=1;
  } else {
    X1=cbind(as.matrix(rep(z[1],m3)),cbind(as.matrix(rep(rel(2*(t4+2*t5)/4),m3)),D3[,1:5]));
    X2=cbind(as.matrix(rep(z[2],m4)),cbind(as.matrix(rep(rel(2*(t4+15*t5)/30),m4)),D4[,1:5]));
    X3=cbind(as.matrix(rep(z[3],m5)),cbind(as.matrix(rep(rel(2*(t4+50*t5)/100),m5)),D5[,1:5]));
    mx.new=c(m.x(X1),m.x(X2),m.x(X3));
    
    
    temp=c(z[1:3],t4,z[5:(7+m+1)],t5);
    
    a=logf(1.3,temp,mx.new,E)-logf(1.3,z,mx.old,E);
    if (a>=log(runif(1,0,1))) {
      z=temp;
      acc[4]=acc[4]+1;
      mx.old=mx.new;
    }
  }
  
  
  temp=c(z[1:4],zn[5]+z[5],z[6:(1+m+7)],z[7+m+2]);
  a=logf(2,temp,mx.old,E)-logf(2,z,mx.old,E)-temp[5]+z[5];
  if (a >=log(runif(1,0,1))) {
    z=temp;
    acc[5]=acc[5]+1;
  }
  
  temp=c(z[1:5],zn[6]+z[6],z[7:(1+m+7)],z[7+m+2]);
  a=logf(3,temp,mx.old,E)-logf(3,z,mx.old,E)-temp[6]+z[6];
  if (a >=log(runif(1,0,1))) {
    z=temp;
    acc[6]=acc[6]+1;
  }
  
  temp=c(z[1:6],z[7]+zn[7],z[8:(7+m+1)],z[7+m+2]);
  a=logf(4,temp,mx.new,E)-logf(4,z,mx.old,E)-temp[7]+z[7];
  if (a >=log(runif(1,0,1))) {
    z=temp;
    acc[7]=acc[7]+1;
  }
  
  temp[7]=z[7];
  temp[(7+1):(7+m3)]=as.vector(E[[1]]$vectors%*%diag(sqrt(abs(E[[1]]$values)))%*%as.matrix(zn[(7+1):(7+m3)]))+z[(7+1):(7+m3)];
  a=logf(4.1,temp,mx.old,E)-logf(4.1,z,mx.old,E);
  if (a>log(runif(1,0,1))) {
    z=temp;
    acc[8]=acc[8]+1;
  }
  
  temp[(7+1):(7+m3)]=z[(7+1):(7+m3)]
  temp[(7+1+m3):(7+m3+m4)]=as.vector(E[[2]]$vectors%*%diag(sqrt(abs(E[[2]]$values)))%*%as.matrix(zn[(7+1+m3):(7+m3+m4)]))+z[(7+1+m3):(7+m3+m4)];
  a=logf(4.2,temp,mx.old,E)-logf(4.2,z,mx.old,E);
  if (a>log(runif(1,0,1))) {
    z=temp;
    acc[9]=acc[9]+1;
  }
  
  temp[(7+1):(7+m3+m4)]=z[(7+1):(7+m3+m4)]
  temp[(7+1+m3+m4):(7+m)]=as.vector(E[[3]]$vectors%*%diag(sqrt(abs(E[[3]]$values)))%*%as.matrix(zn[(7+1+m3+m4):(7+m)]))+z[(7+1+m3+m4):(7+m)];
  a=logf(4.3,temp,mx.old,E)-logf(4.3,z,mx.old,E);
  if (a>log(runif(1,0,1))) {
    z=temp;
    acc[10]=acc[10]+1;
  }
  
  temp[(7+1+m3+m4):(7+m)]=z[(7+1+m3+m4):(7+m)];
  temp[7+m+1]=z[7+m+1]+zn[7+m+1];
  A=R.disc(log.inv(temp[7+m+1]),D3[,5]); E1=eigen(A,symmetric=TRUE);
  A=R.disc(log.inv(temp[7+m+1]),D4[,5]); E2=eigen(A,symmetric=TRUE);
  A=R.disc(log.inv(temp[7+m+1]),D5[,5]); E3=eigen(A,symmetric=TRUE);
  E.new=list(E1,E2,E3);
  a=logf(5,temp,mx.old,E.new)-logf(5,z,mx.old,E)+log(exp(-z[7+m+1])/(1+exp(-z[7+m+1]))^2)-log(exp(-temp[7+m+1])/(1+exp(-temp[7+m+1]))^2);
  if (a>log(runif(1,0,1))) {
    z[7+m+1]=temp[7+m+1];
    acc[11]=acc[11]+1;
    E=E.new;
  }
  
  
           
           if(i %% 1000==0) {
             
             print(acc/1000);
             acc=rep(0,11);
             stop=Sys.time();
             print(stop-start);
             start=stop;
           }  

  MC[i,]=z;
}

X1=cbind(rep(Mode(MC[,1]),m3),cbind(rep(rel((Mode(MC[,4]+2*MC[,7+m+2])/2)),m3),D3[,1:5]));
X2=cbind(rep(Mode(MC[,2]),m4),cbind(rep(rel((Mode(MC[,4]+15*MC[,7+m+2])/15)),m4),D4[,1:5]));
X3=cbind(rep(Mode(MC[,3]),m5),cbind(rep(rel((Mode(MC[,4]+50*MC[,7+m+2])/50)),m5),D5[,1:5]));

par(mfrow=c(2,3),mar=c(3,2,3,2));
disc.m.3=apply(MC[,(7+1):(7+m3)],2,Mode);
disc.m.4=apply(MC[,(7+m3+1):(7+m3+m4)],2,Mode);
disc.m.5=apply(MC[,(7+m3+m4+1):(7+m)],2,Mode);
plot(D3[,5],y.exp[1:m3],type='p',pch='+',cex=2); lines(D3[,5],m.x(X1)+disc.m.3,type='p',pch='*',col='blue',cex=2);
plot(D4[,5],y.exp[(m3+1):(m3+m4)],type='p',pch='+',cex=2); lines(D4[,5],m.x(X2)+disc.m.4,type='p',pch='*',col='blue',cex=2);
plot(D5[,5],y.exp[(m3+m4+1):(m3+m4+m5)],type='p',pch='+',cex=2); lines(D5[,5],m.x(X3)+disc.m.5,type='p',pch='*',col='blue',cex=2);

results.3.un=MC;


MC.trans=MC;
MC.trans[,1]=means[1]+MC[,1]*sds[1];
MC.trans[,2]=means[1]+MC[,2]*sds[1];
MC.trans[,3]=means[1]+MC[,3]*sds[1];

par=rep(0,8);
par[1:3]=apply(MC.trans[,1:3],2,Mode)
par[4:5]=apply(MC[,c(4,7+m+2)],2,Mode);
par[6]=Mode(MC[,4]+2*MC[,7+m+2]);
par[7]=Mode(MC[,4]+15*MC[,7+m+2]);
par[8]=Mode(MC[,4]+50*MC[,7+m+2]);
par_g_pp_cv=par;
names(par_g_pp_cv)=c("Int Mod 4","Int Mod 30","Int Mod 100","Intercept","Slope","Abs Th 4","Abs Th 30","Abs Th 100")

sd.s=rep(0,8);
sd.s[1:3]=apply(MC.trans[which(1:5000%%50==0),1:3],2,sd)
sd.s[4:5]=apply(MC[which(1:5000%%50==0),c(4,7+m+2)],2,sd);
sd.s[6]=sd(MC[which(1:5000%%50==0),4]+2*MC[which(1:5000%%50==0),7+m+2]);
sd.s[7]=sd(MC[which(1:5000%%50==0),4]+15*MC[which(1:5000%%50==0),7+m+2]);
sd.s[8]=sd(MC[which(1:5000%%50==0),4]+50*MC[which(1:5000%%50==0),7+m+2]);
sd.g_pp_cv=sd.s;
names(sd.g_pp_cv)=c("Int Mod 4","Int Mod 30","Int Mod 100","Intercept","Slope","Abs Th 4","Abs Th 30","Abs Th 100")


par(mfrow=c(2,3),mar=c(4,4,2,1))

hist(MC.trans[,1],main="4 Micron",xlab="Interphase Modulus (GPa)",freq=FALSE)
hist(MC.trans[,2],main="30 Micron",xlab="Interphase Modulus(GPa)",freq=FALSE)
hist(MC.trans[,3],main="100 Micron",xlab="Interphase Modulus (GPa)",freq=FALSE)
hist(MC[,4]+2*MC[,7+m+2],main="",xlab="Absolute Thickness (microns)",freq=F);
hist(MC[,4]+15*MC[,7+m+2],main="",xlab="Absolute Thickness (microns)",freq=F);
hist(MC[,4]+50*MC[,7+m+2],main="",xlab="Absolute Thickness (microns)",freq=F);



par(mfrow=c(1,2),mar=c(4,4,2,1))
hist(MC[,4],main="Intercept Parameter",xlab="Intercept Parameter",freq=FALSE)
hist(MC[,7+m+2],main="Slope Parameter",xlab="Slope Parameter",freq=FALSE)




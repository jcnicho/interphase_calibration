library(readxl)
E_avg_Glass_Polypropylene_081720 <- read_excel("E_avg_Glass-Polypropylene_081720.xlsx")
E_avg_Glass_Polypropylene_071320 <- read_excel("E_avg_Glass-Polypropylene_071320.xlsx")
E_avg_Glass_Polystyrene_081720 <- read_excel("E_avg_Glass-Polystyrene_081720.xlsx")
E_avg_Glass_Polycarbonate_071320 <- read_excel("E_avg_Glass-Polycarbonate_071320.xlsx")
E_avg_Glass_Styrene_acrylonitrile_copolymer_071320 <- read_excel("E_avg_Glass-Styrene-acrylonitrile copolymer_071320.xlsx")

data.1=rbind(E_avg_Glass_Polystyrene_081720)
rm(E_avg_Glass_Polystyrene_081720)

data.1=as.matrix(data.1[,4:12]);
D=as.matrix(data.1[,c(1,2,7)]); ycomp=as.vector(data.1[,8]);


data.2=E_avg_Alumina_Polystyrene_101020_1_;
rm(E_avg_Alumina_Polystyrene_101020_1_);

data.2=as.matrix(data.2[,4:12]);
D=as.matrix(data.2[,c(1,2,7)]); y1=as.vector(data.2[,8]);

data.3=E_avg_Glass_Polycarbonate_071320
rm(E_avg_Glass_Polycarbonate_071320);

data.3=as.matrix(data.3[,4:12]);
D=as.matrix(data.3[,c(1,2,7)]); ycomp=as.vector(data.3[,8]);

data.4=rbind(E_avg_Alumina_Polycarbonate_071320,E_avg_Alumina_Polycarbonate_081720)
rm(E_avg_Alumina_Polycarbonate_081720);

data.4=as.matrix(data.4[,4:12]);
D=as.matrix(data.4[,c(1,2,7)]); y1=as.vector(data.4[,8]);

data.5=E_avg_Glass_Styrene_acrylonitrile_copolymer_071320
rm(E_avg_Glass_Styrene_acrylonitrile_copolymer_071320);

data.5=as.matrix(data.5[,4:12]);
D=as.matrix(data.5[-117,c(1,2,7)]); ycomp=as.vector(data.5[-117,8]);

data.6=E_avg_Alumina_Styrene_acrylonitrile_copolymer_101020_1_
rm(E_avg_Alumina_Styrene_acrylonitrile_copolymer_101020_1_);

data.6=as.matrix(data.6[,4:12]);
D=as.matrix(data.6[,c(1,2,7)]); y1=as.vector(data.6[,8]);
D.add=cbind(seq(.01,50,length=20),0,0); y1.add=rep(3.8,20);
D=rbind(D,D.add); y1=c(y1,y1.add);

data.7=rbind(E_avg_Glass_Polypropylene_071320,E_avg_Glass_Polypropylene_081720);
rm(E_avg_Glass_Polypropylene_081720,E_avg_Glass_Polypropylene_071320);

data.7=as.matrix(data.7[,4:12]);
D=as.matrix(data.7[,c(1,2,7)]); ycomp=as.vector(data.7[,8]);


ycomp=scale(ycomp,center=TRUE,scale=TRUE);
D=scale(D,center=TRUE,scale=TRUE);

means=c(attr(D,"scaled:center"),attr(ycomp,"scaled:center"));
sds=c(attr(D,"scaled:scale"),attr(ycomp,"scaled:scale"));


#D=rbind(D,fixed.pts[,1:7]);
#y1=c(y1,fixed.pts[,8]);

D=as.matrix(D);
for (i in 1:3 ) {
  D[,i]=(D[,i]-means[i])/sds[i];
}
ycomp=(ycomp-means[4])/sds[4]



Data.exp=as.matrix(read.csv("exp.data.csv"))

s=c(1,3,2,4,5,6);
D1=Data.exp[1:5,1:6]; D1=D1[,s]; D1[,5]=D1[,5]/100;
D2=Data.exp[1:5,8:13]; D2=D2[,s]; D2[,5]=D2[,5]/100;
D3=Data.exp[11:14,15:20]; D3=D3[,s]; D3[,5]=D3[,5]/100; 
D4=Data.exp[12:15,22:27]; D4=D4[,s]; D4[,5]=D4[,5]/100;
D5=Data.exp[13:16,29:34]; D5=D5[,s]; D5[,5]=D5[,5]/100;

for (i in c(5:6)) {
  D1[,i]=(D1[,i]-means[i-2])/sds[i-2];
  D2[,i]=(D2[,i]-means[i-2])/sds[i-2];
  D3[,i]=(D3[,i]-means[i-2])/sds[i-2];
  D4[,i]=(D4[,i]-means[i-2])/sds[i-2];
  D5[,i]=(D5[,i]-means[i-2])/sds[i-2];
}

Dex = c(D3[, 5], D4[, 5], D5[, 5]);
sizes = c(rep(log(4), dim(D3)[1]), rep(log(30), dim(D4)[1]), rep(log(100), dim(D5)[1]));
Dex = cbind(Dex, sizes);





b.1=(c(0,50)-means[1])/sds[1];
b.2=(c(0,.25)-means[2])/sds[2];


calib=list();

calib[[1]]=as.matrix(calib[[1]][,4:12]);
D.cal=as.matrix(calib[[1]][,-c(8,9)]); y.cal=as.vector(calib[[1]][,8]);

for (i in 1:2 ) {
  D.cal[,i]=(D.cal[,i]-means[i])/sds[i];
}
D.cal[,7]=(D.cal[,7]-means[3])/sds[3];
y.cal=(y.cal-means[4])/sds[4]



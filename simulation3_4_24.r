incir = 1000
outcir = 1000
#生成两个1000 X 2 的矩阵，用于存放分位数
quantile.R1 = matrix(NA,outcir,8)
quantile.R2 = matrix(NA,outcir,8)
quantile.R = matrix(NA,outcir,8)

#分位数的取值,相对于incir次内循环
a=0.00204;b=3;

quantil = c(1/incir,0.025,0.05,0.1,0.9,0.95,0.975,1)

s0=0;s1=0.25;s2=0.5;s3=0.75;s4=1;

mu=c(a+b*s0,a+b*s1,a+b*s2,a+b*s3,a+b*s4)

sigmaa = sqrt(0.0000001764)
sigmab = sqrt(0.0000010816)

#delta
delt = 250;

#epsilon的矩阵行列数
input.row = c(10,20,30,40); input.col =15;

omega = 10
t = 3876

sto = function(mrow,mcol,mu,sigmaa,sigmab,delta){ 
	i= mrow; 
	j = mcol; 
	mu=mu;
	sigma.a = sigmaa;
	sigma.b = sigmab;
	delta.t = delta;
	Y = matrix(NA,i,j)
	mu.matrix = matrix(data = mu , nrow = i,ncol = j );
	alpha = matrix(NA,i,j); alpha.pre = matrix(data = 1,i,j);epsilon = matrix(rnorm(i*j,mean = 0, sd = sigma.b / sqrt(delta.t)),i,j);
	sigma.matrix = matrix(rnorm(i,mean = 0,sd = sigma.a),1,i);
	for (t in 1:i){
		alpha[t,] = alpha.pre[t,]*sigma.matrix[1,t]}
	Y = mu.matrix + alpha + epsilon;
	Y.row.bar = matrix(rowMeans(Y),1,i);
	Y.bar = matrix(rowMeans(Y.row.bar),1,i);
	Y.bar.int = Y.bar[1]
	se.pre = matrix(1,i,j);
	for (r in 1:i){
		se.pre[r,] = se.pre[r,]*Y.row.bar[1,r]};
	SE = sum(rowSums((Y - se.pre)^2));	
	SB = rowSums((Y.row.bar - Y.bar)^2) * j;
	return(c(SE,SB,Y.bar.int,Y.row.bar))
}
	

for(k in 1:outcir){

#h = 1
#while(R1.matrix[1000] < 0){
#	d = sto(input.row,input.col)
#	if (d[1] > 0){
#	R1.matrix[h] = d[1]
#	R2.matrix[h] = d[2]
#	h = h + 1
#	}
#}

R1 = c()
R2 = c()
R = c()

d1 = sto(input.row[1],input.col,mu[2],sigmaa,sigmab,delt)
d2 = sto(input.row[2],input.col,mu[3],sigmaa,sigmab,delt)
d3 = sto(input.row[3],input.col,mu[4],sigmaa,sigmab,delt)
d4 = sto(input.row[4],input.col,mu[5],sigmaa,sigmab,delt)

sse=sum(c(d1[1],d2[1],d3[1],d4[1]))
ssb=sum(c(d1[2],d2[2],d3[2],d4[2]))
F=sum(c(input.row[1]*s1^2,input.row[2]*s2^2,input.row[3]*s3^2,input.row[4]*s4^2))
G=sum(c(input.row[1]*s1,input.row[2]*s2,input.row[3]*s3,input.row[4]*s4))
H=sum(c(input.row[1]*s1*d1[3],input.row[2]*s2*d2[3],input.row[3]*s3*d3[3],input.row[4]*s4*d4[3]))
J=sum(c(input.row[1]*d1[3],input.row[2]*d2[3],input.row[3]*d3[3],input.row[4]*d4[3]))
N=sum(input.row)
a.hat=(J*F-H*G)/(N*F-G^2)
b.hat=(N*H-J*G)/(N*F-G^2)
mu0.hat=a.hat
z1=matrix(d1[4:length(d1)],length(d1)-3,1)
z2=matrix(d2[4:length(d2)],length(d2)-3,1)
z3=matrix(d3[4:length(d3)],length(d3)-3,1)
z4=matrix(d4[4:length(d4)],length(d4)-3,1)
ab1=matrix(a.hat+b.hat*s1,length(d1)-3,1)
ab2=matrix(a.hat+b.hat*s2,length(d2)-3,1)
ab3=matrix(a.hat+b.hat*s3,length(d3)-3,1)
ab4=matrix(a.hat+b.hat*s4,length(d4)-3,1)

A.in=sum(c(sum((z1-ab1)^2),sum((z2-ab2)^2),sum((z3-ab3)^2),sum((z4-ab4)^2)))

A=F/(N*F-G^2)*A.in/(N-2)

while(length(R1)<incir){
	u = rchisq(1,df = N*input.col-N);
	v = rchisq(1,df = N - 4);
	w = rt(1,df=N-2,ncp=0);
	R.1 = ((ssb / v) - (sse / u)) / input.col;
	R.2 = a.hat - (w * sqrt(A));
	R.3 = delt*sse/u
	h1 = (R.2*t - omega)/(sqrt(R.3*t + R.1*t^2))
	h2 = (2*R.2*omega/R.3) + (2*R.1*omega^2 / R.3^2)
	h3 = - (2*R.1*omega*t + R.3*(R.2*t + omega)) / (R.3*sqrt(R.3*t + R.1*t^2))
	R.pre= 1 - pnorm(h1) - exp(h2)*pnorm(h3)
	if(R.1>0  && is.na(R.pre)==0 ){
	R1 = append(R1,R.1)
	R2 = append(R2,R.2)
	R = append(R,R.pre)
	}
	if(R.1>0  && is.na(R.pre)==1 ){
	R1 = append(R1,R.1)
	R2 = append(R2,R.2)
        R.pre= 1 - pnorm(h1)
	R = append(R,R.pre)
	}
	}


R1 = sort(R1)
R2 = sort(R2)
R = sort(R)


for(b in seq(8)){
	quantile.R1[k,b] = R1[quantil[b]*incir]
	quantile.R2[k,b] = R2[quantil[b]*incir]
	quantile.R[k,b] = R[quantil[b]*incir]
	}

}

l = matrix(0,3,6)
jj = 0.9
for(h in 1:outcir){
	if(quantile.R1[h,6] > sigmaa^2 && quantile.R1[h,3] < sigmaa^2){
	l[1,1] = l[1,1]+1
	}
	if(quantile.R1[h,7] > sigmaa^2 && quantile.R1[h,2] < sigmaa^2){
	l[1,2] = l[1,2]+1
	}
	if(quantile.R1[h,8] > sigmaa^2 && quantile.R1[h,4] < sigmaa^2){
	l[1,3] = l[1,3]+1
	}
	if(quantile.R1[h,8] > sigmaa^2 && quantile.R1[h,3] < sigmaa^2){
	l[1,4] = l[1,4]+1
	}
	if(quantile.R1[h,5] > sigmaa^2 && quantile.R1[h,1] < sigmaa^2){
	l[1,5] = l[1,5]+1
	}
	if(quantile.R1[h,6] > sigmaa^2 && quantile.R1[h,1] < sigmaa^2){
	l[1,6] = l[1,6]+1
	}
	
	if(quantile.R2[h,6] > a && quantile.R2[h,3] < a){
	l[2,1] = l[2,1]+1
	}
	if(quantile.R2[h,7] > a && quantile.R2[h,2] < a){
	l[2,2] = l[2,2]+1
	}
	if(quantile.R2[h,8] > a && quantile.R2[h,4] < a){
	l[2,3] = l[2,3]+1
	}
	if(quantile.R2[h,8] > a && quantile.R2[h,3] < a){
	l[2,4] = l[2,4]+1
	}
	if(quantile.R2[h,5] > a && quantile.R2[h,1] < a){
	l[2,5] = l[2,5]+1
	}
	if(quantile.R2[h,6] > a && quantile.R2[h,1] < a){
	l[2,6] = l[2,6]+1
	}

	if(quantile.R[h,6] > jj && quantile.R[h,3] < jj){
	l[3,1] = l[3,1]+1
	}
	if(quantile.R[h,7] > jj && quantile.R[h,2] < jj){
	l[3,2] = l[3,2]+1
	}
	if(quantile.R[h,8] > jj && quantile.R[h,4] < jj){
	l[3,3] = l[3,3]+1
	}
	if(quantile.R[h,8] > jj && quantile.R[h,3] < jj){
	l[3,4] = l[3,4]+1
	}
	if(quantile.R[h,5] > jj && quantile.R[h,1] < jj){
	l[3,5] = l[3,5]+1
	}
	if(quantile.R[h,6] > jj && quantile.R[h,1] < jj){
	l[3,6] = l[3,6]+1
	}
	
       v11=(quantile.R1[h,6]-quantile.R1[h,3] )
       v12=(quantile.R1[h,7]-quantile.R1[h,2])
       v13=(quantile.R1[h,4])
       v14=(quantile.R1[h,3])
       v15=(quantile.R1[h,5])
       v16=(quantile.R1[h,6])
       
       v21=(quantile.R2[h,6]-quantile.R2[h,3] )
       v22=(quantile.R2[h,7]-quantile.R2[h,2])
       v23=(quantile.R2[h,4])
       v24=(quantile.R2[h,3])
       v25=(quantile.R2[h,5])
       v26=(quantile.R2[h,6])
       
       v31=(quantile.R[h,6]-quantile.R[h,3] )
       v32=(quantile.R[h,7]-quantile.R[h,2])
       v33=(quantile.R[h,4])
       v34=(quantile.R[h,3])
       v35=(quantile.R[h,5])
       v36=(quantile.R[h,6])

       
}

rownames(l) = c("R1","R2","R3")
colnames(l) = c("90%intv","95%intv","10%dlim","5%dlim","90%ulim","95%ulim")
e=matrix(0,3,6)
e[1,1]=mean(v11);e[1,2]=mean(v12);e[1,3]=mean(v13);e[1,4]=mean(v14);e[1,5]=mean(v15);e[1,6]=mean(v16);

e[2,1]=mean(v21);e[2,2]=mean(v22);e[2,3]=mean(v23);e[2,4]=mean(v24);e[2,5]=mean(v25);e[2,6]=mean(v26);

e[3,1]=mean(v31);e[3,2]=mean(v32);e[3,3]=mean(v33);e[3,4]=mean(v34);e[3,5]=mean(v35);e[3,6]=mean(v36);













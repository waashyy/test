incir = 10000
outcir = 10000
#生成两个1000 X 2 的矩阵，用于存放分位数
quantile.R1 = matrix(NA,outcir,8)
quantile.R2 = matrix(NA,outcir,8)
quantile.R = matrix(NA,outcir,8)

#分位数的取值,相对于incir次内循环

quantil = c(1/incir,0.025,0.05,0.1,0.9,0.95,0.975,1)

mu =0.00204
sigma = sqrt(0.0000001764)
sigmab = sqrt(0.0000010816)

#delta
delt = 250;

#epsilon的矩阵行列数
input.row = 10; input.col = 15;

omega = 10
t = 3876



sto = function(mrow,mcol,mu,sigmaa,sigmab,delta){ 
	i= mrow; j = mcol; 
	mu=mu;
	sigma.a = sigmaa;
	sigma.b = sigmab;
	delta.t = delta;
	Y = matrix(NA,i,j)
	x = matrix(rnorm(i*j,0,sigma.b),i,j)

	xt=matrix(0,i,1)

	for(o in 1:j){
		xt.temp = matrix(rnorm(i,o*mu*delta.t,sqrt(o*delta.t)*sigma.a),i,1)
		xt=cbind(xt,xt.temp)
	}
	xt=xt[,-1]

	phase1 = x+xt
	phase2 = matrix(phase1[,1],i,1)
	for(o in 1:j-1){
		phase2=cbind(phase2,(phase1[,o+1]-phase1[,o]))
	}
	
	rs = rowSums(phase2)
	rs_bar = mean(rs)
	
	a22 = rs-rs_bar
	
	A22 = sum(a22^2)
	
	phase3=phase2[,-1]
	rs3 = rowSums(phase3)
	rs3_bar = mean(rs3)
	
	a11 = rs3-rs3_bar
	
	A11 = sum(a11^2)
	
	A12 = sum(a22*a11)
	
	return(c(A11,A12,A22))
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
d = sto(input.row,input.col,mu,sigma,sigmab,delt)
while(length(R1)<incir){
	u = rchisq(1,df = input.row*input.col-input.row);
	v = rchisq(1,df = input.row - 1);
	w = rnorm(1,mean = 0, sd = 1);
	R.1 = ((d[2] / v) - (d[1] / u)) / input.col;
	R.2 = d[3] - ((w / sqrt(input.row*input.col))*sqrt(d[2] / v));
	R.3 = delt*d[1]/u
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
	
	if(quantile.R2[h,6] > mu && quantile.R2[h,3] < mu){
	l[2,1] = l[2,1]+1
	}
	if(quantile.R2[h,7] > mu && quantile.R2[h,2] < mu){
	l[2,2] = l[2,2]+1
	}
	if(quantile.R2[h,8] > mu && quantile.R2[h,4] < mu){
	l[2,3] = l[2,3]+1
	}
	if(quantile.R2[h,8] > mu && quantile.R2[h,3] < mu){
	l[2,4] = l[2,4]+1
	}
	if(quantile.R2[h,5] > mu && quantile.R2[h,1] < mu){
	l[2,5] = l[2,5]+1
	}
	if(quantile.R2[h,6] > mu && quantile.R2[h,1] < mu){
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













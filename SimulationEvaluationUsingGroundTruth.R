# check individual time usage part 1
#Rprof(tmp <- tempfile())

time_start=proc.time()[[3]]

#------------configuration begin--------------
n=100 #sample size in simulation
loopNum=200 #repeat x times


## for local
library(EntropyEstimation)
library(entropy)
library(praznik)

#------------configuration end--------------



#MI.test tests the independence between a feature and the outcome. It prints out p.value, the smaller the p.value, the stronger evidence of dependence between them.
MI.test=function(feature, outcome, k1, k2){ #feature is the vector of X, outcome is the vector of Y, k1 and k2 are the corresponding number of categories, X and Y must have the same sample size.
  n=length(feature);
  test.stat=2*n*MI.z(table(feature, outcome))+(k1-1)*(k2-1);
  p.value=pchisq(test.stat,(k1-1)*(k2-1),lower.tail = F);
  return(p.value);
}

CASMI.featureSelection=function(data){ #outcome must be in the last column
  # Step 1: return dependent features, step1Index[]

  alpha=0.1
  
  step1Index=vector()
  count=0
  k2=length(table(data[length(data)])) #outcome categories
  for(fi in 1:(length(data)-1)){
    k1=length(table(data[fi])) #feature categories
    p_value=MI.test(data[[fi]],data[[length(data)]], k1, k2)
    if(p_value<=alpha){
      count=count+1
      step1Index[count]=fi
    }
  }
  
  indexCASMI=vector()

  if(count==0){print("Warning: No dependent feature exists. Step 2 skipped.")
  } else{
    # Step 2: select features by joint SMI with Hz estimator
    # select the best feature
    maxKappaStar=0
    maxIndex=0
    for(index in step1Index){
      #valueKappaStar=kappa.star(data[[index]], data[[length(data)]])
      feature=data[[index]]; outcome=data[[length(data)]];
      
      #valueKappaStar=MI.z(table(feature, outcome))/Entropy.z(table(outcome))*(1-length(which(table(feature)==1))/n)
      valueKappaStar=MI.z(table(feature, outcome))/Entropy.z(table(outcome))*(1-sum(table(feature)==1L)/n)
      
      if(valueKappaStar > maxKappaStar){
        maxKappaStar=valueKappaStar
        maxIndex=index
      }
    }
    indexCASMI=c(indexCASMI, maxIndex)
    
    # select the 2nd, 3rd, ... best features by joint
    maxmaxKappaStar=0
    while(maxKappaStar>maxmaxKappaStar & length(indexCASMI)<count){
      maxmaxKappaStar=maxKappaStar
      
      step1Index=step1Index[!step1Index==maxIndex]
      maxKappaStar=0
      maxIndex=0
      for(index in step1Index){
        tmpIndex1=c(index, indexCASMI, length(data))
        tmpIndex2=c(index, indexCASMI)
        ftable=ftable(table(data[tmpIndex1]))
        ftableOfFeatures=ftable(table(data[tmpIndex2]))
        #valueKappaStar=kappa.star2(ftable, ftableOfFeature, data[[length(data)]])
        outcome=data[[length(data)]]
        
        #valueKappaStar=MI.z(ftable)/Entropy.z(table(outcome))*(1-length(which(ftableOfFeatures==1))/n)
        valueKappaStar=MI.z(ftable)/Entropy.z(table(outcome))*(1-sum(ftableOfFeatures==1L)/n)
        
        if(valueKappaStar > maxKappaStar){
          maxKappaStar=valueKappaStar
          maxIndex=index
        }
      }
      
      if(maxKappaStar>maxmaxKappaStar+10^-14){ # +10^-14 is to solve the problem of precision
        indexCASMI=c(indexCASMI, maxIndex)
      }
    }
  }

  return(indexCASMI)
}


#Warning: this function is model dependent. Check the form of y.
#replace x6 by x4, or remove x6 if x4 is selected
#remove y independent features
processFeatures=function(selectedFeatures){
  if(6 %in% selectedFeatures){
    selectedFeatures[which(selectedFeatures==6)]=4
  }
  if(length(which(selectedFeatures==4))==2){
    selectedFeatures=selectedFeatures[-(which(selectedFeatures==4)[2])]
  }
  unrelevant_x=c(7,8,9,10)
  if(any(selectedFeatures %in% unrelevant_x)){
    selectedFeatures=selectedFeatures[-(which(selectedFeatures %in% unrelevant_x))]
  }
  return(selectedFeatures)
}


#entropy of x
Hx=function(idx_selected_x, p_x){
  p_selected_x=p_x[idx_selected_x];
  xProbTable=Reduce(outer, p_selected_x)
  return(-sum(xProbTable*log(xProbTable)))
}

#combined mutual information of selected features
#p_x: entire features with error; xList: entire features with error; xTerms: calculated terms of features with error term; idx_selected_x: indices of selected features
#Input features (p_selected_x or p_x) must be independent. (xTerms) are calculated terms of x, y=term(x1)+term(x2)+...
#p_x & xTerms must have same length
MIwPlus=function(p_x, xList, xTerms, idx_selected_x){
  xProbTable=Reduce(outer, p_x)
  yValueTable=Reduce(function (x, y) outer(x, y, "+"), xTerms) #y=term(x1)+term(x2)+... + error
  p=as.vector(xProbTable)
  y=as.vector(yValueTable)
  
  uniq_y=unique(y)
  yProbVec=vector()
  for(i in 1:length(uniq_y)){
    yProbVec=c(yProbVec,sum(p[which(y==unique(y)[i])]))
  }
  Hy=-sum(yProbVec*log(yProbVec))
  
  table=cbind(expand.grid(xList),y,p)
  s1=c(idx_selected_x,(ncol(table)-1),ncol(table))
  table2=table[s1]
  c=c(1:(length(s1)-1))
  tableGroup=aggregate(p ~ ., data= table2[, c('p', names(table2)[c])], sum)
  p_xy=tableGroup[length(tableGroup)]
  Hxy=-sum(p_xy*log(p_xy))
  
  return(Hx(idx_selected_x, p_x) + Hy - Hxy)
}

#---------------------------------main---------------------------------

#population
#calculate probability of x
p_x1=vector(); p_x1[1]=pnorm(-3); p_x1[2]=pnorm(-0.5)-pnorm(-3); p_x1[3]=pnorm(0.5)-pnorm(-0.5); p_x1[4]=pnorm(3)-pnorm(0.5); p_x1[5]=1-pnorm(3);
p_x2=vector(); p_x2[1]=ppois(0,2); p_x2[2]=ppois(1,2)-ppois(0,2); p_x2[3]=ppois(2,2)-ppois(1,2); p_x2[4]=ppois(4,2)-ppois(2,2); p_x2[5]=1-ppois(4,2); 
p_x3=vector(); p_x3[1:5]=0.2
p_x4=vector(); p_x4[1]=pbinom(0,4,0.1); for(i in 1:4){p_x4[i+1]=pbinom(i,4,0.1)-pbinom(i-1,4,0.1)};
p_x5=vector(); p_x5[1]=pnorm(-0.5); p_x5[2]=pnorm(-0.2)-pnorm(-0.5); p_x5[3]=pnorm(0.2)-pnorm(-0.2); p_x5[4]=pnorm(0.6)-pnorm(0.2); p_x5[5]=1-pnorm(0.6);
#p_x6=p_x4; #x6=x4
#p_x7=p_x2; p_x8=p_x3; p_x9=p_x4; p_x10=p_x5; #distributions are the same even though the variables are independent
p_error=c(1/3,1/3,1/3)

#save values of x
v_x1=c(-3.5,-1.4,0,1,2.2)
v_x2=c(-5,-3,0,2.4,5.4)
v_x3=c(-2,-1,0,1,2)
v_x4=c(-2,-1,0,1,5)
v_x5=c(-2.5,-2,1.7,2,4)
v_x6=v_x4
v_error=c(-1,0,1)

#experiment 1
#calculate total mutual information H(x1,x2,x3,x4,x5,x6)
p_x_data1=list(p_x1,p_x2,p_x3,p_x4,p_x5,p_error) #x6=x4 so ignored
v_x_data1=list(v_x1,v_x2,v_x3,v_x4,v_x5,v_error)
v_xTerms_data1=list(v_x1, v_x2, v_x3^3, (-0.5*v_x4^2+v_x6), abs(v_x5), v_error) #x6=x4 so x6 term merged to x4 #y1=x1+x2+x3^3-0.5*x4^2+abs(x5)+x6+error
totalFeature1=c(1:6); total1=processFeatures(totalFeature1);

totalMI1=MIwPlus(p_x_data1, v_x_data1, v_xTerms_data1, total1) #input features must be independent, 2nd & 3nd para. must have same length


#MIRtable BEGIN---------------
IIRtableCol1=list(c(1), c(2), c(3), c(4), c(5), c(1,2), c(1,3), c(1,4), c(1,5), c(2,3), c(2,4), c(2,5), c(3,4), c(3,5), c(4,5), c(1,2,3), c(1,2,4), c(1,2,5), c(1,3,4), c(1,3,5), c(1,4,5), c(2,3,4), c(2,3,5), c(2,4,5), c(3,4,5), c(1,2,3,4), c(1,2,3,5), c(1,2,4,5), c(1,3,4,5), c(2,3,4,5), c(1,2,3,4,5))

IIRtableCol2=vector()
for(i in 1:length(IIRtableCol1)){
  #partialMI/totalMI
  IIRtableCol2[i]=MIwPlus(p_x_data1, v_x_data1, v_xTerms_data1, IIRtableCol1[[i]])/totalMI1 #input features must be independent, 2nd & 3nd para. must have same length
}
#MIRtable END---------------


CASMIpct1=vector()
MIMpkGpct1=vector()
JMIpkGpct1=vector()
CMIMpct1=vector()
MRMRpct1=vector()
DISRpct1=vector()
NJMIMpct1=vector()
autoFeatureNum1=vector()

CASMItime1=vector()
MIMpkGtime1=vector()
JMIpkGtime1=vector()
CMIMtime1=vector()
MRMRtime1=vector()
DISRtime1=vector()
NJMIMtime1=vector()

for(i_loop in 1:loopNum){
  
  #generate dataset
  x1=sample(v_x1,n,replace = T,prob = p_x1)
  x2=sample(v_x2,n,replace = T,prob = p_x2)
  x3=sample(v_x3,n,replace = T,prob = p_x3)
  x4=sample(v_x4,n,replace = T,prob = p_x4)
  x5=sample(v_x5,n,replace = T,prob = p_x5)
  x6=x4; #redundent feature
  x7=rpois(n,2)-2; x7[x7>=2]=2
  x8=runif(n,-1,1); x8[x8<=(-0.6)]=-2; x8[(x8>(-0.6)) & (x8<(-0.2))]=-1.2; x8[(x8>=(-0.2)) & (x8<=0.2)]=0; x8[x8>=0.6]=2.5; x8[(x8>0.2) & (x8<0.6)]=1.5
  x9=rbinom(n,6,0.2)-6*0.2
  x10=rnorm(n); x10[x10<(-1.5)]=-2; x10[(x10>=(-1.5)) & (x10<(-0.7))]=-1.5; x10[(x10>=(-0.7)) & (x10<=0.7)]=0; x10[x10>1.5]=2; x10[(x10>0.7) & (x10<=1.5)]=1.5;
  error=runif(n,0,1); error[error<=(1/3)]=-1; error[(error>(1/3)) & (error<(2/3))]=0; error[error>=(2/3)]=1;
  
  #---------------------Experiment 1---------------------
  y1=x1+x2+x3^3-0.5*x4^2+abs(x5)+x6+error #x6=x4, others are independent
  data1=list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,y1)

  #feature selection
  CASMIstart1=proc.time()[[3]]
  CASMIfeature1=CASMI.featureSelection(data1)
  CASMIend1=proc.time()[[3]]
  CASMItime1[i_loop]=CASMIend1-CASMIstart1
  autoFeatureNum1[i_loop]=length(CASMIfeature1)

  if(autoFeatureNum1[i_loop]!=0){
    
    #for the Praznik package
    X=data.frame(matrix(unlist(data1[-length(data1)]), nrow=n, byrow=F))
    Y=as.factor(data1[[length(data1)]])
    
    MIMpkGstart1=proc.time()[[3]]
    MIMpkGfeature1=as.vector(MIM(X, Y, autoFeatureNum1[i_loop])$selection)
    MIMpkGend1=proc.time()[[3]]
    MIMpkGtime1[i_loop]=MIMpkGend1-MIMpkGstart1
    
    JMIpkGstart1=proc.time()[[3]]
    JMIpkGfeature1=as.vector(JMI(X, Y, autoFeatureNum1[i_loop])$selection)
    JMIpkGend1=proc.time()[[3]]
    JMIpkGtime1[i_loop]=JMIpkGend1-JMIpkGstart1
    
    CMIMstart1=proc.time()[[3]]
    CMIMfeature1=as.vector(CMIM(X, Y, autoFeatureNum1[i_loop])$selection)
    CMIMend1=proc.time()[[3]]
    CMIMtime1[i_loop]=CMIMend1-CMIMstart1
    
    MRMRstart1=proc.time()[[3]]
    MRMRfeature1=as.vector(MRMR(X, Y, autoFeatureNum1[i_loop])$selection)
    MRMRend1=proc.time()[[3]]
    MRMRtime1[i_loop]=MRMRend1-MRMRstart1
    
    DISRstart1=proc.time()[[3]]
    DISRfeature1=as.vector(DISR(X, Y, autoFeatureNum1[i_loop])$selection)
    DISRend1=proc.time()[[3]]
    DISRtime1[i_loop]=DISRend1-DISRstart1
    
    NJMIMstart1=proc.time()[[3]]
    NJMIMfeature1=as.vector(NJMIM(X, Y, autoFeatureNum1[i_loop])$selection)
    NJMIMend1=proc.time()[[3]]
    NJMIMtime1[i_loop]=NJMIMend1-NJMIMstart1
    
    #evaluation by feature recovered information using mutual information
    #process x4 x6 redundency
    CASMI1=processFeatures(CASMIfeature1)
    if(length(CASMI1)!=0){
      #Information Recovery Ratio using the IIRtable
      for(index in 1:length(IIRtableCol1)){
        if(all(sort(CASMI1)==IIRtableCol1[[index]])){
          CASMIpct1[i_loop]=IIRtableCol2[[index]]
          break
        }
      }
    }else{
      CASMIpct1[i_loop]=0
    }
    
    #process x4 x6 redundency
    MIMpkG1=processFeatures(MIMpkGfeature1)
    if(length(MIMpkG1)!=0){
      #Information Recovery Ratio using the IIRtable
      for(index in 1:length(IIRtableCol1)){
        if(all(sort(MIMpkG1)==IIRtableCol1[[index]])){
          MIMpkGpct1[i_loop]=IIRtableCol2[[index]]
          break
        }
      }
    }else{
      MIMpkGpct1[i_loop]=0
    }
    
    #process x4 x6 redundency
    JMIpkG1=processFeatures(JMIpkGfeature1)
    if(length(JMIpkG1)!=0){
      #Information Recovery Ratio using the IIRtable
      for(index in 1:length(IIRtableCol1)){
        if(all(sort(JMIpkG1)==IIRtableCol1[[index]])){
          JMIpkGpct1[i_loop]=IIRtableCol2[[index]]
          break
        }
      }
    }else{
      JMIpkGpct1[i_loop]=0
    }
    
    #process x4 x6 redundency
    CMIM1=processFeatures(CMIMfeature1)
    if(length(CMIM1)!=0){
      #Information Recovery Ratio using the IIRtable
      for(index in 1:length(IIRtableCol1)){
        if(all(sort(CMIM1)==IIRtableCol1[[index]])){
          CMIMpct1[i_loop]=IIRtableCol2[[index]]
          break
        }
      }
    }else{
      CMIMpct1[i_loop]=0
    }
    
    #process x4 x6 redundency
    MRMR1=processFeatures(MRMRfeature1)
    if(length(MRMR1)!=0){
      #Information Recovery Ratio using the IIRtable
      for(index in 1:length(IIRtableCol1)){
        if(all(sort(MRMR1)==IIRtableCol1[[index]])){
          MRMRpct1[i_loop]=IIRtableCol2[[index]]
          break
        }
      }
    }else{
      MRMRpct1[i_loop]=0
    }
    
    #process x4 x6 redundency
    DISR1=processFeatures(DISRfeature1)
    if(length(DISR1)!=0){
      #Information Recovery Ratio using the IIRtable
      for(index in 1:length(IIRtableCol1)){
        if(all(sort(DISR1)==IIRtableCol1[[index]])){
          DISRpct1[i_loop]=IIRtableCol2[[index]]
          break
        }
      }
    }else{
      DISRpct1[i_loop]=0
    }
    
    #process x4 x6 redundency
    NJMIM1=processFeatures(NJMIMfeature1)
    if(length(NJMIM1)!=0){
      #Information Recovery Ratio using the IIRtable
      for(index in 1:length(IIRtableCol1)){
        if(all(sort(NJMIM1)==IIRtableCol1[[index]])){
          NJMIMpct1[i_loop]=IIRtableCol2[[index]]
          break
        }
      }
    }else{
      NJMIMpct1[i_loop]=0
    }
    
  }else{
    #percentage
    CASMIpct1[i_loop]=0
    MIMpkGpct1[i_loop]=0
    JMIpkGpct1[i_loop]=0
    CMIMpct1[i_loop]=0
    MRMRpct1[i_loop]=0
    DISRpct1[i_loop]=0
    NJMIMpct1[i_loop]=0
    
    CASMItime1[i_loop]=NULL
    MIMpkGtime1[i_loop]=NULL
    JMIpkGtime1[i_loop]=NULL
    CMIMtime1[i_loop]=NULL
    MRMRtime1[i_loop]=NULL
    DISRtime1[i_loop]=NULL
    NJMIMtime1[i_loop]=NULL
  }

}

#results
print(n)
print(loopNum)

print(mean(autoFeatureNum1))
print(mean(CASMIpct1))
print(mean(MIMpkGpct1))
print(mean(JMIpkGpct1))
print(mean(CMIMpct1))
print(mean(MRMRpct1))
print(mean(DISRpct1))
print(mean(NJMIMpct1))


#95% confidence interval
if(loopNum%%200==0){
  lowbound=sort(CASMIpct1)[(loopNum*0.025+1)]; upbound=sort(CASMIpct1)[(loopNum*0.975)]
  CASMIci1=paste0('[',lowbound,', ',upbound,']')
  
  lowbound=sort(MIMpkGpct1)[(loopNum*0.025+1)]; upbound=sort(MIMpkGpct1)[(loopNum*0.975)]
  MIMpkGci1=paste0('[',lowbound,', ',upbound,']')
  
  lowbound=sort(JMIpkGpct1)[(loopNum*0.025+1)]; upbound=sort(JMIpkGpct1)[(loopNum*0.975)]
  JMIpkGci1=paste0('[',lowbound,', ',upbound,']')
  
  lowbound=sort(CMIMpct1)[(loopNum*0.025+1)]; upbound=sort(CMIMpct1)[(loopNum*0.975)]
  CMIMci1=paste0('[',lowbound,', ',upbound,']')
  
  lowbound=sort(MRMRpct1)[(loopNum*0.025+1)]; upbound=sort(MRMRpct1)[(loopNum*0.975)]
  MRMRci1=paste0('[',lowbound,', ',upbound,']')
  
  lowbound=sort(DISRpct1)[(loopNum*0.025+1)]; upbound=sort(DISRpct1)[(loopNum*0.975)]
  DISRci1=paste0('[',lowbound,', ',upbound,']')
  
  lowbound=sort(NJMIMpct1)[(loopNum*0.025+1)]; upbound=sort(NJMIMpct1)[(loopNum*0.975)]
  NJMIMci1=paste0('[',lowbound,', ',upbound,']')
  
  print(CASMIci1)
  print(MIMpkGci1)
  print(JMIpkGci1)
  print(CMIMci1)
  print(MRMRci1)
  print(DISRci1)
  print(NJMIMci1)
}

# print running time
time_end=proc.time()[[3]]
runningTime=time_end-time_start
print(runningTime)





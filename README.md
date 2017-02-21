# gene-project high dimension data analysis

setwd('~/desktop/')
load('gene.rdata')

now we have three dataset  pheno[include the id,age,sex,sbp]  map[position of those snp]  and geno[snp data ]
#data explore:
dim(geno)  # 1000  4106  high dimension data  hard to deal with
pheno1=pheno

sum(is.na(pheno1)) # check whether include NA value.
attach(pheno)
qplot(id,sbp_v1)   # we get two outliers which sbp_v1= 0 and we all know the blood pressure could not be zero, so it should be missing value.we need to deal with the NA value at first.


pheno1[which(pheno1$stroke=='NA'),]=1
pheno1[which(pheno1$sbp_v1=='0'),]=mean(pheno1$sbp_v1,na.rm=T)
sex=ifelse(pheno1$gender=='F',0,1)

#check the mean of different group, healthy and disease
t.test(stroke0$sbp_v1,stroke1$sbp_v1)

''p-value = 0.2538.so remain the null hypithesis : difference is equal to 0 .so we could say the sbp do not impact the disease.Drop sbp feature''

t.test(stroke0$age_v1,stroke1$age_v1)

''p-value = 0.002035 .so reject the null hypithesis : difference is equal to 0.so we could say the age do impact the disease.''
mean(x)=54.33120  
mean (y)=57.04688
So keep age feature


male=pheno1[which(pheno1[,2]=='M'),]
female=pheno1[which(pheno1[,2]=='F'),]
t.test(female$stroke,male$stroke)

p-value = 0.09728.so remain the null hypithesis : difference is equal to 0 . so we could say the gender do not impact the disease.But the sex is not the t distribution, so not t.test

fm1=glm(pheno1$stroke~pheno1$age_v1+pheno1$sbp_v1+sex,family = 'binomial')

fm2=glm(pheno1$stroke~pheno1$age_v1+sex,family = 'binomial')
summary(fm)

aa=geno

#anova to fliter the features
pvalue=as.numeric()
for(i in 1:4106){
  pvalue[i]=anova(glm(pheno1$stroke~pheno1$age_v1+sex+aa[,i],
                      family = 'binomial'),test = 'Chisq')[5,5]
}
plot(map$base_pair_position,y=-log10(pvalue),type='h')

which(-log10(pvalue)>3)

685 1576 1876 1979  have siginificant effect to control the stroke.

pvalue=as.numeric()
for(i in 1:4106){
  pvalue[i]=anova(glm(pheno1$stroke~pheno1$age_v1+sex+new[,i],
                      family = 'binomial'),test = 'Chisq')[4,5]
}
plot(map$base_pair_position,y=-log10(pvalue),type='h')

which(-log10(pvalue)>2)

176  189  206  241  272  457  599  623  685  864  964 1013 1129 1159 1346 1576 1577 1631 1710 1787 1962 1979 2030 2209 2235 2258 2326 2429 2448 2547 2613 2780 2887 3209 3290 3330

#try to complete the geno with none missing value then use the lasso method and randomforest to do variable selection

new=impute.knn(x)
x=as.matrix(data.frame(x,genocomplete$data)
fm3=cv.glmnet(x,y=as.factor(pheno1$stroke),family = 'binomial')            
fm4=glmnet(x,y=as.factor(pheno1$stroke),family = 'binomial',
          alpha = 1,lambda = fm3$lambda.1se)
aa=as.matrix(coef(fm4))
id=which(aa!=0)
order(aa[id])[1:10]

176 189 206 215 241 258 272 284 358 363
##the result is similar with the anova feature selection (GWAS part) and some conflict may due to our complete method knn.

#random forest

rf=randomForest(x=gene1,y=pheno1$stroke,ntree = 5000,importance = T)
imp=importance(rf)
imp1=imp[,ncol(imp)-1]]
rf.gene=names(imp)[order(imp,decreasing=T)[1:30]]
rf.genes1
 [1]  206 2398  844 4091 3225  136 3370  801 3925 3958 3896 3573 2613   53 1250 3591 4008 3794
[19] 3472  241 2531 4061 2184 2014  253 2905 2671 2918  988 3849

rf.gene1=order(imp,decreasing=T)[1:30]

[1] "SNP_A-2116945_C" "SNP_A-2097288_C" "SNP_A-1893339_A" "SNP_A-8379052_A" "SNP_A-8681544_T"
 [6] "SNP_A-2271733_C" "SNP_A-2066767_A" "SNP_A-8581747_A" "SNP_A-4256643_A" "SNP_A-8282222_C"
[11] "SNP_A-2117385_G" "SNP_A-1944211_A" "SNP_A-4206226_T" "SNP_A-2305747_G" "SNP_A-8634574_G"
[16] "SNP_A-1907102_T" "SNP_A-8707038_G" "SNP_A-8696522_C" "SNP_A-2240872_A" "SNP_A-1796659_G"
[21] "SNP_A-8289795_T" "SNP_A-2031914_G" "SNP_A-1795319_G" "SNP_A-1993437_T" "SNP_A-8476241_G"
[26] "SNP_A-8548519_A" "SNP_A-2013882_A" "SNP_A-4274431_G" "SNP_A-8374102_C" "SNP_A-8579143_A"

Because we just have the snp geno data which is 0,1,2. If we have the geno project as continue number, we can use method like this.
use the loocv method to get the best parameter of those model:

vars=list()
vars$randomForest=list(ntree=c(500,750,1000),mtry=c(5,10,20),
                       fs.meth=list(list('all'),list('rf',30),
                                    list('varclus',30,50)))
vars$knn=list(k=c(3,5,7,11),norm=c(T,F),
              fs.meth=list(list('all'),
              list('rf',30),
              list('varclus',30,50)))

all=join(knn,randomforest,by='variants')
rankSystems(all,top=10,max=T) # then we could get the important gene feature.






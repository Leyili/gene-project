# gene-project
#high dimension

#For the whole project, we first do some data clean and then we build our model.

#in data clean:
#1.  we found stroke 85 == NA   we can delete or complete this case
#2.  and we also found pheno$sbp_v1==0   77 and 260  which is not reliable 
#so we need to delete or complete this two data

#3.also we find snp data set have so many NA values.

#First of all, we need to deal with the Q 1,2 
#For logistics regression model, it can deal with those NA values. But we can also compare the result with NA or not.

#For the Q3, we choose knn method to complete those value.

#In logistics, we do the GWAS add each snp into the model.

#When we get the result, I think the regular method is not so good. So I try to use the random Forest to do the classify regression. To check which snp may be highly related with this disease.


load('~/desktop/FinalExam2.rdata')

pheno1=pheno
sex=ifelse(pheno1$gender=='F',0,1)

t.test(stroke0$sbp_v1,stroke1$sbp_v1)
# p-value = 0.2538
# so remain the null hypithesis : difference is equal to 0 .
# so we could say the sbp do not impact the disease.

Drop sbp feature

t.test(stroke0$age_v1,stroke1$age_v1)
#p-value = 0.002035
# so reject the null hypithesis : difference is equal to 0 .
# so we could say the age do impact the disease.
#mean of x   mean of y 
# 54.33120   57.04688
#Keep age feature


male=pheno1[which(pheno1[,2]=='M'),]
female=pheno1[which(pheno1[,2]=='F'),]
t.test(female$stroke,male$stroke)
#p-value = 0.09728
#so remain the null hypithesis : difference is equal to 0 . 
#so we could say the gender do not impact the disease.
#But the sex is not the t distribution, so not t.test

fm=glm(pheno1$stroke~pheno1$age_v1+pheno1$sbp_v1+sex,family = 'binomial')

fm=glm(pheno1$stroke~pheno1$age_v1+sex,family = 'binomial')
summary(fm)
#Coefficients:
               Estimate Std. Error z value Pr(>|z|)    
(Intercept)   -3.206072   0.322127  -9.953   <2e-16 ***
pheno1$age_v1  0.005371   0.004425   1.214   0.2248    
sex            0.444913   0.261146   1.704    0.0884 .  

aa=geno[-c(77,216),]

pvalue=as.numeric()
for(i in 1:4106){
  pvalue[i]=anova(glm(pheno1$stroke~pheno1$age_v1+sex+aa[,i],
                      family = 'binomial'),test = 'Chisq')[5,5]
}
plot(map$base_pair_position,y=-log10(pvalue),type='h')

which(-log10(pvalue)>3)
#685 1576 1876 1979  have siginificant effect to control the stroke.



pvalue=as.numeric()
for(i in 1:4106){
  pvalue[i]=anova(glm(pheno1$stroke~pheno1$age_v1+sex+new[,i],
                      family = 'binomial'),test = 'Chisq')[4,5]
}
plot(map$base_pair_position,y=-log10(pvalue),type='h')

which(-log10(pvalue)>2)
 [1]  176  189  206  241  272  457  599  623  685  864  964 1013 1129 1159 1346 1576 1577 1631
[19] 1710 1787 1962 1979 2030 2209 2235 2258 2326 2429 2448 2547 2613 2780 2887 3209 3290 3330

rf=randomForest(x=new,y=pheno1$stroke,ntree = 5000,importance = T)
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

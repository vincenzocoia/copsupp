# check of code on hjcondvine.R

library(CopulaModel)
source("hjcondvine.R")
d=5
np=matrix(1,d,d)
ntrunc=4
udat=matrix(c(.1,.2,.3,.4,.5, .6,.7,.8,.9,.5),2,5,byrow=T)
parvec=c(2,2,2,2,1.5,1.5,1.5,1.2,1.2,1.1)

# D-vine
D5=Dvinearray(d)
pcondmat=matrix("pcondgum",d,d)
qcondmat=matrix("qcondgum",d,d)

pmatd=rvinepcond(parvec,udat,D5,ntrunc,pcondmat,np)
print(pmatd)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.7385935 0.8312337 0.8782010
#[2,]  0.6 0.7328919 0.8830648 0.9754268 0.1403602

qmatd=rvineqcond(pmatd,D5,ntrunc,parvec,np,qcondmat,pcondmat,iinv=F)
print(qmatd) # same as udat

q5=rvineqcond(pmatd,D5,ntrunc,parvec,np,qcondmat,pcondmat,iinv=T)
print(q5)
# [1] 0.5 0.5
pnewd=pmatd; pnewd[,d]=0.9
q5b=rvineqcond(pnewd,D5,ntrunc,parvec,np,qcondmat,pcondmat,iinv=T)
print(q5b)
# [1] 0.5227173 0.8994989

# truncated vine
pmatd2=rvinepcond(parvec,udat,D5,2,pcondmat,np)
print(pmatd2)
#      [,1]      [,2]      [,3]      [,4]       [,5]
#[1,]  0.1 0.4938008 0.7385935 0.7430672 0.76377507
#[2,]  0.6 0.7328919 0.8830648 0.9508236 0.08757126
qmatd2=rvineqcond(pmatd2,D5,2,parvec,np,qcondmat,pcondmat,iinv=F)
print(qmatd2) # same as udat
q52=rvineqcond(pmatd2,D5,2,parvec,np,qcondmat,pcondmat,iinv=T)
print(q52)
# [1] 0.5 0.5
pnewd=pmatd2; pnewd[,d]=0.9
q52b=rvineqcond(pnewd,D5,2,parvec,np,qcondmat,pcondmat,iinv=T)
print(q52b)
# [1] 0.6190225 0.9281953

# ============================================================

# C-vine
C5=Cvinearray(d)
pmatc=rvinepcond(parvec,udat,C5,ntrunc,pcondmat,np)
print(pmatc)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.7228406 0.8433683 0.9167192
#[2,]  0.6 0.7328919 0.8700505 0.9776872 0.1262401

qmatc=rvineqcond(pmatc,C5,ntrunc,parvec,np,qcondmat,pcondmat,iinv=F)
print(qmatc) # same as udat

q5=rvineqcond(pmatc,C5,ntrunc,parvec,np,qcondmat,pcondmat,iinv=T)
print(q5)
# [1] 0.5 0.5
pnewc=pmatc; pnewc[,d]=0.9
q5b=rvineqcond(pnewc,C5,ntrunc,parvec,np,qcondmat,pcondmat,iinv=T)
print(q5b)
# [1] 0.4800962 0.8916635

# different truncation levels
pmatc3=rvinepcond(parvec,udat,C5,3,pcondmat,np)
print(pmatc3)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.7228406 0.8433683 0.9319560
#[2,]  0.6 0.7328919 0.8700505 0.9776872 0.1922732

qmatc3=rvineqcond(pmatc3,C5,3,parvec,np,qcondmat,pcondmat,iinv=F)
print(qmatc3) # same as udat
q53=rvineqcond(pmatc3,C5,3,parvec,np,qcondmat,pcondmat,iinv=T)
print(q53)
# [1] 0.5 0.5
pnewc=pmatc3; pnewc[,d]=0.9
q53b=rvineqcond(pnewc,C5,3,parvec,np,qcondmat,pcondmat,iinv=T)
print(q53b)
# [1] 0.4562103 0.8437685
 
# different quantiles of case 1

newpmatc=matrix(pmatc[1,],9,d,byrow=T)
newpmatc[,d]=seq(.1,.9,.1)
print(newpmatc)
#      [,1]      [,2]      [,3]      [,4] [,5]
# [1,]  0.1 0.4938008 0.7228406 0.8433683  0.1
# [2,]  0.1 0.4938008 0.7228406 0.8433683  0.2
# [3,]  0.1 0.4938008 0.7228406 0.8433683  0.3
# [4,]  0.1 0.4938008 0.7228406 0.8433683  0.4
# [5,]  0.1 0.4938008 0.7228406 0.8433683  0.5
# [6,]  0.1 0.4938008 0.7228406 0.8433683  0.6
# [7,]  0.1 0.4938008 0.7228406 0.8433683  0.7
# [8,]  0.1 0.4938008 0.7228406 0.8433683  0.8
# [9,]  0.1 0.4938008 0.7228406 0.8433683  0.9

newqmatc=rvineqcond(newpmatc,C5,ntrunc,parvec,np,qcondmat,pcondmat,iinv=T)
print(newqmatc) 
# [1] 0.06239024 0.10819359 0.15071992 0.19263156 0.23561800 0.28141551 0.33255194
# [8] 0.39398474 0.48009621

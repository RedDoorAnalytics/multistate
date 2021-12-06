
library(mstate)
data(ebmt3)
tmat=trans.illdeath()
tglong <- msprep(time=c(NA,"prtime","rfstime"),status=c(NA,"prstat","rfsstat"),
                 data=ebmt3,trans=tmat)

cx <- coxph(Surv(Tstart/365.24,Tstop/365.24,status)~strata(trans),
            data=tglong,method="breslow")
summary(cx)
# new data, to check whether results are the same for transition 1 as
# those in appendix E.1 of Therneau & Grambsch (2000)
newdata <- data.frame(trans=1:3,strata=1:3)
HvH <- msfit(cx,newdata,trans=tmat)
# probtrans
pt <- probtrans(HvH,predt=0)
plot(pt)
# predictions from state 1
pt[[1]]
library(Epi)
library(mstate)
data(ebmt3)
bmt <- Lexis(exit = list(tft = rfstime/365.25), exit.status = factor(rfsstat, labels = c("Tx", "RD")), data = ebmt3)
bmtr <- cutLexis(bmt, cut = bmt$prtime/365.25, precursor.states = "Tx", new.state = "PR")
summary(bmtr)
boxes(bmtr, boxpos = TRUE)
bmtx <- splitLexis(bmtr, time.scale = "tft", breaks = c(0:19/10, seq(2, 10, 1)))

i.kn <- c(0,0.2, 2,4, 6,10)
b.kn <- c(0, 7)

library( splines )
DM.Ins <- glm( (lex.Xst=="PR") ~ Ns( tft, knots=c(0.1,7)) ,
               family=poisson, offset=log(lex.dur),
               data = subset(bmtx,lex.Cst=="Tx") )
DM.Dead <- glm( (lex.Xst=="RD") ~ Ns( tft, knots=i.kn) ,
                family=poisson, offset=log(lex.dur),
                data = subset(bmtx,lex.Cst=="Tx") )
Ins.Dead <- glm( (lex.Xst=="RD") ~ Ns( tft, knots=i.kn) ,
                 family=poisson, offset=log(lex.dur),
                 data = subset(bmtx,lex.Cst=="PR") )

Tr <- list( "Tx" = list( "PR"       = DM.Ins,
                         "RD"      = DM.Dead  ),
            "PR" = list( "RD" = Ins.Dead ) )
lapply( Tr, names )

ini <- subset(bmtx,select=1:14)[NULL,]
ini[1:1,"lex.Cst"] <- "Tx"
ini[1:1,"tft"] <- 0
str(ini)

system.time(simL <- simLexis( Tr, ini, time.pts=seq(0,7.7,0.385), N=100000 ))
#summary( simL )
system.time(nSt <- nState( simL, at=seq(0,7.7,0.385), from=0,time.scale="tft"))
system.time(pp <- pState( nSt))
plot(pp)

#217s

library(flexsurv)

tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
tmat=trans.illdeath()
tglong <- msprep(time=c(NA,"prtime","rfstime"),status=c(NA,"prstat","rfsstat"),
                 data=ebmt3,trans=tmat)
bexp <- flexsurvreg(Surv(Tstart,Tstop, status) ~ trans, data=tglong, dist="weibull")
bspl <- flexsurvspline(Surv(Tstart,Tstop, status) ~ trans ,
                       data=tglong, k=0)

bexp.list <- vector(3, mode="list")
for (i in 1:3) { 
  bexp.list[[i]] <- flexsurvspline(Surv(Tstart,Tstop, status) ~ 1, subset=(trans==i),
                                data=tglong, k=4)
}

pmatrix.simfs(bexp.list, t=s[i], trans=tmat,M=100000)

bspl <- flexsurvspline(Surv(Tstart,Tstop, status) ~ trans,
                       data=bosms3, k=4)
bexp <- flexsurvreg(Surv(Tstart,Tstop, status) ~ trans, data=tglong, dist="weibull")
test <- function() {
	s <- seq(0,7.7,0.385)
	for(i in length(s)) {
p = 	pmatrix.simfs(bspl, t=s[i], trans=tmat,M=1000)
}
}
system.time(test())

pmatrix.simfs(bexp, t=5, trans=tmat,M=100000)
pmatrix.simfs(bexp, t=10, trans=tmat)




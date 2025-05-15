library(cpt)
library(changepoint.np)
library(mosum)
library(ecp)
library(wbs)
library(cumSeg)
#mix数据 eta检测 显著性水平为0.1和0.5
td1 <- testData(model = "mix",seed = 1234)
ts <- td1$x
rankts <- rank(td1$x)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(4, 1))
plot(ts(ts),col = "black",main="原序列",xlab="",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,41,61,101,121,161,201,251,301,361,421,491),lwd = 3.5, col = "blue",lty=3)
#G = 10 alpha = 0.1 and 0.5
m1 <- mosum(rankts,G = 10,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m1$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 10 alpha = 0.1 and 0.5",ylab = "")
abline(v = m1$cpts, col = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.1), col = 4)
m1 <- mosum(rankts,G = 10,alpha =0.5 ,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m1$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.5), col = 4,lty = 2)
#G = 25 alpha = 0.1 and 0.5
m2 <- mosum(rankts,G = 25,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m2$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 25 alpha = 0.1 and 0.5",ylab = "")
abline(v = m2$cpts, col = 2)
abline(h = mosum.criticalValue(n, 25, 25, 0.1), col = 4)
m2 <- mosum(rankts,G = 25,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m2$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 25, 25, 0.5), col = 4,lty = 2)
#G = 50 alpha = 0.1 and 0.5
m3 <- mosum(rankts,G = 50,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m3$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 50 alpha = 0.1 and 0.5",ylab = "")
abline(v = m3$cpts, col = 2)
abline(h = mosum.criticalValue(n, 50, 50, 0.1), col = 4)
m3 <- mosum(rankts,G = 50,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m3$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 50, 50, 0.5), col = 4,lty = 2)
#G = 60 alpha = 0.1 and 0.5
m4 <- mosum(rankts,G = 60,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m4$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 60 alpha = 0.1 and 0.5",ylab = "")
abline(v = m4$cpts, col = 2)
abline(h = mosum.criticalValue(n, 60, 60, 0.1), col = 4)
m4 <- mosum(rankts,G = 60,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m4$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 60, 60, 0.5), col = 4,lty = 2)
#teeth10数据 eta检测 显著性水平为0.1和0.5
td1 <- testData(model = "teeth10",seed = 1234)
ts <- td1$x
rankts <- rank(td1$x)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(4, 1))
plot(ts(ts),col = "black",main="原序列",xlab="",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,31,41,51,61,71,81,91,101,111,121,131),lwd = 3.5, col = "blue",lty=3)
#G = 8 alpha = 0.1 and 0.5
m1 <- mosum(rankts,G = 8,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m1$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 10 alpha = 0.1 and 0.5",ylab = "")
abline(v = m1$cpts, col = 2)
abline(h = mosum.criticalValue(n, 8, 8, 0.1), col = 4)
m1 <- mosum(rankts,G = 8,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m1$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 8, 8, 0.5), col = 4,lty = 2)
#G = 10 alpha = 0.1 and 0.5
m2 <- mosum(rankts,G = 10,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m2$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",,xlab = "G = 25 alpha = 0.1 and 0.5",ylab = "")
abline(v = m2$cpts, col = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.1), col = 4)
m2 <- mosum(rankts,G = 10,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m2$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.5), col = 4,lty = 2)
#G = 25 alpha = 0.1 and 0.5
m3 <- mosum(rankts,G = 20,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m3$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 50 alpha = 0.1 and 0.5",ylab = "")
abline(v = m3$cpts, col = 2)
abline(h = mosum.criticalValue(n, 20, 20, 0.1), col = 4)
m3 <- mosum(rankts,G = 20,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m3$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 20, 20, 0.5), col = 4,lty = 2)
#Stairs10数据 eta检测 显著性水平为0.1和0.5
td1 <- testData(model = "stairs10",seed = 1234)
ts <- td1$x
rankts <- rank(td1$x)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(4, 1))
plot(ts(ts),col = "black",main="原序列",xlab="",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141),lwd = 3.5, col = "blue",lty=3)
#G = 8 alpha = 0.1 and 0.5
m1 <- mosum(rankts,G = 8,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m1$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 10 alpha = 0.1 and 0.5",ylab = "")
abline(v = m1$cpts, col = 2)
abline(h = mosum.criticalValue(n, 8, 8, 0.1), col = 4)
m1 <- mosum(rankts,G = 8,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m1$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 8, 8, 0.5), col = 4,lty = 2)
#G = 10 alpha = 0.1 and 0.5
m2 <- mosum(rankts,G = 10,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m2$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 25 alpha = 0.1 and 0.5",ylab = "")
abline(v = m2$cpts, col = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.1), col = 4)
m2 <- mosum(rankts,G = 10,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m2$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.5), col = 4,lty = 2)
#G = 20 alpha = 0.1 and 0.5
m3 <- mosum(rankts,G = 20,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m3$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 50 alpha = 0.1 and 0.5",ylab = "")
abline(v = m3$cpts, col = 2)
abline(h = mosum.criticalValue(n, 20, 20, 0.1), col = 4)
m3 <- mosum(rankts,G = 20,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m3$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 20, 20, 0.5), col = 4,lty = 2)

#第二个样本
#mix数据 eta检测 显著性水平为0.1和0.5
td1 <- testData(model = "mix",seed = 2345)
ts <- td1$x
rankts <- rank(td1$x)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(4, 1))
plot(ts(ts),col = "black",main="原序列",xlab="",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,41,61,101,121,161,201,251,301,361,421,491),lwd = 3.5, col = "blue",lty=3)
#G = 10 alpha = 0.1 and 0.5
m1 <- mosum(rankts,G = 10,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m1$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 10 alpha = 0.1 and 0.5",ylab = "")
abline(v = m1$cpts, col = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.1), col = 4)
m1 <- mosum(rankts,G = 10,alpha =0.5 ,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m1$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.5), col = 4,lty = 2)
#G = 25 alpha = 0.1 and 0.5
m2 <- mosum(rankts,G = 25,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m2$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 25 alpha = 0.1 and 0.5",ylab = "")
abline(v = m2$cpts, col = 2)
abline(h = mosum.criticalValue(n, 25, 25, 0.1), col = 4)
m2 <- mosum(rankts,G = 25,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m2$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 25, 25, 0.5), col = 4,lty = 2)
#G = 50 alpha = 0.1 and 0.5
m3 <- mosum(rankts,G = 50,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m3$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 50 alpha = 0.1 and 0.5",ylab = "")
abline(v = m3$cpts, col = 2)
abline(h = mosum.criticalValue(n, 50, 50, 0.1), col = 4)
m3 <- mosum(rankts,G = 50,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m3$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 50, 50, 0.5), col = 4,lty = 2)
#G = 60 alpha = 0.1 and 0.5
m4 <- mosum(rankts,G = 60,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m4$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 60 alpha = 0.1 and 0.5",ylab = "")
abline(v = m4$cpts, col = 2)
abline(h = mosum.criticalValue(n, 60, 60, 0.1), col = 4)
m4 <- mosum(rankts,G = 60,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m4$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 60, 60, 0.5), col = 4,lty = 2)
#teeth10数据 eta检测 显著性水平为0.1和0.5
td1 <- testData(model = "teeth10",seed = 2345)
ts <- td1$x
rankts <- rank(td1$x)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(4, 1))
plot(ts(ts),col = "black",main="原序列",xlab="",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,31,41,51,61,71,81,91,101,111,121,131),lwd = 3.5, col = "blue",lty=3)
#G = 8 alpha = 0.1 and 0.5
m1 <- mosum(rankts,G = 8,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m1$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 10 alpha = 0.1 and 0.5",ylab = "")
abline(v = m1$cpts, col = 2)
abline(h = mosum.criticalValue(n, 8, 8, 0.1), col = 4)
m1 <- mosum(rankts,G = 8,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m1$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 8, 8, 0.5), col = 4,lty = 2)
#G = 10 alpha = 0.1 and 0.5
m2 <- mosum(rankts,G = 10,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m2$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",,xlab = "G = 25 alpha = 0.1 and 0.5",ylab = "")
abline(v = m2$cpts, col = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.1), col = 4)
m2 <- mosum(rankts,G = 10,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m2$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.5), col = 4,lty = 2)
#G = 25 alpha = 0.1 and 0.5
m3 <- mosum(rankts,G = 20,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m3$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 50 alpha = 0.1 and 0.5",ylab = "")
abline(v = m3$cpts, col = 2)
abline(h = mosum.criticalValue(n, 20, 20, 0.1), col = 4)
m3 <- mosum(rankts,G = 20,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m3$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 20, 20, 0.5), col = 4,lty = 2)
#Stairs10数据 eta检测 显著性水平为0.1和0.5
td1 <- testData(model = "stairs10",seed = 2345)
ts <- td1$x
rankts <- rank(td1$x)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(4, 1))
plot(ts(ts),col = "black",main="原序列",xlab="",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141),lwd = 3.5, col = "blue",lty=3)
#G = 8 alpha = 0.1 and 0.5
m1 <- mosum(rankts,G = 8,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m1$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 10 alpha = 0.1 and 0.5",ylab = "")
abline(v = m1$cpts, col = 2)
abline(h = mosum.criticalValue(n, 8, 8, 0.1), col = 4)
m1 <- mosum(rankts,G = 8,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m1$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 8, 8, 0.5), col = 4,lty = 2)
#G = 10 alpha = 0.1 and 0.5
m2 <- mosum(rankts,G = 10,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m2$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 25 alpha = 0.1 and 0.5",ylab = "")
abline(v = m2$cpts, col = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.1), col = 4)
m2 <- mosum(rankts,G = 10,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m2$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 10, 10, 0.5), col = 4,lty = 2)
#G = 20 alpha = 0.1 and 0.5
m3 <- mosum(rankts,G = 20,alpha = 0.1,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
plot(m3$stat,type = "l",lwd = 2,lty = 1,main = "Rmosum",xlab = "G = 50 alpha = 0.1 and 0.5",ylab = "")
abline(v = m3$cpts, col = 2)
abline(h = mosum.criticalValue(n, 20, 20, 0.1), col = 4)
m3 <- mosum(rankts,G = 20,alpha = 0.5,criterion = "eta", eta = 0.15,var.est.method = 'mosum')
abline(v = m3$cpts, col = 2,lty = 2)
abline(h = mosum.criticalValue(n, 20, 20, 0.5), col = 4,lty = 2)

#自下而上合并窗宽
options (warn = - 1)
td1 <- testData(model = "mix",seed = 1234)
ts <- td1$x
rankts <- rank(td1$x)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(2, 1))
plot(ts(ts),col = "black",main="原序列",xlab="",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,41,61,101,121,161,201,251,301,361,421,491),lwd = 3.5, col = "blue",lty=3)
mbu1 <- multiscale.bottomUp(rankts,G = c(10, 25, 50, 60),eta = 2/3,var.est.method = 'mosum')
plot(mbu1,type = "l",lty = 2,main = "mosum",xlab = "G = 10 alpha = 0.1 and 0.5",ylab = "")
plot(mbu1,display = "data", shaded = 'none') 
plot(mbu1, display = "significance", shaded = "CI")
summary(mbu1)

threshold.custom <- function(G, n, alpha) {
  mosum.criticalValue(n, G, G, alpha) * log(n/G)^0.1
}
mbu2 <- multiscale.bottomUp(rankts, G = c(10, 25, 50, 60),var.est.method = 'mosum',threshold = "custom",threshold.function = threshold.custom)
plot(mbu2,type = "l",lty = 2,main = "mosum",xlab = "G = 10 alpha = 0.1 and 0.5",ylab = "")
plot(mbu2,display = "data", shaded = 'none') 
plot(mbu2, display = "significance", shaded = "CI")
summary(mbu2)

mlp <- multiscale.localPrune(rankts, G = c(10, 25, 50, 60),var.est.method = 'mosum',criterion = "epsilon",epsilon =0.15)
plot(mlp,type = "l",lty = 2,main = "mosum",xlab = "G = 10 alpha = 0.1 and 0.5",ylab = "")
plot(mlp,display = "data", shaded = 'none') 
plot(mlp, display = "significance", shaded = "CI")
summary(mlp)



options (warn = - 1)
td1 <- testData(model = "mix")
ts <- td1$x
rankts <- rank(td1$x)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(2, 1))
plot(ts(ts),col = "black",xlab="Primary Sequence",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,41,61,101,121,161,201,251,301,361,421,491),lwd = 3.5, col = "blue",lty=3)
mbu1 <- multiscale.bottomUp(rankts,G = c(10, 25, 50, 60),eta = 2/3,var.est.method = 'mosum')
plot(mbu1,type = "l",lwd = 2,lty = 2,main = "mosum",xlab = "Rmosum",ylab = "")


plot(mbu1, display = "significance", shaded = "CI")
mlp <- multiscale.localPrune(rankts, G = c(10, 25, 50, 60),var.est.method = 'mosum',criterion = "epsilon",epsilon =0.15)
plot(mlp, display = "significance", shaded = "CI")

options (warn = - 1)
td1 <- testData(model = "teeth10")
ts <- td1$x
rankts <- rank(td1$x)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(2, 1))
plot(ts(ts),col = "black",xlab="Primary Sequence",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,31,41,51,61,71,81,91,101,111,121,131),lwd = 3.5, col = "blue",lty=3)
mbu1 <- multiscale.bottomUp(rankts,G = c(10, 25, 50, 60),eta = 2/3,var.est.method = 'mosum')
plot(mbu1,type = "l",lwd = 2,lty = 2,main = "mosum",xlab = "Rmosum",ylab = "")


plot(mbu1, display = "significance", shaded = "CI")
mlp <- multiscale.localPrune(rankts, G = c(10, 25, 50, 60),var.est.method = 'mosum',criterion = "epsilon",epsilon =0.1)
plot(mlp, display = "significance", shaded = "CI")




options (warn = - 1)
td1 <- testData(model = "stairs10")
ts <- td1$x
rankts <- rank(td1$x)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(2, 1))
plot(ts(ts),col = "black",xlab="Primary Sequence",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141),lwd = 3.5, col = "blue",lty=3)
mbu1 <- multiscale.bottomUp(rankts,G = c(10, 25, 50, 60),eta = 2/3,var.est.method = 'mosum')
plot(mbu1,type = "l",lwd = 2,lty = 2,main = "mosum",xlab = "Rmosum",ylab = "")


plot(mbu1, display = "significance", shaded = "CI")
mlp <- multiscale.localPrune(rankts, G = c(10, 25, 50, 60),var.est.method = 'mosum',criterion = "epsilon",epsilon =0.15)
plot(mlp, display = "significance", shaded = "CI")


td1 <- testData(lengths = c(10,10,20,20,30,30,40,40,50,50,60,60,70,70),means = c(7,-7,6,-6,5,-5,4,-4,3,-3,2,-2,1,-1),sds = rep(0,14))
set.seed(1234)
ts = td1$x + 2*rt(560,2)
rankts = rank(ts)
n <- length(signal)
signal <- td1$mu
par(mfrow = c(2, 1))
plot(ts(ts),col = "black",xlab="Primary Sequence",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,41,61,101,121,161,201,251,301,361,421,491),lwd = 3.5, col = "blue",lty=3)
mbu1 <- multiscale.bottomUp(rankts,G = c(10, 25, 50, 60),eta = 2/3,var.est.method = 'mosum')
plot(mbu1,type = "l",lwd = 2,lty = 2,main = "mosum",xlab = "Rmosum",ylab = "")


plot(mbu1, display = "significance", shaded = "CI")
mlp <- multiscale.localPrune(rankts, G = c(10, 25, 50, 60),var.est.method = 'mosum',criterion = "epsilon",epsilon =0.1)
plot(mlp, display = "significance", shaded = "CI")

td1 <- testData(lengths = c(10,10,20,20,30,30,40,40,50,50,60,60,70,70),means = c(7,-7,6,-6,5,-5,4,-4,3,-3,2,-2,1,-1),sds = rep(0,14))
set.seed(1234)
ts = td1$x + 2*rt(560,3)
rankts = rank(ts)
n <- length(signal)
signal <- td1$mu
par(mfrow = c(2, 1))
plot(ts(ts),col = "black",xlab="Primary Sequence",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,41,61,101,121,161,201,251,301,361,421,491),lwd = 3.5, col = "blue",lty=3)
mbu1 <- multiscale.bottomUp(rankts,G = c(10, 25, 50, 60),eta = 2/3,var.est.method = 'mosum')
plot(mbu1,type = "l",lwd = 2,lty = 2,main = "mosum",xlab = "G = 10 alpha = 0.1 and 0.5",ylab = "")


plot(mbu1, display = "significance", shaded = "CI")
mlp <- multiscale.localPrune(rankts, G = c(10, 25, 50, 60),var.est.method = 'mosum',criterion = "epsilon",epsilon =0.1)
plot(mlp, display = "significance", shaded = "CI")

td1 <- testData(lengths = rep(10,14),means = c(0,1,0,1,0,1,0,1,0,1,0,1,0,1),sds = rep(0,14))
set.seed(1234)
ts = td1$x + 0.2*rt(140,2)
rankts = rank(ts)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(2, 1))
plot(ts(ts),col = "black",xlab="Primary Sequence",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,31,41,51,61,71,81,91,101,111,121,131),lwd = 3.5, col = "blue",lty=3)
mbu1 <- multiscale.bottomUp(rankts,G = c(10, 25, 50, 60),eta = 2/3,var.est.method = 'mosum')
plot(mbu1,type = "l",lwd = 2,lty = 2,main = "mosum",xlab = "Rmosum",ylab = "")


plot(mbu1, display = "significance", shaded = "CI")
mlp <- multiscale.localPrune(rankts, G = c(10, 25, 50, 60),var.est.method = 'mosum',criterion = "epsilon",epsilon =0.1)
plot(mlp, display = "significance", shaded = "CI")

td1 <- testData(lengths = rep(10,14),means = c(0,1,0,1,0,1,0,1,0,1,0,1,0,1),sds = rep(0,14))
set.seed(1234)
ts = td1$x + 0.2*rt(140,3)
rankts = rank(ts)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(2, 1))
plot(ts(ts),col = "black",xlab="Primary Sequence",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,31,41,51,61,71,81,91,101,111,121,131),lwd = 3.5, col = "blue",lty=3)
mbu1 <- multiscale.bottomUp(rankts,G = c(10, 25, 50, 60),eta = 2/3,var.est.method = 'mosum')
plot(mbu1,type = "l",lwd = 2,lty = 2,main = "mosum",xlab = "Rmosum",ylab = "")


plot(mbu1, display = "significance", shaded = "CI")
mlp <- multiscale.localPrune(rankts, G = c(10, 25, 50, 60),var.est.method = 'mosum',criterion = "epsilon",epsilon =0.1)
plot(mlp, display = "significance", shaded = "CI")

td1 <- testData(lengths = rep(10,15),means = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),sds = rep(0,15))
ts = td1$x + 0.1*rt(150,2)
rankts = rank(ts)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(2, 1))
plot(ts(ts),col = "black",xlab="Primary Sequence",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141),lwd = 3.5, col = "blue",lty=3)
mbu1 <- multiscale.bottomUp(rankts,G = c(10, 25, 50, 60),eta = 2/3,var.est.method = 'mosum')
plot(mbu1,type = "l",lwd = 2,lty = 2,main = "mosum",xlab = "Rmosum",ylab = "")


plot(mbu1, display = "significance", shaded = "CI")
mlp <- multiscale.localPrune(rankts, G = c(10, 25, 50, 60),var.est.method = 'mosum',criterion = "epsilon",epsilon =0.15)
plot(mlp, display = "significance", shaded = "CI")

td1 <- testData(lengths = rep(10,15),means = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),sds = rep(0,15))
ts = td1$x + 0.1*rt(150,3)
rankts = rank(ts)
signal <- td1$mu
n <- length(signal)
par(mfrow = c(2, 1))
plot(ts(ts),col = "black",xlab="Primary Sequence",ylab="")
lines(signal, col = 2, lwd = 2, lty = 2)
abline(v = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141),lwd = 3.5, col = "blue",lty=3)
mbu1 <- multiscale.bottomUp(rankts,G = c(10, 25, 50, 60),eta = 2/3,var.est.method = 'mosum')
plot(mbu1,type = "l",lwd = 2,lty = 2,main = "mosum",xlab = "Rmosum",ylab = "")


plot(mbu1, display = "significance", shaded = "CI")
mlp <- multiscale.localPrune(rankts, G = c(10, 25, 50, 60),var.est.method = 'mosum',criterion = "epsilon",epsilon =0.1)
plot(mlp, display = "significance", shaded = "CI")


# 20.1.2021, Matti Pirinen (matti.pirinen@helsinki.fi), Heidi Hautakangas
# Create FIGURE 3 in Hautakangas et al. 2021 
# "Genome-wide analysis of 102,084 migraine cases identifies 123 risk loci and subtype-specific risk alleles"

#Includes:
#1. Bayesian model comparison between subtype-specific models and shared effect model and null model.
#   Each model has prior prob of 0.25.
#   Shared effect model is a 50%:50% mixture of fixed effect model and independent effect model.
#
#2. Plot -log10 P-values of each subtype and color according to model if posterior > 95%. (Panel A)
#
#3. Plot subtype effects and color if effects are different after Bonferroni correction. (Panel B)
#

# Explanation of the Bayesian approach: 
#  https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS9.html

# Read in functions:
log.dmvnorm <- function(x, mu = rep(0, length(x)), S = diag(1, length(x)) ){
  #returns log of density of MV-Normal(mean = mu, var = S) at x 
  K = length(mu)
  stopifnot(all(dim(S) == K))
  stopifnot(length(x) == K)
  chol.S = chol(S) #Cholesky decomposition
  log.det = 2*sum(log(diag(chol.S))) #log of det(S)
  inv.chol.S = solve(t(chol.S)) #inverse of cholesky^T
  return(-K/2*log(2*pi) - 0.5*(log.det + crossprod(inv.chol.S %*% as.numeric(x-mu))))
}

abf.mv <- function(b.est, Sigmas, prior = rep(1,length(Sigmas))){
  #Returns posterior probabilities of the models listed in Sigmas by their 
  #   total variance matrix (= sum of prior + likelihood variance matrices)
  #Returns also approximate Bayes factors (ABFs) w.r.t the first model in Sigmas.
  
  M = length(Sigmas) #number of models
  K = length(b.est) #number of studies
  prior = prior/sum(prior)
  log.abf = sapply(Sigmas, function(x){log.dmvnorm(b.est, S = x)})
  abf = exp(log.abf - log.abf[1]) #abf w.r.t the first model
  posterior = prior*abf
  posterior = posterior/sum(posterior)
  return(list(posterior = posterior, abf = abf))
}

#Input data has effect estimates for 123 migraine risk variants.
# Estimates are given for 3 phenotypes: migraine, 
# migraine with aura (MA) and migraine without aura (MO).
x <- read.table("MIG_MA_MO_123loci.txt", as.is = T, header = T)

r = 0.148 #MA-MO correlation for the likelihood; was computed empirically from GWAS data (see paper)
tau2 = 0.2^2 #prior variance of effect size -- same for all non-null models
S.null = matrix(0, 2, 2) #null model has 0 effects
S.fix = tau2 * matrix(1, 2, 2) #fixed-effect model has 1s in the correlation matrix
S.ind = tau2 * diag(2) #diagonal matrix for independent effects model, off-diagonals = 0
S.ma = tau2 * matrix(c(1, 0, 0, 0), 2, 2) #Effect only in MA but not in MO
S.mo = tau2 * matrix(c(0, 0, 0, 1), 2, 2) #Effect only in MO but not in MA
priors = c(1, 0.5, 0.5, 1, 1) #NULL, FIX, IND, MA, MO
#Note that model BOTH will be a union of FIX and IND and hence they 
# have only half the prior of other models --> BOTH has equal prior with other models

res = c()
for(ii in 1:nrow(x)){
  #Sigma is the variance of likelihood function 
  Sigma = diag(x[ii,c("SE_MA","SE_MO")]) %*% matrix(c(1,r,r,1),2,2) %*% diag(x[ii,c("SE_MA","SE_MO")])
  Var.matrices = list(Sigma + S.null, Sigma + S.fix, Sigma + S.ind, Sigma + S.ma, Sigma + S.mo)
  res = rbind(res,abf.mv(b.est = x[ii,c("Beta_MA","Beta_MO")], Sigmas = Var.matrices, prior = priors)$posterior)
}
colnames(res) = c("NULL","FIXED","IND","MA","MO")

res2 = cbind(res[,1], res[,2] + res[,3], res[,4], res[,5])
colnames(res2) = c("NULL","BOTH","MA","MO")


#test for same effects
v = x$SE_MA^2 + x$SE_MO^2 - 2 * x$SE_MA * x$SE_MO * r #variance of difference (b.ma - b.mo)
pval.same = pchisq( (x$Beta_MA - x$Beta_MO)^2 / v, df = 1, lower = F) #testing whether effect are equal
x$P_same <- pval.same


#For plotting, remove variants that have MAF < 5%
#because these have high SEs in subtype analyses and spoil the effect size plot 
maf.thr = 0.05
ii.keep = (x$EAF_mig > maf.thr & x$EAF_mig < (1 - maf.thr) )
y = x[ii.keep,]
pr = res2[ii.keep,] #pr = (posterior) probability
n = nrow(y)

#Make Figure 3

#png(paste0("Fig3.png"), width = 800, height = 400)
pdf(paste0("Fig3.pdf"), width = 10, height = 5)

par(mfrow = c(1,2))
par(mar = c(4.5,4.5,1,1))

cols = rep("gray90",n)
pr.thr = 0.95
ii.ma = (pr[,"MA"] > pr.thr)
ii.mo = (pr[,"MO"] > pr.thr)
ii.both = (pr[,"BOTH"] > pr.thr)

cols[ii.ma] = "red"
cols[ii.mo] = "blue"
cols[ii.both] = "orange"

texts = rep(NA, n)
texts[ii.ma | ii.mo | ii.both] = gsub("_"," ",y[ ii.ma | ii.mo | ii.both, "Locus_ID"])

insy = rep(NA, n)
insy[ ii.ma | ii.mo | ii.both ] <- 0.35
insy[1] <- -0.4 #PRDM16
insx = rep(NA, n)
insx[ ii.ma | ii.mo | ii.both ] <- 0.35
insx[63] <- 1.2 # PLCE1
insx[c(107)] <- 2.3
insx[c(94)] <- 2.0
insx[c(73)] <- 2.1
insy[c(107,94,73)] <- 0.1

# 3a
plot(NULL, xlim = c(0,30), ylim = c(0,15), xlab = "-log10 P-value MO", ylab = "-log10 P-value MA")
abline(h = -log10(5e-8), lwd = 0.5, lty = 2)
abline(v = -log10(5e-8), lwd = 0.5, lty = 2)

ii = (ii.ma | ii.mo)
points(-log10(y$P_MO), -log10(y$P_MA), col = cols, pch = 19, cex = 1)
text(-log10(y$P_MO)+insx, -log10(y$P_MA)+insy, labels = texts, cex = 0.6)

legend("topright", legend = c(paste0(c("Pr(MA) > ", "Pr(MO) > ", "Pr(BOTH) > "),pr.thr), "OTHER"),
       col = c("red","blue","orange","gray90"), pch = 19, cex = 1.1)

# 3b

ii = (y$P_same < 0.05/123) #bonferroni for # lead SNPs
cols = rep("indianred1", n) 
cols[abs(y$Beta_MO) > abs(y$Beta_MA)] = "dodgerblue"
cols[!ii] = "gray90"

texts = rep(NA, n)
texts[ii] = gsub("_"," ",y[ii, "Locus_ID"])

insy = rep(NA, n)
insy[ii] = .008
insy[31] = -.005 #near SPINK2
insy[c(107,94,73,106)] = .005 
insx = rep(NA, n)
insx[ii] = .005
insx[c(107,94,73)] = .015 # CACNA1A, HMOX2, MPPED2

ran = range(y$Beta_MA, y$Beta_MO, 0.2)
plot(NULL, xlim = ran, ylim = ran, xlab = "logOR MO", ylab = "logOR MA")
abline(0,1, lwd = 0.5, lty = 2)
abline(h = 0, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)

arrows(y$Beta_MO[ii], y$Beta_MA[ii]-1.96*y$SE_MA[ii], 
       y$Beta_MO[ii], y$Beta_MA[ii]+1.96*y$SE_MA[ii], 
       col = "indianred1", lwd = 0.5, code = 3, length = 0)

arrows(y$Beta_MO[ii]-1.96*y$SE_MO[ii], y$Beta_MA[ii], 
       y$Beta_MO[ii]+1.96*y$SE_MO[ii], y$Beta_MA[ii], 
       col = "dodgerblue", lwd = 0.5, code = 3, length = 0)

points(y$Beta_MO, y$Beta_MA, col = cols, pch = 18, cex = 1.3)

text(y$Beta_MO + insx, y$Beta_MA + insy, labels = texts, cex = 0.6)

legend("topright", legend = c("MA > MO", "MO > MA", "OTHER"),
       col = c("indianred1","dodgerblue","gray90"), pch = 18, cex = 1.1, bg="white")

dev.off()


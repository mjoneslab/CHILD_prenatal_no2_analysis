# run the below in the order as it appears

######################################
# irwsva required internal functions #
######################################

# edge.lfdr function is internal to irwsva function below
# mono function is internal to edge.lfdr
# they have to be initialized before it will run

mono <- function(lfdr){
  .Call("monotone", as.numeric(lfdr), PACKAGE="sva")
}

edge.lfdr <- function(p, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, lambda=0.8, ...) {
  pi0 <- mean(p >= lambda)/(1 - lambda)
  pi0 <- min(pi0, 1)
  n <- length(p)
  transf <- match.arg(transf)
  if(transf=="probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1-eps)
    x <- qnorm(p)
    myd <- density(x, adjust=adj)
    mys <- smooth.spline(x=myd$x, y=myd$y)
    y <- predict(mys, x)$y
    lfdr <- pi0*dnorm(x)/y
  }
  if(transf=="logit") {
    x <- log((p+eps)/(1-p+eps))
    myd <- density(x, adjust=adj)
    mys <- smooth.spline(x=myd$x, y=myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x) / (1+exp(x))^2
    lfdr <- pi0 * dx/y
  }
  if(trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if(monotone) {	
    lfdr <- lfdr[order(p)]
    lfdr <- mono(lfdr)
    lfdr <- lfdr[rank(p)]
  }
  return(lfdr)
}


##################################
# numsv function internal to sva #
##################################

#determine the number of req SVs
#this is numsv function code pulled out

#this is what is fed into this function from the sva function
dat = betas
mod = mod
#method="be" # this is the default and what is called by numsv
#vfilter=NULL # we dont apply a filter and its not used for be method
B=20
# seed=1234 # this is for a separate part of the function
# also we set seed before runnign sva

warn <- NULL
n <- ncol(dat)
m <- nrow(dat)
H <- mod %*% solve(t(mod) %*% mod) %*% t(mod)
res <- dat - t(H %*% t(dat))
uu <- svd(res)
ndf <- min(m, n) - ceiling(sum(diag(H)))
dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
dstat0 <- matrix(0, nrow = B, ncol = ndf)
for (i in 1:B) {
  res0 <- t(apply(res, 1, sample, replace = FALSE))
  res0 <- res0 - t(H %*% t(res0))
  uu0 <- svd(res0)
  dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
  print(i)
}
psv <- rep(1, n)
for (i in 1:ndf) {
  psv[i] <- mean(dstat0[, i] >= dstat[i])
}
for (i in 2:ndf) {
  psv[i] <- max(psv[(i - 1)], psv[i])
}

# this is the variable holding the number of SVs required 
nsv <- sum(psv <= 0.1)

  
####################################
# code from the main SVA function  #
####################################
  
dat =  betas
mod = mod
mod0 = mod0 
n.sv = nsv # from code above (run internally when sva called)
controls = NULL # we dont include controls - only for supervised method
method = "irw" # this isnt actually used but would be called by sva
vfilter = NULL # we dont set this when running sva
B = 5 
numSVmethod = "be" # this is the default called by sva (code above)
  

###################################
# irwsva function internal to sva #
###################################

# this code is really want does all the work in sva
# we use the irw method so this is the code for that method
# there are also supervised and two step methods (not used by us)

# irw function takes input from sva
# these were also defined above
dat = dat 
mod = mod
mod0 = mod0
n.sv = n.sv
B = B

# this is the rest of the code from the irwsva function 
n <- ncol(dat)
m <- nrow(dat)

Id <- diag(n)
resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% 
                    t(mod))
uu <- eigen(t(resid) %*% resid)
vv <- uu$vectors
ndf <- n - dim(mod)[2]
pprob <- rep(1, m)
one <- rep(1, n)
Id <- diag(n)
df1 <- dim(mod)[2] + n.sv
df0 <- dim(mod0)[2] + n.sv
rm(resid)
cat(paste("Iteration (out of", B, "):"))

for (i in 1:B) {
  mod.b <- cbind(mod, uu$vectors[, 1:n.sv])
  mod0.b <- cbind(mod0, uu$vectors[, 1:n.sv])
  ptmp <- f.pvalue(dat, mod.b, mod0.b)
  pprob.b <- (1 - edge.lfdr(ptmp))
  mod.gam <- cbind(mod0, uu$vectors[, 1:n.sv])
  mod0.gam <- cbind(mod0)
  ptmp <- f.pvalue(dat, mod.gam, mod0.gam)
  pprob.gam <- (1 - edge.lfdr(ptmp))
  pprob <- pprob.gam * (1 - pprob.b)
  dats <- dat * pprob
  dats <- dats - rowMeans(dats)
  uu <- eigen(t(dats) %*% dats)
  cat(paste(i, " "))
}

# this is what is normally returned by sva
# problem is that it doesnt include the diagonal matrix we want
sv = svd(dats)$v[, 1:n.sv, drop = FALSE] 
svobj <- list(sv = sv, pprob.gam = pprob.gam, pprob.b = pprob.b, 
               n.sv = n.sv)

# to get all info (including diagonal)
svobj_all <- svd(dats)

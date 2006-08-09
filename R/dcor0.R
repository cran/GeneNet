### dcor0.R  (2004-03-15)
###
###    Distribution of the Correlation Coefficient (rho=0) 
###    and Related Functions
###    
###
### Copyright 2003-04 Juliane Schaefer and Korbinian Strimmer
###
### This file is part of the `GeneNet' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA



# exact likelihood for rho=0
dcor0 <- function (x, kappa, log = FALSE)
{
    r <- x
    
    log.f <- log(kappa - 1) - kappa * log(2) +
             (kappa - 3) * log(1 - r * r)/2 +
             lgamma(kappa) - 2 * lgamma((1 + kappa)/2)
    
    for (i in which(abs(r) == 1) ) # special treatment for abs(r)=1
    {
      if (kappa > 3)  log.f[i] <- -Inf
      if (kappa == 3) log.f[i] <- -log(2)      
      if (kappa < 3)  log.f[i] <- Inf
    }
       
    if (log)
    {
        return(log.f)
    }
    else
    {
        return(exp(log.f))
    }
}


# random number generator
rcor0 <- function(n, kappa)
{
  #If n is a vector, length(n) values will be generated
  if (length(n) > 1)
    num.values <- length(n)
  else
    num.values <- n
  
  cvec <- rep(NA, num.values)
  for (i in 1:num.values)
  {
     x1 <- rnorm(kappa+1)
     x2 <- rnorm(kappa+1)
     cvec[i] <- cor(x1,x2)
  }  
    
  return(cvec)
}


# incomplete beta function
ibeta <- function(z, a, b)
{
  pbeta(z, a, b)*beta(a,b)
}


# distribution function for rho=0
pcor0 <- function(q, kappa, lower.tail=TRUE, log.p=FALSE)
{
  if (kappa == 1)
  {
    value <- rep(0.5, length(q))  
    value[q==-1.0] <- 0.0
    value[q== 1.0] <- 1.0
  }
  else
  { 
    z <- q*q
    a <- 0.5
    b <- (-1 + kappa)/2
    gamma.ratio <- exp(lgamma(b+0.5)-lgamma(b))
  
    value <-  0.5*(1 + (q*ibeta(z, a, b)*gamma.ratio)/
                       (sqrt(pi) * sqrt(z)))
  
    value[q==0] <- 0.5 # special treatment for q=0  
  }      
 
  if (lower.tail == FALSE)
     value <- 1-value
     
  if (log.p == TRUE)
    value <- log(value)
  
  return(value)
}


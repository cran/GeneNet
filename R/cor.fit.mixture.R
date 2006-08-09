### cor.fit.mixture.R  (2006-03-10)
###
###    Fit mixture model to empirical distribution of (partial)
###    correlation coefficients.
###    
###
### Copyright 2003-06 Juliane Schaefer and Korbinian Strimmer
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



# fit mixture to empirical (partial) correlations
cor.fit.mixture <- function(r, df=7, plot.locfdr=0)
{
  if ( any(r > 1) || any(r < -1) )
      stop("Data out of range: input correlations must be in [-1; 1]")
  
     require("locfdr")
     
     out <- locfdr(z.transform(r), df=df, plot=plot.locfdr)
     eta0 <- as.double( out$fp0[3] ) # p0
     sigma <- as.double( out$fp0[2] ) # sig
     kappa <- 1/(sigma*sigma) + 2 # Fisher's rule
     prob.nonzero <- 1-out$fdr
  

  return( list(kappa=kappa, eta0=eta0, prob.nonzero=prob.nonzero) )
}





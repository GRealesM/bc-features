##' this is a sandpit, to create a function that finds the slopes of tangents to
##' an ellipse from the origin
library("plotrix")
library(ellipse)
library(magrittr)
## mu = c(5,5)
## V=matrix(c(0.5,.1,.1,.3),2,2)
do_plot=TRUE

## plot(mu[1],mu[2], asp=1,xlim=c(-5,10),ylim=c(-5,10),pch="X")
## lines(ellipse(V,centre=mu))
## points(0,0)


find_tangent_slopes=function(mu, V, do_plot=TRUE) {
                                        # solve numerically, because the algebra gets ugly
  if(do_plot)
    par(mfrow=c(2,3))
  ## first note that slope from origin to ellipse tangents where ellipse centre is
  ## mu === slope from -mu to ellipse tangents where ellipse is centred on origin.
  ## So solve that simpler problem

## ## check V
##   if(V[1,1] < V[2,2]) {
##     diag(V)=rev(diag(V))
##     mu=rev(mu)
##     Vrev=TRUE
## } else {
##   Vrev = FALSE
##   }

  ## essential quantities
  P=-as.vector(mu)
  x0=mu[1]; y0=mu[2]
  ## these next two eqns from https://cookierobotics.com/007/
  evals=sum(diag(V))/2 + c(1,-1) * sqrt((V[1,1]-V[2,2])^2/4 + V[1,2]^2) # eigen values, ordered major, minor axes
  theta=atan2(evals[1]-V[1,1], V[1,2]) # angle of major axis to x axis
  RR=matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2) # rotate -theta
  R=matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2) # rotate theta
  scale_factor=sqrt(qchisq(0.95,2,lower=TRUE) * evals)

  if(do_plot) {
    E=ellipse(V)
    limits=c(min(c(as.vector(E),mu)),max(c(as.vector(E),mu)))*1.1
    Et=t(R %*% t(E)) # %*% diag(1/scale_factor)
    Ets=t(R %*% t(E)) %*% diag(1/(scale_factor))
    plot(0,0,pch="o",asp=1,xlim=limits,ylim=limits,main="forward transformations",sub="black -> red -> magenta")
    points(E,col="black")
    points(x0,y0,pch="x")
    points(Et,col="red")
    points(Ets,col="magenta")
  }

  P_rotated = R %*% P
  ## shrink by eigenvalues so that we have a unit circle
  P_shrunk = P_rotated / scale_factor
  ## now find tangents to the unit circle
  a=P_shrunk[1,1]; b=P_shrunk[2,1]
  yt1=(b-sqrt(a^4+a^2*b^2-a^2))/(a^2+b^2) # from algebra
  yt2=(b+sqrt(a^4+a^2*b^2-a^2))/(a^2+b^2)
  xt1=c(1,-1)*sqrt(1-yt1^2)
  xt2=c(1,-1)*sqrt(1-yt2^2)
  ## xt1, xt2 can be +/- ... pick that which satisfies slope condition (b-y)/(a-x)=-x/y
  slope=function(x,y) (b-y)/(a-x) + x/y
  w1=which.min(abs(slope(xt1,yt1))) # have to use < small rather than == 0 because of numerical errors? not sure, but close to zero results "look right"
  w2=which.min(abs(slope(xt2,yt2)))
  if(slope(xt1,yt1)[w1] > 1e-8 | slope(xt2,yt2)[w2] > 1e-8)
    stop("cannot find tangent points")
  T1=c(xt1[w1],yt1)
  T2=c(xt2[w2],yt2)
  limits=c(min(c(T1,T2,-1,1,a,b)),max(c(T1,T2,-1,1,a,b)))*1.1
  ## viz
  if(do_plot) {
    ## viz destination
    plot(0,0, asp = 1,ylim=limits,xlim=limits,pch="o", main="fully transformed")
    draw.circle(0, 0, 1, nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)
    points(a,b,pch="X")
    segments(a,b,T1[1],T1[2])
    segments(a,b,T2[1],T2[2])
    points(T1[1],T1[2],pch="+")
    points(T2[1],T2[2],pch="+")
  }

  ## transform back - first scaling
  T1_unshrunk=T1 * scale_factor
  T2_unshrunk=T2 * scale_factor
  a_unshrunk=a * scale_factor[1]
  b_unshrunk=b * scale_factor[2]
  ##viz
  if(do_plot) {
    limits=c(min(c(T1_unshrunk,T2_unshrunk,a_unshrunk,b_unshrunk)),
             max(c(T1_unshrunk,T2_unshrunk,a_unshrunk,b_unshrunk)))*1.1
    plot(0,0,asp=1,xlim=limits,ylim=limits,pch="o", main="transformed and unshrunk")
    lines(ellipse(0,scale=sqrt(evals)),col="red")
    points(a_unshrunk,b_unshrunk,pch="x")
    segments(a_unshrunk,b_unshrunk,T1_unshrunk[1],T1_unshrunk[2])
    segments(a_unshrunk,b_unshrunk,T2_unshrunk[1],T2_unshrunk[2])
  }

  ## transform back - second rotate
  RR= R * matrix(c(1,-1,-1,1),2,2)
  T1_rotated=RR %*% T1_unshrunk
  T2_rotated=RR %*% T2_unshrunk
  mu_rotated=RR %*% c(a_unshrunk, b_unshrunk)
  a_rotated=mu_rotated[1]
  b_rotated=mu_rotated[2]
  e_rotated=t(RR %*% t(ellipse(0,scale=sqrt(evals))))
  ## viz
  if(do_plot) {
    limits=c(min(c(T1_rotated,T2_rotated,a_rotated,b_rotated)),
             max(c(T1_rotated,T2_rotated,a_rotated,b_rotated)))*1.1
    plot(0,0,asp=1,xlim=limits,ylim=limits,pch="o",main="fully transformed and back again")
    lines(ellipse(V))
    lines(e_rotated,col="red")
    points(a_rotated,b_rotated,pch="x")
    segments(a_rotated,b_rotated,T1_rotated[1],T1_rotated[2])
    segments(a_rotated,b_rotated,T2_rotated[1],T2_rotated[2])
  }
  ## viz to check slope values
  s1=(T1_rotated[2]-b_rotated) / (T1_rotated[1] - a_rotated)
  s2=(T2_rotated[2]-b_rotated) / (T2_rotated[1] - a_rotated)
 if(do_plot) {
    limits=c(-.1,.1)
    ## min(c(T1_unshrunk,T2_unshrunk,a_unshrunk,b_unshrunk)),
    ##          max(c(T1_unshrunk,T2_unshrunk,a_unshrunk,b_unshrunk)))*1.1
    plot(0,0,asp=1,xlim=limits,ylim=limits,pch="o", main="original scale")
    lines(ellipse(V,centre=mu))
    abline(a=0, b=s1)
    abline(a=0, b=s2)
  }

  ## if(Vrev)
  ##  return(1/c(s1,s2))
  ## else
    return(c(s1,s2))
}

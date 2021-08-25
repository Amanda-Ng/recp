######### Required packages ########
require(clue, quietly=T)
require(energy, quietly = T)
require(randtoolbox, quietly = T)
require(pracma, quietly = T)
require(kernlab, quietly = T)
require(crossmatch, quietly = T)
require(HHG, quietly = T)
require(gTests, quietly = T)
require(ramify)

######### Data generation ##########
m=200
n=200
data1=cbind(rcauchy(m,0,1),rcauchy(m,0,1))
data2=cbind(rcauchy(n,0.5,1),rcauchy(n,0,1))

######### Computing Rank Energy Statistic #########
computestatistic=function(x,y,m=nrow(x),n=nrow(y),dim=ncol(x),gridch=torus(m+n,dim))
{
  comdata=rbind(x,y)
  distmat=matrix(0,nrow=m+n,ncol=m+n)
  for(i in 1:(m+n))
    distmat[i,]=apply((comdata[i,]-t(gridch)),2,Norm)^2
  assignmentFUN=solve_LSAP(distmat)
  assignmentSOL=cbind(seq_along(assignmentFUN),assignmentFUN)
  randenergySTAT=eqdist.etest(gridch[assignmentSOL[,2],],sizes = c(m,n), R=1)
  return(randenergySTAT$statistic)
            
####### scaled sample measure of divergence #########
EqDisttest=((m+n)/(m*n))(randenergySTAT)

######## estimate change point location ##########
changept=argmax(EqDisttest)}

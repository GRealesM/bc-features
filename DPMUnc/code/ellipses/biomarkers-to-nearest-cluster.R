# Once you're in ~/rds/rds-cew54-basis/03-Bases/cell_basis_v3_varimax/, you have the projections in Projections/ and the clustering data in DPMUnc/psm_results/.
d="~/rds/rds-cew54-basis/03-Bases/cell_basis_v3_varimax"
d1=file.path(d,"Projections")
d2=file.path(d,"DPMUnc/psm_results")
list.files(d1)
f1=file.path(d1,"Projection_cell_basis_v3_FinnGenLD_20221005-v1.tsv")
list.files(d2)
f2=file.path(d2,  "UKBB_PC58_IMD_psm_data.rds")
f3=file.path(d2,  "UKBB_PCgrac_IMD_psm_data.rds")

library(pheatmap)
library(randomFunctions)
library(data.table)
library(magrittr)
library(Matrix)
source("ellipse-tangents.R")
## data=fread(f1)
pc58=readRDS(f2)
pc58b=sub("IMD","BMK",f2) %>% readRDS()
pcgr=readRDS(f3)
pcgrb=sub("IMD","BMK",f3) %>% readRDS()

data=fread("~/rds/rds-cew54-basis/03-Bases/cell_basis_v3_varimax/Projections/PT_sig_202211_UKBB.tsv")
data[,fdr:=p.adjust(P),by=c("Collection","PC")]
data=data[ !(Trait %in% c("MS_IMSGC_24076602_1","ph471_PanUKBB_PanUKBBR1_1"))][ !grepl("Chen_32888493_6",Trait)]
data=data[!grepl("ph401_PanUKBB_PanUKBBR2_1",Trait)]

## what are the assigned classes for diseases?
dt1=data.table(Label=names(pc58$calls$cl),cl58=unname(pc58$calls$cl))
dt2=data.table(Label=names(pcgr$calls$cl),clgr=unname(pcgr$calls$cl))
dt=merge(dt1,dt2,by="Label",all=TRUE)
dt[is.na(cl58), cl58:=0]
dt[is.na(clgr), clgr:=0]
library(ggalluvial)
dtsum=dt[,.(y=.N),by=c("cl58","clgr")]
theme_set(theme_minimal())
ggplot(dtsum, aes(y=y, axis1=cl58,axis2=clgr)) + geom_alluvium(aes(fill=factor(cl58))) +
    geom_stratum(width = 1/5,  fill="grey", color = "white") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("PC58", "PCgran"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  scale_y_continuous("Count") +
  ggtitle("Relationship between disease clusterings")
################################################################################

## 2D

wanted=c(names(pc58$calls$cl), names(pc58b$calls$cl))
mat=dcast(data[PC %in% c("PC5","PC8") & Label %in% wanted , .(Label,PC,value=Delta)],#ifelse(fdr < 0.01, sign(z), 0))],
           Label ~ PC, value.var="value")
mat$cl=pc58$calls$cl[ mat$Label ]
avg=mat[!is.na(cl), .(PC5=mean(PC5),PC8=mean(PC8)),by="cl"]
##    cl         PC5        PC8
## 1:  3 -0.12102896 0.05396896
## 2:  1 -0.05428827 0.01880943
## 3:  2 -0.01460196 0.03082377

## make them unit vectors
avg[,sum:=sqrt(PC5^2+PC8^2)]
avg[,PC5:=(PC5/sum)][,PC8:=(PC8/sum)]
avg
##    cl        PC5       PC8       sum
## 1:  3 -0.9133114 0.4072618  7.546223
## 2:  1 -0.9448928 0.3273800 17.405101
## 3:  2 -0.4281158 0.9037239 29.319060

## distance
f=function(x,y,a,b) sqrt(a^2 + b^2 - (a*x+b*y)^2 )
d1=f(avg[cl==1]$PC5, avg[cl==1]$PC8, mat$PC5, mat$PC8)
d2=f(avg[cl==2]$PC5, avg[cl==2]$PC8, mat$PC5, mat$PC8)
d3=f(avg[cl==3]$PC5, avg[cl==3]$PC8, mat$PC5, mat$PC8)

mat[grep("EGPA",Label)]
##                                    Label        PC5        PC8 cl
## 1: ANCA negative EGPA (European) / Lyons -0.1291607 0.05819574  3
## 2:               EGPA (European) / Lyons -0.1273751 0.05738087  3

## avg$Label=avg$cl
## ref=mat[!is.na(cl)]
## mat=rbind(mat[is.na(cl)], avg)

## library(pheatmap)
## pheatmap(data.frame(row.names=mat$Label,PC5=mat$PC5,PC8=mat$PC8,d1=d1,d2=d2,d3=d3),
##          annotation_row=data.frame(row.names=mat$Label,cl=factor(mat$cl)))

mat[,inferred_group:=apply(cbind(d1,d2,d3),1,which.min)]
mat=merge(mat, avg[,.(inferred_group=cl, mean_PC5=PC5, mean_PC8=PC8)], by="inferred_group")
mat[ is.na(cl) & (sign(PC5)!=sign(mean_PC5) | sign(PC8)!=sign(mean_PC8)), inferred_group:=NA]
head(mat)
with(mat, table(cl,inferred_group))
mat[cl==2 & inferred_group==3]
mat[is.na(inferred_group)]
mat[!is.na(cl), inferred_group:=cl]
mat[inferred_group==3, inferred_group:=1]

paletteLength=50
myColor <- colorRampPalette(c("yellow", "white", "blue"))(paletteLength)
x=c(mat$PC5,mat$PC8)
myBreaks <- c(seq(min(x), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(x)/paletteLength, max(x), length.out=floor(paletteLength/2)))
## pheatmap(test, color=myColor, breaks=myBreaks)
mat2 = as.data.frame(mat)
rownames(mat2)=mat2$Label
mat2=mat2[!is.na(mat2$inferred_group),]
mat2=mat2[order(mat2$inferred_group),]
mat2$cl  %<>%  as.factor()
mat2$inferred_group  %<>%  as.factor()
pheatmap(mat2[,c("PC5","PC8")],
         annotation_row=mat2[,c("cl","inferred_group")],
         breaks=myBreaks, color=myColor, cluster_rows=FALSE)


pheatmap(data.frame(row.names=mat$Label,PC5=mat$PC5,PC8=mat$PC8),
         annotation_row=data.frame(row.names=mat$Label,cl=factor(mat$cl), infgr=factor(mat$inferred_group)))
mat[inferred_group==2]

options(width=200)
mat[,.(Label,inferred_group)][order(inferred_group),]

################################################################################

## 4D

## try without dcast to write neater code

wanted=c(names(pcgr$calls$cl), names(pcgrb$calls$cl))
mat=data[ Label %in% wanted & PC %in% c("PC5","PC8","PC2","PC11"), .(Label, PC, Delta)]
mat$cl=pcgr$calls$cl[ mat$Label ]
avg=mat[!is.na(cl), .(Delta=mean(Delta)),by=c("PC","cl")]
## make them unit vectors
avg[,sum:=sqrt(sum(Delta^2)),by="cl"]
avg[,Delta:=(Delta/sum)]
avg
##       PC cl       Delta        sum
##  1:  PC2  1  0.16280136 0.06283447
##  2:  PC5  1 -0.88309695 0.06283447
##  3:  PC8  1  0.43570228 0.06283447
##  4: PC11  1  0.06163608 0.06283447
##  5:  PC2  2 -0.40906160 0.04250321
##  6:  PC5  2 -0.42654758 0.04250321
##  7:  PC8  2  0.66466547 0.04250321
##  8: PC11  2  0.45710567 0.04250321
##  9:  PC2  3  0.78861894 0.05313605
## 10:  PC5  3 -0.29747238 0.05313605
## 11:  PC8  3  0.48778997 0.05313605
## 12: PC11  3 -0.22726922 0.0531360

## distance
f=function(u,v,w,x,a,b,c,d) sqrt(a^2 + b^2 + c^2 + d^2 - (a*u+b*v+c*w+d*x)^2 )
g=function(vec,point) sqrt( sum(point^2) - sum(vec*point)^2)

v=rnorm(4); p=rnorm(4)
all.equal(g(v,p), f(v[1],v[2],v[3],v[4], p[1],p[2],p[3],p[4]))

mat=mat[order(PC)]
avg=avg[order(PC)]
mat[,d1:=g(avg[cl==1]$Delta, Delta), by="Label"]
mat[,d2:=g(avg[cl==2]$Delta, Delta), by="Label"]
mat[,d3:=g(avg[cl==3]$Delta, Delta), by="Label"]
mat[,inferred_group:=apply(cbind(d1,d2,d3),1,which.min)]

mat=merge(mat, avg[,.(PC,inferred_group=cl,mean_Delta=Delta)],by=c("PC","inferred_group"))
mat[!(PC=="PC11" & inferred_group==1), flag:=any(sign(Delta) != sign(mean_Delta), na.rm=TRUE),by="Label"] # don't infer if wrong side of axis
mat[,flag:=any(flag==TRUE,na.rm=TRUE),by="Label"] # except cl=1, PC11 which is almost along axis
mat[grep("glutamyl trans",Label)]
mat[grep("Eos",Label)]
mat[is.na(cl) & flag==TRUE, inferred_group:=NA]
mat[!is.na(cl), inferred_group:=cl]

library(randomFunctions)
mat2=dcast(mat, cl + inferred_group + d1 + d2 + d3 + Label ~ PC, value.var="Delta") %>% as.data.frame()
rownames(mat2)=mat2$Label
mat2$cl  %<>%  as.factor()
mat2$inferred_group  %<>%  as.factor()
library(pheatmap)
paletteLength=50
myColor <- colorRampPalette(c("yellow", "white", "blue"))(paletteLength)
myBreaks <- c(seq(min(mat$Delta), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(mat$Delta)/paletteLength, max(mat$Delta), length.out=floor(paletteLength/2)))
## pheatmap(test, color=myColor, breaks=myBreaks)
mat2=mat2[!is.na(mat2$inferred_group),]
mat2=mat2[order(mat2$inferred_group),]
pheatmap(mat2[,c("PC2","PC5","PC8","PC11")],
         annotation_row=mat2[,c("cl","inferred_group")],breaks=myBreaks, color=myColor, cluster_rows=FALSE)

################################################################################

## one function to rule them all

library(rmeta)
library(randomFunctions)
library(pheatmap)
labels_wanted=c(names(pcgr$calls$cl), names(pcgrb$calls$cl))
pcs_wanted=c("PC5","PC8","PC2","PC11")
calls=pcgr$calls

labels_wanted=c(names(pc58$calls$cl), names(pc58b$calls$cl))
pcs_wanted=c("PC5","PC8")
calls=pc58$calls

align_biomarkers=function(pcs_wanted, labels_wanted, calls, show_vectors=FALSE, debug=FALSE) {
  mat=data[ Label %in% labels_wanted & PC %in% pcs_wanted, .(Label, PC, Delta, SE=sqrt(Var.Delta))]
  mat$cl=calls$cl[ mat$Label ]
  tmp=split(mat,mat[,c("cl","PC")])
  for(i in seq_along(tmp)) {
    d=meta.summaries(tmp[[i]]$Delta, tmp[[i]]$SE)
    tmp[[i]][,mean.Delta:=d$summary][,mean.se:=d$se.summary]
  }
  avg=lapply(tmp, head, 1) %>% rbindlist()
  avg[,Label:=NULL][,Delta:=NULL][,SE:=NULL]
  ## make them unit vectors
  avg[,sum:=sqrt(sum(mean.Delta^2)),by="cl"]
  avg[,mean.Delta:=(mean.Delta/sum)]
  avg[,mean.se:=(mean.se/sum)]
  ## avg_old=mat[!is.na(cl), .(Delta=mean(Delta)),by=c("PC","cl")]
  ## avg_old[,sum:=sqrt(sum(Delta^2)),by="cl"]
  ## avg_old[,Delta:=(Delta/sum)]
  ## distance
  g=function(vec,point) sqrt( sum(point^2) - sum(vec*point)^2)
  mat=mat[order(PC)]
  avg=avg[order(PC)]
  mat[,d1:=g(avg[cl==1]$mean.Delta, Delta), by="Label"]
  mat[,d2:=g(avg[cl==2]$mean.Delta, Delta), by="Label"]
  mat[,d3:=g(avg[cl==3]$mean.Delta, Delta), by="Label"]
  mat[,inferred_group:=apply(cbind(d1,d2,d3),1,which.min)]

  mat=merge(mat, avg[,.(PC,inferred_group=cl,mean.Delta,mean.se)],by=c("PC","inferred_group"))
  mat[!(PC=="PC11" & inferred_group==1), flag:=sign(Delta) != sign(mean.Delta),by="Label"] # don't infer if wrong side of axis
  mat[,flag:=any(flag==TRUE,na.rm=TRUE),by="Label"] # except cl=1, PC11 which is almost along axis
  mat[is.na(cl) & flag==TRUE, inferred_group:=NA] # remove inference if flagged
  mat[!is.na(cl), inferred_group:=cl] #keep original group if disease
  if(length(pcs_wanted)==2)
    mat[inferred_group==3, inferred_group:=1]

  ## pheatmap likes data.frames
  mat2=dcast(mat, cl + inferred_group + d1 + d2 + d3 + Label ~ PC, value.var="Delta") %>% as.data.frame()
  rownames(mat2)=mat2$Label
  mat2$cl  %<>%  as.factor()
  mat2$inferred_group  %<>%  as.factor()
  paletteLength=50
  myColor <- colorRampPalette(c("yellow", "white", "blue"))(paletteLength)
  myBreaks <- c(seq(min(mat$Delta), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(mat$Delta)/paletteLength, max(mat$Delta), length.out=floor(paletteLength/2)))
  ## pheatmap(test, color=myColor, breaks=myBreaks)
  mat2=mat2[!is.na(mat2$inferred_group),]
  mat2=mat2[order(mat2$inferred_group),]
  if(show_vectors) {
    rs=sqrt(rowSums(mat2[,pcs_wanted]^2))
    for(nm in pcs_wanted)
      mat2[[nm]]=mat2[[nm]]/rs
  }
  pheatmap(mat2[,pcs_wanted],
           annotation_row=mat2[,c("cl","inferred_group")],breaks=myBreaks, color=myColor, cluster_rows=FALSE)
  invisible(mat)
}

align_biomarkers(c("PC5","PC8"), labels_wanted=c(names(pc58$calls$cl), names(pc58b$calls$cl)), calls=pc58$calls)
align_biomarkers(c("PC5","PC8","PC2","PC11"), labels_wanted=c(names(pcgr$calls$cl), names(pcgrb$calls$cl)), calls=pcgr$calls)

################################################################################

## ellipses

library(mvmeta) # multi-variate meta analysis
(load("~/rds/rds-cew54-basis/03-Bases/cell_basis_v3_varimax/basis_building/cell-basis-sparse-3.0.RData"))
## need full var.proj

mpath <- "~/rds/rds-cew54-basis/03-Bases/cell_basis_v2/reduced_datasets"

## project them, and return var.proj for pcs_wanted
calc_var.proj= function (beta, seb, pids, pcs_wanted) {
    if (length(beta) != length(seb) || length(beta) != length(pids) ||
        !length(beta))
        stop("arguments must be equal length vectors > 0")
    if (!all(pids %in% SNP.manifest$pid))
        stop("all pids must be members of sparse basis (SNP.manifest$pid)")
    if (length(pids) < 0.95 * nrow(rot.pca))
        warning("more than 5% sparse basis snps missing")
    b <- beta * shrinkage[pids] - beta.centers[pids]
    proj <- b %*% rot.pca[pids, ]
    v <- seb * shrinkage[pids] * rot.pca[pids, ]
    var.proj <- t(v) %*% LD[pids, pids] %*% v
    var.proj[pcs_wanted, pcs_wanted]
}

get_var.proj=function(f, ...) {
  sm <- fread(file.path(mpath,f))
                                        # Some checks
  sm <- unique(sm)
  sm[sm == ""] <- NA # Some missing data might pass as empty string. This will fix that
  sm <- na.omit(sm, cols = c("pid", "BETA", "SE", "P"))
  dups <- sm$pid[duplicated(sm$pid)]

  if(length(dups) > 0){
    dupmessage= "This file has duplicated pids. I removed them prior to projection. You might want to check it."
    message(dupmessage)
    sm <- sm[!pid %in% dups] # Remove all duplicated instances, to be safe
  }
  calc_var.proj(beta=sm$BETA, seb=sm$SE, pids=sm$pid, ...)
}

## PREPARE INPUTS
INPUT58=list(pcs_wanted=c("PC5","PC8"), labels_wanted=c(names(pc58$calls$cl), names(pc58b$calls$cl)), calls=pc58$calls)
INPUTgr=list(pcs_wanted=c("PC5","PC8","PC2","PC11"), labels_wanted=c(names(pcgr$calls$cl), names(pcgrb$calls$cl)), calls=pcgr$calls)

files_to_project  <- list.files(mpath, pattern = ".tsv")
lookup=with(unique(data[,.(Label,Trait)]), structure(Label,names=Trait))

traits_wanted=unique(data$Trait)
  files_wanted=paste0(traits_wanted, "-ft.tsv")
VARS=lapply(files_wanted, get_var.proj, pcs_wanted=INPUTgr$pcs_wanted)
names(VARS)=lookup[ traits_wanted ]

INPUT58$VARS=lapply(VARS, function(v) v[ INPUT58$pcs_wanted, INPUT58$pcs_wanted ])
INPUTgr$VARS=VARS


calc_meta_v=function(ss, pcs_wanted) {
  tmp=dt2mat(ss[PC %in% pcs_wanted], Label ~ PC, value.var="Delta")
  d=mvmeta(tmp, lapply(VARS[rownames(tmp)], function(v) v[pcs_wanted, pcs_wanted]))
  conf_slopes=find_tangent_slopes(mu=d$coefficients, V=d$vcov)
  ## d[c("coefficients","vcov")]
  tmp=dcast(ss, Label ~ PC, value.var="Delta")
  sc=sqrt(sum(d$coefficients^2))
  for(j in seq_along(pcs_wanted)) {
    tmp[[ paste0(pcs_wanted[j],".mu") ]] = d$coefficients[[j]]/sc
    tmp[[ paste0(pcs_wanted[j],".se") ]] = sqrt(diag(d$vcov))[j]/sc
  }
  tmp[,lower_slope:=conf_slopes[1]]
  tmp[,upper_slope:=conf_slopes[2]]
  attr(tmp, "meta") <- d[c("coefficients","vcov")]
  tmp
}

## TODO
choose_vector=function(label, ss, AVG, pcs_wanted=INPUT$pcs_wanted,do_plot=FALSE, fdr_threshold=0.01) {
  ## test data
  sstmp=ss[[label]][ match(pcs_wanted, PC) ]
  if(!any(sstmp$fdr < fdr_threshold))
    return(rep(FALSE,length(AVG)))
  Vtest=as.matrix(VARS[[label]][pcs_wanted, pcs_wanted])
  mutest=sstmp$Delta
  Etest=ellipse( Vtest, centre=mutest)
  ## disease cluster vectors
  overlap=sapply(seq_along(AVG), function(j) {
    E=cbind(Etest, yl=AVG[[j]]$lower_slope[1] * Etest[,1], yu=AVG[[j]]$upper_slope[1] * Etest[,1])
    ## check for overlap
    !(all(E[,2] < E[,"yl"] & E[,2] < E[,"yu"]) || all(E[,2] > E[,"yl"] & E[,2] > E[,"yu"]))
  })
  angle_ok=sapply(seq_along(AVG), function(j) {
    cf=attr(AVG[[j]],"meta")$coefficients
    angle_disease=atan(cf[2]/cf[1])
    angle_test=atan(mutest[2]/mutest[1])
    abs(angle_disease - angle_test) < 2 * pi/3
  })
  if(do_plot) {
    Elist=lapply(AVG, function(A)
      ellipse(attr(A,"meta")$vcov, centre=attr(A,"meta")$coefficients))
    limits=c(min(c(Etest,unlist(Elist))), max(c(Etest,unlist(Elist))))*1.1
    plot(Etest, type="l", asp=1, xlim=limits, ylim=limits)
    points(0,0,pch="o")
    for(j in seq_along(Elist)) {
      lines(Elist[[j]], col=j+1)
      abline(a=0, b=AVG[[j]]$lower_slope[1], col=j+1)
      abline(a=0, b=AVG[[j]]$upper_slope[1], col=j+1)
     }
  }
  overlap & angle_ok
  }

  ## if(debug) {
  ##   par(mfrow=c(1,1))
  ##   plot(E,type="l",xlim=c(-.1,.1),asp=1);
  ##   abline(a=0,b=ss[[i]]$lower_slope[1]); abline(a=0,b=ss[[i]]$upper_slope[1]);
  ## for(j in seq_along(AVG)) {
  ##   points(ellipse(attr(AVG[[j]],"meta")$vcov, centre=attr(AVG[[j]],"meta")$coefficients),col=j+1)
  ##   abline(a=0,b=AVG[[j]]$lower_slope[1],col=j+1,lty="dashed"); abline(a=0,b=AVG[[j]]$upper_slope[1],col=j+1,lty="dashed");
  ##   abline(a=0,b=attr(AVG[[j]],"meta")$coefficients[2] / attr(AVG[[j]],"meta")$coefficients[1],col=j+1)
  ##          }


## version that works for PC58 but not 4 PCs
align_biomarkers=function(INPUT,debug=FALSE) {
  mat=data[ Label %in% INPUT$labels_wanted & PC %in% INPUT$pcs_wanted, .(Label, PC, Delta, SE=sqrt(Var.Delta))]
  mat$cl=INPUT$calls$cl[ mat$Label ]
  if(length(INPUT$pcs_wanted)==2)
    mat[cl==3,cl:=1]
  ss=split(mat,mat[,c("cl")])
  if(length(INPUT$pcs_wanted==2)) {
    AVG=lapply(ss, calc_meta_v, INPUT$pcs_wanted)
  } else {
    AVG1=lapply(ss, calc_meta_v, INPUT$pcs_wanted[1:2])
    AVG2=lapply(ss, calc_meta_v, INPUT$pcs_wanted[3:4])
  }

  ## this bit for finding nearest by distance, may drop
  avg=lapply(AVG, head, 1) %>% rbindlist()
  avg[,Label:=NULL]
  ## distance
  g=function(vec,point) sqrt( sum(point^2) - sum(vec*point)^2)
  mat=mat[order(PC)]
  avg[,cl:=1:.N]
  avg_long=melt(avg,id.vars=c("cl","lower_slope","upper_slope"),
                measure.vars=patterns(mu=".mu$",se=".se$"))
  avg_long[,PC:=INPUT$pcs_wanted[ variable ]]
  avg_long=avg_long[order(PC)]
  mat[,d1:=g(avg_long[cl==1]$mu, Delta), by="Label"]
  mat[,d2:=g(avg_long[cl==2]$mu, Delta), by="Label"]
  if(nrow(avg)==3) {
    mat[,d3:=g(avg_long[cl==3]$mu, Delta), by="Label"]
  } else {
    mat$d3=Inf
  }
  ## infer nearest
  mat[,inferred_group:=apply(cbind(d1,d2,d3),1,which.min)]
  mat=merge(mat, avg_long[,.(PC,inferred_group=cl,mu,se,lower_slope, upper_slope)],by=c("PC","inferred_group"))
  ## simple flags, based on sign
  mat[!(PC=="PC11" & inferred_group==1), flag_old:=sign(Delta) != sign(mu),by="Label"] # don't infer if wrong side of axis
  mat[,flag_old:=any(flag_old==TRUE,na.rm=TRUE),by="Label"] # except cl=1, PC11 which is almost along axis

  ## ? better flags, based on overlap
  ss=split(mat, mat$Label)
  choices=lapply(names(ss), choose_vector, ss, AVG)
  exact1=sapply(choices,sum)==1
  final_choice=unlist(ifelse(exact1==TRUE, sapply(choices,which), NA))
  names(final_choice)=names(ss)
  mat[,inferred_group_old:=inferred_group]
  mat[,inferred_group:=final_choice[Label]]
  ## FLAGS=sapply(seq_along(ss), function(i) {
  ##   V=as.matrix(VARS[[  names(ss)[i] ]][INPUT$pcs_wanted, INPUT$pcs_wanted])
  ##   E=ellipse( V, centre=ss[[i]]$Delta)
  ##   E=cbind(E, yl=ss[[i]]$lower_slope[1] * E[,1], yu=ss[[i]]$upper_slope[1] * E[,1])
  ##   flag=all(E[,2] < E[,"yl"] & E[,2] < E[,"yu"]) ||
  ##     all(E[,2] > E[,"yl"] & E[,2] > E[,"yu"])
  ## if(debug) {
  ##   par(mfrow=c(1,1))
  ##   plot(E,type="l",xlim=c(-.1,.1),asp=1);
  ##   abline(a=0,b=ss[[i]]$lower_slope[1]); abline(a=0,b=ss[[i]]$upper_slope[1]);
  ## for(j in seq_along(AVG)) {
  ##   points(ellipse(attr(AVG[[j]],"meta")$vcov, centre=attr(AVG[[j]],"meta")$coefficients),col=j+1)
  ##   abline(a=0,b=AVG[[j]]$lower_slope[1],col=j+1,lty="dashed"); abline(a=0,b=AVG[[j]]$upper_slope[1],col=j+1,lty="dashed");
  ##   abline(a=0,b=attr(AVG[[j]],"meta")$coefficients[2] / attr(AVG[[j]],"meta")$coefficients[1],col=j+1)
  ##          }
  ## }
  ## flag
  ## })
  ## names(FLAGS)=names(ss)
  ## mat[,flag:=FLAGS[ Label ]]
  ## with(mat[PC=="PC5"], table(flag, flag_old))
  ## mat[is.na(cl) & flag==TRUE, inferred_group:=NA] # remove inference if flagged
  mat[!is.na(cl), inferred_group:=cl] #keep original group if disease
  ## if(length(INPUT$pcs_wanted)==2)
  ##   mat[inferred_group==3, inferred_group:=1]
  ## mat[,flag:=is.na(inferred_group)]

  ## pheatmap likes data.frames
  mat2=dcast(mat[!is.na(inferred_group)], cl + inferred_group + d1 + d2 + d3 + Label ~ PC, value.var="Delta") %>% as.data.frame()
  rownames(mat2)=mat2$Label
  mat2$cl  %<>%  as.factor()
  mat2$inferred_group  %<>%  as.factor()
  paletteLength=50
  myColor <- colorRampPalette(c("yellow", "white", "blue"))(paletteLength)
  myBreaks <- c(seq(min(mat$Delta), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(mat$Delta)/paletteLength, max(mat$Delta), length.out=floor(paletteLength/2)))
  ## pheatmap(test, color=myColor, breaks=myBreaks)
  ## mat2=mat2[!is.na(mat2$inferred_group),]
  mat2=mat2[order(mat2$inferred_group),]
  pheatmap(mat2[,INPUT$pcs_wanted],
           annotation_row=mat2[,c("cl","inferred_group")],breaks=myBreaks, color=myColor, cluster_rows=FALSE)
  invisible(mat)
  }

plot_heatmap=function(mat,INPUT,...) {
   ## pheatmap likes data.frames
  mat2=dcast(mat[!is.na(inferred_group)], cl + inferred_group + Label ~ PC, value.var="Delta") %>% as.data.frame()
  rownames(mat2)=mat2$Label
  mat2$cl  %<>%  as.factor()
  mat2$inferred_group  %<>%  as.factor()
  paletteLength=50
  myColor <- colorRampPalette(c("yellow", "white", "blue"))(paletteLength)
  myBreaks <- c(seq(min(mat$Delta), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(mat$Delta)/paletteLength, max(mat$Delta), length.out=floor(paletteLength/2)))
  mat2=mat2[order(mat2$inferred_group),]
  pheatmap(mat2[,INPUT$pcs_wanted],
           annotation_row=mat2[,c("cl","inferred_group")],
           breaks=myBreaks, color=myColor, cluster_rows=FALSE, ...)
 }

INPUT=INPUTgr

align_biomarkers=function(INPUT,debug=FALSE) {
  mat=data[ Label %in% INPUT$labels_wanted & PC %in% INPUT$pcs_wanted, .(Label, PC, fdr, Delta, SE=sqrt(Var.Delta))]
  mat$cl=INPUT$calls$cl[ mat$Label ]
  if(length(INPUT$pcs_wanted)==2)
    mat[cl==3,cl:=1]
  ss=split(mat,mat[,c("cl")])
  if(length(INPUT$pcs_wanted)==2) {
    AVG=lapply(ss, calc_meta_v, INPUT$pcs_wanted)
  } else {
    AVG1=lapply(ss, calc_meta_v, INPUT$pcs_wanted[1:2])
    AVG2=lapply(ss, calc_meta_v, INPUT$pcs_wanted[3:4])
  }
  ## mat=merge(mat, avg_long[,.(PC,inferred_group=cl,mu,se,lower_slope, upper_slope)],by=c("PC","inferred_group"))

  ## find matches based on overlap
  ss=split(mat, mat$Label)
  choices=if(length(INPUT$pcs_wanted)==2) {
            lapply(names(ss), choose_vector, ss, AVG, INPUT$pcs_wanted)
          } else {
            lapply(names(ss), function(nm) {
              c(FALSE,choose_vector(nm, ss, AVG2[-1], INPUT$pcs_wanted[3:4])) #, do_plot=TRUE)
            })
          }
  exact1=sapply(choices,sum)==1
  final_choice=unlist(ifelse(exact1==TRUE, sapply(choices,which), NA))
  names(final_choice)=names(ss)
  mat[,inferred_group:=final_choice[Label]]
  mat[!is.na(cl), inferred_group:=cl] #keep original group if disease

 plot_heatmap(mat, INPUT)
  invisible(mat)
  }

result58=align_biomarkers(INPUT58)
resultgr=align_biomarkers(INPUTgr)

unique(data[fdr < 0.01 & Label %in% INPUTgr$labels_wanted & PC %in% c("PC5","PC8")]$Label) %>%
  setdiff(., names(INPUTgr$calls$cl))
unique(data[fdr < 0.01 & Label %in% INPUTgr$labels_wanted & PC %in% c("PC2","PC11")]$Label) %>%
  setdiff(., names(INPUTgr$calls$cl))

plot_heatmap(result58,INPUT58,file="pc58.png",
             cluster_cols=FALSE)
plot_heatmap(resultgr,INPUTgr,file="pcgr.png", cluster_cols=FALSE)

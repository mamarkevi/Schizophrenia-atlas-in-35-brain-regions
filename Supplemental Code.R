library(openxlsx)
library(SAJR)
library('DESeq2')
library('dplyr')
library('org.Hs.eg.db')
library(biomaRt)

# functions
plotGeomxExplVar = function(varf,sgn,fstat=FALSE,df,diags,features='proteins',selected='SZ-related',title=''){
  
  d = list()
  
  if(fstat){
    v = varf[,4:15]
    c1 = paste('resid',colnames(v),diags[1])
    c2 = paste('resid',colnames(v),diags[2])
    v = v/(varf[,c1]+varf[,c2])
    c1 = paste(colnames(v),diags[1])
    c2 = paste(colnames(v),diags[2])
    v = sweep(v,2,df[c1]+df[c2]-2,'*')
  }else
    v = varf[,4:15]*100
  sd = apply(v,2,sd)/sqrt(nrow(varf))
  m = apply(v,2,mean)
  d$all.mean = list(m = m,
                    l = m - sd,
                    h = m + sd,
                    main=paste0('All ',features,', mean +- SD'))
  d$all.median= list(m = apply(v,2,median),
                     l = apply(v,2,quantile,prob=0.25),
                     h = apply(v,2,quantile,prob=0.75),
                     main=paste0('All ',features,', median 25-75% quantiles'))
  sd = apply(v[sgn,],2,sd)/sqrt(sum(sgn))
  m  = apply(v[sgn,],2,mean)
  d$sgn.mean = list(m = m,
                    l = m - sd,
                    h = m + sd,
                    main=paste0(selected,' ',features,', mean +- SD'))
  d$sgn.median = list(m = apply(v,2,median),
                      l = apply(v,2,quantile,prob=0.25),
                      h = apply(v,2,quantile,prob=0.75),
                      main=paste0(selected,' ',features,', median 25-75% quantiles'))
  
  layout(rbind(c(1,3),
               c(2,4),
               c(5,6)))
  par(tck=-0.02,mgp=c(1.3,0.3,0),mar=c(3.1,2.5,1.5,0),oma=c(0,0,1.5,1))
  ylab = ifelse(fstat,'F','% of explained variance')
  for(i in d){
    m = i$m
    b = barplot(cbind(m[1:6],m[7:12]),beside = T,border=NA,ylim=range(0,i$h),ylab=ylab,main=i$main)
    arrows(b,i$l,b,i$h,length = 0.02,angle = 90,code=3)
    axis(1,b[,1],paste0('L',1:6))
    axis(1,b[,2],paste0('L',1:6))
    text(apply(b,2,mean),grconvertY(-0.11,'npc','user'),c('Region 22p','Region 9'),xpd=T,adj=c(0.5,0.5))
  }
  boxplot(at=c(1:6,8:13),v,main=paste0('All ',features),outline=F,frame=F,xlab='',ylab=ylab,xaxt='n')
  axis(1,1:6,paste0('L',1:6))
  axis(1,8:13,paste0('L',1:6))
  text(apply(cbind(1:6,8:13),2,mean),grconvertY(-0.11,'npc','user'),c('Region 22p','Region 9'),xpd=T,adj=c(0.5,0.5))
  boxplot(at=c(1:6,8:13),v[sgn,],main=paste0(selected,' ',features),outline=F,frame=F,xlab='',ylab=ylab,xaxt='n')
  axis(1,1:6,paste0('L',1:6))
  axis(1,8:13,paste0('L',1:6))
  text(apply(cbind(1:6,8:13),2,mean),grconvertY(-0.11,'npc','user'),c('Region 22p','Region 9'),xpd=T,adj=c(0.5,0.5))
  mtext(title,3,outer = T)
  invisible(v)
}


calcGeomxExplVar = function(v,r,l,d,norm=T){
  rs = sort(unique(r))
  ls = sort(unique(l))
  ds = sort(unique(d))
  
  
  m = lm(terms(v ~ r+l+r:l,keep.order = TRUE))
  res = anova(m)[1:3,2]
  names(res) = rownames(anova(m))[1:3]
  resid = residuals(m)
  for(rr in rs){
    for(ll in ls){
      c = resid[r==rr & l==ll & d==ds[1]]
      s = resid[r==rr & l==ll & d==ds[2]]
      total = sum((c(c,s)-mean(c(c,s)))^2)
      c = resid[r==rr & l==ll & d==ds[1]] = c - mean(c)
      s = resid[r==rr & l==ll & d==ds[2]] = s - mean(s)
      res[paste(rr,ll)] = total - sum(c(c,s)^2)
      res[paste('resid',rr,ll,ds[1])] = sum(c^2)
      res[paste('resid',rr,ll,ds[2])] = sum(s^2)
    }
  }
  rls = paste(rep(rs,each=length(ls)),rep(ls,times=length(rs)))
  res = res[c(names(res)[1:3],rls,paste('resid',rep(rls,times=length(ds)),rep(ds,each=length(rls))))]
  res/sum(res)
}

plotBarplot = function(d,o1,o2,cols=RColorBrewer::brewer.pal(3,'Set1')[2:1],pch=NULL,pch.cex=1,pv.srt=0,...){
  m = sapply(d,mean)
  mx = sapply(d,max)
  mn = sapply(d,min)
  sd = sapply(d,function(x)sd(x)/sqrt(length(x)))
  y = rbind(m[o1],m[o2])
  ysd = rbind(sd[o1],sd[o2])
  b=barplot(y,names.arg = paste0('L',1:6),col=cols,beside=T,border=NA,ylim=range(0,unlist(d),y+ysd)*1.15,...)
  arrows(b,y-ysd,b,y+ysd,angle = 90,length = 0.02,code = 3)
  pv = sapply(1:length(o1),function(i)t.test(d[[o1[i]]],d[[o2[i]]])$p.value)
  #browser()
  for(i in 1:length(pv)){
    pvt = as.numeric(strsplit(format(pv[i],digits = 1,scientific = TRUE),'e')[[1]])
    pvt = bquote(.(pvt[1]) %*% "10"^.(pvt[2]))
    text(mean(b[,i]),grconvertY(1,'npc','user'),pvt,adj = c(0.7,0.5),xpd=T,cex=0.7,srt=pv.srt)
  }
  
  
  #arrows(b,rbind(mn[o1],mn[o2]),b,rbind(mx[o1],mx[o2]),angle = 90,length = 0.05,code = 3)
  for(r in 1:nrow(b)){
    for(c in 1:ncol(b)){
      if(r==1)
        y = d[[o1[c]]]
      else
        y = d[[o2[c]]]
      if(is.null(pch))
        p = 19
      else
        p = pch[names(y)]
      points(b[r,c]+runif(length(y),-0.2,0.2),y,pch=p,cex=pch.cex)
    }
  }
}

plotVolcano = function(l2fc,pv=NULL,gnames=NULL,pch=19,cex=c(0.5,2),xlab='l2FC',bty='n',fdr.thr=0.05,topN2name=10,topHigherThan='pv',glabcol='red',plot.pv.thr=TRUE,gncex=1,gnfont=2,...){
  if(is.null(pv)){
    pv = l2fc$pv
    gnames = rownames(l2fc)
    l2fc = l2fc$l2fc
  }
  fdr = p.adjust(pv,m='BH')
  if(!is.null(topHigherThan) & topHigherThan=='pv')
    topN2name = min(topN2name,sum(pv<0.05))
  if(!is.null(topHigherThan) & topHigherThan=='fdr')
    topN2name = min(topN2name,sum(fdr<0.05))
  pv = log10(pv)
  if(sum(fdr<fdr.thr) == 0){
    fdr = log10(fdr.thr/length(pv))
  }else{
    fdr = max(pv[fdr<0.05])
  }
  pvcex = -pv
  pvcex = (pvcex-min(pvcex))/(max(pvcex)-min(pvcex))*(cex[2]-cex[1]) + cex[1]
  plot(l2fc,pv,bty=bty,pch=pch,ylim=rev(range(pv,fdr,log10(0.05))),xlab=xlab,ylab='log10(pv)',cex=pvcex,...)
  
  inx = order(pv)[1:topN2name]
  par(xpd=NA)
  text(l2fc[inx],pv[inx],gnames[inx],col=glabcol,adj=c(-0.1,-0.1),cex=gncex[inx],font=gnfont[inx])
  par(xpd=FALSE)
  if(plot.pv.thr){
    hor = c(log10(0.05),fdr)
    txt = c('P = 0.05','FDR = 0.05')
  }else{
    hor = c(fdr)
    txt = c('FDR = 0.05')
  }
  abline(h=hor,col='gray',lty=2)
  text(grconvertX(1,'npc','user'),hor,txt,adj=c(1,1.1))
  abline(v=0,col='gray',lty=2)
}


plot2LMs = function(m1,m2,...){
  t = table(m1,m2,useNA='always')
  t = t[,ncol(t):1]
  rownames(t)[is.na(rownames(t))] = 'NA'
  colnames(t)[is.na(colnames(t))] = 'NA'
  c = t
  c['NA',] = c[,'NA'] = NA
  t['NA','NA'] = NA
  imageWithText(c,t,...)
}

calcMeanL2Lexp = function(g,m,mlayers,agglayers,agg.fun=mean,test=t.test){
  agglayers = as.character(agglayers)
  g = g[!is.na(mlayers),rownames(m)]
  mlayers = mlayers[!is.na(mlayers)]
  ourlayers = sort(unique(m$Layer))
  means = pvs = matrix(0,ncol=length(ourlayers),nrow=length(agglayers),dimnames = list(agglayers,ourlayers))
  for(al in agglayers){
    g. = g[grep(al,mlayers,fixed = T),,drop=F]
    for(l in ourlayers){
      #if(al=='layer 2/3/4' & l=='L6') browser()
      means[al,l] = agg.fun(g.[,m$Layer==l])
      pvs[al,l] = test(g.[,m$Layer==l],g.[,m$Layer!=l],alternative='greater')$p.value
    }
  }
  list(means=means,pvs=pvs)	
}

countGeneWithMaxExp = function(g,m,mlayers,agglayers){
  g = g[!is.na(mlayers),rownames(m)]
  mlayers = mlayers[!is.na(mlayers)]
  ourlayers = sort(unique(m$Layer))
  res = matrix(0,ncol=length(ourlayers),nrow=length(agglayers),dimnames = list(agglayers,ourlayers))
  irs = unique(m[,c('Region','Individual')])
  
  uml = unique(mlayers)
  uml = revList( lapply(setNames(agglayers,agglayers),function(x)uml[grep(x,uml,fixed = T)]))
  
  for(i in 1:nrow(irs)){
    m. = m[m$Region==irs$Region[i] & m$Individual == irs$Individual[i],]
    g. = g[,rownames(m.)]
    for(j in 1:nrow(g)){
      res[uml[[mlayers[j]]],m.$Layer[which.max(g.[j,])]] = res[uml[[mlayers[j]]],m.$Layer[which.max(g.[j,])]] + 1
    }
  }
  res
}

compGandP = function(g,p,mg,mp,l,r,log=TRUE){
  gs = g[,mg$Segment.tags=='SZ' & mg$Layer == paste0('L',l) & mg$Region == paste0('Region ',r)]
  gc = g[,mg$Segment.tags=='HC' & mg$Layer == paste0('L',l) & mg$Region == paste0('Region ',r)]
  ps = p[,mp$Segment.tags=='SZ' & mp$Layer == paste0('L',l) & mp$Region == paste0('Region ',r)]
  pc = p[,mp$Segment.tags=='HC' & mp$Layer == paste0('L',l) & mp$Region == paste0('Region ',r)]
  if(log){
    gs = log2(gs)
    gc = log2(gc)
    ps = log2(ps)
    pc = log2(pc)
  }
  data.frame(ge = apply(gs,1,mean) - apply(gc,1,mean),
             pe = apply(ps,1,mean) - apply(pc,1,mean))
}

compGandPall = function(g,p,m,rs=c('9','22P'),ls=1:6,log=TRUE,cor.m='p'){
  res = NULL
  for(gn in names(g)[-1])
    for(pn in names(p)[-1])
      for(l in ls)
        for(r in rs){
          t = compGandP(g[[gn]],p[[pn]],m,l,r,log=log)
          ct = cor.test(t[,1],t[,2],m=cor.m)
          res = rbind(res,data.frame(gnorm=gn,pnorm=pn,region=r,layer=l,p=ct$p.value,cc=ct$estimate,gmean=mean(t[,1]),pmean=mean(t[,2])))
        }
  res
}




# RNA 35 regions and 10 regions analysis

# load data
hc = read.xlsx('~/HC_counts.xlsx',sheet = 1, skip=4)
sz = read.table('~/SZ_counts_an.xlsx',sheet = 1, skip=4)
counts = cbind(hc,sz)
rownames(counts) = counts[,1]
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                  dataset = "hsapiens_gene_ensembl", 
                  host = "https://ensembl.org")

protein_genes = getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','ensembl_peptide_id'),
                      mart=ensembl)
protein_gene_id = unique(protein_genes$ensembl_gene_id[protein_genes$ensembl_peptide_id!=''])

counts = counts[protein_gene_id,] # only protein coding
# ST 2 - sample info
meta = read.xlsx('~/metadata.xlsx',sheet = 1, skip=2) # ST 2 - sample info
rownames(meta) = meta$Sample 

# Deseq normalization data
dds = DESeqDataSetFromMatrix(countData = counts,
                                colData = meta,
                                design = ~ diagnosis)
dds = estimateSizeFactors(dds)

# extract the counts
counts.d = counts(dds, normalized=TRUE)

# log-transformatio

# log 
log2_counts = t(apply(counts.d,1, function(x){
  ifelse(x>0,log2(x),log2(x+1))
}))

# reduce RIN effect
rin_counts = t(apply(log2_counts,1, function(x){
  model = lm(x ~ meta[colnames(log2_counts),]$rin)
  residuals(model)
}))

# z-centering
rin_counts_norm_by_donor = lapply(unique(meta$donor), function(d){
    rin_df = rin_counts[,rownames(meta[meta$donor==d,])]
    m = t(apply(rin_df,1,function(x){
      x - mean(x)
    }))
    m
  })
rin_counts_znorm = do.call(cbind,rin_counts_norm_by_donor)

# ANOVA
pv = lapply(1:nrow(rin_counts_znorm),function(i){
  meta$y = rin_counts_znorm[i,]
  res_aov = aov(y ~ Diagnosis*Brain_region,
                data = meta)
  list(summary(res_aov)[[1]]['Pr(>F)'][1,],
       summary(res_aov)[[1]]['Pr(>F)'][2,],
       summary(res_aov)[[1]]['Pr(>F)'][3,])
})

pv_diag_bh = p.adjust(sapply(pv, function(x) x[[1]]), 'BH') 
pv_reg_bh = p.adjust(sapply(pv, function(x) x[[2]]), 'BH') 
pv_diag_reg_bh = p.adjust(sapply(pv, function(x) x[[3]]), 'BH')

# Deconvolution
ref = readRDS('C:/Mary/AS/references/pseudobulk_ba32_cbl.rds') # Siletti reference
counts_genes_info = getGeneLengthAndGCContent(gen_c_list,org ='hsa',mode=c("biomart", "org.db"))
genelength = counts_genes_info[,1]
counts_tpm = convertCounts(countsMatrix=as.matrix(counts),unit='TPM',
                           geneLength=genelength.c)

com_genes = intersect(rownames(counts_tpm),rownames(ref))

ref.f = ref[com_genes,]
counts_tpm = counts_tpm[com_genes,]
genelength1 = genelength[com_genes]

ref.f.tpm = get_TPM(ref.com,effLen=as.vector(genelength1)) 

decon = deconvolute(m = counts_tpm, sigMatrix = ref.f.tpm, use_cores = 1,
                    methods = c('ols','nnls','qprog','qprogwc','dtangle'))

prop = decon$proportions$dtangle_sig1 # the best


# FC by region with cell types adjustment
counts1_by_reg = lapply(unique(meta$Brain_region), function(x){
  counts[,rownames(meta[meta$Brain_region==x & meta[,diag_colm]=='Schizophrenia',])]
})

counts2_by_reg = lapply(unique(meta$Brain_region), function(x){
  counts[,rownames(meta[meta$Brain_region==x & meta[,diag_colm]=='Healthy control',])]
})


names(counts1_by_reg) = unique(meta$Brain_region)
names(counts2_by_reg) = unique(meta$Brain_region)

fc_by_reg = lapply(names(counts2_by_reg), function(r){
  if(is_deconv){
    mean1 = apply(counts1_by_reg[[r]],1,mean)*
      sum(prop[rownames(prop) %in% cell_types,r])/100 
    mean2 = apply(counts2_by_reg[[r]],1,mean)*
      sum(prop[rownames(prop) %in% cell_types,r])/100 
    mean1-mean2  
  }else{
    mean1 = apply(counts1_by_reg[[r]],1,mean)
    mean2 = apply(counts2_by_reg[[r]],1,mean)
    mean1-mean2  
  }
})
names(fc_by_reg) =names(counts1_by_reg)


# GeoMX analysis ####
# load preprocessed data
gmx.ge = readRDS('~/Rdata/gmx.ge.Rdata') # RNA
gmx.pe =  readRDS('~/Rdata/gmx.pe.Rdata') # proteins
gmx.pe.sum = readRDS('~/Rdata/gmx.pe.per-gene.sum.Rdata')
# load metadata
mge =  readRDS('~/Rdata/mge.Rdata') # RNA
mpe =  readRDS('~Rdata/mpe.Rdata') # protein

#normalize RNA per donor and region 
gmx.ge$`sn.Q3 Norm` = log2(gmx.ge$`Q3 Norm`)
for(r in unique(mge$Region))
  for(i in unique(mge$Individual)){
    f = mge$Region == r & mge$Individual == i
    gmx.ge$`sn.Q3 Norm`[,f] = t(scale(t(gmx.ge$`sn.Q3 Norm`[,f]),scale = F))
  }


# RNA vs protein
mpe$group = paste0(substr(mpe$Donor,1,1),'.BA',sapply(strsplit(mpe$Region,' '),'[',2),'.',mpe$Layer)
mge$group = paste0(ifelse(mge$Segment.tags=='HC','C','S'),'.BA',sapply(strsplit(mge$Region,' '),'[',2),'.',mge$Layer)

cmn = intersect(rownames(gmx.ge$`Q3 Norm`),rownames(gmx.pe.sum$`SNR without Ms IgG2a`))
cmn = cmn[-grep('IGHG',cmn)] # remove Isotype Negative Controls

pe = calcMeanCols(log2(gmx.pe.sum$`SNR without Ms IgG2a`[cmn,]),mpe$group)
ge = calcMeanCols(log2(gmx.ge$`Q3 Norm`[cmn,]),mge$group)

table(colnames(pe) == colnames(ge))
table(rownames(pe) == rownames(ge))
inx = c(1:6,13:18)
colnames(pe)[-inx]
hist(sapply(1:nrow(pe),function(i)cor(pe[i,inx],ge[i,inx])),-10:10/10)
hist(sapply(1:nrow(pe),function(i)cor(pe[i,inx]-pe[i,-inx],ge[i,inx]-ge[i,-inx])),-10:10/10)

aov.snr = apply(gmx.pe.sum$`SNR without Ms IgG2a` ,1,function(x)anova(lm(x ~ Region * Segment.tags * Layer,mpe)))
apply(sapply(aov.snr,rownames),1,function(x)length(unique(x)))

pvs.snr = t(sapply(aov.snr,function(x)x[,5]))
colnames(pvs.snr) = rownames(aov.snr[[1]])
qvs.snr = apply(pvs.snr,2,p.adjust,m='BH')
apply(qvs.snr<0.05,2,sum)

# PE test ######
aov.snr = apply(gmx.pe$`SNR without Ms IgG2a` ,1,function(x)anova(lm(x ~ Region * Segment.tags * Layer,mpe)))
apply(sapply(aov.snr,rownames),1,function(x)length(unique(x)))

pvs.snr = t(sapply(aov.snr,function(x)x[,5]))
colnames(pvs.snr) = rownames(aov.snr[[1]])
qvs.snr = apply(pvs.snr,2,p.adjust,m='BH')
apply(qvs.snr[-grep('IGHG',gmx.pe$genes$V2),]<0.05,2,sum)

gmx.pe$genes[grep('IGHG',gmx.pe$genes$V2),]

manova(t(gmx.pe$`SNR without Ms IgG2a`) ~  Region * Segment.tags * Layer,mpe)
aov(t(gmx.pe$`SNR without Ms IgG2a`)[,1] ~  Region * Segment.tags * Layer,mpe)

# layer as numerical
mpe$layer=as.numeric(substr(mpe$Layer,2,3))
aov.snr.lr = apply(gmx.pe$`SNR without Ms IgG2a` ,1,function(x)anova(lm(x ~ Layer + Layer*Region + Region * Segment.tags * layer,mpe)))
apply(sapply(aov.snr.lr,rownames),1,function(x)length(unique(x)))

pvs.snr.lr = t(sapply(aov.snr.lr,function(x)x[,5]))
colnames(pvs.snr.lr) = rownames(aov.snr.lr[[1]])
qvs.snr.lr = apply(pvs.snr.lr,2,p.adjust,m='BH')
apply(qvs.snr.lr<0.05,2,sum)
qvs.snr.lr[qvs.snr.lr[,6]<1,,drop=F]
qvs.snr.lr[order(qvs.snr.lr[,6]),]

pe.tt = t(apply(gmx.pe$`SNR without Ms IgG2a`,1,function(x){
  p = c()
  for(r in unique(mpe$Region))
    for(l in unique(mpe$Layer)){
      f = mpe$Region==r & mpe$Layer==l
      p[paste(r,l)] = t.test(x[f & mpe$Segment.tags=='Control'],x[f & mpe$Segment.tags=='SZ'])$p.value
    }
  p
}))

pe.tt.qv = apply(pe.tt,2,p.adjust,m='BH')
apply(pe.tt.qv<0.05,2,sum)
pe.tt[apply(pe.tt.qv<0.05,1,sum)>0,]


# per layer test ####
pv = NULL
cn = c()
for(r in unique(mpe$Region)){
  for(l in sort(unique(mpe$Layer))){
    f = mpe$Region == r & mpe$Layer == l
    f1 = mpe$Segment.tags[f]==mpe$Segment.tags[1]
    pv = cbind(pv,apply(gmx.pe$`SNR without Ms IgG2a`[,f],1,function(x){t.test(x[f1],x[!f1])$p.value}))
    cn = c(cn,paste0(r,'.',l))
  }
}
colnames(pv) = gsub('Region','',cn)
dim(pv)
pv = pv[-grep('IGHG',gmx.pe$genes$V2),]
pv[1:2,]
qv = apply(pv,2,p.adjust,m='BH')
pdf('figures/sz/geomx/prot.perreg-layer.ttest.pdf',w=9,h=6)
par(mfrow=c(2,3),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(5,2.5,2.5,0),oma=c(0,0,1.5,1))
barplot(apply(qv<0.05,2,sum),las=3,ylab='# of sign. proteins',main='FDR < 0.05')
barplot(apply(qv<0.2,2,sum),las=3,ylab='# of sign. proteins',main='FDR < 0.2')
barplot(apply(qv<0.5,2,sum),las=3,ylab='# of sign. proteins',main='FDR < 0.5')

barplot(apply(pv<0.001,2,sum),las=3,ylab='# of sign. proteins',main='pv < 0.001')
barplot(apply(pv<0.05,2,sum),las=3,ylab='# of sign. proteins',main='pv < 0.05')
barplot(apply(pv<0.1,2,sum),las=3,ylab='# of sign. proteins',main='pv < 0.1')
dev.off()

table(ttest=apply(pv<0.05,1,sum)>0,qvs.snr[-grep('IGHG',gmx.pe$genes$V2),2]<0.05)

pv[order(pv[,7])[1:5],]
paste(rownames(qv[apply(qv<0.05,1,sum)>0,]),collapse = ', ')






library(openxlsx)
library(SAJR)
source('code/r.functions/geomx.F.R')
source('~/skoltech/r.code/util.R')

gmx.ge = readRDS('Rdata/sz/geomx/gmx.ge.Rdata')
mge = readRDS('Rdata/sz/geomx/mge.Rdata')
gmx.pe = readRDS('Rdata/sz/geomx/gmx.pe.Rdata')
mpe = readRDS('Rdata/sz/geomx/mpe.Rdata')
#gmx.pe.sum = readRDS('Rdata/sz/geomx/gmx.pe.per-gene.sum.Rdata')

# load ge #####
# n = getSheetNames('input/geomx/All Data WTA (1).xlsx')
# gmx.ge = lapply(setNames(n,n),function(s)read.xlsx('input/geomx/All Data WTA (1).xlsx',s))
# mge = gmx.ge$SegmentProperties
# gmx.ge = gmx.ge[-1]
# for(i in 1:length(gmx.ge)){
# 	rownames(gmx.ge[[i]]) = gmx.ge[[i]][,1]
# 	gmx.ge[[i]] =  gmx.ge[[i]][,-1]
# }
# 
# for(i in 2:4){
# 	gmx.ge[[i]] =  as.matrix(gmx.ge[[i]])
# }
# 
# class(gmx.ge) = c('sajr','list')
# table(apply(sapply(gmx.ge,rownames),1,function(x)length(unique(x))))
# 
# table(rownames(gmx.ge$TargetProperties) %in% rownames(gmx.ge$`Probe QC`))
# 
# gmx.ge$TargetProperties = gmx.ge$TargetProperties[rownames(gmx.ge$`Probe QC`),]
# table(apply(sapply(gmx.ge,rownames),1,function(x)length(unique(x))))
# table(apply(do.call(cbind,lapply(gmx.ge,colnames)[-1]),1,function(x)length(unique(x))))
# 
# mge[1:2,]
# ann = read.xlsx('input/geomx/Annotation template file WTA.xlsx')
# ann[1:2,]
# table(paste(ann$Scan.name,'|',ann$`ROI.(label)`,'|',ann$`Segment.(Name/.Label)`) == mge$SegmentDisplayName)
# table(gsub(' ','.',paste0(ann$Scan.name,'|',ann$`ROI.(label)`,'|',ann$`Segment.(Name/.Label)`),fixed = T) == colnames(gmx.ge$`Probe QC`))
# rownames(ann) = colnames(gmx.ge$`Probe QC`)
# mge = cbind(ann,mge)
# mge[1:2,]
# sapply(gmx.ge,dim)
# table(apply(gmx.ge$`Q3 Norm`,2,quantile,prob=0.75))
# 
# 
# # rawge = read.xlsx('input/geomx/Initial Dataset.xlsx',2)
# # rawge[1:10,]
# # setdiff(rawge$TargetName,rownames(gmx.ge$`Probe QC`))
# # gmx.ge['NegProbe-WTX',] # probably it is mean?
# #
# # colnames(gmx.ge$`Probe QC`) %in% gsub('.|.','|',colnames(rawge),fixed = T)
# # colnames(rawge) = gsub('.|.','|',colnames(rawge),fixed = T)
# # colnames(gmx.ge$`Probe QC`) %in% colnames(rawge)
# #
# #
# # negprobes = as.matrix(rawge[rawge$TargetName=='NegProbe-WTX',colnames(gmx.ge$`Probe QC`)])
# # rawge = rawge[rawge$TargetName!='NegProbe-WTX',]
# # rownames(rawge) = rawge$TargetName
# # rawge = as.matrix(rawge[rownames(gmx.ge$`Probe QC`),colnames(gmx.ge$`Probe QC`)])
# # table(rawge == gmx.ge$`Probe QC`,useNA='always')
# #
# # # reasonable variants are bit different
# # plot(apply(negprobes,2,mean),gmx.ge$`Probe QC`['NegProbe-WTX',],log='xy')
# # plot(apply(negprobes,2,median),gmx.ge$`Probe QC`['NegProbe-WTX',],log='xy')
# # plot(apply(negprobes,2,function(x)exp(mean(log(x)))),gmx.ge$`Probe QC`['NegProbe-WTX',],log='xy')
# # plot(apply(negprobes,2,function(x)exp(mean(log(x))))-gmx.ge$`Probe QC`['NegProbe-WTX',])
# # abline(a=0,b=1,col='red')
# #
# # plot(apply(negprobes,2,function(x){mean(x)+2*sd(x)}),mge$LOQ,log='xy')
# # plot(apply(negprobes,2,function(x){x=log(x);exp(mean(x))+2*(exp(sd(x)))}),mge$LOQ,log='xy')
# # plot(apply(negprobes,2,function(x){x=log(x);exp(mean(x)+2*sd(x))}),mge$LOQ,log='xy')
# # plot(apply(negprobes,2,function(x){x=log(x);exp(mean(x))*exp(sd(x))^2}),mge$LOQ,log='xy')
# #
# # sapply(2:4,function(j)sapply(1:ncol(gmx.ge$`Probe QC`),function(i)cor(gmx.ge$`Probe QC`[,i],gmx.ge[[j]][,i])))
# #
# gmx.ge$passed.loq = matrix(FALSE,ncol=ncol(gmx.ge$`Probe QC`),nrow=nrow(gmx.ge$`Probe QC`),dimnames = dimnames(gmx.ge$`Probe QC`))
# for(i in 1:ncol(gmx.ge$`Probe QC`))
# 	gmx.ge$passed.loq[,i] = gmx.ge$`Probe QC`[,i] >= mge$LOQ[i]
# hist(apply(gmx.ge$passed.loq,2,sum),50)
# 
# # _load layer markers ####
# # NN2017
# lm2017 = read.xlsx('input/geomx/NN2017Genes.xlsx',1)
# lm2017[1:2,]
# boxplot(lm2017$`Section.specificity.in.human.(Log2FC)` ~ is.na(lm2017$Layer.marker.in.human))
# t=table(lm2017$Gene.symbol[!is.na(lm2017$Layer.marker.in.human)])
# t[t>1]
# lm2017[lm2017$Gene.symbol=='CRYBG3',]
# # ENSG00000233280.2 was Retired
# lm2017 = lm2017[lm2017$EnsemblID!='ENSG00000233280.2',]
# 
# markers = unique(lm2017[!is.na(lm2017$Layer.marker.in.human) & lm2017$Gene.symbol %in% rownames(gmx.ge$`Probe QC`),c('Gene.symbol','Layer.marker.in.human','Section.specificity.in.human.(Log2FC)')])
# dim(markers)
# table(markers$Layer.marker.in.human)
# table(table(markers$Gene.symbol))
# rownames(markers) = markers$Gene.symbol
# #markers = markers[order(markers$Layer.marker.in.human,markers$`Section.specificity.in.human.(Log2FC)`),]
# gmx.ge$TargetProperties$lm.nn2017 = NA
# cmn = intersect(rownames(gmx.ge$TargetProperties),rownames(markers))
# gmx.ge$TargetProperties[cmn,'lm.nn2017'] = markers[cmn,'Layer.marker.in.human']
# 
# # Nowakowski
# nowakowsk.markers = openxlsx::read.xlsx('input/geomx/aap8809_Nowakowski_SM-Tables-S9-S11.xlsx',1)
# rownames(nowakowsk.markers) = nowakowsk.markers$Marker
# nowakowsk.markers[1:2,]
# nowakowsk.markers
# table(nowakowsk.markers$Marker %in% rownames(gmx.ge$`Q3 Norm`))
# cmn = intersect(nowakowsk.markers$Marker ,rownames(gmx.ge$`Q3 Norm`))
# 
# gmx.ge$TargetProperties$lm.nowakowski = NA
# gmx.ge$TargetProperties[cmn,'lm.nowakowski'] = nowakowsk.markers[cmn,'Annotation']
# 
# # Maynard
# maynard.markers = read.xlsx('input/geomx/Maynard et al. visium layer markers 41593_2020_787_MOESM3_ESM.xlsx',8)
# maynard.markers[1:2,]
# maynard.markers = sapply(split(maynard.markers$models,maynard.markers$gene_name),function(x)paste0('layer ',paste0(sort(unique(unlist(strsplit(gsub('Layer','',x,fixed = T),'+',fixed = T)))),collapse = '/')))
# 
# cmn = intersect(names(maynard.markers),rownames(gmx.ge$`Q3 Norm`))
# 
# gmx.ge$TargetProperties$lm.maynard = NA
# gmx.ge$TargetProperties[cmn,'lm.maynard'] = maynard.markers[cmn]
# gmx.ge$TargetProperties[1:2,]
# mge$Region = sub('P','p',mge$Region)
# saveRDS(gmx.ge,'Rdata/sz/geomx/gmx.ge.Rdata')
# saveRDS(mge,'Rdata/sz/geomx/mge.Rdata')
# load pe #####

# n = getSheetNames('input/geomx/proteins/All Data Human Neuro Protein.20210422.olgas.annotation.xlsx')
# gmx.pe = lapply(setNames(n,n),function(s)read.xlsx('input/geomx/proteins/All Data Human Neuro Protein.20210422.olgas.annotation.xlsx',s))
# sapply(gmx.pe,dim)
# gmx.pe$`Probe Groups`
# 
# gene.names = read.table('input/geomx/proteins/genesSel (1).txt',sep='\t')
# gene.names$V1 = gsub(' ','.',gene.names$V1,fixed = T)
# setdiff(gene.names$V2, rownames(gmx.ge$`Probe QC`))
# setdiff(gene.names$V1 , colnames(gmx.pe$QC))
# setdiff(colnames(gmx.pe$QC),gene.names$V1)
# 
# t = table(gene.names$V2)
# gene.names[gene.names$V2 %in% names(t)[t>1],]
# 
# sapply(gmx.pe,dim)
# probe.groups = gmx.pe$`Probe Groups`
# gmx.pe = gmx.pe[-1]
# 
# gmx.pe$QC[1:10,1:10]
# table(apply(sapply(gmx.pe,colnames),1,function(x)length(unique(x))))
# 
# roicn = !(colnames(gmx.pe$QC) %in% gene.names$V1)
# mpe = gmx.pe$QC[,roicn]
# for(i in 1:length(gmx.pe))
# 	print(table(gmx.pe$QC[,roicn] == gmx.pe[[i]][,roicn]))
# 
# colnames(gmx.pe$QC)
# 
# gmx.pe. = gmx.pe
# 
# myroi = read.csv('input/geomx/proteins/my.roi.ann.csv')
# rownames(mpe)  = paste0('s',substr(mpe$Scan_ID,1,1),'.r',mpe$ROI_ID)
# rownames(myroi) = paste0('s',myroi$slide,'.r',myroi$roi)
# table(rownames(mpe) %in% rownames(myroi))
# myroi[setdiff(rownames(myroi),rownames(mpe)),] # as far as i remember this one was used as negative control
# myroi = myroi[rownames(mpe),]
# table(myroi$layer==mpe$Layer)
# table(myroi$region,mpe$Region)
# table(myroi$donor,mpe$Segment.tags)
# mpe$Donor = myroi$donor
# 
# mge[1:2,]
# mpe[1:2,]
# 
# for(i in names(gmx.pe)){
# 	gmx.pe[[i]] = t(as.matrix(gmx.pe[[i]][,gene.names$V1]))
# 	rownames(gmx.pe[[i]]) = gene.names$V2
# 	colnames(gmx.pe[[i]]) = rownames(mpe)
# }
# 
# gmx.pe$QC[1:10,1:10]
# range(sapply(2:length(gmx.pe),function(j)sapply(1:ncol(gmx.pe$QC),function(i)cor(gmx.pe$QC[,i],gmx.pe[[j]][,i])))) # all normalizations are proportional
# 
# gmx.pe = c(genes=list(gene.names),gmx.pe)
# class(gmx.pe) = c('sajr','list')
# 
# # mult probes per protein
# t = table(gmx.pe$genes$V2)
# t[t>1]
# t = gmx.pe$genes[gmx.pe$genes$V2 %in% names(t)[t>1],]
# t = t[order(t$V2),]
# t
# 
# pdf('figures/sz/geomx/comp.proteins.with.multiple.probes.pdf',w=7,h=7)
# for(pid in unique(t$V2))
# 	pairs(t(gmx.pe$`HK Norm without GAPDH`[gmx.pe$genes$V2==pid,]),labels = gmx.pe$genes$V1[gmx.pe$genes$V2==pid],panel = function(x,y,...)
# 		{plotLine(x,y,new = F,pch=19,line.lwd = 2,plot.ci = T,...);abline(a=0,b=1,col='gray',lty=2)})
# dev.off()
# # deduplicate; main variant should have highest expression
# gmx.pe$genes$pid = gmx.pe$genes$V2
# gmx.pe$genes$pid[gmx.pe$genes$V1=='Amyloid-Beta.1-40'] = 'APP40'
# gmx.pe$genes$pid[gmx.pe$genes$V1=='Amyloid-Beta.1-42'] = 'APP42'
# gmx.pe$genes$pid[gmx.pe$genes$V1=='Ms.IgG1'] = 'IGHG1m'
# gmx.pe$genes$pid[gmx.pe$genes$V1=='Phospho-Tau.(S404)'] = 'pMAPT'
# gmx.pe$genes$pid[gmx.pe$genes$V1=='Phospho-Alpha-synuclein.(S129)'] = 'pSNCA'
# gmx.pe$genes$pid[gmx.pe$genes$V1=='Phospho-Tdp-43.(S409/S410)'] = 'pTARDBP'
# table(table(gmx.pe$genes$pid))
# for(i in names(gmx.pe))
# 	rownames(gmx.pe[[i]]) = gmx.pe$genes$pid
# gmx.pe['pSNCA',]
# 
# 
# p = gmx.pe
# for(i in names(p)[-1]){
# 	p[[i]] = t(sapply(split.data.frame(p[[i]],p$genes$V2),function(x)apply(x,2,sum)))
# }
# p$genes=data.frame(V2=rownames(p$QC))
# rownames(p$genes) = p$genes$V2
# apply(sapply(p, rownames),1,function(x)length(unique(x)))
# 
# saveRDS(gmx.pe,'Rdata/sz/geomx/gmx.pe.Rdata')
# saveRDS(p,'Rdata/sz/geomx/gmx.pe.per-gene.sum.Rdata')
# saveRDS(mpe,'Rdata/sz/geomx/mpe.Rdata') 

# check layer markers #####
pdf('figures/sz/geomx/comp.markers.from.papers.pdf',w=7,h=7)
par(tck=-0.01,mgp=c(6,0.7,0),mar=c(7,7,1.5,0),oma=c(0,0,0,1),las=2)
plot2LMs(gmx.ge$TargetProperties$lm.nowakowski,gmx.ge$TargetProperties$lm.maynard,xlab='Nowakowski',ylab='Maynard')
plot2LMs(gmx.ge$TargetProperties$lm.nn2017,gmx.ge$TargetProperties$lm.nowakowski,ylab='Nowakowski',xlab='NN2017')
plot2LMs(gmx.ge$TargetProperties$lm.nn2017,gmx.ge$TargetProperties$lm.maynard,ylab='Nowakowski',xlab='NN2017')
dev.off()

# norm per induvidual/region
gmx.ge$`sn.Q3 Norm` = log2(gmx.ge$`Q3 Norm`)
for(r in unique(mge$Region))
	for(i in unique(mge$Individual)){
		f = mge$Region == r & mge$Individual == i
		gmx.ge$`sn.Q3 Norm`[,f] = t(scale(t(gmx.ge$`sn.Q3 Norm`[,f]),scale = F))
	}

# gmx.ge$`na.Q3 Norm` = gmx.ge$`Q3 Norm`
# gmx.ge$`na.Q3 Norm`[!gmx.ge$passed.loq] = NA


pdf('figures/sz/geomx/layer.markers.in.our.data.pdf',w=7,h=7)
layout(matrix(1:6,ncol=3),heights = c(2,1))

exp.matrixes = c('Q3 norm + splice norm' = 'sn.Q3 Norm','Q3 norm'='Q3 Norm')
filters = list(Control=mge$Segment.tags=='HC',SZ=mge$Segment.tags=='SZ',all=rep(T,nrow(mge)))
for(filter in names(filters)){
	flt=filters[[filter]]
	for(en in names(exp.matrixes)){
		exp.name = exp.matrixes[en]
		par(oma=c(0,0,1,1),bty='n')
		for(lmsource in c('lm.nowakowski','lm.maynard')){
			par(tck=-0.01,mgp=c(6,0.7,0),mar=c(3.5,7,1.5,0),las=2)
			ml =  sort(unique(gmx.ge$TargetProperties[[lmsource]][!is.na(gmx.ge$TargetProperties[[lmsource]])]),decreasing = T)
			z=calcMeanL2Lexp(gmx.ge[[exp.name]][,flt],mge[flt,],gmx.ge$TargetProperties[[lmsource]],ml)
			pvt=apply(z$pvs,2,function(x){z=x;z[]='';z[x<0.05]='*';z[x<0.01]='**';z[x<0.001]='***';z})
			title=first2Upper(substr(lmsource,4,100))
			imageWithText(scale(t(z$means)),t(pvt),col=rev(RColorBrewer::brewer.pal(name = 'RdYlBu',n=10)),xlab='',ylab=title,main=title)
			mtext('Out Geomx',1,2.3,las=1,cex=0.7)
			par(mgp=c(2,0.7,0),mar=c(3.5,3.5,1.5,0))
			z=calcMeanL2Lexp(gmx.ge[[exp.name]][,flt],mge[flt,],gmx.ge$TargetProperties[[lmsource]],6:1)
			rownames(z$means) = 	rownames(z$pvs) = paste0('L',rownames(z$means))
			pvt=apply(z$pvs,2,function(x){z=x;z[]='';z[x<0.05]='*';z[x<0.01]='**';z[x<0.001]='***';z})
			imageWithText(scale(t(z$means)),t(pvt),col=rev(RColorBrewer::brewer.pal(name = 'RdYlBu',n=10)),xlab='',ylab=title,main=title)
			mtext('Out Geomx',1,2.3,las=1,cex=0.7)
		}
		plot.new()
		legend('topleft',bty='n',c('*   <0.05','**  <0.01','*** <0.001'))
		
		lmsource = 'lm.nn2017'
		z=calcMeanL2Lexp(gmx.ge[[exp.name]][,flt],mge[flt,],gmx.ge$TargetProperties[[lmsource]],6:1)
		rownames(z$means) = 	rownames(z$pvs) = paste0('L',rownames(z$means))
		pvt=apply(z$pvs,2,function(x){z=x;z[]='';z[x<0.05]='*';z[x<0.01]='**';z[x<0.001]='***';z})
		title=substr(lmsource,4,100)
		imageWithText(scale(t(z$means)),t(pvt),col=rev(RColorBrewer::brewer.pal(name = 'RdYlBu',n=10)),xlab='',ylab=title,main=title)
		mtext('Out Geomx',1,2.3,las=1,cex=0.7)
		mtext(paste0(filter,', ',en),3,las=1,outer = T)
	}
}
dev.off()

ml =  sort(unique(gmx.ge$TargetProperties$lm.maynard[!is.na(gmx.ge$TargetProperties$lm.maynard)]),decreasing = T)
par(mar=c(4,6,1,1),las=2)

t = t(countGeneWithMaxExp(gmx.ge$`Q3 Norm`,mge[mge$Region=='Region 22P',],gmx.ge$TargetProperties$lm.maynard,1:6))
imageWithText(scale(t),t)

plm = gmx.ge$TargetProperties[rownames(gmx.ge$TargetProperties) %in% gmx.pe$genes$V2,]
plm = plm[!is.na(plm$lm.nn2017) | !is.na(plm$lm.maynard) | !is.na(plm$lm.nowakowski),]

boxplot(gmx.pe$`SNR without Ms IgG2a`['CALB1',] ~ mpe$Layer_1)


# P-vs-G ####
mge[1:2,]
mpe[1:2,]

mpe$group = paste0(substr(mpe$Donor,1,1),'.BA',sapply(strsplit(mpe$Region,' '),'[',2),'.',mpe$Layer)
mge$group = paste0(ifelse(mge$Segment.tags=='HC','C','S'),'.BA',sapply(strsplit(mge$Region,' '),'[',2),'.',mge$Layer)
table(mpe$group %in% mge$group)
table(mge$group %in% mpe$group)


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



pdf('figures/sz/geomx/ge.vs.pe.sz.log2fc.pdf',w=15,h=5)
par(mfrow=c(2,6),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,1.5,1))
sgn = rownames(pe) %in% rownames(qvs.snr)[qvs.snr[,2]<0.05]
for(r in c('9','22p'))
	for(l in 1:6){
		s = paste0('S.BA',r,'.L',l)
		c = paste0('C.BA',r,'.L',l)
		plotLine(ge[,s]-ge[,c],pe[,s]-pe[,c],xlab='GE, Q3 Norm',ylab='PE, SNR',main=paste0('R',r,'-L',l),bty='n',pch=19,line.lwd = 3,plot.ci = T,col=ifelse(sgn,'blue','gray'))
		abline(v=0,col='gray',lty=2)
		abline(h=0,col='gray',lty=2)
	}
dev.off()

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


pdf('figures/sz/geomx/protein.barplots.pdf',w=20,h=10)
e = gmx.pe$`SNR without Ms IgG2a`
o22 = paste('Region 22p',rep(c('Control','SZ'),each=6),rep(c(paste0('L',1:6)),times=2))
o9  = paste('Region 9',rep(c('Control','SZ'),each=6),rep(c(paste0('L',1:6)),times=2))
par(mfrow=c(4,8),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,2.5,2.5,0),oma=c(0,0,1.5,1))
donors = unique(mpe$Donor)
donpch = match(mpe$Donor,donors)
names(donpch) = rownames(mpe)
for(p in order(pvs.snr[,2])){
	pe = split(e[p,],paste(mpe$Region,mpe$Segment.tags,mpe$Layer))[c(o22,o9)]
	par(mar=c(3,2.5,1.5,0))
	plotBarplot(pe,o22[1:6],o22[7:12],main='',ylab='Expression',xlab='BA22p',pch=donpch)
	text(grconvertX(1,'nfc','user'),grconvertY(1,'nfc','user'),
			 paste0(gsub('.',' ',gmx.pe$genes[p,1],fixed=T),' (',gmx.pe$genes[p,2],', p=',format(pvs.snr[p,2],digits = 2,scientific = T),', padj=',format(qvs.snr[p,2],digits = 2,scientific = T),')'),xpd=NA,adj=c(0.5,1))
	par(mar=c(3,1.5,1.5,1))
	plotBarplot(pe,o9[1:6],o9[7:12],main='',ylab='',xlab='BA9',pch=donpch)
}
plot.new()
legend('topleft',bty='n',fill=RColorBrewer::brewer.pal(3,'Set1')[2:1],legend=c('Control','Schizophrenia'))
legend('bottomleft',bty='n',pch=1:length(donors),legend=donors)
dev.off()

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

# P.G E expl var #####
aov.snr = apply(gmx.pe$`SNR without Ms IgG2a` ,1,function(x)anova(lm(x ~ Region * Segment.tags * Layer,mpe)))
apply(sapply(aov.snr,rownames),1,function(x)length(unique(x)))

pvs.snr = t(sapply(aov.snr,function(x)x[,5]))
colnames(pvs.snr) = rownames(aov.snr[[1]])
qvs.snr = apply(pvs.snr,2,p.adjust,m='BH')
apply(qvs.snr<0.05,2,sum)
apply(qvs.snr[-grep('IGHG',gmx.pe$genes$V2),]<0.05,2,sum)


aov.ge = apply(gmx.ge$`Q3 Norm` ,1,function(x)anova(lm(x ~ Region * Segment.tags * Layer,mge)))
apply(sapply(aov.ge,rownames),1,function(x)length(unique(x)))

gmxf = readRDS('input/geomx/symbols_detected_genes_sample500tr_union_BA9_and_BA22p_unions.rds')

pvs.ge = t(sapply(aov.ge,function(x)x[,5]))
colnames(pvs.ge) = rownames(aov.ge[[1]])
qvs.ge = apply(pvs.ge,2,p.adjust,m='BH')
apply(qvs.ge<0.05,2,sum)

qvs.gef = apply(pvs.ge[gmxf,],2,p.adjust,m='BH')
apply(qvs.gef<0.05,2,sum)
rownames(qvs.gef)[qvs.gef[,2]<0.05] %in% gmx.pe$genes$V2
paste(rownames(qvs.gef)[qvs.gef[,2]<0.05],collapse=', ')


#layer as covar
mge$layer=as.numeric(substr(mge$Layer,2,3))
aov.ge.lr = apply(gmx.ge$`Q3 Norm` ,1,function(x)anova(lm(x ~ Layer + Layer*Region + Region * Segment.tags * layer,mge)))
apply(sapply(aov.ge.lr,rownames),1,function(x)length(unique(x)))

pvs.ge.lr = t(sapply(aov.ge.lr,function(x)x[,5]))
colnames(pvs.ge.lr) = rownames(aov.ge.lr[[1]])
qvs.ge.lr = apply(pvs.ge.lr,2,p.adjust,m='BH')
apply(qvs.ge.lr<0.05,2,sum)

qvs.gef.lr = apply(pvs.ge.lr[gmxf,],2,p.adjust,m='BH')
apply(qvs.gef.lr<0.05,2,sum)


varp = t(apply((gmx.pe$`SNR without Ms IgG2a`),1,function(v)calcGeomxExplVar(v,mpe$Region,mpe$Layer,mpe$Segment.tags)))
varpl = t(apply(log(gmx.pe$`SNR without Ms IgG2a`),1,function(v)calcGeomxExplVar(v,mpe$Region,mpe$Layer,mpe$Segment.tags)))
varg = t(apply((gmx.ge$`Q3 Norm`),1,function(v)calcGeomxExplVar(v,mge$Region,mge$Layer,mge$Segment.tags)))
vargl = t(apply(log(gmx.ge$`Q3 Norm`),1,function(v)calcGeomxExplVar(v,mge$Region,mge$Layer,mge$Segment.tags)))

pdf('figures/sz/geomx/proteins.explained.variance.pdf',w=8,h=8)
plotGeomxExplVar(varp,qvs.snr[,2]<0.05,title='Proteins, % var')
plotGeomxExplVar(varp,qvs.snr[,2]<0.05,T,table(paste(mpe$Region,mpe$Layer,mpe$Segment.tags)),unique(mpe$Segment.tags),title='Proteins, F')

plotGeomxExplVar(varpl,qvs.snr[,2]<0.05,title='Proteins log, % var')
plotGeomxExplVar(varpl,qvs.snr[,2]<0.05,T,table(paste(mpe$Region,mpe$Layer,mpe$Segment.tags)),unique(mpe$Segment.tags),title='Proteins log, F')
dev.off()

pdf('figures/sz/geomx/mRNA.explained.variance.pdf',w=8,h=8)
plotGeomxExplVar(varg,rownames(varg) %in% gmxf,title='mRNA, % var',features = 'genes',selected = 'filtered')
plotGeomxExplVar(varg,rownames(varg) %in% gmxf,title='mRNA, F',features = 'genes',selected = 'filtered',T,table(paste(mge$Region,mge$Layer,mge$Segment.tags)),unique(mge$Segment.tags))

plotGeomxExplVar(vargl,rownames(vargl) %in% gmxf,title='mRNA log, % var',features = 'genes',selected = 'filtered')
plotGeomxExplVar(vargl,rownames(vargl) %in% gmxf,title='mRNA log, F',features = 'genes',selected = 'filtered',T,table(paste(mge$Region,mge$Layer,mge$Segment.tags)),unique(mge$Segment.tags))
dev.off()

library(agricolae)
v = varp[qvs.snr[,2]<0.05,4:15]
v = varp[,4:15]
g = rep(colnames(varp)[4:15],each=nrow(v))
a = aov(as.numeric(v) ~ g)
anova(a)

hsd = HSD.test(a,'g')
plot(hsd,las=3)

# r = residuals(lm(log(gmx.pe$`SNR without Ms IgG2a`[1,]) ~ mpe$Region+mpe$Layer+mpe$Region:mpe$Layer))
# f = mpe$Region=='Region 9' & mpe$Layer=='L3'
# anova(lm(r[f] ~ mpe$Segment.tags[f]))[1,4]



sz.genes = readRDS('Rdata/bulk_clusters_DEGs_treshold.rds')
hgtf = loadEnsGTF('processed/index/Homo_sapiens.GRCh38.101.gtf.gz','gene')
rownames(hgtf) = hgtf$gene_id
sz.genes$gene.name = hgtf[sz.genes$gene,'gene_name']

# 22p-9 corrr ######
artem.det.genes = readRDS('input/geomx/symbols_detected_genes_sample500tr_union_BA9_and_BA22p_unions.rds')

mpe$group = paste0(substr(mpe$Donor,1,1),'.BA',sapply(strsplit(mpe$Region,' '),'[',2),'.',mpe$Layer)
mge$group = paste0(ifelse(mge$Segment.tags=='HC','C','S'),'.BA',sapply(strsplit(mge$Region,' '),'[',2),'.',mge$Layer)

pe = calcMeanCols(log2(gmx.pe.sum$`SNR without Ms IgG2a`),mpe$group)
ge = log2(gmx.ge$`Q3 Norm`)
ge[!gmx.ge$passed.loq] = NA
ge = calcMeanCols(ge,mge$group)

pba9  = pe[,paste0('S.BA9.L',1:6)]   - pe[,paste0('C.BA9.L',1:6)]
pba22 = pe[,paste0('S.BA22p.L',1:6)] - pe[,paste0('C.BA22p.L',1:6)]

gba9  = ge[,paste0('S.BA9.L',1:6)]   - ge[,paste0('C.BA9.L',1:6)]
gba22 = ge[,paste0('S.BA22p.L',1:6)] - ge[,paste0('C.BA22p.L',1:6)]


ge. = calcMeanCols(log2(gmx.ge$`Q3 Norm`),mge$group)
gba9.  = ge.[,paste0('S.BA9.L',1:6)]   - ge.[,paste0('C.BA9.L',1:6)]
gba22. = ge.[,paste0('S.BA22p.L',1:6)] - ge.[,paste0('C.BA22p.L',1:6)]



pdf('figures/sz/geomx/ba9-ba22p.sz-cor.pdf',w=7,h=7)
par(las=2,mar=c(6,6,2,2),bty='n')
c = cor(cbind(pba9,pba22))
rownames(c)=colnames(c) = substr(colnames(c),3,20)
imageWithText(c,xlab='',ylab='',col=getPal(),main='Proteins')

c = cor(cbind(gba9,gba22),u='p')
rownames(c)=colnames(c) = substr(colnames(c),3,20)
imageWithText(c,xlab='',ylab='',col=getPal(),main='mRNA')

par(mfrow=c(2,2),tck=-0.02,mgp=c(1.3,0.3,0),mar=c(3.1,2.5,1.5,0),oma=c(0,0,1.5,1),bty='n')
imageWithText(cor(pba9,pba22),xlab='BA9',ylab='BA22p',xaxlab = paste0('L',1:6),yaxlab =  paste0('L',1:6),col=getPal(),main='Proteins')
imageWithText(cor(gba9,gba22,u='p'),xlab='BA9',ylab='BA22p',xaxlab = paste0('L',1:6),yaxlab =  paste0('L',1:6),col=getPal(),main='mRNA')

dev.off()

pdf('figures/sz/geomx/ba9-ba22p.sz-cor.layers.scatterplots.pdf',w=9,h=6)
par(mfrow=c(2,3),tck=-0.02,mgp=c(1.3,0.3,0),mar=c(3.1,2.5,1.5,0),oma=c(0,0,1.5,1),bty='n')
lim = max(abs(c(pba9,pba22)))
lim=c(-lim,lim)

for(i in 1:6){
	plotLine(pba9[,i],pba22[,i],pch=19,line.lwd = 3,plot.ci = T,xlab='BA09',ylab='BA22p',xlim=lim,ylim=lim,main=paste0('L',i))
	abline(a=0,b=1,lty=2)
}
mtext('Proteins',side = 3,outer = T,las=1)


# f = apply(gmx.ge$passed.loq,1,sum) > 19
# table(f)


lim = max(abs(c(gba9.[artem.det.genes,],gba22.[artem.det.genes,])),na.rm=T)
lim=c(-lim,lim)
for(i in 1:6){
	plotLine(gba9.[artem.det.genes,i],gba22.[artem.det.genes,i],pch=19,line.lwd = 3,plot.ci = T,xlab='BA09',ylab='BA22p',xlim=lim,ylim=lim,main=paste0('L',i))
	abline(a=0,b=1,lty=2)
}
mtext(paste0('WTA, ',length(artem.det.genes),' detected'),side = 3,outer = T,las=1)


f=T
lim = max(abs(c(gba9[f,],gba22[f,])),na.rm=T)
lim=c(-lim,lim)

for(i in 1:6){
	plotLine(gba9[f,i],gba22[f,i],pch=19,line.lwd = 3,plot.ci = T,xlab='BA09',ylab='BA22p',xlim=lim,ylim=lim,main=paste0('L',i))
	abline(a=0,b=1,lty=2)
}
mtext('WTA, > loq',side = 3,outer = T,las=1)

dev.off()

f = apply(gmx.ge$passed.loq,1,sum)
barplot(table(f))
table(f>29)
imageWithText(cor(cbind(gba9[f>29,],gba22[f>29,])))
imageWithText(cor(gba9[f>29,],gba22[f>29,]))


# pseudobulk
ge = readRDS('Rdata/sz/ge.Rdata')
meta = readRDS('Rdata/sz/meta35.Rdata')
areas = readRDS('Rdata/sz/areas35.Rdata')

pepb = calcMeanCols(log2(gmx.pe.sum$`SNR without Ms IgG2a`),paste0(substr(mpe$Donor,1,1),'.',gsub('Region ','',mpe$Region)))
gepb = calcMeanCols(log2(gmx.ge$`Q3 Norm`),paste0(ifelse(mge$Segment.tags=='HC','C','S'),'.BA',sapply(strsplit(mge$Region,' '),'[',2)))


cmn = intersect(rownames(gmx.ge$`Q3 Norm`),rownames(gmx.pe.sum$`SNR without Ms IgG2a`))
cmn = cmn[-grep('IGHG',cmn)]

g = ge[ge$genes$gene_name %in% rownames(pepb),]
length(g)
gn = g$genes$gene_name
g = calcMeanCols(g$norm.ilia,paste0(meta$condition,'.',meta$Acronym))
rownames(g) = gn


pdf('figures/sz/geomx/pseudobulk.comp.pearson.pdf',w=6,h=6)
corm = 'p'
par(mfrow=c(2,2),tck=-0.02,mgp=c(1.3,0.3,0),mar=c(3.1,2.5,1.5,0),oma=c(0,0,1.5,1),bty='n')

plotLine(pepb[cmn,"S.22p"]-pepb[cmn,"C.22p"],gepb[cmn,"S.BA22p"]-gepb[cmn,"C.BA22p"],cor.method = corm,xlab='Protein l2FC',ylab='geomx mRNA l2FC',main='BA22p',line.lwd = 3,plot.ci = T,pch=19)
plotLine(pepb[cmn,"S.9"]-pepb[cmn,"C.9"],gepb[cmn,"S.BA9"]-gepb[cmn,"C.BA9"],cor.method = corm,xlab='Protein l2FC',ylab='geomx mRNA l2FC',main='BA9',line.lwd = 3,plot.ci = T,pch=19)

plotLine(pepb[gn,"S.22p"]-pepb[gn,"C.22p"],g[,"X.BA22p"]-g[,"H.BA22p"],cor.method = corm,xlab='Protein l2FC',ylab='bulk mRNA l2FC',main='BA22p',line.lwd = 3,plot.ci = T,pch=19)
plotLine(pepb[gn,"S.9"]-pepb[gn,"C.9"],g[,"X.BA9"]-g[,"H.BA9"],cor.method = corm,xlab='Protein l2FC',ylab='bulk mRNA l2FC',main='BA9',line.lwd = 3,plot.ci = T,pch=19)
dev.off()

ge[1:2,1:2]

table(table(ge$genes$gene_name))


# look for examples ######
cl = readRDS('Rdata/bulk_clusters_DEGs_treshold.rds')
ab = read.csv('input/names_genes_antibodies.csv')
ab[1:2,]
table(ab$gene %in% gmxf)

apply(qvs.ge[ab$gene,]<0.05,2,sum)
table(ab$gene %in% rownames(qvs.snr))
qvs.snr[intersect(ab$gene , rownames(qvs.snr)),]
g = intersect(ab$gene , rownames(qvs.ge)[qvs.ge[,2]<0.05])
ab[ab$gene %in% g,]

table(ab$ens_id %in% cl$gene)
cl[cl$gene %in% ab[ab$gene %in% g,'ens_id'],]

g %in% gmxf
par(mar=c(10,4,2,1))
boxplot(gmx.ge$`Q3 Norm`['F3',] ~  mge$Segment.tags + mge$Layer + mge$Region,las=3,col=2:3,xlab='')
boxplot(gmx.ge$`Q3 Norm`['CSRP1',] ~  mge$Segment.tags + mge$Layer + mge$Region,las=3,col=2:3,xlab='')
boxplot(gmx.ge$`Q3 Norm`['SLC1A3',] ~  mge$Segment.tags + mge$Layer + mge$Region,las=3,col=2:3,xlab='')

# number of detected genes ####
expn = apply(gmx.ge$passed.loq,2,sum)
mge[1:2,]

pdf('figures/sz/geomx/numer.of.detected.genes.pdf',w=5,h=4)
par(mfrow=c(1,1),tck=-0.02,mgp=c(1.3,0.3,0),mar=c(3.1,2.5,1.5,0),oma=c(0,0,1.5,1),bty='n')
plot(as.numeric(substr(mge$Layer,2,3))+runif(nrow(mge),-0.2,0.2),expn,col=(mge$Region=='Region 9')+2,pch=19,log='y',xlab='Layer',ylab='# of detected',bty='n')
legend('topright',col=2:3,legend = c('Region 22p','Region 9'),pch=19,bty='n')
dev.off()

boxplot(log(gmx.pe$`SNR without Ms IgG2a`))
o = order(apply(gmx.pe$`SNR without Ms IgG2a`,1,mean))
col=(grepl('IGHG',gmx.pe$genes$V2[o]))+1
boxplot(log(t(gmx.pe$`SNR without Ms IgG2a`[o,])),col=col,border=col)
r = apply(gmx.pe$`SNR without Ms IgG2a`,2,rank)
o = order(apply(r,1,mean),decreasing = T)
col=(grepl('IGHG',gmx.pe$genes$V2[o]))+1

pdf('figures/sz/geomx/protein.detection.pdf',w=8,h=15)
par(mfrow=c(1,1),tck=-0.01,mgp=c(1.6,0.6,0),mar=c(3.1,6.5,1.5,0),oma=c(0,0,1.5,1),bty='n')
boxplot(t(r[o,]),col=col,border=col,horizontal=T,frame=F,xlab='Rank',las=1)
dev.off()

# detected proteins #####
d=gmx.pe$`SNR without Ms IgG2a`
f = rep(F,nrow(d))
for(l in unique(mpe$Layer))
	for(r in unique(mpe$Region)){
		t = apply(d[,mpe$Layer==l & mpe$Region==r & mpe$Segment.tags=='SZ'     ]>=3,1,sum) > 1 &
				apply(d[,mpe$Layer==l & mpe$Region==r & mpe$Segment.tags=='Control']>=3,1,sum) > 1
		f = f | t
	}
table(f,grepl('IGHG',gmx.pe$genes$V2))
at = readRDS('input/detected_proteins_union_BA9_and_BA22p_unions.rds')
table(f,gmx.pe$genes$V1 %in% at)
gmx.pe$genes$good = f
saveRDS(gmx.pe,'Rdata/sz/geomx/gmx.pe.Rdata')

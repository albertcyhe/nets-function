tab<-merge(read.delim('data/processed/GSE125989/GSE125989.scores.tsv'), read.delim('data/processed/GSE125989/GSE125989.covariates.tsv'), by='sample_id')
ph<-read.delim('data/processed/GSE125989/GSE125989.pheno.tsv')
df<-ph; df$group<-ifelse(grepl('brain metast', tolower(df$title)),'met','primary')
pid<-apply(df,1,function(r){m<-grep('paired sample:', tolower(as.character(r)), value=TRUE); if(length(m)>0) sub('.*paired sample: ','',m[1]) else NA})
df$pair_id<-pid
tab<-merge(tab, df[,c('sample_id','group','pair_id')], by='sample_id')
agg<-aggregate(deltaFM_singscore ~ pair_id + group, data=tab, FUN=function(x) mean(x, na.rm=TRUE))
wide<-reshape(agg, idvar='pair_id', timevar='group', direction='wide')
print(wide[1:10,])
d1 <- wide$`deltaFM_singscore.met` - wide$`deltaFM_singscore.primary`
print(d1)
print(summary(d1))
print(wilcox.test(d1, mu=0))

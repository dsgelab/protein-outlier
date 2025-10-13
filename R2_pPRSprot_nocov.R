library(reshape2)
args=commandArgs(T)

dir= "/scratch/project_2007428/projects/prj_100_pprs_discordance/data/"

files=list.files(paste0(dir,"pred_real_protvalues/"))

alldata= data.frame()
r2tab= data.frame()
for (f in files) {
  r2table=data.frame(ncol=7)
  r2table$protname= strsplit(f, "_OP")[[1]][1]
  data= read.table(paste0(dir,"pred_real_protvalues/",f), header=F)
  colnames(data)= c("sample","predicted","real")
  model=lm(real ~ predicted, data = data)
  r2table$sample= nrow(data)
  r2table$R2=summary(model)$r.squared
  r2table$adj_R2=summary(model)$adj.r.squared
  r2table$beta= summary(model)$coefficients[2, "Estimate"]
  r2table$se= summary(model)$coefficients[2, "Std. Error"]
  r2table$pvalue= summary(model)$coefficients[2, "Pr(>|t|)"]
  r2tab=rbind(r2tab, r2table)
  
}
r2tab$ncol=NULL

write.table(r2tab, file=args[1], quote=F, sep="\t", 
            col.names = T, row.names = F)


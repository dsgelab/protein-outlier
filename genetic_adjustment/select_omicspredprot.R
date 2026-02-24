library(data.table)
library(dplyr)
args=commandArgs(T)

omicspred= fread(args[1])
omicspred$`UniProt ID`= gsub(";","_", omicspred$`UniProt ID`)
omicspred$Gene= gsub(";","_", omicspred$Gene)
omicspred$`UniProt ID`=ifelse(omicspred$`UniProt ID`=="Q8NEV9_Q14213",
                              "Q14213_Q8NEV9", omicspred$`UniProt ID`) # manual correction since in ukbb is inverted
omicspred$`UniProt ID`=ifelse(omicspred$`UniProt ID`=="P29460_P29459",
                              "P29459_P29460", omicspred$`UniProt ID`) # manual correction since in ukbb is inverted

olink=fread(args[2])

all= merge(omicspred, olink, by.x="UniProt ID",by.y="uniprot")
genes= all %>% select(`OMICSPRED ID`, Gene, `UniProt ID`)
write.table(genes, args[3], sep="\t", row.names = F,quote=F,col.names = T)
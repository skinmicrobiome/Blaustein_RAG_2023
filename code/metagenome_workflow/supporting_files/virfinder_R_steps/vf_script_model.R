library(VirFinder)

## specify the directory of the new model file "VF.modEPV_k8.rda", and load the new model to the work space
modFile <- "VF.modEPV_k8.rda"
load(modFile)

## predict the contigs using the new model
predResult <- VF.pred.user("filt_scaffolds.fasta", modEPV)

## estimate q-values (false discovery rates) based on p-values
predResult$qvalue <- VF.qvalue(predResult$pvalue)
predResult$BH_adj <- p.adjust(predResult$pvalue, method="BH")

## write tables
write.table(predResult, file="table.txt", sep = "\t",  row.names = TRUE, col.names = NA)
write.table(predResult[order(predResult$qvalue),], file="table_qval.txt", sep = "\t",  row.names = TRUE, col.names = NA)
write.table(predResult[order(predResult$BH_adj),], file="table_bh.txt", sep = "\t",  row.names = TRUE, col.names = NA)

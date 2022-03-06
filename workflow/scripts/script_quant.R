library(tximport)
library(EnsDb.Hsapiens.v86)
k <- keys(EnsDb.Hsapiens.v86,"TXID")#capture transcript IDs from Ensembl annotation
tx2gene <- ensembldb::select(EnsDb.Hsapiens.v86, k, "GENENAME","TXID")#create object to collapse transcript IDs into gene IDs
samples=read.table("./metadata/samples.tsv",sep = "\t",header = T)
files <- file.path("./results/salmon", paste0(samples$Run, "_quant"),"quant.sf")#create vector with the path to each sample
names(files) <- samples$Sample_Name
all(file.exists(files))#control that all files are reachable 
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,ignoreAfterBar = T,ignoreTxVersion = T,countsFromAbundance = "lengthScaledTPM")#read files and collapse to TPM per gene
head(txi.salmon$counts)
data=as.data.frame(txi.salmon$counts)#create table with TPM 
data$GeneName=row.names(data)
log_data=log2(data[,1:24]+1)
log_data$GeneName=row.names(log_data)
write.table(data,"./results/salmon/TPM_quant_table.tsv",row.names=F,header=T,sep="\t",quote=F)

data_transcripts=tximport(files, type = "salmon",ignoreAfterBar = T,ignoreTxVersion = T,countsFromAbundance = "lengthScaledTPM",txOut = T)
table_transcripts=as.data.frame(data_transcripts$counts)
table_transcripts$code=row.names(table_transcripts)

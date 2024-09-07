# Functional enrichment analysis
# Chun Shen, 2024

# Load R packages
library(gprofiler2)
library(data.table)

# Read PWAS results
Result_SI_m2 <- fread('Result_SI_protein_assoc_M2.csv',data.table = F)
Result_LO_m2 <- fread('Result_LO_protein_assoc_M2.csv',data.table = F)

#fdr correction performed together
Result_com <- rbind(Result_SI_m2,Result_LO_m2)
Result_com$pheno <- rep(c('SI','LO'),each=nrow(Result_LO_m2))
Result_com$fdr_com <- p.adjust(Result_com$pval,method='fdr')

fdr_SI <- Result_com[Result_com$fdr_com<0.05 & Result_com$pheno=='SI',]
fdr_LO <- Result_com[Result_com$fdr_com<0.05 & Result_com$pheno=='LO',]

#version information
get_version_info(organism = "hsapiens")

SI_gene <- fdr_SI$V1
SI <- gost(query = SI_gene,
           organism = 'hsapiens', ordered_query = FALSE,
           domain_scope = "custom",custom_bg = Result_SI_m2$V1,
           multi_query = FALSE, exclude_iea = FALSE,
           measure_underrepresentation = FALSE, evcodes = FALSE,
           sources = c('GO:BP','GO:MF','KEGG','REAC','WP'),
           significant = TRUE, user_threshold = 0.05,correction_method = "fdr")

LO_gene <- fdr_LO$V1
LO <- gost(query = LO_gene,
           organism = 'hsapiens', ordered_query = FALSE,
           domain_scope = "custom",custom_bg = Result_LO_m2$V1,
           multi_query = FALSE, exclude_iea = FALSE,
           measure_underrepresentation = FALSE, evcodes = FALSE,
           sources = c('GO:BP','GO:MF','KEGG','REAC','WP'),
           significant = TRUE, user_threshold = 0.05,correction_method = "fdr")

##save results
SI_result <- SI$result
write.csv(SI_result[,1:13],file = 'pathway_SI_fdr05.csv',row.names = F)

LO_result <- LO$result
write.csv(LO_result[,1:13],file = 'pathway_LO_fdr05.csv',row.names = F)


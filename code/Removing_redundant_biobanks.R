p.table <- fread("../data/Projection_cell_basis_v3_20220804-v1.tsv")
qc.table <- fread("../data/QC_cell_basis_v3_20220804-v1.tsv")
metadata <- fread("../data/Metadata_20230307-v1.tsv")
qc.table <- merge.data.table(qc.table, metadata, all.x = TRUE)
qc.table[, c("Collection", "Chip", "File_ID"):=NULL]


fwrite(qc.table, "../tables/Table_S3_Full_QC_table.tsv", sep = "\t")
fwrite(p.table,  "../tables/Table_S4_Full_Projection_table.tsv", sep = "\t")


# Remove redundant biobank GWAS efforts
qc.table <- qc.table[!Reference %in% c("PanUKBBR1", "UKBB", "FinnGenR5", "FinnGenR7")]
p.table <- p.table[Trait %in% unique(qc.table$Trait)]

fwrite(p.table, "../data/ptable_20230313-v1.tsv", sep = "\t")
fwrite(qc.table, "../data/qctable_20230313-v1.tsv", sep = "\t")


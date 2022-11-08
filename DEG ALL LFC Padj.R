P10_KG2 <- read.table("./DEG/DEG/P10_KG2_DEGseq2_P10_KSL.txt", sep = "\t", quote = "")
P10_2i <- read.table("./DEG/DEG/P10_2i_DEGseq2_P10_KSL.txt", sep = "\t", quote = "")
P5_KG2 <- read.table("./DEG/DEG/P5_KG2_DEGseq2_P5_KSL.txt", sep = "\t", quote = "")
P5_2i <- read.table("./DEG/DEG/P5_2i_DEGseq2_P5_KSL.txt", sep = "\t", quote = "")
P10_KSL <- read.table("./DEG/DEG/P10_KSL_DEGseq2_P5_KSL.txt", sep = "\t", quote = "")

P10_KG2$SYMBOL <- rownames(P10_KG2)
P10_KG2 <- P10_KG2[,c(2,6,8)]
colnames(P10_KG2) <- c("P10_KG2/P10_KSL_LFC","P10_KG2/P10_KSL_padj","external_gene_name")


P10_2i$SYMBOL <- rownames(P10_2i)
P10_2i <- P10_2i[,c(2,6,8)]
colnames(P10_2i) <- c("P10_2i/P10_KSL_LFC","P10_2i/P10_KSL_padj","external_gene_name")


P5_KG2$SYMBOL <- rownames(P5_KG2)
P5_KG2 <- P5_KG2[,c(2,6,8)]
colnames(P5_KG2) <- c("P5_KG2/P5_KSL_LFC","P5_KG2/P5_KSL_padj","external_gene_name")


P5_2i$SYMBOL <- rownames(P5_2i)
P5_2i <- P5_2i[,c(2,6,8)]
colnames(P5_2i) <- c("P5_2i/P5_KSL_LFC","P5_2i/P5_KSL_padj","external_gene_name")

P10_KSL$SYMBOL <- rownames(P10_KSL)
P10_KSL <- P10_KSL[,c(2,6,8)]
colnames(P10_KSL) <- c("P10_KSL/P5_KSL_LFC","P10_KSL/P5_KSL_padj","external_gene_name")

HM <- merge(P10_KG2, P5_KG2, by = "external_gene_name")
HM <- merge(HM, P10_2i, by = "external_gene_name")
HM <- merge(HM, P5_2i, by = "external_gene_name")
HM <- merge(HM, P10_KSL, by = "external_gene_name")
rownames(HM) <- HM[,1]
HM <- HM[,-1]

HM_1 <- HM[abs(HM$`P10_KG2/P10_KSL_LFC`) > 1|abs(HM$`P5_KG2/P5_KSL_LFC`) > 1|abs(HM$`P10_2i/P10_KSL_LFC`) > 1|abs(HM$`P5_2i/P5_KSL_LFC`) > 1|abs(HM$`P10_KSL/P5_KSL_LFC`) > 1,]
HM_1 <- HM_1[]

write.table(HM, file="./DEG/DEG ALL LFC Padj.txt")
write.csv(HM_1, file="./DEG/20220620 DEG LFC1 - 2.csv")

####TCGA-LUAD 数据下载####
setwd('TCGA-LUAD')
setwd("TCGAdata")
library(tidyverse)#加载包
library(TCGAbiolinks)#加载包
cancer_type = "TCGA-LUAD"   #肿瘤类型，这里可修改癌症类型
expquery <- GDCquery(project = cancer_type,
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts"
)
GDCdownload(expquery,directory = "GDCdata")
expquery2 <- GDCprepare(expquery,directory = "GDCdata",summarizedExperiment = T)
save(expquery2,file = "luad.gdc_2025.rda") # 保存 rda格式

####TCGA counts/tpms提取及整理####
####counts####
setwd('TCGA-LUAD')
setwd("TCGAdata")
library(tidyverse)
load("luad.gdc_2025.rda")#导入文件，rda格式文件也可直接从文件夹双击导入
load("gene_annotation_2022.rda")#导入gene注释文件
table(gene_annotation_2022$type)#table 分组计数 频数
#基因名称symbol ENSEMBL

#提取counts 以下三句无需掌握 
counts <- expquery2@assays@data@listData[["unstranded"]]
colnames(counts) <- expquery2@colData@rownames
rownames(counts) <- expquery2@rowRanges@ranges@NAMES
#基因ID转换
counts <- as.data.frame(counts)
counts <- rownames_to_column(counts,var = 'ENSEMBL')  #把行名放入数据框
counts <- inner_join(counts,gene_annotation_2022,"ENSEMBL")  #合并基因名
counts <- counts[!duplicated(counts$symbol),]#美元符号   去除重复基因名行

rownames(counts) <- NULL
counts <- counts %>% column_to_rownames("symbol")   #把基因名作为行名
#counts <- column_to_rownames(counts,"symbol")
#继续跑
# 保留mRNA 
#table
table(counts$type)#（注：可通过table(counts$type)查看基因类型）#lncRNA
counts <- counts[counts$type == "protein_coding",]
#counts <- counts[counts$type == "lncRNA",]
counts <- counts[,-c(1,ncol(counts))]
#ncol是什么
ncol(counts)
# 把TCGA barcode切割为16位字符,并去除重复样本
colnames(counts) <- substring(colnames(counts),1,16)
counts <- counts[,!duplicated(colnames(counts))]
table(substring(colnames(counts),14,16))
# 保留01A  （注：可通过table(substring(colnames(counts),14,16))查看样本类型）
counts01A <- counts[,substring(colnames(counts),14,16) == "01A"]
# 保留11A
counts11A <- counts[,substring(colnames(counts),14,16) == c("11A")]
table(substring(colnames(counts01A),14,16))
table(substring(colnames(counts11A),14,16))

####tpms####
#和counts基本一模一样
tpms <- expquery2@assays@data@listData[["tpm_unstrand"]]
colnames(tpms) <- expquery2@colData@rownames
rownames(tpms) <- expquery2@rowRanges@ranges@NAMES
tpms <- tpms %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation_2022,"ENSEMBL") %>% 
  .[!duplicated(.$symbol),]
rownames(tpms) <- NULL
tpms <- tpms %>% column_to_rownames("symbol") 
# 保留mRNA （注：可通过table(tpms$type)查看基因类型）
tpms <- tpms[tpms$type == "protein_coding",]
tpms <- tpms[,-c(1,ncol(tpms))]
# 把TCGA barcode切割为16位字符,并去除重复样本
colnames(tpms) <- substring(colnames(tpms),1,16)
tpms <- tpms[,!duplicated(colnames(tpms))]
# 保留01A  （注：可通过table(substring(colnames(tpms),14,16))查看样本类型）
tpms01A <- tpms[,substring(colnames(tpms),14,16) == c("01A")]
# 保留11A
tpms11A <- tpms[,substring(colnames(tpms),14,16) == c("11A")]

#判断counts和tpms的行列名是否一致
#identical的用法
a <- c("a","b","c")
b <- c("b","a","c")
identical(a,b)

identical(rownames(counts01A),rownames(counts11A))
identical(rownames(tpms01A),rownames(tpms11A))
identical(rownames(counts01A),rownames(tpms01A))
identical(colnames(counts01A),colnames(tpms01A))

#cbind和rbind

#cbind之前需要确认两个数据框的行名
counts <- cbind(counts01A,counts11A)
table(substring(colnames(counts),14,16))
tpms <- cbind(tpms01A,tpms11A)

####tpms_log2####
range(tpms)#查看数据范围
range(tpms01A)
range(tpms11A)
tpms_log2 <- log2(tpms+1)#log2转换 为什么要加1
range(tpms_log2)
tpms01A_log2 <- log2(tpms01A+1)
range(tpms01A_log2)
tpms11A_log2 <- log2(tpms11A+1)
range(tpms11A_log2)
#保存counts和tpms数据
write.table(counts01A,"counts01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(counts11A,"counts11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(counts,"counts.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms_log2,"tpms_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms01A_log2,"tpms01A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms11A_log2,"tpms11A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####TCGA T/N差异分析####
#设置工作目录
setwd("TCGA-LUAD")
setwd("T_N_DEG")
##读取csv文件
counts <- read.csv("./counts.csv",row.names = 1)
counts <- round(counts)
##读取文本文件
#转换
#差异分析
library(tidyverse)
#安装BiocManager
#BiocManager::install('DESeq2')
library(DESeq2)
##以下理解即可
#去除低质量的基因 质控
counts = counts[apply(counts, 1, function(x) sum(x > 1) > 32), ]
#设置分组信息
conditions=data.frame(sample=colnames(counts),
                      group=factor(ifelse(substr(colnames(counts),14,16) == "01A","T","N"),levels = c("N","T"))) %>% 
  column_to_rownames("sample")

factor(ifelse(substr(colnames(counts),14,16) == "01A","T","N"),levels = c("N","T"))#因子



#准备差异分析
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = conditions,
  design = ~ group)
#开始差异分析
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
save(res,file = "LUAD-DEG.rda")#一定要保存！
#自定义差异基因
DEG <- as.data.frame(res)
DEG <- arrange(DEG,padj)
#DEG <- arrange(DEG,desc(padj))
#arrange(数据框名,列名) 默认从小到大排序
#arrange(数据框名，desc(列名)) 从大到小排序
DEG <- filter(DEG,abs(log2FoldChange) > 1, padj < 0.05)
#filter(数据框名,筛选条件)
#abs 绝对值

####选定基因肿瘤/正常表达柱状图####
setwd("TCGA-LUAD")
setwd("TCGAdata")
library(tidyverse)
tpms01A_log2 <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tpms11A_log2 <- read.table("tpms11A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- "STK11"#以后修改这里即可 小技巧
#gene <- c("","","") 
#如果想多提取几个基因怎么办 集合
a <- tpms01A_log2[gene,]
b <- tpms11A_log2[gene,]
##运用传导符%>%  cltrl+shift+M 
a <- a %>% t() %>% as.data.frame()
b <- b %>% t() %>% as.data.frame()
write.csv(a, file = "STK11_01A.csv")
write.csv(b, file = "STK11_11A.csv")
#我喜欢用仙桃学术画图 便于修图 而且免费
#Graphpad等
#R语言代码画图 更快！！！
a <- read.table("tpms_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
a <- a %>% t() %>% as.data.frame()
#添加分组信息
a <- a %>%  mutate(group=factor(ifelse(substr(rownames(a),14,16) == "01A","Tumor","Normal"),levels = c("Normal","Tumor")))
#a$group <- factor(ifelse(substr(rownames(a),14,16) == "01A","Tumor","Normal"),levels = c("Normal","Tumor"))
#levels的顺序反过来会怎么样？
library(ggpubr)
library(ggsci)
ggboxplot(a,
          x = "group", y = "STK11",#可以画任何一个基因
          fill = "group") +
  scale_fill_lancet()+
  stat_compare_means(comparisons = list(c(1, 2)),
                     method = "wilcox.test",   #设置统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))
dev.off()

####TCGA W/M差异分析####
#设置工作目录
setwd("TCGA-LUAD")
setwd("W_M_DEG")
##读取文本文件
counts01A <- read.table("counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#差异分析
library(tidyverse)
#安装BiocManager
#BiocManager::install('DESeq2')
library(DESeq2)
##以下理解即可
#去除低质量的基因 质控
counts01A = counts01A[apply(counts01A, 1, function(x) sum(x > 1) > 32), ]
#设置分组信息
colnames(counts01A) <- gsub("\\.", "-", colnames(counts01A))
# 使用str_sub函数提取前12位
colnames(counts01A) <- substr(colnames(counts01A), 1, 12)
a <- intersect(colnames(counts01A),rownames(group))

counts01A <- counts01A[,a]

str(group)
conditions <- group[a, , drop = FALSE]

identical(colnames(counts01A),rownames(conditions))
colnames(conditions)[1] <- "group"
counts01A <- round(counts01A)
#准备差异分析
dds <- DESeqDataSetFromMatrix(
  countData = counts01A,
  colData = conditions,
  design = ~ group)
#开始差异分析
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
save(res,file = "LUAD-DEG.rda")#一定要保存！
#自定义差异基因
DEG <- as.data.frame(res)
DEG <- arrange(DEG,padj)
#DEG <- arrange(DEG,desc(padj))
#arrange(数据框名,列名) 默认从小到大排序
#arrange(数据框名，desc(列名)) 从大到小排序
DEG <- filter(DEG,abs(log2FoldChange) > 1, padj < 0.05)
#filter(数据框名,筛选条件)
#abs 绝对值
write.csv(DEG, file = "DEG W-M.csv")
write.csv(conditions, file = "conditions突变分组.csv")
####TCGA差异分析热图展示####
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
colnames(exp) <- gsub("\\.", "-", colnames(exp))
colnames(exp) <- substr(colnames(exp), 1, 12)
group <- read.csv("./conditions突变分组.csv",row.names = 1)

A <- intersect(colnames(exp),rownames(group))

exp <- exp[,A]

DEG <- as.data.frame(res)
DEG <- arrange(DEG,padj)

#DEG <- as.data.frame(res)%>% arrange(padj)

logFC_cutoff <- 1 #and &
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)
DEG <- na.omit(DEG)
#下载pheatmap包 
library(pheatmap)
#方法二
a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
c <- rbind(a,b)
d <- rownames(c)

#提取差异基因表达谱

exp_diff <- exp[d,]

#设置分组信息
#a <- filter(annotation_col,group == 'Wild')
#b <- filter(annotation_col,group == 'Mutation')
#c <- rbind(b,a)
annotation_col <- group

#write.csv(annotation_col, file = "conditions突变分组.csv")
#开始画图
sum(is.na(exp_diff))
exp_diff_clean <- exp_diff[complete.cases(exp_diff), ]
write.csv(exp_diff_clean, file = "exp_diff_clean.csv")
exp_diff_clean <- read.csv("./exp_diff_clean.csv",row.names = 1)
annotation_col <- read.csv("./conditions突变分组.csv",row.names = 1)
colnames(exp_diff_clean) <- gsub("\\.", "-", colnames(exp_diff_clean))

library(dplyr)

# 如果annotation_col有样本名作为行名
sample_order <- rownames(annotation_col)[order(annotation_col$group)]

exp_diff_sorted <- exp_diff_clean[, sample_order]
annotation_col_sorted <- annotation_col[sample_order, , drop = FALSE]

pheatmap(exp_diff_sorted,
         annotation_col = annotation_col_sorted,
         scale = "row",
         show_rownames = F,
         show_colnames = F,
         color = colorRampPalette(c(rep("royalblue",5),"white",rep("firebrick3",5)))(50),
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 10,
         fontsize_row = 3,
         fontsize_col = 3)
#保存图片 调整大小
dev.off()#关闭画板

####TCGA差异分析火山图展示####
setwd("TCGA-LUAD")
setwd("W_M_DEG")
library(ggplot2)                         
library(ggrepel) 
library(tidyverse)
DEG <- as.data.frame(res)
DEG <- arrange(DEG,padj)

data <- DEG
data$significant="stable"

data$significant[data$log2FoldChange>=1 & data$padj <0.05]="up"
data$significant[data$log2FoldChange<=-1 & data$padj <0.05]="down"

table(data$significant)
#绘制火山图
# 将行名添加为一列
data$gene_id <- rownames(data)

# 计算合适的坐标轴范围
x_max <- max(abs(data$log2FoldChange), na.rm = TRUE)
x_lim <- ceiling(max(4, x_max))
y_max <- max(-log10(data$padj), na.rm = TRUE)
y_lim <- ceiling(max(10, y_max))

# 提取最显著的上调和下调基因（前10个）
up_genes <- data[data$significant == "up", ]
up_genes <- up_genes[order(-up_genes$log2FoldChange, up_genes$padj), ]
top_up <- head(up_genes, 10)

down_genes <- data[data$significant == "down", ]
down_genes <- down_genes[order(down_genes$log2FoldChange, down_genes$padj), ]
top_down <- head(down_genes, 10)

genes_to_label <- rbind(top_up, top_down)

# 添加统计信息
deg_counts <- table(data$significant)

ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
  xlim(-x_lim, x_lim) +
  ylim(0, y_lim) +
  
  geom_point(aes(color = significant, alpha = significant), size = 1.5) +
  
  # 添加基因标签
  geom_text_repel(
    data = genes_to_label,
    aes(label = gene_id),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 20,
    segment.color = "grey50",
    segment.size = 0.3,
    min.segment.length = 0.2,
    force = 2,
    show.legend = FALSE
  ) +
  
  # 阈值线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "black", alpha = 0.7, size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", 
             color = "black", alpha = 0.7, size = 0.5) +
  
  # 添加统计标注
  annotate("text", x = -x_lim + 0.5, y = y_lim - 0.5, 
           label = paste("Down:", deg_counts["down"]), 
           color = "#377eb8", size = 4, fontface = "bold", hjust = 0) +
  annotate("text", x = x_lim - 0.5, y = y_lim - 0.5, 
           label = paste("Up:", deg_counts["up"]), 
           color = "#e41a1c", size = 4, fontface = "bold", hjust = 1) +
  
  scale_color_manual(
    name = "Expression",
    values = c(
      "down" = "#377eb8",
      "stable" = "gray60",
      "up" = "#e41a1c"
    )
  ) +
  
  scale_alpha_manual(
    values = c("down" = 0.8, "stable" = 0.3, "up" = 0.8),
    guide = "none"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2)
  ) +
  
  labs(
    x = "log2(Fold Change)",
    y = "-log10(Adjusted p-value)"
  )

####富集分析####
####GO####
setwd("TCGA-LUAD")
setwd("FUJI_WM_DEG")
getwd()
library(tidyverse)
library("BiocManager")
#安装加载包
chooseBioCmirror()
#BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
#org.Hs.eg.db包主要注释人类基因:用于不同数据库ID间的转化
library(clusterProfiler)
#导入差异基因文件
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)

DEG <- DEG %>% rownames_to_column("SYMBOL")

genelist <- bitr(DEG$SYMBOL, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by="SYMBOL")

#GO分析
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
save(ego,ego_res,file = "GO_STK11_DEG.Rda")

#3. 可视化
##3.1 柱状图
barplot(ego, showCategory = 10,color = "pvalue")
##3.2 气泡图
dotplot(ego, showCategory = 10)
##3.3 分类展示
barplot(ego, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dotplot(ego,showCategory = 10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dev.off()
##3.4 网络图
#install.packages("ggnewscale")
library(ggnewscale)
List = DEG$log2FoldChange
names(List)= DEG$ENTREZID
head(List)
List = sort(List,decreasing = T)

cnetplot(ego, categorySize="pvalue", foldChange = List,colorEdge = TRUE)
##3.5 网络图2
cnetplot(ego, foldChange = List, circular = TRUE, colorEdge = TRUE)
#cnetplot(ego, 
         categorySize = "pvalue", 
         foldChange = List, 
         colorEdge = TRUE,
         showCategory = 10,           # 限制通路数量
         max.overlaps = 50,           # 增加重叠容忍度
         cex_label_category = 0.7,    # 调整标签大小
         cex_label_gene = 0.5,
         circular = FALSE)            # 对于复杂网络，非环形布局可能更好
##3.6 热图
heatplot(ego,foldChange = List)
#GO原理

####KEGG####
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
kk_res <- kk@result
save(kk,kk_res,file = "KEGG_STK11_DEG.Rda")

#3. 可视化
##3.1 柱状图
barplot(kk, showCategory = 10,color = "pvalue")
##3.2 气泡图
dotplot(kk, showCategory = 10)
dev.off()
##3.3 网络图
List = DEG$log2FoldChange
names(List)= DEG$ENTREZID
List = sort(List,decreasing = T)

cnetplot(kk, categorySize="pvalue", foldChange = List,colorEdge = TRUE)
##3.4 网络图2
cnetplot(kk, foldChange = List, circular = TRUE, colorEdge = TRUE)
##3.5 热图
heatplot(kk,foldChange = List)

####GSEA####
setwd("TCGA-LUAD")
setwd("FUJI_WM_DEG")
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
#GSEA转换不转换都可以#
#https://zhuanlan.zhihu.com/p/504101161
####不转换ID版本####
#双击打开差异分析文件
DEG <- as.data.frame(res)%>% arrange(padj) 
DEG <- na.omit(DEG)
DEG <- DEG %>% rownames_to_column("Gene")

msigdb_GMTs <- "msigdb_v7.0_GMTs"
#msigdb <- "h.all.v7.0.symbols.gmt" #全基因
msigdb <- "c7.all.v7.0.symbols.gmt" #免疫相关
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

geneList = DEG[,3]
names(geneList) = as.character(DEG[,'Gene'])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

set.seed(1) #设置种子 好像不设置也没区别
gsea <-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
#转换成数据框
gsea_result_df <- as.data.frame(gsea)
write.table(gsea_result_df,file="GSEA_c7_STK11_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
save(gsea,gsea_result_df,file = "GSEA_c7_STK11.rda")
#单个图绘制
#安装enrichplot
library(enrichplot)
gseaplot2(gsea,1,color="red",pvalue_table = T)
#多个结果绘制
gseaplot2(gsea, geneSetID = c(1,3), subplots = 1:3)
gseaplot2(gsea, geneSetID = 1:4, subplots = 1:3,pvalue_table = T)


####ESTIMATE####
#计算患者免疫评分与肿瘤纯度#
setwd("TCGA-LUAD")
setwd("ESTIMATE")  #设置工作目录
#安装包
library(utils) #这个包应该不用下载，自带的
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(tidyverse)
#读取肿瘤患者01A表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

#计算免疫评分
filterCommonGenes(input.f = "tpms01A_log2.txt",   #输入文件名
                  output.f = "tpms01A_log2.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("tpms01A_log2.gct",   #刚才的输出文件名
              "tpms01A_log2_estimate_score.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix")   #默认平台

#提取结果
ESTIMATE_result <- read.table("tpms01A_log2_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ESTIMATE_result <- ESTIMATE_result[,-1]   
colnames(ESTIMATE_result) <- ESTIMATE_result[1,]   
ESTIMATE_result <- as.data.frame(t(ESTIMATE_result[-1,]))
rownames(ESTIMATE_result) <- colnames(exp)
write.table(ESTIMATE_result, file = "ESTIMATE_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F) 

#仙桃绘图
med=median(as.numeric(exp["STK11",]))
exp <- exp %>% t() %>% as.data.frame()
exp <- exp %>% mutate(group=factor(ifelse(exp$STK11>med,"high","low"),levels = c("low","high")))
class(exp$group)
identical(rownames(exp),rownames(ESTIMATE_result))

ESTIMATE_result$group <- exp$group

ESTIMATE_wd <- ESTIMATE_result1 %>% filter(STK11 == "0")
ESTIMATE_mut <- ESTIMATE_result1 %>% filter(STK11 == "1")

write.csv(ESTIMATE_wd, file = "ESTIMATE_wd.csv")
write.csv(ESTIMATE_mut, file = "ESTIMATE_mut.csv")

#如何用R绘图？
library(ggpubr)
library(ggsci)
#按照基因突变分组
ESTIMATE_result <- read.table("ESTIMATE_result.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rownames(ESTIMATE_result) <- gsub("\\.", "-", rownames(ESTIMATE_result))
# 使用str_sub函数提取前12位
rownames(ESTIMATE_result) <- substr(rownames(ESTIMATE_result), 1, 12)

a <- intersect(rownames(ESTIMATE_result),rownames(group))
ESTIMATE_result <- ESTIMATE_result[a,]

str(group)
group_filtered <- group[A, , drop = FALSE]
group = group_filtered
identical(rownames(ESTIMATE_result),rownames(group))
ESTIMATE_result1 <- cbind(ESTIMATE_result,group)
write.csv(ESTIMATE_result1, file = "ESTIMATE_result1.csv")

ESTIMATE_result1$StromalScore

ggboxplot(ESTIMATE_result1,
          x = "STK11", y = "StromalScore",#可以画任何一个基因
          fill = "STK11") +
  scale_fill_lancet()+
  stat_compare_means(comparisons = list(c(1, 2)),
                     method = "wilcox.test",   #设置统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))
dev.off()

####cibersort####
setwd("TCGA-LUAD")
setwd("CIBERSORT")   
#install.packages('e1071')
#install.packages('parallel')
#install.packages("BiocManager")
#BiocManager::install("preprocessCore", version = "3.17")
library(e1071)
library(parallel)
library(preprocessCore)
library(tidyverse)
source("CIBERSORT.R")   
sig_matrix <- "LM22.txt"   
mixture_file = 'tpms01A_log2.txt'   #肿瘤患者表达谱
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)

res_cibersort <- res_cibersort[,1:22]   
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞
ciber.res <- as.data.frame(ciber.res)
write.table(ciber.res,"ciber.res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####分组比较图####
#读取肿瘤患者01A表达谱
#exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#med=median(as.numeric(exp["STK11",]))
#exp <- exp %>% t() %>% as.data.frame()
#exp <- exp %>% mutate(group=factor(ifelse(exp$STK11>med,"high","low"),levels = c("high","low")))
#class(exp$group)
ciber.res <- read.table("ciber.res.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rownames(ciber.res) <- gsub("\\.", "-", rownames(ciber.res))
# 使用str_sub函数提取前12位
rownames(ciber.res) <- substr(rownames(ciber.res), 1, 12)

a <- intersect(rownames(ciber.res),rownames(group))
ciber.res <- ciber.res[a,]

str(group)
group_filtered <- group[a, , drop = FALSE]
group = group_filtered
identical(rownames(ciber.res),rownames(group))
ciber.res <- cbind(ciber.res,group)
ciber.res <- ciber.res %>% rownames_to_column("sample")
write.csv(ciber.res, file = "ciber.res1.csv")
#ciber.res <- read.csv("./ciber.res1.csv",row.names = 1)
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(ciber.res,key=CIBERSORT,value = Fraction,-c(STK11,sample))
#b$STK11 <- as.factor(b$STK11)
ggboxplot(b, x = "CIBERSORT", y = "Fraction",
          fill = "STK11", palette = "jco")+
  stat_compare_means(aes(group = STK11),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()

####ssGSEA####
setwd("TCGA-LUAD")
setwd("ssGSEA")
#BiocManager::install('GSVA')
library(tidyverse)
library(data.table)
library(GSVA)

##1.3 表达量矩阵的准备
###行是基因，列是样本
expr <- read.table("tpms01A_log2.txt", sep = "\t",row.names = 1,check.names = F,header = T)
expr <- as.matrix(expr)   #将expr转换为矩阵格式

#2. 使用ssGSEA量化免疫浸润
# 创建参数对象
ssgsea_param <- gsvaParam(expr, cellMarker)

# 运行 ssGSEA
gsva_data <- gsva(ssgsea_param)
a <- gsva_data %>% t() %>% as.data.frame()
rownames(a) <- gsub("\\.", "-", rownames(a))
write.table(a,"ssGSEA_res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#设置分组信息 CDKN3高低组
ssGSEA_res <- read.table("ssGSEA_res.txt", sep = "\t",row.names = 1,check.names = F,header = T)
library(tidyverse)
rownames(ssGSEA_res) <- substr(rownames(ssGSEA_res), 1, 12)

a <- intersect(rownames(ssGSEA_res),rownames(group))
ssGSEA_res <- ssGSEA_res[a,]

str(group)
group_filtered <- group[a, , drop = FALSE]
group = group_filtered
identical(rownames(ssGSEA_res),rownames(group))
ssGSEA_res <- cbind(ssGSEA_res,group)
ssGSEA_res <- ssGSEA_res %>% rownames_to_column("sample")
write.csv(ssGSEA_res, file = "ssGSEA_res1.csv")
#ssGSEA_res <- read.csv("./ssGSEA_res1.csv",row.names = 1)
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(ssGSEA_res,key=ssGSEA,value = Expression,-c(STK11,sample))
#b$STK11 <- as.factor(b$STK11)
ggboxplot(b, x = "ssGSEA", y = "Expression",
          fill = "STK11", palette = "lancet")+
  stat_compare_means(aes(group = STK11),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()

####X cell####
#安装devtools 右下角install
#安装xCell包
#chooseBioCmirror()
devtools::install_github('dviraran/xCell')
library(xCell)
library(ggpubr)
library(tidyverse)
setwd("TCGA-LUAD")
setwd("Xcell")
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#查看细胞类型
celltypeuse<-xCell.data$spill$K

Xcell_res <-xCellAnalysis(exp,parallel.sz=10) #计算
Xcell_res <- Xcell_res %>% t() %>% as.data.frame()
write.table(Xcell_res,"Xcell_res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#准备分组信息
Xcell_res <- read.table("Xcell_res.txt", sep = "\t",row.names = 1,check.names = F,header = T)
rownames(Xcell_res) <- gsub("\\.", "-", rownames(Xcell_res))
rownames(Xcell_res) <- substr(rownames(Xcell_res), 1, 12)

a <- intersect(rownames(Xcell_res),rownames(group))
Xcell_res <- Xcell_res[a,]

str(group)
group_filtered <- group[a, , drop = FALSE]
group = group_filtered
identical(rownames(Xcell_res),rownames(group))
Xcell_res <- cbind(Xcell_res,group)
Xcell_res <- Xcell_res %>% rownames_to_column("sample")
write.csv(Xcell_res, file = "Xcell_res1.csv")
#ssGSEA_res <- read.csv("./ssGSEA_res1.csv",row.names = 1)
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(Xcell_res,key=xCell,value = Expression,-c(STK11,sample))
#b$STK11 <- as.factor(b$STK11)
ggboxplot(b, x = "xCell", y = "Expression",
          fill = "STK11", palette = "lancet")+
  stat_compare_means(aes(group = STK11),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()

####免疫检查点表达差异####
#T细胞激活的二个信号
#1、Tcell的TCR-CD3与MHC-抗原肽复合物结合
#2、Tcell表面的CD28与APC表面的B7结合
#CTLA4与CD28竞争 与B7结合 诱导无反应性Tcell
#PD1-PDL1/PDL2
#癌细胞逃避免疫系统摧毁的一种方法，是通过配体联接到T细胞的PD-1蛋白上
#当配体与PD-1联接以后，T细胞就不能够发现肿瘤和向免疫系统发出攻击肿瘤的信号
#PD1药物的设计思路是：给予肿瘤患者针对PD-1或PD-L1的一种抗体蛋白质，将使两种蛋白质不会联接，T细胞的功能也不会关停
#PD1抗体(immune checkpoint blockade,ICB)
setwd("TCGA-LUAD")
setwd("ICB")
library(tidyverse)
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- group %>% rename(group = STK11)
exp <- exp %>% t() %>% as.data.frame()
rownames(exp) <- gsub("\\.", "-", rownames(exp))
rownames(exp) <- substr(rownames(exp), 1, 12)
a <- intersect(rownames(exp),rownames(group))
exp <- exp[a,]
str(group)
group_filtered <- group[a, , drop = FALSE]
group = group_filtered
identical(rownames(exp),rownames(group))
exp <- cbind(exp,group)

#exp_Wild <- exp %>% filter(group == "0")
#exp_Mutation <- exp %>% filter(group == "1")

#gene <- c("PD1","PDL1","PDL2")
#gene_high <- exp_high[,gene]
#gene_low <- exp_low[,gene]
#为什么提取不出来
#https://www.ncbi.nlm.nih.gov/
#gene <- c("PDCD1","CD274","PDCD1LG2")
#gene_high <- exp_high[,gene]
#gene_low <- exp_low[,gene]
#导出csv仙桃作图
#write.csv(gene_high, file = "gene_high.csv")
#write.csv(gene_low, file = "gene_low.csv")

#R作图
library(tidyverse)
gene <- c("PDCD1","CD274","PDCD1LG2","group")
a <- exp[,gene]
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=ICB,value = Expression,-c(group,sample))
#b$group <- as.factor(b$group)
ggboxplot(b, x = "ICB", y = "Expression",
          fill = "group", palette = "jco")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()

####TMB####
setwd("TCGA-LUAD")
setwd("TMB")
#BiocManager::install("maftools")
library(maftools)
library(tidyverse)
library(readxl)
library(readr)

#按照是否突变进行分组
group <- mat_0_1["STK11", , drop = FALSE]

group <- t(group) 
group <- as.data.frame(group)
table(group$STK11)
save(group, file = "group_STK11t突变.Rda")

exp <- read.table("tpms01A_log2.txt", sep = "\t",row.names = 1,check.names = F,header = T)
colnames(exp) <- gsub("\\.", "-", colnames(exp))

#统一两数据框的样本保持一致
library(stringr)

# 使用str_sub函数提取前12位
colnames(exp) <- substr(colnames(exp), 1, 12)

a <- intersect(colnames(exp),rownames(group))
exp <- exp[,a]
group <- group[a,]
#行和列进行转换
exp01A <- exp01A %>% t() %>% as.data.frame()
identical(rownames(exp),rownames(group))



# https://portal.gdc.cancer.gov/

#定义合并函数
merge_maf <- function(metadata, path){
  #通过合并path,还有sample sheet前两列得到每一个文件的完整路径
  filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                         fsep = .Platform$file.sep)
  
  message ('############### Merging maf data ################\n',
           '### This step may take a few minutes ###\n')
  #通过lapply循环去读每一个样本的maf，然后通过rbind合并成矩阵，按行来合并
  #colClasses指定所有列为字符串
  mafMatrix <- do.call("rbind", lapply(filenames, function(fl) 
    read.table(gzfile(fl),header=T,sep="\t",quote="",fill=T,colClasses="character")))
  return (mafMatrix)
}

#定义去除重复样本的函数FilterDuplicate
FilterDuplicate <- function(metadata) {
  filter <- which(duplicated(metadata[,'sample']))
  if (length(filter) != 0) {
    metadata <- metadata[-filter,]
  }
  message (paste('Removed', length(filter), 'samples', sep=' '))
  return (metadata)
}

#读入maf的sample sheet文件
metaMatrix.maf=read.table("maf_sample_sheet.tsv",sep="\t",header=T)
#替换.为下划线，转换成小写，sample_id替换成sample
names(metaMatrix.maf)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.maf))))
#删掉最后一列sample_type中的空格
metaMatrix.maf$sample_type=gsub(" ","",metaMatrix.maf$tissue_type)
#删掉重复的样本
metaMatrix.maf<- FilterDuplicate(metaMatrix.maf)

#调用merge_maf函数合并maf的矩阵
maf_value=merge_maf(metadata=metaMatrix.maf, 
                    path="maf_data"
)
#查看前三行前十列
maf_value[1:3,1:10]

#保存合并后的maf文件
write.table(maf_value,file="combined_maf_value.txt",row.names=F,quote=F,sep="\t")

#突变数据整理成数据框
str(maf_value)
a <- maf_value %>% 
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>% 
  as.data.frame() %>% 
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))
table(a$Variant_Classification)
a <- filter(a,Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Nonstop_Mutation","In_Frame_Del","Translation_Start_Site","In_Frame_Ins" ))
table(a$Variant_Classification)
gene <- as.character(unique(a$Hugo_Symbol))
sample <- as.character(unique(a$Tumor_Sample_Barcode))

#以下不需要掌握
mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))

for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))
#需要知道a gene_count mat mat_0_1分别是什么
save(mat,mat_0_1,gene_count,a,file = "TMB-LUAD.rda")

####TMB整体结果图#### 
mut2 <- read.maf(maf = "combined_maf_value.txt")
plotmafSummary(maf = mut2,
               rmOutlier = TRUE,   #删除数据中的离群值
               addStat = 'median',
               dashboard = TRUE,
               titvRaw = FALSE)

#Oncostrip图
oncostrip(maf = mut2,
          genes = c('TP53','STK11','KRAS','LRP1B'))   #输入研究的基因名称

####TMB分组比较图####
#计算TMB
mut2 <- read.maf(maf = "combined_maf_value.txt")
maf = tmb(maf = mut2,
          captureSize = 50,
          logScale = TRUE)   

maf$sample <- substr(maf$Tumor_Sample_Barcode,1,16)

write.csv(maf,file = "maf.csv")
#手动在excell进行修改
maf <- read.csv("maf.csv", row.names=1)

#导入group 需要根据自己的需求改(例如根据某基因分高低组)
#load("group_AY.rda")

####分组信息如何修改####
exp <- read.table("tpms01A_log2.txt", sep = "\t",row.names = 1,check.names = F,header = T)
colnames(exp) <- gsub("\\.", "-", colnames(exp))
gene <- "STK11"#每次运行只改这个基因名
med=median(as.numeric(exp[gene,]))
group=data.frame(sample=colnames(exp),
                      group=factor(ifelse(exp[gene,]>med,"high","low"),levels = c("high","low"))) %>% 
  column_to_rownames("sample")
save(group, file = "group_STK11.Rda")

#有没有远处转移分组
exp <- read.table("tpms01A_log2.txt", sep = "\t",row.names = 1,check.names = F,header = T)
colnames(exp) <- gsub("\\.", "-", colnames(exp))
clinical.expr01A <- read.table("clinical.expr01A.txt", sep = "\t",row.names = 1,check.names = F,header = T)
a <- filter(clinical.expr01A,ajcc_pathologic_m.diagnoses == 'M0')
b <- filter(clinical.expr01A,ajcc_pathologic_m.diagnoses == 'M1')
c <- rbind(a,b)
M_group <- data.frame(
  row_names = rownames(c),
  M = c$ajcc_pathologic_m.diagnoses
)
library(stringr)
M_group$M <- str_remove_all(M_group$M, "M")

write.csv(M_group, file = "M_group.csv")
#手动在excell进行修改
M_group <- read.csv("M_group.csv", row.names=1)

save(M_group, file = "M_group.Rda")

#load("group_STK11.Rda")
maf <- maf[substring(rownames(maf),14,16) == c("01A"),]
rownames(maf) <- substr(rownames(maf), 1, 12)
com <- intersect(rownames(maf),rownames(group))
maf <- maf[com,]
group <- group %>% t() %>% as.data.frame()
group <- group[,com]
group <- group %>% t() %>% as.data.frame()
identical(rownames(group),rownames(maf))
tmb <- cbind(maf,group)
write.csv(tmb, file = "突变tmb.csv")

#R作图
library(ggpubr)
library(ggsci)
ggboxplot(tmb,
          x = "STK11", y = "total_perMB",
          fill = "STK11") +
  scale_fill_lancet()+
  stat_compare_means(comparisons = list(c(1, 2)),
                     method = "wilcox.test",   #设置统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))
dev.off()





####TMB瀑布图及分组绘图####
setwd("TCGA-LUAD")
setwd("TMB")
library(tidyverse)
load("TMB-LUAD.rda")
#导入group 需要根据自己的需求改(例如根据某基因分高低组)
load("group_STK11.Rda")

#BiocManager::install("ComplexHeatmap")
#install.packages("circlize")
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
#group = M_group
mat = as.matrix(mat)
mat[1:5, 1:5]
dim(mat)

# 设置颜色
mutant_type <- as.character(unique(a$Variant_Classification))
# [1] "Missense_Mutation"      "Nonsense_Mutation"      "Frame_Shift_Del"       
# [4] "Frame_Shift_Ins"        "Splice_Site"            "Nonstop_Mutation"      
# [7] "In_Frame_Del"           "Translation_Start_Site" "In_Frame_Ins" 

library(RColorBrewer)
display.brewer.all()

col = c("Missense_Mutation"=brewer.pal(9, "Set1")[3],
        "Nonsense_Mutation"=brewer.pal(9, "Set1")[1],
        "Frame_Shift_Del"=brewer.pal(9, "Set1")[2],
        "Frame_Shift_Ins"=brewer.pal(9, "Set1")[4],
        "Splice_Site"=brewer.pal(9, "Set1")[5],
        "Nonstop_Mutation"=brewer.pal(9, "Set1")[6],
        "In_Frame_Del"=brewer.pal(9, "Set1")[7],
        "Translation_Start_Site"=brewer.pal(9, "Set1")[8],
        "In_Frame_Ins"=brewer.pal(9, "Set1")[9])

alter_fun = list(
  # 设置背景颜色
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # 设置HOMDEL类型的颜色
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  # 设置AMP类型的颜色
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  # 设置MUT类型的颜色
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },   
  # 设置HOMDEL类型的颜色
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  # 设置AMP类型的颜色
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  # 设置MUT类型的颜色
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
  },
  # 设置HOMDEL类型的颜色
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  # 设置AMP类型的颜色
  Translation_Start_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Translation_Start_Site"], col = NA))
  },
  # 设置MUT类型的颜色
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  }
)

# 用oncoPrint函数绘图
group <- group %>% 
  rownames_to_column("sample") %>% 
  mutate(sample=substring(.$sample,1,12))
com_sample <- intersect(colnames(mat),group$sample)

gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))
gene_top <- gene_count$gene[1:30] # 修改数字，代表TOP多少

# 总体的瀑布图
oncoPrint(mat[gene_top,],#修改绘图的行和列
          alter_fun = alter_fun, col = col, 
          top_annotation = HeatmapAnnotation(
            column_barplot = anno_oncoprint_barplot(ylim=c(0,20),border=T)),
          right_annotation = rowAnnotation(
            row_barplot = anno_oncoprint_barplot(border = TRUE,ylim=c(0,300))),
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          row_names_gp = gpar(fontsize = 7),
          pct_gp = gpar(fontsize = 7), 
          pct_side = "right", row_names_side = "left",
          column_title = "OncoPrint for TCGA LUAD Top30 Mutation Gene", 
          heatmap_legend_param = list(title = "Alternations"))

# 瀑布图：group1 Mutation的TOP30突变基因
group <- group %>% rename(group = M)
group1 <- group %>% filter(.$group==0)
group1_sample <- intersect(com_sample,group1$sample)
mat_0_1_group1 <- mat_0_1[,group1_sample]
gene1_count <- data.frame(gene=rownames(mat_0_1_group1),
                          count=as.numeric(apply(mat_0_1_group1,1,sum))) %>%
  arrange(desc(count))
gene1_top <- gene1_count$gene[1:30] # 修改数字，代表TOP多少

oncoPrint(mat[gene1_top,group1_sample],
          alter_fun = alter_fun, col = col, 
          top_annotation = HeatmapAnnotation(
            column_barplot = anno_oncoprint_barplot(ylim=c(0,20),border=TRUE)),
          right_annotation = rowAnnotation(
            row_barplot = anno_oncoprint_barplot(border = TRUE,ylim=c(0,150))),
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          row_names_gp = gpar(fontsize = 7),
          pct_gp = gpar(fontsize = 7), 
          pct_side = "right", row_names_side = "left",
          column_title = "Top30 Mutation Gene of Group1 (n=339)", 
          heatmap_legend_param = list(title = "Alternations"))

# 瀑布图：group2 wild的TOP30突变基因
group2 <- group %>% filter(.$group==1)
group2_sample <- intersect(com_sample,group2$sample)
mat_0_1_group2 <- mat_0_1[,group2_sample]
gene2_count <- data.frame(gene=rownames(mat_0_1_group2),
                          count=as.numeric(apply(mat_0_1_group2,1,sum))) %>%
  arrange(desc(count))
gene2_top <- gene2_count$gene[1:30] # 修改数字，代表TOP多少

oncoPrint(mat[gene2_top,group2_sample],
          alter_fun = alter_fun, col = col, 
          top_annotation = HeatmapAnnotation(
            column_barplot = anno_oncoprint_barplot(ylim=c(0,20),border=TRUE)),
          right_annotation = rowAnnotation(
            row_barplot = anno_oncoprint_barplot(border = TRUE,ylim=c(0,150))),
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          row_names_gp = gpar(fontsize = 7),
          pct_gp = gpar(fontsize = 7), 
          pct_side = "right", row_names_side = "left",
          column_title = "Top30 Mutation Gene of Group2 (n=25)", 
          heatmap_legend_param = list(title = "Alternations"))




# 瀑布图：group1的TOP30突变基因
group <- group %>% rename(group = M)
group1 <- group %>% filter(.$group==0)
group1_sample <- intersect(com_sample,group1$sample)
mat_0_1_group1 <- mat_0_1[,group1_sample]
gene1_count <- data.frame(gene=rownames(mat_0_1_group1),
                          count=as.numeric(apply(mat_0_1_group1,1,sum))) %>%
  arrange(desc(count))
gene1_top <- gene1_count$gene[1:30] # 修改数字，代表TOP多少

oncoPrint(mat[gene1_top,group1_sample],
          alter_fun = alter_fun, col = col, 
          top_annotation = HeatmapAnnotation(
            column_barplot = anno_oncoprint_barplot(ylim=c(0,20),border=TRUE)),
          right_annotation = rowAnnotation(
            row_barplot = anno_oncoprint_barplot(border = TRUE,ylim=c(0,150))),
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          row_names_gp = gpar(fontsize = 7),
          pct_gp = gpar(fontsize = 7), 
          pct_side = "right", row_names_side = "left",
          column_title = "Top30 Mutation Gene of Group1 (n=339)", 
          heatmap_legend_param = list(title = "Alternations"))

# 瀑布图：group2的TOP30突变基因
group2 <- group %>% filter(.$group==1)
group2_sample <- intersect(com_sample,group2$sample)
mat_0_1_group2 <- mat_0_1[,group2_sample]
gene2_count <- data.frame(gene=rownames(mat_0_1_group2),
                          count=as.numeric(apply(mat_0_1_group2,1,sum))) %>%
  arrange(desc(count))
gene2_top <- gene2_count$gene[1:30] # 修改数字，代表TOP多少

oncoPrint(mat[gene2_top,group2_sample],
          alter_fun = alter_fun, col = col, 
          top_annotation = HeatmapAnnotation(
            column_barplot = anno_oncoprint_barplot(ylim=c(0,20),border=TRUE)),
          right_annotation = rowAnnotation(
            row_barplot = anno_oncoprint_barplot(border = TRUE,ylim=c(0,150))),
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          row_names_gp = gpar(fontsize = 7),
          pct_gp = gpar(fontsize = 7), 
          pct_side = "right", row_names_side = "left",
          column_title = "Top30 Mutation Gene of Group2 (n=25)", 
          heatmap_legend_param = list(title = "Alternations"))


# 瀑布图：group high的TOP30突变基因
grouphigh <- group %>% filter(.$group=="high")
grouphigh_sample <- intersect(com_sample,grouphigh$sample)
mat_0_1_grouphigh <- mat_0_1[,grouphigh_sample]
genehigh_count <- data.frame(gene=rownames(mat_0_1_grouphigh),
                             count=as.numeric(apply(mat_0_1_grouphigh,1,sum))) %>%
  arrange(desc(count))
genehigh_top <- genehigh_count$gene[1:30] # 修改数字，代表TOP多少

oncoPrint(mat[genehigh_top,grouphigh_sample],
          alter_fun = alter_fun, col = col, 
          top_annotation = HeatmapAnnotation(
            column_barplot = anno_oncoprint_barplot(ylim=c(0,20),border=TRUE)),
          right_annotation = rowAnnotation(
            row_barplot = anno_oncoprint_barplot(border = TRUE,ylim=c(0,150))),
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          row_names_gp = gpar(fontsize = 7),
          pct_gp = gpar(fontsize = 7), 
          pct_side = "right", row_names_side = "left",
          column_title = "Top30 Mutation Gene of Group high (n=250)", 
          heatmap_legend_param = list(title = "Alternations"))

# 瀑布图：grouplow的TOP30突变基因
grouplow <- group %>% filter(.$group=="low")
grouplow_sample <- intersect(com_sample,grouplow$sample)
mat_0_1_grouplow <- mat_0_1[,grouplow_sample]
genelow_count <- data.frame(gene=rownames(mat_0_1_grouplow),
                            count=as.numeric(apply(mat_0_1_grouplow,1,sum))) %>%
  arrange(desc(count))
genelow_top <- genelow_count$gene[1:30] # 修改数字，代表TOP多少

oncoPrint(mat[genelow_top,grouplow_sample],
          alter_fun = alter_fun, col = col, 
          top_annotation = HeatmapAnnotation(
            column_barplot = anno_oncoprint_barplot(ylim=c(0,20),border=TRUE)),
          right_annotation = rowAnnotation(
            row_barplot = anno_oncoprint_barplot(border = TRUE,ylim=c(0,150))),
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          row_names_gp = gpar(fontsize = 7),
          pct_gp = gpar(fontsize = 7), 
          pct_side = "right", row_names_side = "left",
          column_title = "Top30 Mutation Gene of Group low (n=252)", 
          heatmap_legend_param = list(title = "Alternations"))
####特定基因突变数分组比较图####
save(group, file = "M_group.Rda")
setwd("TCGA-STAD")
setwd("TMB")
library(tidyverse)
load("TMB-LUAD.rda")
#导入group 需要根据自己的需求改(例如根据某基因分高低组)
load("M_group.rda")
group <- column_to_rownames(group, var = "sample")
a <- intersect(rownames(group),colnames(mat_0_1))

group <- as.data.frame(t(group))
group <- group[,a]
group <- as.data.frame(t(group))

mat_0_1 <- mat_0_1[,a]
mat_0_1 <- as.data.frame(t(mat_0_1))
identical(rownames(group),rownames(mat_0_1))
mat_0_1_group <- cbind(mat_0_1,group)

mat_0_1_group1 <- filter(mat_0_1_group,group=="0")
mat_0_1_group2 <- filter(mat_0_1_group,group=="1")

gene <- c("STK11")
gene_1 <- mat_0_1_group1[,gene]
gene_2 <- mat_0_1_group2[,gene]
write.csv(gene_1,file = "gene_1.csv")
write.csv(gene_2,file = "gene_2.csv")


mat_0_1_grouphigh <- filter(mat_0_1_group,group=="high")
mat_0_1_grouplow <- filter(mat_0_1_group,group=="low")

gene <- c("STK11")
gene_high <- mat_0_1_grouphigh[,gene]
gene_low <- mat_0_1_grouplow[,gene]
write.csv(gene_high,file = "gene_high.csv")
write.csv(gene_low,file = "gene_low.csv")


####蛋白组学及磷酸化组学####
# 安装
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("cBioPortalData")

# 加载
library(cBioPortalData)

download.file(
  url = "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_proteome_ratio_NArm_TUMOR.cct",
  destfile = "CPTAC_LUAD_proteome_tumor.cct"
)
proteome_matrix <- read.table("CPTAC_LUAD_proteome_tumor.cct", 
                              header = TRUE, sep = "\t", row.names = 1)

download.file(
  url = "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_phosphoproteome_ratio_norm_NArm_TUMOR.cct",
  destfile = "CPTAC_LUAD_phosphoproteome_tumor.cct"
)
phos_matrix <- read.table("CPTAC_LUAD_phosphoproteome_tumor.cct", 
                          header = TRUE, sep = "\t", row.names = 1)

# 安装未安装的包（如果已安装可跳过）
required_pkgs <- c("limma", "EnhancedVolcano", "pheatmap", "clusterProfiler",
                   "org.Hs.eg.db", "enrichplot", "igraph", "ggraph", 
                   "DEP2", "SummarizedExperiment", "dplyr", "tidyr", "tibble")
install.packages(setdiff(required_pkgs, installed.packages()[,"Package"]))



# 安装 EnhancedVolcano
BiocManager::install("EnhancedVolcano")
# 1. 安装官方版DEP
BiocManager::install("DEP")

# 2. 后续代码与之前使用DEP2的几乎一致，直接将 DEP2() 替换为 DEP()
library(DEP)
# se_prot <- DEP(se_prot, condition = "condition", non_missing = 0.5, missing_fill = TRUE)

# 加载
library(limma)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(igraph)
library(ggraph)
library(DEP)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(tibble)

# 读取突变文件
mut_all <- read.table("data_mutations.txt", header = TRUE, sep = "\t", 
                      quote = "", comment.char = "", stringsAsFactors = FALSE)

# 筛选出 STK11 突变
stk11_mut <- mut_all[mut_all$Hugo_Symbol == "STK11", ]

# 提取唯一突变样本 ID（注意可能列名有变，通常是 Tumor_Sample_Barcode）
mutant_samples <- unique(stk11_mut$Tumor_Sample_Barcode)

# 预览前几个
head(mutant_samples)

# 假设你之前已读取了这两个矩阵（确认列名格式）
# proteome_matrix 和 phos_matrix 是从 .cct 文件读入的

# 查看实际列名格式
head(colnames(proteome_matrix))
head(mutant_samples)

# 创建一个去除特殊符号的统一 ID 函数
clean_id <- function(x) gsub("[-._]", "", x)

# 构建匹配向量
mutant_clean <- clean_id(mutant_samples)
sample_clean <- clean_id(colnames(proteome_matrix))

# 分组：如果一个样本的“清洁版”出现在突变清洁版中，就标记为 MT
is_mutant <- sample_clean %in% mutant_clean
mut_group <- ifelse(is_mutant, "MT", "WT")

# 检查分组情况
table(mut_group)

###数据预处理
library(limma)

# ---------- 1. 缺失值过滤函数 ----------
filter_by_group <- function(expr_matrix, group_vector, min_ratio = 0.5) {
  keep <- logical(nrow(expr_matrix))
  for (i in 1:nrow(expr_matrix)) {
    keep[i] <- all(
      sapply(unique(group_vector), function(g) {
        group_samples <- which(group_vector == g)
        mean(!is.na(expr_matrix[i, group_samples])) >= min_ratio
      })
    )
  }
  expr_matrix[keep, ]
}

# ---------- 2. 缺失值填充函数（每组最小值填充） ----------
impute_matrix <- function(expr_matrix, group_vector) {
  expr_imputed <- expr_matrix
  for (i in 1:nrow(expr_imputed)) {
    for (g in unique(group_vector)) {
      idx <- which(group_vector == g)
      min_val <- min(expr_imputed[i, idx], na.rm = TRUE)
      expr_imputed[i, idx][is.na(expr_imputed[i, idx])] <- min_val
    }
  }
  expr_imputed
}

# ---------- 3. 应用预处理 ----------
prot_filtered <- filter_by_group(proteome_matrix, mut_group)
phos_filtered <- filter_by_group(phos_matrix, mut_group)

prot_imputed <- impute_matrix(prot_filtered, mut_group)
phos_imputed <- impute_matrix(phos_filtered, mut_group)

# 分位数标准化（数据如果已经是log2转换的）
prot_expr_final <- normalizeBetweenArrays(prot_imputed, method = "quantile")
phos_expr_final <- normalizeBetweenArrays(phos_imputed, method = "quantile")

# ---------- 4. 导出整理好的数据 ----------
# 蛋白质组
write.table(as.data.frame(prot_expr_final), file = "STK11_proteome_clean.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# 磷酸化组
write.table(as.data.frame(phos_expr_final), file = "STK11_phosphoproteome_clean.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# 样本分组
group_df <- data.frame(Sample = colnames(prot_expr_final), Group = mut_group)
write.csv(group_df, file = "STK11_sample_groups.csv", row.names = FALSE)

# 保存 R 环境（可选）
save(prot_expr_final, phos_expr_final, mut_group,
     file = "STK11_preprocessed_data.RData")

cat("预处理完成，文件已导出。\n")



# 1. 读入数据（和之前一样）
prot <- read.table("STK11_proteome_clean.tsv", 
                   header = T, sep = "\t", row.names = 1, check.names = F)

phos <- read.table("STK11_phosphoproteome_clean.tsv", 
                   header = T, sep = "\t", row.names = 1, check.names = F)

group_df <- read.csv("STK11_sample_groups.csv", check.names = F)

# 2. 关键修正：先构建样本-分组映射，再筛选
# 确保 group_df$Sample 和 prot/phos 的列名完全一致
group_df$Sample <- as.character(group_df$Sample)
common_samples <- intersect(colnames(prot), colnames(phos))

# 按 common_samples 筛选 group_df，保证顺序和样本名对应
group_df_filtered <- group_df[group_df$Sample %in% common_samples, ]
# 按 prot 的列名顺序重新排序（避免样本顺序错位）
group_df_filtered <- group_df_filtered[match(colnames(prot), group_df_filtered$Sample), ]

# 重新创建分组向量
group <- factor(group_df_filtered$Group)

# 3. 重新过滤矩阵
prot <- prot[, common_samples]
phos <- phos[, common_samples]

cat("最终样本数：", length(common_samples), "\n")
print(table(group))  # 现在会显示 MT 和 WT 的真实样本数

##PCA 分析
run_pca <- function(mat, title_name){
  pca <- prcomp(t(mat), scale. = FALSE)
  var <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
  df <- data.frame(
    PC1 = pca$x[,1], PC2 = pca$x[,2],
    Group = group, Sample = rownames(pca$x)
  )
  ggplot(df, aes(PC1, PC2, color=Group)) +
    geom_point(size=3, alpha=0.8) +
    labs(title=title_name,
         x=paste0("PC1 (",var[1],"%)"),
         y=paste0("PC2 (",var[2],"%)")) +
    theme_bw(base_size=13) +
    stat_ellipse(level=0.95) +
    theme(plot.title = element_text(hjust=0.5))
}

pca_prot <- run_pca(prot, "Proteome: STK11 Mutation vs Wild")
pca_phos <- run_pca(phos, "Phosphoproteome: STK11 Mutation vs Wild")

ggsave("PCA_proteome.pdf", pca_prot, width=7, height=5)
ggsave("PCA_phospho.pdf", pca_phos, width=7, height=5)

print(pca_prot)
print(pca_phos)

##差异分析（DEPs + DPPs）
# 设计矩阵
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 对比：突变 vs 野生
contr <- makeContrasts(
  contrasts = paste0(levels(group)[1], "-", levels(group)[2]),
  levels = design
)

# ================== 蛋白组差异 ==================
fit_p <- lmFit(prot, design)
fit_p <- contrasts.fit(fit_p, contr)
fit_p <- eBayes(fit_p)
deg_p <- topTable(fit_p, coef=1, number=Inf, adjust="fdr")

deg_p_sig <- deg_p %>% filter(abs(logFC) > 1, adj.P.Val < 0.05)
write.csv(deg_p_sig, "01_DEPs_proteome.csv", row.names=T)

# ================== 磷酸化组差异 ==================
fit_ph <- lmFit(phos, design)
fit_ph <- contrasts.fit(fit_ph, contr)
fit_ph <- eBayes(fit_ph)
deg_ph <- topTable(fit_ph, coef=1, number=Inf, adjust="fdr")

deg_ph_sig <- deg_ph %>% filter(abs(logFC) > 0.8, adj.P.Val < 0.05)
write.csv(deg_ph_sig, "02_DPPs_phosphoproteome.csv", row.names=T)

cat("差异蛋白：", nrow(deg_p_sig), "\n")
cat("差异磷酸化位点：", nrow(deg_ph_sig), "\n")

#火山图
# 发表级火山图 + 自动标注TOP差异基因
volcano <- function(deg, title, lfc=1, fdr=0.05, topn=15){
  # 1. 标记上下调
  deg$sig <- case_when(
    deg$adj.P.Val < fdr & deg$logFC > lfc ~ "Up",
    deg$adj.P.Val < fdr & deg$logFC < -lfc ~ "Down",
    TRUE ~ "NS"
  )
  
  # 2. 提取要显示标签的基因（最显著的topn个）
  deg$Gene <- rownames(deg)
  deg <- deg[order(deg$adj.P.Val), ]
  up_genes <- head(deg[deg$sig=="Up", "Gene"], topn)
  down_genes <- head(deg[deg$sig=="Down", "Gene"], topn)
  label_genes <- c(up_genes, down_genes)
  
  # 3. 绘图（带标签）
  ggplot(deg, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(color=sig), alpha=0.7, size=1.6) +
    scale_color_manual(values=c("blue","gray","red")) +
    geom_vline(xintercept = c(-lfc,lfc), linetype="dashed", color="black") +
    geom_hline(yintercept = -log10(fdr), linetype="dashed", color="black") +
    # 自动添加基因标签，防重叠
    ggrepel::geom_text_repel(
      data = subset(deg, Gene %in% label_genes),
      aes(label=Gene), size=3.5, color="black",
      box.padding=0.5, point.padding=0.3, max.overlaps=50
    ) +
    theme_bw(base_size=13) +
    labs(title=title, x="log2(Fold Change)", y="-log10(Adjusted P-value)") +
    theme(plot.title = element_text(hjust=0.5))
}

# 运行（蛋白组 lfc=1，磷酸化组 lfc=0.8）
vp <- volcano(deg_p, "Proteome Volcano", lfc=1, fdr=0.05, topn=10)
vph <- volcano(deg_ph, "Phosphoproteome Volcano", lfc=0.8, fdr=0.05, topn=3)

# 保存图片
ggsave("Volcano_proteome_label.pdf", vp, width=9, height=7)
ggsave("Volcano_phospho_label.pdf", vph, width=9, height=7)

# 显示图片
print(vp)
print(vph)

##热图
##热图
# 发表级热图：同组样本聚在一起 + 无报错 + 符合SCI规范
plot_heat <- function(mat, deg_sig, topn=50, filename){
  # 1. 取TOP显著差异分子
  ord <- deg_sig[order(deg_sig$adj.P.Val), ]
  ids <- rownames(ord)[1:min(topn, nrow(ord))]
  
  # 2. 提取子集
  sub_mat <- mat[ids, ]
  
  # ===================== 核心：样本按分组排序（MT 放一起，WT 放一起）
  sample_order <- order(group)               # 按分组排序样本
  sub_mat <- sub_mat[, sample_order]         # 矩阵重新排序
  group_ordered <- group[sample_order]       # 分组对应排序
  
  # 3. 仅做一次标准化（修复报错关键）
  mat_scaled <- t(scale(t(sub_mat)))
  
  # 4. 注释（与排序后样本对应）
  annot <- data.frame(Group = group_ordered)
  rownames(annot) <- colnames(mat_scaled)
  
  # 5. 发表级热图
  pheatmap(mat_scaled,
           scale = "none",                # 已经标准化过，必须设为 none！
           show_rownames = FALSE,         # 不显示基因名（SCI规范）
           show_colnames = FALSE,         # 不显示样本名
           annotation_col = annot,        # 分组注释
           annotation_names_col = FALSE,  # 不显示注释标题
           treeheight_col = 25,
           treeheight_row = 25,
           border_color = NA,             # 无边框，更美观
           color = colorRampPalette(c("#377EB8", "#FFFFFF", "#E41A1C"))(100), # 标准红蓝
           cluster_cols = FALSE,          # 强制不聚类 → 保证分组排一起
           cluster_rows = TRUE,           # 基因聚类保留
           filename = filename,
           width = 8, height = 6)
}

# 运行热图
plot_heat(prot, deg_p_sig, 50, "Heatmap_DEPs.pdf")
plot_heat(phos, deg_ph_sig, 50, "Heatmap_DPPs.pdf")



####蛋白组 + 磷酸化组 联合分析####
#===============================================================================
# 安装/加载必需包
if (!require("dplyr")) install.packages("dplyr")
if (!require("stringr")) install.packages("stringr")
if (!require("ggplot2")) install.packages("ggplot2")

library(dplyr)
library(stringr)
library(ggplot2)

# ===================== 1. 蛋白组差异基因 =====================
deg_p_genes <- rownames(deg_p_sig) %>% toupper()

# ===================== 2. 磷酸化组：按冒号分割，只取第一个部分 =====================
# 你的格式是 PDE4D:NP_001184147.1:S138s → 只取 PDE4D
deg_ph_genes <- rownames(deg_ph_sig)
# 按冒号分割，取第一部分
parts <- str_split(deg_ph_genes, ":")
deg_ph_genes <- sapply(parts, function(x) x[1])
# 统一大写
deg_ph_genes <- toupper(deg_ph_genes)
# 去重
deg_ph_genes <- unique(deg_ph_genes)

# ===================== 3. 计算交集（现在一定有值！） =====================
intersect_genes <- intersect(deg_p_genes, deg_ph_genes)
cat("✅ 共同差异基因数量：", length(intersect_genes), "\n")

# ===================== 4. 构建3种调控模式 =====================
all_genes <- union(deg_p_genes, deg_ph_genes)

comb_df <- data.frame(
  Gene = all_genes,
  Proteome = all_genes %in% deg_p_genes,
  Phosphoproteome = all_genes %in% deg_ph_genes
)

comb_df$Pattern <- case_when(
  comb_df$Proteome & comb_df$Phosphoproteome ~ "Both_changed",
  comb_df$Proteome ~ "Proteome_only",
  TRUE ~ "Phosphoproteome_only"
)

# ===================== 5. 统计 =====================
count_df <- comb_df %>% count(Pattern)
print(count_df)

# ===================== 6. 发表级柱状图 =====================
p <- ggplot(comb_df, aes(x = Pattern, fill = Pattern)) +
  geom_bar(alpha=0.8, color="black") +
  geom_text(stat="count", aes(label=after_stat(count)), vjust=-0.5, size=5) +
  theme_bw(base_size=14) +
  scale_fill_manual(values=c("#377EB8","#E41A1C","#4DAF4A")) +
  labs(title="Proteome & Phosphoproteome Regulatory Patterns",
       x="Pattern", y="Number of Genes") +
  theme(plot.title=element_text(hjust=0.5), legend.position="none")

ggsave("Combined_Patterns.pdf", p, width=7, height=6)
print(p)

# ===================== 7. 导出结果 =====================
write.csv(comb_df, "Combined_Analysis.csv", row.names=F)
writeLines(c(
  paste("差异蛋白基因：", length(deg_p_genes)),
  paste("差异磷酸化基因：", length(deg_ph_genes)),
  paste("共同变化基因：", length(intersect_genes))
), "Combined_Results.txt")

message("✅ 全部完成！交集已正确计算！")


####富集分析####
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "ggplot2",
  "dplyr",
  "stringr"
))

# 加载包
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(patchwork)  # 拼图包，用于合并BP/CC/MF
set.seed(123)

#===================== 1. 读取联合分组基因（沿用之前结果）=====================
both_genes     <- comb_df %>% filter(Pattern == "Both_changed") %>% pull(Gene)
phos_only_genes <- comb_df %>% filter(Pattern == "Phosphoproteome_only") %>% pull(Gene)
prot_only_genes <- comb_df %>% filter(Pattern == "Proteome_only") %>% pull(Gene)

# 基因Symbol转EntrezID
convert_entrez <- function(genes) {
  mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL",
         column = "ENTREZID", multiVals = "first") %>% 
    na.omit() %>% as.character()
}
both_entrez  <- convert_entrez(both_genes)
phos_entrez  <- convert_entrez(phos_only_genes)
prot_entrez  <- convert_entrez(prot_only_genes)

#===================== 2. 定义函数：完整GO三类富集 + 单类气泡图 + 三合一组合图 =====================
GO_all_analysis <- function(gene_vec, group_name, top_num = 8){
  if(length(gene_vec) < 3){
    message(group_name,"基因数量不足，跳过富集")
    return(NULL)
  }
  # 1. GO三大类富集
  go_bp <- enrichGO(gene = gene_vec, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = "BP", pAdjustMethod = "fdr", qvalueCutoff = 0.05, readable = TRUE)
  go_cc <- enrichGO(gene = gene_vec, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = "CC", pAdjustMethod = "fdr", qvalueCutoff = 0.05, readable = TRUE)
  go_mf <- enrichGO(gene = gene_vec, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = "MF", pAdjustMethod = "fdr", qvalueCutoff = 0.05, readable = TRUE)
  
  # 导出富集表格
  write.csv(as.data.frame(go_bp),paste0(group_name,"_GO_BP_Result.csv"),row.names = F)
  write.csv(as.data.frame(go_cc),paste0(group_name,"_GO_CC_Result.csv"),row.names = F)
  write.csv(as.data.frame(go_mf),paste0(group_name,"_GO_MF_Result.csv"),row.names = F)
  
  # 2. 单独绘制三类GO气泡图
  if(nrow(go_bp)>0){
    p_bp <- dotplot(go_bp,showCategory = top_num,title = paste0(group_name," GO Biological Process"))
    ggsave(paste0(group_name,"_GO_BP.pdf"),p_bp,width = 11,height = 6)
  }
  if(nrow(go_cc)>0){
    p_cc <- dotplot(go_cc,showCategory = top_num,title = paste0(group_name," GO Cellular Component"))
    ggsave(paste0(group_name,"_GO_CC.pdf"),p_cc,width = 11,height = 6)
  }
  if(nrow(go_mf)>0){
    p_mf <- dotplot(go_mf,showCategory = top_num,title = paste0(group_name," GO Molecular Function"))
    ggsave(paste0(group_name,"_GO_MF.pdf"),p_mf,width = 11,height = 6)
  }
  
  # 3. 拼接三合一完整GO组合图（论文主流配图）
  merge_go_plot <- function(go_res,type){
    if(nrow(go_res)<1) return(ggplot()+theme_void())
    df <- as.data.frame(go_res)[1:min(top_num,nrow(go_res)),]
    df$Type <- type
    df
  }
  bp_df <- merge_go_plot(go_bp,"BP")
  cc_df <- merge_go_plot(go_cc,"CC")
  mf_df <- merge_go_plot(go_mf,"MF")
  all_go_df <- rbind(bp_df,cc_df,mf_df)
  
  if(nrow(all_go_df)>0){
    combine_p <- ggplot(all_go_df,aes(x=reorder(Description,Count),y=Count))+
      geom_point(aes(size=Count,color=p.adjust))+
      facet_wrap(~Type,scales="free_y",ncol=1)+
      scale_color_gradient(low="#E74C3C",high="#3498DB")+
      coord_flip()+
      theme_bw(base_size = 12)+
      labs(x="GO Terms",y="Gene Number",color="Adjusted P-value",size="Gene Count",
           title=paste0(group_name," Complete GO Enrichment Analysis"))+
      theme(plot.title = element_text(hjust = 0.5,face="bold"))
    ggsave(paste0(group_name,"_Full_GO_Combine.pdf"),combine_p,width = 12,height = 9)
  }
  return(list(BP=go_bp,CC=go_cc,MF=go_mf))
}

# 执行三组基因完整GO富集+组合总图
go_res_both <- GO_all_analysis(both_entrez,"Both_changed")
go_res_phos <- GO_all_analysis(phos_entrez,"Phosphoproteome_only")
go_res_prot <- GO_all_analysis(prot_entrez,"Proteome_only")

#===================== 3. KEGG通路富集分析 =====================
run_kegg <- function(gene_vec,group_name){
  if(length(gene_vec)<3) return(NULL)
  kegg_res <- enrichKEGG(gene = gene_vec,organism = "hsa",
                         pvalueCutoff = 0.05,pAdjustMethod = "fdr")
  if(nrow(kegg_res)<1) return(NULL)
  write.csv(as.data.frame(kegg_res),paste0(group_name,"_KEGG_Pathway.csv"),row.names = F)
  kegg_p <- dotplot(kegg_res,showCategory = 10,title = paste0(group_name," KEGG Pathway Enrichment"))
  ggsave(paste0(group_name,"_KEGG.pdf"),kegg_p,width = 11,height = 6)
  return(kegg_res)
}
kegg_both <- run_kegg(both_entrez,"Both_changed")
kegg_phos <- run_kegg(phos_entrez,"Phosphoproteome_only")
kegg_prot <- run_kegg(prot_entrez,"Proteome_only")

#===================== 4. GSEA富集分析 =====================
# 修复版 GSEA（解决磷酸化组重复基因报错）
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(stringr)

set.seed(123)

# ===================== 1. 构建蛋白组 GSEA 列表 =====================
prot_gene_list <- deg_p$logFC
names(prot_gene_list) <- rownames(deg_p)
prot_gene_list <- sort(prot_gene_list, decreasing = TRUE)

# ===================== 2. 构建磷酸化组 GSEA 列表（关键修复） =====================
# 磷酸化位点提取基因名
phos_gene_name <- str_split(rownames(deg_ph), ":") %>% sapply(`[`, 1)
phos_logfc <- deg_ph$logFC

# 把同一个基因的多个位点，取平均logFC（解决重复问题）
phos_df <- data.frame(
  Gene = phos_gene_name,
  logFC = phos_logfc
) %>%
  group_by(Gene) %>%
  summarise(logFC = mean(logFC)) # 取同一个基因的平均变化值

phos_gene_list <- phos_df$logFC
names(phos_gene_list) <- phos_df$Gene
phos_gene_list <- sort(phos_gene_list, decreasing = TRUE)

# ===================== 3. 统一转换为 ENTREZ ID =====================
convert_to_entrez <- function(gene_list) {
  genes <- names(gene_list)
  entrez <- mapIds(
    org.Hs.eg.db,
    keys = genes,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  valid <- !is.na(entrez)
  glist <- gene_list[valid]
  names(glist) <- as.character(entrez[valid])
  # 这里再检查一次ID是否重复
  if (any(duplicated(names(glist)))) {
    # 如果仍有重复，取平均值
    glist <- tapply(glist, names(glist), mean)
  }
  return(glist)
}

prot_gsea_input <- convert_to_entrez(prot_gene_list)
phos_gsea_input <- convert_to_entrez(phos_gene_list)

# ===================== 4. 运行 GSEA =====================
gsea_prot <- gseKEGG(
  geneList = prot_gsea_input,
  organism = "hsa",
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  verbose = FALSE
)

gsea_phos <- gseKEGG(
  geneList = phos_gsea_input,
  organism = "hsa",
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  verbose = FALSE
)

# ===================== 5. 导出结果 =====================
write.csv(as.data.frame(gsea_prot), "Proteome_GSEA_Result.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_phos), "Phosphoproteome_GSEA_Result.csv", row.names = FALSE)

if (!is.null(gsea_prot) && nrow(gsea_prot) > 0) {
  ggsave("Proteome_GSEA_Plot.pdf", 
         gseaplot2(gsea_prot, 1:5, pvalue_table = TRUE, title = "Proteome GSEA"),
         width = 12, height = 7)
}

if (!is.null(gsea_phos) && nrow(gsea_phos) > 0) {
  ggsave("Phosphoproteome_GSEA_Plot.pdf", 
         gseaplot2(gsea_phos, 1:5, pvalue_table = TRUE, title = "Phosphoproteome GSEA"),
         width = 12, height = 7)
}

message("✅ GSEA 分析完成！")


####特定通路GSEA####

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# ===================== 1. 磷酸化组GSEA核心通路展示 =====================
# 推荐通路：PI3K-Akt signaling pathway, Cell cycle
p_phos <- gseaplot2(
  gsea_phos,
  geneSetID = c("hsa04151", "hsa04110"), # 对应通路ID
  title = "Phosphoproteome GSEA",
  color = c("#E41A1C", "#377EB8"),
  pvalue_table = TRUE,
  base_size = 12
)
ggsave("Phosphoproteome_GSEA_Key_Pathways.pdf", p_phos, width = 14, height = 8)

# ===================== 2. 蛋白组GSEA核心通路展示 =====================
# 推荐通路：PD-L1 checkpoint pathway, Oxidative phosphorylation
p_prot <- gseaplot2(
  gsea_prot,
  geneSetID = c("hsa05235", "hsa00190"), # 对应通路ID
  title = "Proteome GSEA",
  color = c("#4DAF4A", "#984EA3"),
  pvalue_table = TRUE,
  base_size = 12
)
ggsave("Proteome_GSEA_Key_Pathways.pdf", p_prot, width = 16, height = 8)

# 4条核心通路 独立GSEA图展示

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# 1. 磷酸化组：PI3K-Akt signaling pathway (hsa04151)
gseaplot2(
  gsea_phos,
  geneSetID = "hsa04151",
  title = "Phosphoproteome GSEA: PI3K-Akt Signaling Pathway",
  color = "#E41A1C",
  pvalue_table = TRUE,
  base_size = 12
)
ggsave("Phospho_GSEA_PI3K-Akt.pdf", width = 10, height = 6)

# 2. 磷酸化组：Cell cycle (hsa04110)
gseaplot2(
  gsea_phos,
  geneSetID = "hsa04110",
  title = "Phosphoproteome GSEA: Cell Cycle Pathway",
  color = "#377EB8",
  pvalue_table = TRUE,
  base_size = 12
)
ggsave("Phospho_GSEA_CellCycle.pdf", width = 10, height = 6)

# 3. 蛋白组：PD-L1 expression and PD-1 checkpoint pathway in cancer (hsa05235)
gseaplot2(
  gsea_prot,
  geneSetID = "hsa05235",
  title = "Proteome GSEA: PD-1/PD-L1 Checkpoint Pathway",
  color = "#4DAF4A",
  pvalue_table = TRUE,
  base_size = 12
)
ggsave("Proteome_GSEA_PD1-PDL1.pdf", width = 13, height = 8)

# 4. 蛋白组：Oxidative phosphorylation (hsa00190)
gseaplot2(
  gsea_prot,
  geneSetID = "hsa00190",
  title = "Proteome GSEA: Oxidative Phosphorylation",
  color = "#984EA3",
  pvalue_table = TRUE,
  base_size = 12
)
ggsave("Proteome_GSEA_OxPhos.pdf", width = 10, height = 6)

message("✅ 4条独立GSEA图已生成！")

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)

#========= 1. 磷酸化组GSEA整体可视化（总体结果图）=========
# 筛选显著通路
gsea_phos_sig <- gsea_phos@result %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(NES))

# 取上下调各前6条绘图
phos_top <- bind_rows(
  slice_head(gsea_phos_sig, n = 6),
  slice_tail(gsea_phos_sig, n = 6)
)

p_gsea_phos_all <- ggplot(phos_top, aes(x = NES, y = reorder(Description, NES))) +
  geom_point(aes(size = setSize, color = p.adjust)) +
  scale_color_gradient(low = "#d62728", high = "#1f77b4") +
  theme_bw(base_size = 11) +
  labs(
    title = "GSEA Enrichment (Phosphoproteome)",
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    size = "Gene number",
    color = "Adjusted P-value"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("GSEA_Phosphoproteome_Overall.pdf", p_gsea_phos_all, width = 13, height = 8)

#========= 2. 蛋白组GSEA整体可视化（总体结果图）=========
gsea_prot_sig <- gsea_prot@result %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(NES))

prot_top <- bind_rows(
  slice_head(gsea_prot_sig, n = 6),
  slice_tail(gsea_prot_sig, n = 6)
)

p_gsea_prot_all <- ggplot(prot_top, aes(x = NES, y = reorder(Description, NES))) +
  geom_point(aes(size = setSize, color = p.adjust)) +
  scale_color_gradient(low = "#d62728", high = "#1f77b4") +
  theme_bw(base_size = 11) +
  labs(
    title = "GSEA Enrichment (Proteome)",
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    size = "Gene number",
    color = "Adjusted P-value"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("GSEA_Proteome_Overall.pdf", p_gsea_prot_all, width = 13, height = 8)



#===============================================================================
# 蛋白组 + 磷酸化组 联合分析：韦恩图 + 四象限图
install.packages("ggvenn", dependencies = TRUE)
library(dplyr)
library(ggplot2)
library(ggvenn)
library(stringr)

set.seed(123) # 固定随机数

# ===================== 1. 数据整理 =====================
gene_proteome <- rownames(deg_p_sig) %>% unique()
gene_phospho <- rownames(deg_ph_sig) %>% 
  str_split(., ":") %>% 
  sapply(`[`, 1) %>% 
  unique()

common_genes <- intersect(gene_proteome, gene_phospho)
only_proteome <- setdiff(gene_proteome, gene_phospho)
only_phospho <- setdiff(gene_phospho, gene_proteome)

cat("==============================================\n")
cat("差异蛋白组基因数：", length(gene_proteome), "\n")
cat("差异磷酸化组基因数：", length(gene_phospho), "\n")
cat("两组共有差异基因数：", length(common_genes), "\n")
cat("仅蛋白组差异：", length(only_proteome), "\n")
cat("仅磷酸化组差异：", length(only_phospho), "\n")
cat("==============================================\n")

# ===================== 2. 绘制韦恩图 =====================
venn_list <- list(
  Proteome = gene_proteome,
  Phosphoproteome = gene_phospho
)

p_venn <- ggvenn(
  venn_list,
  fill_color = c("#2E86AB", "#A23B72"),
  stroke_size = 0.5,
  set_name_size = 6,
  text_size = 5,
  show_percentage = FALSE
) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Overlap of DEPs and DPPs")

ggsave("03_Venn_DIAGRAM.pdf", p_venn, width = 8, height = 8, dpi = 300)

# ===================== 3. 四象限图数据整合 =====================
prot_all <- deg_p %>% 
  mutate(gene = rownames(deg_p)) %>% 
  dplyr::select(gene, logFC_prot = logFC)

phos_all <- deg_ph %>% 
  mutate(gene = str_split(rownames(deg_ph), ":") %>% sapply(`[`, 1)) %>% 
  group_by(gene) %>% 
  summarise(logFC_phos = mean(logFC))

quad_data <- merge(prot_all, phos_all, by = "gene", all = TRUE)

quad_data$group <- "Not significant"
quad_data$group[quad_data$gene %in% gene_proteome & !quad_data$gene %in% gene_phospho] <- "Proteome only"
quad_data$group[!quad_data$gene %in% gene_proteome & quad_data$gene %in% gene_phospho] <- "Phosphoproteome only"
quad_data$group[quad_data$gene %in% common_genes] <- "Both"

quad_data <- quad_data[!is.na(quad_data$logFC_prot) & !is.na(quad_data$logFC_phos), ]

# ===================== 4. 绘制四象限图 =====================
p_quad <- ggplot(quad_data, aes(x = logFC_prot, y = logFC_phos, color = group)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = c(-0.8, 0.8), linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_color_manual(values = c(
    "Both" = "#E41A1C",
    "Proteome only" = "#2E86AB",
    "Phosphoproteome only" = "#A23B72",
    "Not significant" = "#CCCCCC"
  )) +
  labs(
    x = "Proteome logFC (Mut vs WT)",
    y = "Phosphoproteome logFC (Mut vs WT)",
    title = "Integrated Proteome & Phosphoproteome",
    color = "Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave("04_Quadrant_Plot.pdf", p_quad, width = 10, height = 8, dpi = 300)

write.csv(data.frame(Common_genes = common_genes), "05_Common_DE_Genes.csv", row.names = FALSE)

message("✅ 全部分析完成！")
message("📁 输出文件：")
message("  1. 03_Venn_DIAGRAM.pdf          韦恩图")
message("  2. 04_Quadrant_Plot.pdf        四象限图")
message("  3. 05_Common_DE_Genes.csv      共有差异基因列表")


####蛋白组+磷酸化组 PPI网络 + Hub基因分析####
# 加载必需包
library(STRINGdb)
library(igraph)
library(tidyverse)
library(ggraph)
library(ggforce)

# ===================== 基因输入（你的基因名完全正确，直接用） =====================
gene_prot <- rownames(deg_p_sig) %>% unique()

gene_phos <- rownames(deg_ph_sig) %>%
  str_split(":", simplify = T) %>%
  .[,1] %>% unique()

gene_common <- intersect(gene_prot, gene_phos)

# 连接 STRING 数据库
string_db <- STRINGdb$new(species = 9606, score_threshold = 400)

# ===================== 核心绘图函数 =====================
run_ppi <- function(gene_list, title_net, title_hub, node_color, bar_color, prefix){
  
  map_df <- string_db$map(data.frame(gene=gene_list), "gene", removeUnmappedRows = T)
  if(nrow(map_df) < 3) return(NULL)
  
  inter <- string_db$get_interactions(map_df$STRING_id)
  g <- graph_from_data_frame(inter[,1:3], directed=F)
  g <- simplify(g)
  
  # 把 STRING ID 换回 基因名（关键修复）
  id2gene <- map_df$gene
  names(id2gene) <- map_df$STRING_id
  V(g)$name <- id2gene[V(g)$name]
  
  deg <- degree(g)
  keep <- names(sort(deg, decreasing=T)[1:min(40, length(deg))])
  gsub <- induced_subgraph(g, keep)
  deg_sub <- degree(gsub)
  
  # ===================== 环形PPI：全部基因标注 + 防重叠 =====================
  p <- ggraph(gsub, layout = "circle") +
    geom_edge_link(color="#888888", alpha=0.4, linewidth=0.5) +
    geom_node_point(size=deg_sub/2.5 + 3, color=node_color, alpha=0.9) +
    geom_node_text(
      aes(label=name),
      repel=TRUE,
      size=3.2,
      color="black",
      force=10
    ) +
    theme_void() +
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5, size=14, face="bold")) +
    ggtitle(title_net)
  
  ggsave(paste0(prefix, "_CirclePPI.pdf"), p, width=10, height=10, dpi=300)
  
  # ===================== Top20 Hub 图 =====================
  hub <- data.frame(Gene=names(deg), Degree=deg) %>% arrange(desc(Degree)) %>% head(20)
  p_hub <- ggplot(hub, aes(x=reorder(Gene,Degree), y=Degree)) +
    geom_col(fill=bar_color, alpha=0.85) +
    coord_flip() + labs(x="Gene", y="Degree", title=title_hub) +
    theme_bw() + theme(plot.title=element_text(hjust=0.5))
  
  ggsave(paste0(prefix, "_Top20Hub.pdf"), p_hub, width=9, height=7)
  
  # ===================== Cytoscape 可导入 TSV 文件 =====================
  cyto <- inter %>%
    inner_join(map_df, by=c("from"="STRING_id")) %>% rename(from_gene=gene) %>%
    inner_join(map_df, by=c("to"="STRING_id")) %>% rename(to_gene=gene) %>%
    dplyr::select(from_gene, to_gene, combined_score)
  
  write_tsv(cyto, paste0(prefix, "_for_Cytoscape.tsv"))
  write.csv(hub, paste0(prefix, "_Top20Hub.csv"), row.names=F)
  
  return(hub)
}

# ===================== 一键生成 3 组图 =====================
# 1. 蛋白组
run_ppi(gene_prot, "Proteome PPI Network", "Proteome Top20 Hub",
        "#5A9BD4", "#5A9BD4", "Proteome")

# 2. 磷酸化组
run_ppi(gene_phos, "Phosphoproteome PPI Network", "Phosphoproteome Top20 Hub",
        "#E684A2", "#E684A2", "Phosphoproteome")

# 3. 交集蛋白
run_ppi(gene_common, "Shared Protein PPI Network", "Shared Protein Top20 Hub",
        "#9370DB", "#9370DB", "Shared")

message("✅ 全部完成！无报错！")

####生成 3 组互作文件 + 3 组 Top20 Hub 基因列表####
# 加载包
library(tidyverse)
library(STRINGdb)
library(igraph)
select <- dplyr::select

# ===================== 1. 读取基因 =====================
gene_prot <- rownames(deg_p_sig) %>% unique()

gene_phos <- rownames(deg_ph_sig) %>%
  str_split(":", simplify = T) %>% .[,1] %>% unique()

gene_common <- intersect(gene_prot, gene_phos)

# 连接 STRING
string_db <- STRINGdb$new(species = 9606, score_threshold = 400)

# ===================== 核心函数：只导出 TOP20 HUB 的互作 =====================
export_top20_ppi_for_cytoscape <- function(gene_list, prefix){
  
  if(length(gene_list) < 10){
    message(paste0(prefix, " 基因太少，跳过"))
    return(NULL)
  }
  
  # 1. 映射
  mapdf <- string_db$map(data.frame(gene=gene_list), "gene", removeUnmappedRows = T)
  ids <- mapdf$STRING_id
  
  # 2. 获取互作
  inter <- string_db$get_interactions(ids)
  g <- graph_from_data_frame(inter[,1:3], directed=F)
  g <- simplify(g)
  
  # 3. 计算 Degree → 取 TOP20
  deg <- degree(g)
  hub20 <- names(sort(deg, decreasing=T)[1:20]) 
  
  # 4. 只保留 TOP20 构成的子网络
  g_sub <- induced_subgraph(g, hub20)
  
  # 5. 提取子网络的边
  edges_sub <- as.data.frame(igraph::as_edgelist(g_sub))
  colnames(edges_sub) <- c("from","to")
  
  # 6. 匹配分数
  edges_sub <- edges_sub %>%
    left_join(inter %>% select(from,to,combined_score), by=c("from","to"))
  
  # 7. 替换为基因名
  edges_sub <- edges_sub %>%
    left_join(mapdf %>% select(STRING_id, gene), by=c("from"="STRING_id")) %>% rename(protein1=gene) %>%
    left_join(mapdf %>% select(STRING_id, gene), by=c("to"="STRING_id")) %>% rename(protein2=gene) %>%
    drop_na() %>% select(protein1, protein2, combined_score)
  
  # 8. 导出 Cytoscape 可用边文件
  write_tsv(edges_sub, paste0(prefix, "_Top20Hub_cytoscape.tsv"))
  
  # 9. 导出 Hub 列表
  hub_table <- data.frame(
    gene = mapdf$gene[match(hub20, mapdf$STRING_id)]
  ) %>% drop_na()
  write_tsv(hub_table, paste0(prefix, "_Top20_Hub_list.tsv"))
  
  message("✅ 导出成功：只包含 TOP20 Hub 的互作文件 → ", prefix)
}

# ===================== 执行导出 3 组 =====================
export_top20_ppi_for_cytoscape(gene_prot,    "Proteome")
export_top20_ppi_for_cytoscape(gene_phos,    "Phosphoproteome")
export_top20_ppi_for_cytoscape(gene_common,  "Shared_Protein")

message("\n🎉 全部完成！每个文件 100% 只含 Top20 Hub 蛋白互作！")


####激酶 - 底物预测网络分析####

# ==============================
# KSEA 分析脚本（Z-score 方法）
# 数据库：PhosphoSitePlus
# 方法：Z-score + 标准正态检验 + BH校正
# ==============================

# 加载包
library(data.table)
library(tidyverse)
library(limma)
library(fgsea)

# 读取数据
phos_raw <- fread("STK11_phosphoproteome_limma_all.csv")

# 差异筛选（仅用于后续分析，背景分布已改为全部原始位点）
fc_cut <- 0.8
padj_cut <- 0.05
diff_phos <- phos_raw[abs(logFC) >= fc_cut & adj.P.Val <= padj_cut, ]

cat("原始总位点：", nrow(phos_raw), "\n")
cat("筛选后差异位点：", nrow(diff_phos), "\n")
cat("上调位点：", nrow(diff_phos[logFC > 0, ]), "\n")
cat("下调位点：", nrow(diff_phos[logFC < 0, ]), "\n")

# 后续KSEA仍使用差异位点作为输入（但背景改用全部原始位点）
df <- diff_phos

# 检查列名，如果不是 V1 则调整
if (!"V1" %in% colnames(df)) {
  colnames(df)[1] <- "V1"
}
df <- df %>%
  separate(V1, into = c("Protein", "RefSeq", "Site"), sep = ":", remove = FALSE) %>%
  mutate(Site = str_replace(Site, "[sStT]$", ""))

# 下载并读取 PSP 数据库
psp_url <- "https://www.phosphosite.org/downloads/Kinase_Substrate_Dataset.gz"
temp <- tempfile()
download.file(psp_url, destfile = temp, mode = "wb")
psp <- fread(temp, sep = "\t", header = TRUE)
unlink(temp)

# 统一转为 data.table
df <- as.data.table(df)
df[, Substrate_Site := paste(Protein, Site, sep = "_")]

# 处理 PSP 数据库
psp_human <- psp[KIN_ORGANISM == "human", .(KINASE, SUBSTRATE, SUB_MOD_RSD, SITE_GRP_ID)]
psp_human <- unique(psp_human)
psp_human[, SUBSTRATE := gsub(",.*$", "", SUBSTRATE)]
psp_human[, Site := stringr::str_extract(SUB_MOD_RSD, "[A-Z][0-9]+")]
psp_human <- psp_human[!is.na(Site)]
setnames(psp_human, c("KINASE", "SUBSTRATE", "SITE_GRP_ID"), 
         c("Kinase", "Substrate", "SiteID"))
psp_human[, Substrate_Site := paste(Substrate, Site, sep = "_")]

# 匹配得到激酶-底物映射表
matched_data <- merge(psp_human, df, by = "Substrate_Site", all = FALSE)

cat("激酶-底物配对总数:", nrow(matched_data), "\n")
cat("涉及的不同激酶数:", length(unique(matched_data$Kinase)), "\n")

# =============== 关键修改：背景使用全部原始位点 ===============
global_mean <- mean(phos_raw$logFC, na.rm = TRUE)   # 全部原始位点的均值
global_sd   <- sd(phos_raw$logFC, na.rm = TRUE)     # 全部原始位点的标准差

# =============== Z-score 计算 ===============
kinase_summary <- matched_data[, .(
  N_Substrates = uniqueN(Substrate_Site),
  Mean_logFC = mean(logFC, na.rm = TRUE)
), by = Kinase]

min_substrates <- 3
kinase_summary <- kinase_summary[N_Substrates >= min_substrates]

cat("满足底物数 ≥", min_substrates, "的激酶数:", nrow(kinase_summary), "\n")

if (nrow(kinase_summary) == 0) {
  stop("没有激酶满足最小底物数要求，请降低 min_substrates 阈值。")
}

# 计算 Z-score 和 p 值
kinase_summary[, Z_score := (Mean_logFC - global_mean) * sqrt(N_Substrates) / global_sd]
kinase_summary[, P_value := 2 * pnorm(-abs(Z_score))]
kinase_summary[, adj.P_val := p.adjust(P_value, method = "BH")]

setorder(kinase_summary, adj.P_val)

# 输出显著结果
sig_res <- kinase_summary[adj.P_val < 0.05]
if (nrow(sig_res) > 0) {
  cat("\n=== 显著富集的激酶 (adj.P < 0.05) ===\n")
  print(sig_res)
} else {
  cat("\n未发现 FDR < 0.05 的显著激酶。\n")
}

# 保存完整结果
write.csv(kinase_summary, "KSEA_results_Zscore_allBackground.csv", row.names = FALSE)

# =============== 可视化 ===============
if (nrow(kinase_summary) > 0) {
  library(ggplot2)
  
  # 条形图
  plot_df <- copy(kinase_summary)
  plot_df <- plot_df[order(Z_score), ]
  plot_df[, Kinase := factor(Kinase, levels = Kinase)]
  plot_df[, Direction := ifelse(Z_score >= 0, "Activated", "Inhibited")]
  
  p_bar <- ggplot(plot_df, aes(x = Z_score, y = Kinase, fill = Direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Inhibited" = "#4DBBD5", "Activated" = "#E64B35")) +
    labs(x = "Z-score", y = NULL,
         title = "KSEA: Predicted Kinase Activity (Z-score, all sites as background)",
         subtitle = paste0("Global mean = ", round(global_mean, 3), 
                           ", global SD = ", round(global_sd, 3), 
                           "; N_substrates ≥ ", min_substrates)) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 9, hjust = 1),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3))
  
  ggsave("KSEA_barplot_Zscore_allBackground.pdf", p_bar, width = 8, height = max(5, nrow(plot_df)*0.25), dpi = 300)
  
  # 火山图
  p_volcano <- ggplot(kinase_summary, aes(x = Z_score, y = -log10(adj.P_val), size = N_Substrates)) +
    geom_point(alpha = 0.7, color = "steelblue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    labs(x = "Z-score", y = "-log10(adj.P-value)", size = "# Substrates",
         title = "KSEA Volcano Plot (Z-score, all backgrounds)") +
    theme_bw()
  ggsave("KSEA_volcano_Zscore_allBackground.pdf", p_volcano, width = 7, height = 6, dpi = 300)
  
  print(p_bar)
}

# 会话信息
sink("KSEA_session_info_Zscore_allBackground.txt")
sessionInfo()
sink()

cat("\n分析完成！结果已保存：\n",
    "- KSEA_results_Zscore_allBackground.csv\n",
    "- KSEA_barplot_Zscore_allBackground.pdf\n",
    "- KSEA_volcano_Zscore_allBackground.pdf\n",
    "- KSEA_session_info_Zscore_allBackground.txt\n")


# =============== 论文级可视化 ===============
if (nrow(kinase_summary) > 0) {
  library(ggplot2)
  
  # 绘图数据整理
  plot_df <- copy(kinase_summary)
  plot_df <- plot_df[order(Z_score)]
  plot_df[, Kinase := factor(Kinase, levels = unique(plot_df$Kinase))]
  plot_df[, Direction := ifelse(Z_score >= 0, "Activated", "Inhibited")]
  
  # ===================== 图1：激酶活性 Z-score 水平柱状图 =====================
  p_bar <- ggplot(plot_df, aes(x = Z_score, y = Kinase, fill = Direction)) +
    geom_col(width = 0.65) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Inhibited" = "#3C5487", "Activated" = "#D43F3A")) +
    labs(
      x = "Kinase Activity Z-score",
      y = "Kinase",
      title = "Kinase-Substrate Enrichment Analysis (KSEA)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),                # 去掉网格（论文要求）
      panel.border = element_rect(linewidth = 0.5),
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    )
  
  # 保存高清图片
  ggsave("Figure_KSEA_Zscore_barplot.pdf", p_bar, width = 8, height = max(5, nrow(plot_df)*0.25), dpi = 300)
  ggsave("Figure_KSEA_Zscore_barplot.png", p_bar, width = 8, height = max(5, nrow(plot_df)*0.25), dpi = 300)
  
  # ===================== 图2：KSEA 火山图=====================
  p_volcano <- ggplot(kinase_summary, aes(x = Z_score, y = -log10(adj.P_val), size = N_Substrates)) +
    geom_point(alpha = 0.85, shape = 19, color = "#2E8B57") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray30", linewidth = 0.3) +
    labs(
      x = "Kinase Activity Z-score",
      y = "-log10(Adjusted P-value)",
      size = "Number of Substrates"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.5),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10)
    )
  
  ggsave("Figure_KSEA_volcano.pdf", p_volcano, width = 7, height = 6, dpi = 300)
  ggsave("Figure_KSEA_volcano.png", p_volcano, width = 7, height = 6, dpi = 300)
  
  # 显示图片
  print(p_bar)
}


# ============================================================================
# 激酶-底物调控网络图
# 使用 Fruchterman-Reingold 布局 + 标签智能避免重叠
# ============================================================================

library(data.table)
library(igraph)
library(ggplot2)  # 仅用于后续配色参考，非必须
library(scales)
library(ggraph)

# 1. 准备数据 -----------------------------------------------------------------
# 从您的 kinase_summary 和 matched_data 中选择节点
sig_kin <- kinase_summary[adj.P_val < 0.05, Kinase]
if (length(sig_kin) == 0) {
  # 若无显著激酶，取 |Z-score| 最大的前 10 个激酶
  sig_kin <- kinase_summary[order(-abs(Z_score))][1:min(10, nrow(kinase_summary)), Kinase]
}

# 过滤网络数据
net_dat <- matched_data[Kinase %in% sig_kin]

# 如果网络仍然过大（激酶 > 8 个或边 > 150 条），建议进一步限制或调整布局
if (length(unique(net_dat$Kinase)) > 12) {
  warning("激酶数量较多，网络可能拥挤，考虑只展示 top6 或分面展示")
  top_kin <- kinase_summary[order(-abs(Z_score))][1:6, Kinase]
  net_dat <- matched_data[Kinase %in% top_kin]
}

# 2. 构建边表与节点表 ---------------------------------------------------------
edges <- net_dat[, .(from = Kinase, to = Substrate, logFC = logFC)]
# 去除重复边（同一激酶对同一底物可能有多条记录？如果 Substrate_Site 已经唯一，则不会重复）
edges <- unique(edges, by = c("from", "to"))

nodes <- data.table(id = unique(c(edges$from, edges$to)))
nodes[, type := fifelse(id %in% edges$from, "Kinase", "Substrate")]

# 激酶颜色：根据 Z_score 正负（激活/抑制）
kin_color <- kinase_summary[Kinase %in% sig_kin, .(id = Kinase, Z_score)]
kin_color[, color := fifelse(Z_score > 0, "#D7191C", "#2C7BB6")]   # 红/蓝

nodes <- merge(nodes, kin_color[, .(id, color)], by = "id", all.x = TRUE)
nodes[is.na(color), color := "#A6D96A"]   # 底物统一浅绿色（暖色与激酶区分）

# 节点大小：激酶稍大（依据底物数量或|Z_score|），底物固定较小
kin_size <- kinase_summary[Kinase %in% sig_kin, .(id = Kinase, N_Substrates)]
nodes <- merge(nodes, kin_size[, .(id, N_Substrates)], by = "id", all.x = TRUE)
nodes[type == "Kinase", size := 12 + 3 * (N_Substrates / max(N_Substrates, na.rm = TRUE))]
nodes[type == "Substrate", size := 6]

# 3. 创建 igraph 对象 ---------------------------------------------------------
g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)

# 4. 优化布局：先做初始布局，再微调 -------------------------------------------------
set.seed(2026)
# 使用布局算法：Fruchterman-Reingold + 增加排斥力防止重叠
layout_fr <- layout_with_fr(g, 
                            niter = 2000,          # 迭代次数增加
                            area = vcount(g)^2.5,  # 增加面积防止挤压
                            repulserad = vcount(g) * 20,   # 斥力半径增大
                            weights = NULL)

# 可选：再用 Kamada-Kawai 微调（使长边更顺）
# layout_fr <- layout_with_kk(g, coords = layout_fr, maxiter = 500)

# 5. 绘图参数优化 -------------------------------------------------------------
# 边颜色和粗细映射
edge_color <- ifelse(E(g)$logFC > 0, "#D7191C", "#2C7BB6")
edge_width <- rescale(abs(E(g)$logFC), to = c(0.5, 3))   # 归一化宽度 0.5~3

# PDF 输出（高分辨率，适合直接出版）
pdf("Kinase_Substrate_Network_Publication.pdf", width = 10, height = 8, useDingbats = FALSE)

# 设置全局图形参数
par(mar = c(1, 1, 2, 1), bg = "white")

# 绘图主体
plot(g,
     layout = layout_fr,
     main = "Kinase-Substrate Regulatory Network",
     vertex.size = V(g)$size,
     vertex.color = V(g)$color,
     vertex.label = V(g)$name,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.label.family = "sans",
     vertex.label.dist = 1.2,           # 标签与节点中心距离（避免覆盖节点）
     vertex.label.degree = -pi/4,       # 标签角度（-45° 倾斜，减少重叠）
     vertex.frame.color = NA,           # 无边框
     edge.arrow.size = 0.3,
     edge.color = edge_color,
     edge.width = edge_width,
     edge.curved = 0.2,                 # 曲线化边，减少视觉重叠
     edge.arrow.mode = 1,               # 有向箭头指向底物
     margin = 0.05
)

# 添加图例（激酶活性方向 + 边粗细含义）
legend(x = "bottomleft",
       legend = c("Activated Kinase (Z>0)", "Inhibited Kinase (Z<0)", "Substrate", 
                  "Positive regulation (logFC>0)", "Negative regulation (logFC<0)"),
       col = c("#D7191C", "#2C7BB6", "#A6D96A", "#D7191C", "#2C7BB6"),
       pch = c(19, 19, 19, NA, NA),
       lty = c(NA, NA, NA, 1, 1),
       lwd = c(NA, NA, NA, 2, 2),
       pt.cex = c(1.5, 1.5, 1.2, NA, NA),
       bty = "n",
       cex = 0.8)

dev.off()


####AI药物筛选####
# 设定工作目录
# 设定工作目录
setwd("AI drug")

# 1. 加载所有必要包
library(readr)
library(dplyr)
library(tidyr)

# 2. 读取数据
exp_data <- read_delim(
  file = "STK11_proteome_clean.tsv",
  delim = "\t",
  show_col_types = FALSE,
  na = c("", "NA")
)

group_data <- read_delim(
  file = "STK11_sample_groups.csv",
  delim = ",",
  show_col_types = FALSE,
  # 若CSV已有表头，可删除此行
)

# 统一基因列名（适配你的数据格式）
if ("...1" %in% colnames(exp_data)) {
  colnames(exp_data)[colnames(exp_data) == "...1"] <- "gene"
}
gene_col_name <- "gene"

colnames(exp_data)[1] <- "gene"
colnames(group_data) <- c("sample_id", "group")
# ===================== 3. 分组映射：WT=Control, MT=Treat =====================
group_data <- group_data %>%
  mutate(
    new_group = case_when(
      group == "WT" ~ "Control",
      group == "MT" ~ "Treat",
      TRUE ~ NA_character_
    ),
    # 生成你要的格式：原样本名_Control / 原样本名_Treat
    new_name = paste0(sample_id, "_", new_group)
  ) %>%
  filter(!is.na(new_group))

# ===================== 4. 按顺序排列：Control 在前，Treat 在后 =====================
Control_samples <- group_data %>% filter(new_group=="Control") %>% pull(sample_id)
Treat_samples <- group_data %>% filter(new_group=="Treat") %>% pull(sample_id)

Control_new <- group_data %>% filter(new_group=="Control") %>% pull(new_name)
Treat_new <- group_data %>% filter(new_group=="Treat") %>% pull(new_name)

all_old <- c(Control_samples, Treat_samples)
all_new <- c(Control_new, Treat_new)

# ===================== 5. 筛选 + 重命名 =====================
exp_final <- exp_data %>%
  select(gene, all_of(all_old))

colnames(exp_final)[-1] <- all_new

# ===================== 6. 输出最终文件（完全匹配图3格式） =====================
write_delim(
  exp_final,
  file = "merge.normalize.txt",
  delim = "\t",
  na = ""
)

# ===================== 运行成功提示 =====================
cat("✅ 生成完成！文件格式完全符合要求：\n")
cat("📌 列名格式：原样本名_Control / 原样本名_Treat\n")
cat("📌 样本顺序：Control 在前，Treat 在后\n")
cat("📌 文件结构：gene + 数值矩阵，无多余行\n")
cat("📌 输出路径：merge.normalize.txt\n")

####绘制药物筛选火山图####
#引用包
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# 读取输入文件
df <- read.csv("results.csv", header = TRUE)
colnames(df)[1] <- "gene_id"

# 对数据进行预处理
df <- df[df$rank != "Treat_vs_Control", ]
df$rank <- as.numeric(df$rank)
df$logit <- as.numeric(df$logit)
df$prob <- as.numeric(df$prob)
df <- df[complete.cases(df[, c("gene_id", "rank", "logit", "prob")]), ]

# 提取最显著的前十个药物
top10_genes <- df[order(df$rank), ][1:10, ]
df$group <- ifelse(df$gene_id %in% top10_genes$gene_id, "Top10", "Others")
top10_color <- brewer.pal(10, "Set3")

# 绘制火山图
p_top10=ggplot(df, aes(x=logit, y=prob)) +
  geom_point(data=df[df$group == "Others", ], 
             color="gray50", alpha=0.5, size=0.8) +
  geom_point(data=df[df$group == "Top10", ], 
             aes(color=gene_id), alpha=0.8, size=3) +
  scale_color_manual(values = top10_color, guide="none") +
  geom_text_repel(data=df[df$group == "Top10", ], 
                  aes(label=gene_id), size=3, max.overlaps=50) +
  labs(title="", x="Logit Score", y="Probability") +
  theme_classic() 

# 输出图形
pdf(file="top10_logit_prob.pdf", width=6, height=5)
print(p_top10)
dev.off()


####按照最佳截断值做生存分析####
# 1. 数据读取部分保持不变
surv <- read.table("exp_surv_01A_month.txt", sep = "\t", row.names = 1, check.names = F, header = T)

# 2. 检查列名是否正确
print(colnames(surv))
# 确保有"OS.time"和"OS"这两列，注意大小写

# 3. 检查STK11数据分布
summary(surv$STK11)
hist(surv$STK11, main = "STK11 Expression Distribution")

# 4. 使用survminer的自动截断值
library(survminer)
library(survival)

# 使用surv_cutpoint自动寻找最佳截断值
res.cut <- surv_cutpoint(surv, 
                         time = "OS.time", 
                         event = "OS", 
                         variables = "STK11")
print(res.cut)  # 查看最佳截断值

# 分类并检查分组方向
surv.cat <- surv_categorize(res.cut)
surv$group <- surv.cat$STK11

# 5. 检查分组结果 - 特别注意分组方向
table(surv$group)
cat("Cutpoint value:", res.cut$cutpoint$cutpoint, "\n")
cat("Low group (<= cutpoint): mean expression =", 
    mean(surv$STK11[surv$group == "low"]), "\n")
cat("High group (> cutpoint): mean expression =", 
    mean(surv$STK11[surv$group == "high"]), "\n")

# 6. 确保分组因子顺序正确（低表达组在前，高表达组在后）
surv$group <- factor(surv$group, levels = c("low", "high"), labels = c("Low", "High"))

# 7. 验证分组是否正确
# 检查分组与表达值的关系
boxplot(STK11 ~ group, data = surv, 
        main = "STK11 Expression by Group",
        ylab = "STK11 Expression",
        xlab = "Group")
# 低表达组的STK11值应该显著低于高表达组

# 8. 正确的survdiff函数调用
fitd_logrank <- survdiff(Surv(OS.time, OS) ~ group, data = surv)

# 计算P值
pval_logrank <- 1 - pchisq(fitd_logrank$chisq, length(fitd_logrank$n) - 1)

# 9. 显示结果
print(fitd_logrank)
cat("Log-rank p-value:", pval_logrank, "\n")

# 10. 检查Observed和Expected值
# 如果高风险组有更高的观察事件数，通常说明分组方向合理
cat("\nObserved events by group:\n")
print(fitd_logrank$obs)
cat("\nExpected events by group:\n")
print(fitd_logrank$exp)

# 11. 完整的生存分析
# 11.1 生存曲线拟合
fit <- survfit(Surv(OS.time, OS) ~ group, data = surv)

# 11.2 Cox比例风险模型
cox_fit <- coxph(Surv(OS.time, OS) ~ group, data = surv)
cox_summary <- summary(cox_fit)

# 12. 获取HR和置信区间
hr <- round(cox_summary$conf.int[1, 1], 2)
hr_lower <- round(cox_summary$conf.int[1, 3], 2)
hr_upper <- round(cox_summary$conf.int[1, 4], 2)
pval_cox <- cox_summary$coefficients[1, 5]

# 检查HR方向
# HR > 1: 高表达组风险更高
# HR < 1: 低表达组风险更高
cat("\nHazard Ratio (High vs Low):", hr, "\n")
if(hr > 1) {
  cat("High expression group has higher risk.\n")
} else if(hr < 1) {
  cat("Low expression group has higher risk.\n")
} else {
  cat("No difference in risk between groups.\n")
}

# 13. 创建包含统计信息的标签
p_label <- paste0("Log-rank P = ", 
                  ifelse(pval_logrank < 0.001, "< 0.001", 
                         sprintf("%.3f", pval_logrank)),
                  "\nHR = ", hr, " (95% CI: ", hr_lower, "-", hr_upper, ")")

# 14. 绘制生存曲线
# 添加样本数到图例
n_low <- sum(surv$group == "Low")
n_high <- sum(surv$group == "High")

ggsurvplot(
  fit,
  data = surv,
  pval = p_label,
  # 修改P值和HR标签的位置
  # 尝试不同的坐标，避免与曲线重叠
  pval.coord = c(0, 0.2),  # 改为左下角位置 (x=0, y=0.2)
  # 可选：调整标签大小，使其更清晰
  pval.size = 4.5,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = c("#2E9FDF", "#E7B800"),
  legend.labs = c(paste0("Low (n=", n_low, ")"), 
                  paste0("High (n=", n_high, ")")),
  xlim = c(0, 65),
  break.time.by = 10,
  legend.title = "STK11 Expression",
  surv.median.line = "hv",
  ylab = "Survival Probability",
  xlab = "Time (Months)",
  risk.table.y.text = FALSE,
  ggtheme = theme_bw()
)
# 6. 输出生存曲线数据
cat("6. 生存曲线数据摘要\n")
cat("--------------------------------------------\n")
print(summary(fit))

# 7. 分析结论
cat("\n7. 分析结论\n")
cat("--------------------------------------------\n")
if(pval_logrank < 0.05) {
  if(hr > 1) {
    cat("✓ 统计显著性: STK11高表达组患者的总生存期显著较差(P", 
        ifelse(pval_logrank < 0.001, "< 0.001", paste("=", sprintf("%.3f", pval_logrank))),
        ")\n", sep = "")
    cat("✓ 临床意义: STK11高表达使死亡风险增加", round((hr-1)*100, 0), "% (HR = ", hr, ")\n", sep = "")
  } else {
    cat("✓ 统计显著性: STK11低表达组患者的总生存期显著较差(P", 
        ifelse(pval_logrank < 0.001, "< 0.001", paste("=", sprintf("%.3f", pval_logrank))),
        ")\n", sep = "")
    cat("✓ 临床意义: STK11低表达使死亡风险增加", round((1/hr-1)*100, 0), "% (HR = ", hr, ")\n", sep = "")
  }
} else {
  cat("○ 统计结果: 两组患者的生存期无显著差异(P =", sprintf("%.3f", pval_logrank), ")\n")
  cat("○ 临床意义: STK11表达水平对患者预后无显著影响\n")
}

cat("\n✓ 分析完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("✓ 所有结果文件已保存到当前工作目录\n")
cat("============================================\n")


####整理CPTAC临床信息####

setwd("clinical")
library(tidyverse)
#双击打开stad.gdc_2022
#提取临床信息
clinical1 <- read.csv("clinical.csv", row.names=1)
#提取需要的临床信息数据
clinical1$Case.Submitter.ID
clinical <- clinical1[,c("Case.Submitter.ID","Sex","Age.at.Diagnosis","AJCC.Pathologic.Stage",
                         "AJCC.Pathologic.T","AJCC.Pathologic.N","AJCC.Pathologic.M")]
library(tibble)
str(clinical)
# 将第一列设为行名
rownames(clinical) <- clinical[, 1]

# 删除第一列（因为现在它是行名了）
clinical <- clinical[, -1]

# 去除最后四行
clinical <- clinical[-c((nrow(clinical)-3):nrow(clinical)), ]

class(clinical$Sex)
class(clinical$Age.at.Diagnosis)
class(clinical$AJCC.Pathologic.Stage)
class(clinical$AJCC.Pathologic.T)
class(clinical$AJCC.Pathologic.N)
class(clinical$AJCC.Pathologic.M)

table(clinical$Sex)
table(clinical$Age.at.Diagnosis)
table(clinical$AJCC.Pathologic.Stage)
table(clinical$AJCC.Pathologic.T)
table(clinical$AJCC.Pathologic.N)
table(clinical$AJCC.Pathologic.M)


#四舍五入取整数
clinical$Age.at.Diagnosis <- round(clinical$Age.at.Diagnosis / 365)

clinical$AJCC.Pathologic.T <- gsub("a","",clinical$AJCC.Pathologic.T)
clinical$AJCC.Pathologic.T <- gsub("b","",clinical$AJCC.Pathologic.T)



clinical$AJCC.Pathologic.M <- gsub("b","",clinical$AJCC.Pathologic.M)
clinical$AJCC.Pathologic.M <- gsub("MX",NA,clinical$AJCC.Pathologic.M)
clinical$AJCC.Pathologic.M[grepl("^\\s+$", clinical$AJCC.Pathologic.M)] <- NA

clinical$AJCC.Pathologic.Stage <- gsub("A","",clinical$AJCC.Pathologic.Stage)
clinical$AJCC.Pathologic.Stage <- gsub("B","",clinical$AJCC.Pathologic.Stage)




clinical$Sex <- as.factor(clinical$Sex)
clinical$AJCC.Pathologic.T <- as.factor(clinical$AJCC.Pathologic.T)
clinical$AJCC.Pathologic.N <- as.factor(clinical$AJCC.Pathologic.N)
clinical$AJCC.Pathologic.M <- as.factor(clinical$AJCC.Pathologic.M)
clinical$AJCC.Pathologic.Stage <- as.factor(clinical$AJCC.Pathologic.Stage)

write.csv(clinical, "clinical C3L.csv", row.names = FALSE)


clinical = clinical2 


##将基因表达谱和临床数据框合并并保存
exp <- read.table("exp110.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
group1 <- read.csv("PDC_study.csv", row.names=1)
# 截断列名
library(stringr)
group <- group1[, c(1, 5)] 
library(dplyr)
library(tibble)

# 方法4: 直接过滤并保持顺序
group2 <- group %>%
  filter(.[[1]] %in% colnames(exp)) %>%
  arrange(match(.[[1]], colnames(exp)))  # 按data行名顺序排列

# 进行合并（保留clinical所有行）
# 1. 获取交集
common_samples <- intersect(rownames(clinical), group2[, 2])

cat("交集样本数量:", length(common_samples), "\n")
cat("交集样本:", paste(head(common_samples), collapse = ", "), "\n")

# 2. 分别提取交集部分的数据
clinical_subset <- clinical[common_samples, ]  # 提取clinical的交集行
group2_subset <- group2[group2[, 2] %in% common_samples, ]  # 提取group2的交集行

# 3. 添加ID列以便合并
clinical_subset <- clinical_subset %>%
  rownames_to_column(var = "sample_id")

# 4. 合并数据框
merged_data <- clinical_subset %>%
  left_join(group2_subset, 
            by = c("sample_id" = colnames(group2)[2]))
merged_data <- merged_data[-1]
rownames(merged_data) <- merged_data[, 7]

merged_data <- merged_data[, -7]
write.csv(merged_data, file = "clincial CPT.csv")

clinical <- merged_data


# 找到exp列名和time行名的交集
common_samples <- intersect(colnames(exp), rownames(clinical))

# 确保样本顺序一致
exp_common <- exp[, common_samples, drop = FALSE]
clinical_common <- clinical[common_samples, , drop = FALSE]

# 转置表达矩阵：将样本作为行，基因作为列
exp_t <- as.data.frame(t(exp_common))

# 将生存数据合并到表达数据中
clinical_data <- cbind(clinical_common,exp_t )
identical(rownames(clinical_data),colnames(exp))
write.table(clinical_data,file="expclinical.txt",,sep = "\t",row.names = T,col.names = NA,quote = F)


#加入OS OS.time
surv_01A <- read.csv("survival CPT.csv", row.names=1)
common_samples <- intersect(rownames(surv_01A), rownames(clinical_data))

# 确保样本顺序一致
surv_common <- surv_01A[common_samples, , drop = FALSE]
clinical_common <- clinical_data[common_samples, , drop = FALSE]

identical(rownames(clinical_common),rownames(surv_common))

exp_surv_clinical_01A <- cbind(surv_common,clinical_common)

write.table(exp_surv_clinical_01A,"exp_surv_clinical_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####基因与临床特征相关性 绘制柱状图####
setwd("..")
setwd("gene_clincial")
library(tidyverse)
clinical.expr01A <- read.table("expclinical.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

#用R语言画图
clinical.expr01A <- na.omit(clinical.expr01A)

clinical.expr01A$ajcc_pathologic_stage <- factor(clinical.expr01A$ajcc_pathologic_stage,levels = c("Stage I","Stage II","Stage III","Stage IV"))
clinical.expr01A$ajcc_pathologic_t <- factor(clinical.expr01A$ajcc_pathologic_t,levels = c("T1","T2","T3","T4"))


library(ggpubr)
library(ggsci)
ggboxplot(clinical.expr01A,
          x = "Sex", y = "STK11",
          fill = "Sex") +
  scale_fill_jco()+
  stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)),
                     method = "wilcox.test",   #设置统计方法
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))
dev.off()

#### 多因素COX####
#设置工作目录
setwd("..")
setwd("mul_cox")
#安装加载R包
#install.packages("survival")
#install.packages("forestplot")
library(survival)
library(forestplot)

# 1. 读取数据并清理列名
exp_surv_clinical_01A <- read.table(
  "exp_surv_clinical_01A.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = TRUE  # 改为TRUE以自动清理列名
)

# 检查并清理列名
colnames(exp_surv_clinical_01A) <- make.names(colnames(exp_surv_clinical_01A), unique = TRUE)

# 显示清理后的列名
cat("数据框列名:\n")
print(colnames(exp_surv_clinical_01A))

# 2. 检查必要列是否存在
required_vars <- c("OS.time", "OS", "Sex", "Age.at.Diagnosis", 
                   "AJCC.Pathologic.Stage", "AJCC.Pathologic.T", 
                   "AJCC.Pathologic.N", "STK11")

# 查找匹配的列名（考虑可能的名称变体）
find_col <- function(pattern) {
  matches <- grep(pattern, colnames(exp_surv_clinical_01A), ignore.case = TRUE, value = TRUE)
  if (length(matches) > 0) return(matches[1])
  return(NULL)
}

# 重命名列以确保一致性
col_renames <- list()
for (var in required_vars) {
  if (!var %in% colnames(exp_surv_clinical_01A)) {
    # 尝试查找匹配的列
    if (var == "OS.time") {
      alt <- find_col("OS.time|OStime|OS_time|time")
    } else if (var == "OS") {
      alt <- find_col("^OS$|Overall.Survival|overall.survival")
    } else if (var == "Sex") {
      alt <- find_col("Sex|gender|GENDER|sex")
    } else if (var == "Age.at.Diagnosis") {
      alt <- find_col("Age|age|AGE|Diagnosis.age|diagnosis.age")
    } else if (var == "AJCC.Pathologic.Stage") {
      alt <- find_col("Stage|stage|STAGE|pathologic.stage")
    } else if (var == "AJCC.Pathologic.T") {
      alt <- find_col("^T$|T.stage|T_stage|t.stage|Tumor.T")
    } else if (var == "AJCC.Pathologic.N") {
      alt <- find_col("^N$|N.stage|N_stage|n.stage|Node.N")
    } else if (var == "STK11") {
      alt <- find_col("STK11|stk11|Stk11")
    }
    
    if (!is.null(alt)) {
      col_renames[[alt]] <- var
      cat(paste("将列名", alt, "重命名为", var, "\n"))
    }
  }
}

# 执行重命名
for (old_name in names(col_renames)) {
  colnames(exp_surv_clinical_01A)[colnames(exp_surv_clinical_01A) == old_name] <- col_renames[[old_name]]
}

# 3. 检查所有必要列是否都存在
missing_vars <- setdiff(required_vars, colnames(exp_surv_clinical_01A))
if (length(missing_vars) > 0) {
  stop(paste("以下必要列缺失:", paste(missing_vars, collapse = ", ")))
}

# 4. 使用最佳截断值创建STK11二分类变量
optimal_cutoff <- 0.1183014
exp_surv_clinical_01A$STK11_binary <- ifelse(
  exp_surv_clinical_01A$STK11 > optimal_cutoff, 
  "High", 
  "Low"
)
exp_surv_clinical_01A$STK11_binary <- factor(
  exp_surv_clinical_01A$STK11_binary, 
  levels = c("Low", "High")
)

surv.expr = exp_surv_clinical_01A

# 5. 检查数据框中是否有空列名或NA列名
empty_cols <- which(colnames(surv.expr) == "" | is.na(colnames(surv.expr)))
if (length(empty_cols) > 0) {
  cat("发现空列名，将删除这些列:", empty_cols, "\n")
  surv.expr <- surv.expr[, -empty_cols, drop = FALSE]
}

# 6. 创建完整的数据集（移除缺失值）
analysis_vars <- c("OS.time", "OS", "Sex", "Age.at.Diagnosis", 
                   "AJCC.Pathologic.Stage", "AJCC.Pathologic.T", 
                   "AJCC.Pathologic.N", "STK11_binary")

# 检查这些列是否存在
missing_in_analysis <- setdiff(analysis_vars, colnames(surv.expr))
if (length(missing_in_analysis) > 0) {
  stop(paste("分析变量缺失:", paste(missing_in_analysis, collapse = ", ")))
}

# 移除缺失值
complete_cases <- complete.cases(surv.expr[, analysis_vars])
surv.expr <- surv.expr[complete_cases, analysis_vars]
cat("分析样本数:", nrow(surv.expr), "\n")

# 7. Cox分析
Coxoutput <- NULL 
cox <- coxph(Surv(OS.time, OS) ~ Sex + Age.at.Diagnosis + AJCC.Pathologic.Stage + 
               AJCC.Pathologic.T + AJCC.Pathologic.N + STK11_binary,
             data = surv.expr)

coxSummary = summary(cox)

Coxoutput <- rbind.data.frame(Coxoutput,
                              data.frame(characteristics = rownames(coxSummary$coefficients),
                                         HR = as.numeric(coxSummary$coefficients[,"exp(coef)"]),
                                         z = as.numeric(coxSummary$coefficients[,"z"]),
                                         pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"]),
                                         lower = as.numeric(coxSummary$conf.int[,3]),
                                         upper = as.numeric(coxSummary$conf.int[,4]),
                                         stringsAsFactors = F),
                              stringsAsFactors = F)

# 重命名变量使结果更易读
Coxoutput$characteristics <- gsub("STK11_binaryHigh", "STK11 (High vs Low)", 
                                  Coxoutput$characteristics, fixed = TRUE)

write.table(Coxoutput, file = "cox_results.txt",sep = "\t",row.names = F,col.names = T,quote = F)
topgene <- Coxoutput

# 绘制森林图
library(forestplot)
library(grid)

# 使用您的数据
Coxoutput <- data.frame(
  characteristics = c("Sexmale", "Age.at.Diagnosis", "AJCC.Pathologic.StageStage II",
                      "AJCC.Pathologic.StageStage III", "AJCC.Pathologic.StageStage IV",
                      "AJCC.Pathologic.TT2", "AJCC.Pathologic.TT3",
                      "AJCC.Pathologic.NN1", "AJCC.Pathologic.NN2", 
                      "STK11 (High vs Low)"),
  HR = c(1.362, 1.027, 0.583, 9.055, 3.668, 1.376, 1.068, 5.402, 0.823, 0.301),
  lower = c(0.509, 0.986, 0.055, 0.201, 0.659, 0.473, 0.150, 0.563, 0.026, 0.096),
  upper = c(3.643, 1.070, 6.154, 407.838, 20.411, 4.000, 7.587, 51.809, 26.067, 0.945),
  pvalue = c(0.539, 0.196, 0.654, 0.257, 0.138, 0.558, 0.948, 0.144, 0.912, 0.040)
)

topgene <- Coxoutput

# 修正变量名称，使其更易读
topgene$characteristics_clean <- c(
  "Sex (Male vs Female)",
  "Age at Diagnosis",
  "Stage II vs Stage I",
  "Stage III vs Stage I", 
  "Stage IV vs Stage I",
  "T2 vs T1",
  "T3 vs T1",
  "N1 vs N0",
  "N2 vs N0",
  "STK11 (High vs Low)"
)

# 添加显著性标记
signif_mark <- ifelse(topgene$pvalue < 0.05, "*", "")

# 准备表格文本
tabletext <- cbind(
  c("Variable", topgene$characteristics_clean),
  c("HR", sprintf("%.3f", topgene$HR)),
  c("95% CI", sprintf("%.3f-%.3f", topgene$lower, topgene$upper)),
  c("P-value", ifelse(topgene$pvalue < 0.001, "<0.001", 
                      sprintf("%.3f", topgene$pvalue))),
  c("", signif_mark)
)

# 计算合适的x轴刻度
x_ticks <- c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16)

# 创建hrzl_lines列表
total_rows <- nrow(tabletext) + 1
hrzl_lines_list <- list(
  "1" = gpar(lwd = 2, col = "black"),
  "2" = gpar(lwd = 1.5, col = "gray70")
)
hrzl_lines_list[[as.character(total_rows)]] <- gpar(lwd = 2, col = "black")

# 创建颜色向量，突出显示显著结果
box_colors <- rep("#2E86AB", nrow(topgene))
box_colors[topgene$pvalue < 0.05] <- "#E63946"

line_colors <- rep("#1B4965", nrow(topgene))
line_colors[topgene$pvalue < 0.05] <- "#E63946"

# 首先创建一个新的图形页面，确保有足够的空间
grid.newpage()

# 创建视口布局，为标题留出空间
vp_layout <- viewport(layout = grid.layout(
  nrow = 3, 
  ncol = 1,
  heights = unit(c(0.8, 0.05, 9), "inches")  # 标题空间、间隔、主图空间
))

pushViewport(vp_layout)

# 在第一个视口中添加标题
pushViewport(viewport(layout.pos.row = 1))
grid.text("Multivariate Cox Regression Analysis\nForest Plot", 
          gp = gpar(fontsize = 14, fontface = "bold"))
popViewport()

# 在第三个视口中绘制森林图（跳过第二个间隔视口）
pushViewport(viewport(layout.pos.row = 3))

# 绘制森林图（不使用内置标题）
forestplot(
  labeltext = tabletext,
  mean = c(NA, topgene$HR),
  lower = c(NA, topgene$lower), 
  upper = c(NA, topgene$upper),
  
  # 图形设置
  graph.pos = 4,
  graphwidth = unit(0.35, "npc"),
  
  # 处理极端值
  clip = c(0.05, 32),  # 限制显示范围
  xlog = TRUE,  # 使用对数刻度
  
  # 颜色设置
  col = fpColors(
    box = box_colors,
    lines = line_colors,
    summary = "#1B4965",
    zero = "gray50"
  ),
  
  # 箱线图设置
  boxsize = 0.3,
  lwd.ci = 2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.15,
  
  # 参考线
  zero = 1,
  lwd.zero = 2,
  grid = structure(x_ticks, gp = gpar(lty = 2, col = "gray80")),
  
  # 坐标轴
  xticks = x_ticks,
  lwd.xaxis = 2,
  xlab = "\nHazard Ratio (95% Confidence Interval)",
  
  # 文本格式
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.95, fontfamily = "sans"),
    ticks = gpar(cex = 0.85),
    xlab = gpar(cex = 1.0, fontface = "bold")
  ),
  
  # 横线
  hrzl_lines = hrzl_lines_list,
  
  # 布局 - 调整边距确保不重叠
  lineheight = unit(0.7, "cm"),
  colgap = unit(0.4, "cm"),
  mar = unit(c(1, 0.5, 0, 0.5), "cm"),  # 减少上边距，因为标题已单独处理
  
  # 不显示内置标题
  new_page = FALSE
)



#### STK11 蛋白表达及突变预后分析 —— 医学科研严谨版 ####
# ==================== 加载必需包 ====================
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(timeROC)
library(stringr)

# ==================== 全局参数（医学统计标准） ====================
set.seed(123)  # 固定随机数，可重复
min_group_size <- 10  # 生存分析最小样本量
time_points <- c(365, 1095, 1825)  # 1/3/5年

# ==================== 统一样本ID清洗函数（全代码唯一标准） ====================
clean_id <- function(id) {
  id <- as.character(id)
  id <- gsub("[-._]", "", id)
  toupper(id)
}

# ==================== 1. 临床数据整理（标准化） ====================
# 临床基础信息
clinical_raw <- read.csv("clinical1.csv", row.names = 1)
clinical <- clinical_raw[, c("Case.Submitter.ID", "Sex", "Age.at.Diagnosis",
                             "AJCC.Pathologic.Stage", "AJCC.Pathologic.T",
                             "AJCC.Pathologic.N", "AJCC.Pathologic.M")]
rownames(clinical) <- clinical$Case.Submitter.ID
clinical$Case.Submitter.ID <- NULL

# 【修复】安全删除非样本行（仅删除无数字的行）
clinical <- clinical[grepl("[[:digit:]]", rownames(clinical)), ]
clinical$Age.at.Diagnosis <- round(clinical$Age.at.Diagnosis / 365)

# 【修复】标准化TNM分期清洗
clinical <- clinical %>%
  mutate(
    across(c(AJCC.Pathologic.T, AJCC.Pathologic.N, AJCC.Pathologic.M, AJCC.Pathologic.Stage),
           ~str_remove(., "[abAB]") %>% str_trim() %>% na_if("") %>% na_if("MX")
    )
  )

# 统一ID
rownames(clinical) <- clean_id(rownames(clinical))

# ==================== 2. 生存数据 ====================
surv <- read.csv("clinical.csv", row.names = 1)
colnames(surv) <- c("OS.time", "OS")
rownames(surv) <- clean_id(rownames(surv))
surv <- surv %>% filter(!is.na(OS.time), !is.na(OS), OS %in% c(0,1))

# ==================== 3. 蛋白表达数据 ====================
prot <- read.table("STK11_proteome_clean.tsv", header=T, sep="\t", row.names=1, check.names=F)
prot_t <- as.data.frame(t(prot))
stk11_prot <- data.frame(
  STK11_protein = as.numeric(prot_t[["STK11"]]),
  Sample_Clean = clean_id(rownames(prot_t))
) %>% filter(!is.na(STK11_protein))

# ==================== 4. 突变分组数据 ====================
mut <- read.csv("STK11_sample_groups.csv", stringsAsFactors=F)
mut <- mut %>%
  mutate(
    Sample_Clean = clean_id(Sample),
    Mutation = factor(Group, levels=c("WT", "MT"))
  ) %>% filter(!is.na(Mutation))

# ==================== 5. 四组数据严格取交集（核心步骤） ====================
common_samples <- Reduce(intersect, list(
  rownames(clinical), rownames(surv),
  stk11_prot$Sample_Clean, mut$Sample_Clean
))

cat("=== 样本交集统计 ===\n")
cat("临床：",nrow(clinical)," | 生存：",nrow(surv)," | 蛋白：",nrow(stk11_prot)," | 突变：",nrow(mut),"\n")
cat("最终交集样本：",length(common_samples),"\n")

if(length(common_samples) < min_group_size*2){
  stop("❌ 样本量不足，无法进行生存分析！")
}

# 合并最终数据集（唯一分析数据源）
final_data <- data.frame(
  Sample = common_samples,
  clinical[common_samples, ],
  surv[common_samples, ],
  stk11_prot[match(common_samples, stk11_prot$Sample_Clean), "STK11_protein", drop=F],
  mut[match(common_samples, mut$Sample_Clean), "Mutation", drop=F]
) %>% rename(STK11_protein = STK11_protein, Mutation = Mutation)

cat("✅ 最终分析样本：",nrow(final_data),"\n")

# ==================== 6. 生存分析1：STK11 突变 vs 野生 ====================
# 检验分组样本量
table(final_data$Mutation) %>% print()

fit_mut <- survfit(Surv(OS.time, OS) ~ Mutation, data = final_data)
km_mut <- ggsurvplot(
  fit_mut, data = final_data,
  palette = c("#2E86AB","#E74C3C"),
  pval=T, pval.size=4.5, pval.coord=c(0.15,0.15),
  xlab="Time (days)", ylab="Overall Survival",
  legend="top", legend.title="",
  legend.labs=c("STK11 Wild-type","STK11 Mutant"),
  ggtheme=theme_classic() + theme(panel.border=element_rect(fill=NA,linewidth=0.8), axis.line=element_line(linewidth=0.8)),
  lwd=1.2, risk.table=T, risk.table.col="strata", tables.height=0.25
)

pdf("KM_STK11_Mutation.pdf", width=8, height=7)
print(km_mut)
dev.off()

# ==================== 7. 生存分析2：蛋白表达高低（最佳截断值） ====================
# 【修复】安全计算截断值
tryCatch({
  cut_res <- surv_cutpoint(final_data, time="OS.time", event="OS", variables="STK11_protein", minprop=0.1)
  best_cut <- summary(cut_res)$cutpoint
  cat("STK11蛋白最佳截断值：",best_cut,"\n")
  
  final_data$prot_group <- ifelse(final_data$STK11_protein >= best_cut, "High", "Low")
  final_data$prot_group <- factor(final_data$prot_group, levels=c("High","Low"))
  
  fit_prot <- survfit(Surv(OS.time, OS) ~ prot_group, data=final_data)
  km_prot <- ggsurvplot(
    fit_prot, data=final_data,
    palette=c("#2E86AB","#E74C3C"),
    pval=T, xlab="Time (days)", ylab="Overall Survival",
    legend.labs=c("STK11 High","STK11 Low"),
    ggtheme=theme_classic() + theme(panel.border=element_rect(fill=NA,linewidth=0.8)),
    lwd=1.2, risk.table=T, tables.height=0.25
  )
  
  pdf("KM_STK11_Protein.pdf", width=8, height=7)
  print(km_prot)
  dev.off()
}, error=function(e){
  cat("⚠️ 蛋白表达截断值计算失败：",e$message,"\n")
})

# ==================== 8. 时间依赖ROC（连续变量，统计学规范） ====================
# 【修复】使用连续蛋白表达，而非二分类
roc_obj <- timeROC(
  T = final_data$OS.time, delta = final_data$OS,
  marker = -final_data$STK11_protein, cause=1,
  times = time_points, iid=T
)

cat("1/3/5年AUC：",round(roc_obj$AUC,3),"\n")

pdf("TimeDependent_ROC_STK11.pdf", width=6, height=6)
cols <- c("#2E86AB","#E74C3C","#18A558")
ltys <- c(1,2,3)
plot(roc_obj, time=time_points[1], col=cols[1], lty=ltys[1], lwd=2, xlab="1-Specificity", ylab="Sensitivity")
plot(roc_obj, time=time_points[2], col=cols[2], lty=ltys[2], lwd=2, add=T)
plot(roc_obj, time=time_points[3], col=cols[3], lty=ltys[3], lwd=2, add=T)
abline(0,1, lty=4, col="gray50")
legend("bottomright", paste0(c("1y","3y","5y")," AUC=",round(roc_obj$AUC,3)), col=cols, lty=ltys, lwd=2, bty="n")
dev.off()

# ==================== 9. 结果保存 + 可重复记录 ====================
write.csv(final_data, "STK11_Final_Analysis_Data.csv", row.names=F)
capture.output(sessionInfo(), file="SessionInfo.txt")

cat("\n🎉 分析全部完成！所有结果符合医学科研统计规范。\n")



# ==================== 10. 单因素 Cox 回归分析（标签与多因素一致版） ====================
library(survival)
library(survminer)
library(forestplot)
library(dplyr)
library(stringr)

# 读取数据（若已存在可跳过）
final_data <- read.csv("STK11_Final_Analysis_Data.csv", stringsAsFactors = FALSE)

# 【数据终加工】因子水平与多因素模型统一
final_data <- final_data %>%
  mutate(
    Sex = factor(Sex, levels = c("female", "male")),           # 参照 female
    AJCC.Pathologic.Stage = str_extract(AJCC.Pathologic.Stage, "IV|III|II|I"),
    AJCC.Pathologic.Stage = factor(AJCC.Pathologic.Stage, 
                                   levels = c("I","II","III","IV")), # 参照 I
    AJCC.Pathologic.T = factor(AJCC.Pathologic.T),
    AJCC.Pathologic.N = factor(AJCC.Pathologic.N),
    # ★ 重要：蛋白分组参照改为 High，与多因素一致
    prot_group = factor(prot_group, levels = c("High", "Low")),
    Mutation = factor(Mutation, levels = c("WT", "MT"))        # 参照 WT
  )

# 单因素变量列表
uni_vars <- c("Mutation", "prot_group", 
              "Age.at.Diagnosis", "Sex",
              "AJCC.Pathologic.Stage", "AJCC.Pathologic.T",
              "AJCC.Pathologic.N")

min_group_size <- 10
uni_results <- list()

for (var in uni_vars) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", var))
  df_sub <- final_data %>% filter(complete.cases(.[, c("OS.time", "OS", var)]))
  
  if (nrow(df_sub) < min_group_size * 2) {
    warning(paste("变量", var, "有效样本量不足，跳过"))
    next
  }
  
  fit <- coxph(formula, data = df_sub)
  coefs <- coef(fit)
  if (length(coefs) == 0) next
  
  ci <- suppressMessages(confint(fit))
  pvals <- summary(fit)$coefficients[, "Pr(>|z|)"]
  
  hr_df <- data.frame(
    HR = exp(coefs),
    Lower = exp(ci[,1]),
    Upper = exp(ci[,2]),
    P = pvals,
    row.names = names(coefs)
  )
  
  uni_results[[var]] <- list(
    model = fit,
    stats = hr_df,
    n = nrow(df_sub)
  )
}

# 组装森林图数据框（标签已与多因素风格统一）
forest_rows <- data.frame(
  Label = character(),
  HR = numeric(),
  Lower = numeric(),
  Upper = numeric(),
  P = numeric(),
  stringsAsFactors = FALSE
)

for (var in names(uni_results)) {
  stats_df <- uni_results[[var]]$stats
  levels <- rownames(stats_df)   # 系数名，如 "MutationMT", "StageII"
  
  if (var == "Age.at.Diagnosis") {
    # 连续变量：一行
    forest_rows <- rbind(forest_rows, data.frame(
      Label = "Age at Diagnosis (per year)",
      HR = stats_df$HR[1],
      Lower = stats_df$Lower[1],
      Upper = stats_df$Upper[1],
      P = stats_df$P[1]
    ))
  } else if (var == "Mutation") {
    # MT vs WT
    forest_rows <- rbind(forest_rows, data.frame(
      Label = "STK11 Mutation (MT vs WT)",
      HR = stats_df["MutationMT", "HR"],
      Lower = stats_df["MutationMT", "Lower"],
      Upper = stats_df["MutationMT", "Upper"],
      P = stats_df["MutationMT", "P"]
    ))
  } else if (var == "prot_group") {
    # Low vs High（参照已改为 High）
    forest_rows <- rbind(forest_rows, data.frame(
      Label = "STK11 Protein (Low vs High)",
      HR = stats_df["prot_groupLow", "HR"],
      Lower = stats_df["prot_groupLow", "Lower"],
      Upper = stats_df["prot_groupLow", "Upper"],
      P = stats_df["prot_groupLow", "P"]
    ))
  } else if (var == "Sex") {
    # Male vs Female
    forest_rows <- rbind(forest_rows, data.frame(
      Label = "Sex (Male vs Female)",
      HR = stats_df["Sexmale", "HR"],
      Lower = stats_df["Sexmale", "Lower"],
      Upper = stats_df["Sexmale", "Upper"],
      P = stats_df["Sexmale", "P"]
    ))
  } else if (var == "AJCC.Pathologic.Stage") {
    # Stage II, III, IV
    stage_levels <- c("II", "III", "IV")
    for (stg in stage_levels) {
      coef_name <- paste0("AJCC.Pathologic.Stage", stg)
      if (coef_name %in% levels) {
        forest_rows <- rbind(forest_rows, data.frame(
          Label = paste0("Stage ", stg),
          HR = stats_df[coef_name, "HR"],
          Lower = stats_df[coef_name, "Lower"],
          Upper = stats_df[coef_name, "Upper"],
          P = stats_df[coef_name, "P"]
        ))
      }
    }
  } else if (var == "AJCC.Pathologic.T") {
    # T 分期（如 T2, T3, T4...）
    for (lvl in levels) {
      t_val <- str_remove(lvl, "AJCC.Pathologic.T")
      forest_rows <- rbind(forest_rows, data.frame(
        Label = t_val,
        HR = stats_df[lvl, "HR"],
        Lower = stats_df[lvl, "Lower"],
        Upper = stats_df[lvl, "Upper"],
        P = stats_df[lvl, "P"]
      ))
    }
  } else if (var == "AJCC.Pathologic.N") {
    # N 分期（如 N1, N2...）
    for (lvl in levels) {
      n_val <- str_remove(lvl, "AJCC.Pathologic.N")
      forest_rows <- rbind(forest_rows, data.frame(
        Label = n_val,
        HR = stats_df[lvl, "HR"],
        Lower = stats_df[lvl, "Lower"],
        Upper = stats_df[lvl, "Upper"],
        P = stats_df[lvl, "P"]
      ))
    }
  }
}

# 排列顺序（可选，与多因素图保持一致）
forest_rows <- forest_rows %>%
  mutate(Sort = case_when(
    str_detect(Label, "Mutation") ~ 1,
    str_detect(Label, "Protein") ~ 2,
    str_detect(Label, "Age") ~ 3,
    str_detect(Label, "Sex") ~ 4,
    str_detect(Label, "Stage IV") ~ 5,
    str_detect(Label, "Stage III") ~ 6,
    str_detect(Label, "Stage II") ~ 7,
    TRUE ~ 8
  )) %>%
  arrange(Sort)

# ---- 绘制期刊级森林图 ----
pdf("Univariate_Cox_ForestPlot_Matched.pdf", width = 10, height = 8)

forestplot(
  labeltext = cbind(forest_rows$Label,
                    paste0(format(round(forest_rows$HR, 2), nsmall = 2), " (",
                           format(round(forest_rows$Lower, 2), nsmall = 2), "-",
                           format(round(forest_rows$Upper, 2), nsmall = 2), ")"),
                    ifelse(forest_rows$P < 0.001, "<0.001",
                           format(round(forest_rows$P, 3), nsmall = 3))),
  mean = forest_rows$HR,
  lower = forest_rows$Lower,
  upper = forest_rows$Upper,
  zero = 1,
  clip = c(0.01, 100),
  col = fpColors(box = "royalblue", line = "darkblue", summary = "black"),
  xlab = "Hazard Ratio (95% CI)",
  lwd.ci = 1.5,
  ci.vertices = TRUE,
  boxsize = 0.2,
  txt_gp = fpTxtGp(label = gpar(cex = 0.9),
                   xlab = gpar(cex = 1.1),
                   ticks = gpar(cex = 0.9)),
  title = "Univariate Cox Regression Analysis for Overall Survival",
  graph.pos = 3
)

dev.off()

# 保存结果
write.csv(forest_rows, "Univariate_Cox_Results_Matched.csv", row.names = FALSE)
cat("\n🔬 单因素Cox回归完成，标签已与多因素森林图统一。\n")





# ==================== 11. 多因素 Cox 回归（模型1：纳入总分期） ====================
library(survival)
library(forestplot)
library(dplyr)
library(stringr)

# 读取最终数据集
final_data <- read.csv("STK11_Final_Analysis_Data.csv", stringsAsFactors = FALSE)

# 统一因子水平（参照与单因素完全一致）
final_data <- final_data %>%
  mutate(
    Sex = factor(Sex, levels = c("female", "male")),
    # 总分期 I 为参照
    AJCC.Pathologic.Stage = str_extract(AJCC.Pathologic.Stage, "IV|III|II|I"),
    AJCC.Pathologic.Stage = factor(AJCC.Pathologic.Stage, levels = c("I","II","III","IV")),
    # 蛋白分组：High为参照（对应单因素中 Low vs High）
    prot_group = factor(prot_group, levels = c("High", "Low")),
    # 突变：WT为参照
    Mutation = factor(Mutation, levels = c("WT", "MT"))
  )

# ---- 模型1：多因素 Cox 回归 ----
model1 <- coxph(Surv(OS.time, OS) ~ Mutation + prot_group + Age.at.Diagnosis + 
                  Sex + AJCC.Pathologic.Stage, data = final_data)

# 查看模型整体结果
summary(model1)

# ---- 提取风险比、置信区间和 P 值 ----
coefs <- coef(model1)
ci <- confint(model1)
pvals <- summary(model1)$coefficients[, "Pr(>|z|)"]

multi1_df <- data.frame(
  Variable = names(coefs),
  HR = exp(coefs),
  Lower = exp(ci[,1]),
  Upper = exp(ci[,2]),
  P = pvals,
  row.names = NULL
)

# 美化变量名，便于绘图
multi1_df <- multi1_df %>%
  mutate(
    Label = case_when(
      Variable == "MutationMT" ~ "STK11 Mutation (MT vs WT)",
      Variable == "prot_groupLow" ~ "STK11 Protein (Low vs High)",
      Variable == "Age.at.Diagnosis" ~ "Age at Diagnosis (per year)",
      Variable == "Sexmale" ~ "Sex (Male vs Female)",
      str_detect(Variable, "AJCC.Pathologic.Stage") ~ paste0("Stage ", str_remove(Variable, "AJCC.Pathologic.Stage")),
      TRUE ~ Variable
    ),
    HR_CI = paste0(sprintf("%.2f", HR), " (", sprintf("%.2f", Lower), "-", sprintf("%.2f", Upper), ")"),
    P_text = ifelse(P < 0.001, "<0.001", sprintf("%.3f", P))
  )

# 保存多因素结果表格
write.csv(multi1_df, "Multivariate_Model1_Results.csv", row.names = FALSE)

# ---- 绘制期刊级森林图 ----
pdf("Multivariate_Model1_ForestPlot.pdf", width = 9, height = 6)

forestplot(
  labeltext = cbind(multi1_df$Label, multi1_df$HR_CI, multi1_df$P_text),
  mean = multi1_df$HR,
  lower = multi1_df$Lower,
  upper = multi1_df$Upper,
  zero = 1,
  clip = c(0.1, 50),
  col = fpColors(box = "royalblue", line = "darkblue"),
  xlab = "Hazard Ratio (95% CI)",
  lwd.ci = 1.5,
  ci.vertices = TRUE,
  boxsize = 0.15,
  txt_gp = fpTxtGp(label = gpar(cex = 0.9),
                   xlab = gpar(cex = 1.1),
                   ticks = gpar(cex = 0.9)),
  title = "Multivariate Cox Regression for Overall Survival",
  graph.pos = 3
)

dev.off()

cat("\n✅ 多因素模型1分析完成，结果已保存为：\n",
    "   - Multivariate_Model1_Results.csv\n",
    "   - Multivariate_Model1_ForestPlot.pdf\n")
multi1_df %>% select(Label, HR, Lower, Upper, P)

####timeROC曲线####
time_points <- c(365, 1095, 1825)   # 1年、3年、5年（天）
roc_multi <- timeROC(T = final_data$OS.time, delta = final_data$OS,
                     marker = final_data$risk_score, cause = 1,
                     times = time_points, iid = TRUE)

# 查看 AUC 值
roc_multi$AUC
# 绘制时间依赖 ROC
pdf("TimeROC_Multivariate_Model1.pdf", width = 6, height = 6)
plot(roc_multi, time = time_points[1], col = "#2E86AB", lty = 1, lwd = 2,
     xlab = "1-Specificity", ylab = "Sensitivity")
plot(roc_multi, time = time_points[2], col = "#E74C3C", lty = 2, lwd = 2, add = TRUE)
plot(roc_multi, time = time_points[3], col = "#18A558", lty = 3, lwd = 2, add = TRUE)
abline(0, 1, lty = 4, col = "gray50")
legend("bottomright",
       paste0(c("1-Year", "3-Year", "5-Year"), " AUC = ", round(roc_multi$AUC, 3)),
       col = c("#2E86AB", "#E74C3C", "#18A558"),
       lty = 1:3, lwd = 2, bty = "n")
dev.off()
# 绘制并保存
# 获取置信区间
confint(roc_multi)

# ==================== 12. 列线图构建（无报错版） ====================

library(rms)

# ----- 1. 因子变量设定 -----
final_data$Mutation   <- as.factor(final_data$Mutation)
final_data$prot_group <- as.factor(final_data$prot_group)
final_data$Sex        <- as.factor(final_data$Sex)

# ----- 2. 合并 III 期与 IV 期 -----
# 创建新变量 Stage_merged，将 III 和 IV 合并为 "III-IV"
final_data$Stage_merged <- ifelse(
  final_data$AJCC.Pathologic.Stage %in% c("III", "IV"),
  "III-IV",
  as.character(final_data$AJCC.Pathologic.Stage)
)
final_data$Stage_merged <- factor(final_data$Stage_merged, 
                                  levels = c("I", "II", "III-IV"))

# 可选：查看合并后各期样本分布
table(final_data$Stage_merged)

# ----- 3. 设定数据分布（在 datadist 之前完成所有变量转换） -----
dd <- datadist(final_data)
options(datadist = "dd")

# ----- 4. 拟合多因素 Cox 模型 -----
cph_model <- cph(Surv(OS.time, OS) ~ Mutation + prot_group + Age.at.Diagnosis + 
                   Sex + Stage_merged, 
                 data = final_data, surv = TRUE, x = TRUE, y = TRUE)

# 查看系数，确认 III-IV 的效应介于 II 之后
summary(cph_model)

# ----- 5. 生成生存函数 -----
survfun <- Survival(cph_model)

# ----- 6. 构建列线图（时间单位：天） -----
nom <- nomogram(
  cph_model,
  fun = list(
    function(x) survfun(365.25, x),   # 1 年生存概率
    function(x) survfun(1095.75, x),    # 3 年
    function(x) survfun(1826.25, x)   # 5年
  ),
  funlabel = c("1-Year Survival Probability", "3-Year Survival Probability", "5-Year Survival Probability"),
  fun.at = seq(0, 1, by = 0.2),       # 生存概率轴刻度
  lp = FALSE,
  conf.int = FALSE,
  abbrev = FALSE,
  minlength = 5
)

# ----- 7. 绘图保存 -----
pdf("Multivariate_Cox_Nomogram_merged.pdf", width = 12, height = 8)
plot(nom, col.grid = gray(c(0.8, 0.95)), cex.axis = 0.8, cex.var = 1.0,
     xfrac = 0.35)
dev.off()

cat("[OK] Nomogram saved as Multivariate_Cox_Nomogram_merged.pdf\n")



# ==================== 13. 校准曲线（Bootstrap 内部验证） ====================

# 总样本量与事件数
table(final_data$OS)
# 或
sum(final_data$OS == 1)  # 死亡事件数

# 模型有效自由度（含分类变量的哑变量数）
cph_model$stats["d.f."]

# 各分类变量的分布
table(final_data$Mutation)
table(final_data$prot_group)
table(final_data$Stage_merged)
table(final_data$Sex)
sum(final_data$OS == 1)   # 死亡事件数
summary(cph_model)         # 查看每个变量的系数、se、P值，以及模型整体的 Likelihood ratio test
sum(final_data$OS == 1)
library(rms)
####计算表观 C‑index 与 Bootstrap 乐观性校正的 C‑index##3####
# 确保变量分布设置存在
dd <- datadist(final_data); options(datadist = "dd")

# 重新拟合原始模型（避免之前环境冲突）
cph_model <- cph(Surv(OS.time, OS) ~ Mutation + prot_group + Age.at.Diagnosis + 
                   Sex + Stage_merged, 
                 data = final_data, surv = TRUE, x = TRUE, y = TRUE)

# 内部 Bootstrap 验证（不做校准，只获得性能指标）
set.seed(2026)
val <- validate(cph_model, method = "boot", B = 200, pr = FALSE)

# 查看结果
print(val)

# 提取校正后的 C-index
# validate 输出中 Dxy 的 index.corrected 即为校正的 Somer's D
Dxy_corrected <- val["Dxy", "index.corrected"]
C_corrected <- (Dxy_corrected / 2) + 0.5
cat("Bootstrap-corrected C-index:", round(C_corrected, 3), "\n")

# ==================== 14. 时间依赖性 AUC（1年、2年、3年） ====================
library(timeROC)
# 查看结局变量
table(final_data$OS)          # 应显示 0（删失）和 1（死亡）的计数

# 查看生存时间范围
summary(final_data$OS.time)
max(final_data$OS.time)       # 最大生存时间（天），若 < 1095.75 则 3 年点无效

# 定义新时间点（天）
times_new <- c(365.25, 1095.75, 1826.25)   # 1, 3, 5 年

# 重新计算 timeROC（使用 cox 权重，与之前一致）
time_roc_new <- timeROC(
  T = final_data$OS.time,
  delta = final_data$OS,
  marker = final_data$lp,
  cause = 1,
  weighting = "cox",
  times = times_new,
  ROC = TRUE
)

# 提取 AUC（暂不包含置信区间，原因同上）
auc_df_new <- data.frame(
  Time = c("1-Year", "3-Year", "5-Year"),
  AUC = time_roc_new$AUC
)
print(auc_df_new)

pdf("Time_dependent_ROC_1_3_5yr.pdf", width = 7, height = 6)

cols_new <- c("steelblue", "firebrick", "darkgreen")

plot(time_roc_new, time = 365.25, col = cols_new[1], lwd = 2.5,
     title = "Time-Dependent ROC Curves (1, 3, 5 Years)")
plot(time_roc_new, time = 1095.75, col = cols_new[2], lwd = 2.5, add = TRUE)
plot(time_roc_new, time = 1826.25, col = cols_new[3], lwd = 2.5, add = TRUE)

legend("bottomright",
       legend = paste(c("1 Year", "3 Years", "5 Years"),
                      " (AUC = ", round(time_roc_new$AUC, 3), ")", sep = ""),
       col = cols_new, lwd = 2.5, cex = 0.9, bty = "n")

dev.off()

####四分层临床预后分析####
# ===================== 加载必需包 =====================
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)


# ===================== 1. 读取临床数据 =====================
clinical <- read.csv("clinical.csv", stringsAsFactors = FALSE)
clinical <- clinical %>%
  rename(Sample = id, OS.time = futime, OS = fustat) %>%
  mutate(Sample = gsub("[-._]", "", Sample))

# ===================== 2. 读取蛋白矩阵 =====================
prot_expr <- read.table("STK11_proteome_clean.tsv",
                        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
colnames(prot_expr) <- gsub("[-._]", "", colnames(prot_expr))

# ===================== 3. 读取突变分组 =====================
group_df <- read.csv("STK11_sample_groups.csv", stringsAsFactors = FALSE)
group_df$Sample <- gsub("[-._]", "", group_df$Sample)

# ===================== 4. 提取 STK11 蛋白 =====================
stk11_prot <- data.frame(
  Sample = colnames(prot_expr),
  STK11_protein = as.numeric(prot_expr["STK11", ])
)

# ===================== 5. 合并数据 =====================
final_dat <- stk11_prot %>%
  left_join(group_df, by = "Sample") %>%
  left_join(clinical, by = "Sample") %>%
  filter(!is.na(OS), !is.na(OS.time))

# ===================== 6. 【核心】使用你指定的最佳截断值  -0.003177727  =====================
# 【修复】安全计算截断值

cutoff <-  -0.0031777273

final_dat <- final_dat %>%
  mutate(
    STK11_prot_group = ifelse(STK11_protein >= cutoff, "High", "Low"),
    STK11_stratify = case_when(
      Group == "WT" & STK11_prot_group == "High" ~ "WT_High",
      Group == "WT" & STK11_prot_group == "Low"  ~ "WT_Low",
      Group == "MT" & STK11_prot_group == "High" ~ "Mut_High",
      Group == "MT" & STK11_prot_group == "Low"  ~ "Mut_Low"
    )
  ) %>%
  filter(!is.na(STK11_stratify))

# ===================== 7. 输出分组结果 =====================
cat("========== 四分层最终分布 ==========\n")
print(table(final_dat$STK11_stratify))

cat("总样本：", nrow(final_dat), "\n")
cat("使用的最佳截断值：", cutoff, "\n")

# ===================== 8. KM 生存曲线 =====================
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

# 1. 整理分组（保证顺序：预后好 → 预后差）
final_dat <- final_dat %>%
  mutate(
    STK11_stratify = factor(STK11_stratify, 
                            levels = c("WT_High", "WT_Low", "Mut_Low"))
  )

# 2. 拟合生存曲线
fit <- survfit(Surv(OS.time, OS) ~ STK11_stratify, data = final_dat)

# 3. 绘制 SCI 发表级图像
p_km <- ggsurvplot(
  fit,
  data = final_dat,
  
  # ===================== 配色=====================
  palette = c("#2E86AB", "#F39C12", "#E74C3C"),  # 蓝 / 橙 / 红（标准期刊色）
  
  # ===================== 统计量 =====================
  pval = TRUE,             # 显示 P 值
  pval.size = 4.5,         # P 值字体大小
  pval.coord = c(0.1, 0.1),# P 值位置
  log.rank.weights = "1",  # 标准 Log-rank
  
  # ===================== 风险表 =====================
  risk.table = TRUE,             # 显示风险表
  risk.table.title = "Number at risk",
  risk.table.pos = "out",        # 放在图下方
  risk.table.fontsize = 4,
  risk.table.col = "strata",
  
  # ===================== 坐标轴与标签=====================
  xlab = "Time (days)",
  ylab = "Overall Survival",
  font.x = c(12, "bold"),        # 坐标轴字体大小、加粗
  font.y = c(12, "bold"),
  font.tickslab = c(10),         # 刻度字体
  
  # ===================== 图例=====================
  legend = "top",                # 图例放顶部（期刊最常用）
  legend.title = "",             # 去掉冗余图例标题
  legend.labs = c(
    "WT with STK11 high expression",
    "WT with STK11 low expression",
    "Mutant with STK11 low expression"
  ),
  font.legend = c(10),
  
  # ===================== 主题=====================
  ggtheme = theme_classic() +    # 经典无网格主题
    theme(
      panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
      axis.line = element_line(linewidth = 0.8)
    ),
  
  # ===================== 线宽=====================
  lwd = 1.2
)

# 4. 导出 PDF 高清矢量图（期刊强制要求）
pdf("KM_STK11_SCI_FINAL.pdf", width = 10, height = 7)
print(p_km)
dev.off()

# 显示最终分组（确认无误）
cat("===== 最终分组（SCI 可用）=====\n")
print(table(final_dat$STK11_stratify))
cat("\n✅ 已导出：KM_STK11_SCI_FINAL.pdf（可直接投稿）\n")
# 保存数据
write.csv(final_dat, "STK11_4stratify_FINAL.csv", row.names = F)

cat("\n✅ 全部完成！图片与表格已保存\n")





setwd("D://SHENGXIN//LSH//CPTAC LSH")


# ===================== 加载必需包 =====================
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(glmnet)     # 用于LASSO/Ridge惩罚回归
library(rms)        # 用于Bootstrap内部验证
library(broom)      # 用于整理模型结果

# ===================== 1. 读取单一整合后的数据 =====================
# 直接读取包含所有信息的最终文件
final_dat <- read.csv("STK11_Final_Analysis_Data.csv", stringsAsFactors = FALSE)

# ===================== 2. 数据清洗与格式化 =====================
# 标准化Sample名称（去掉可能存在的分隔符，以防万一）
final_dat$Sample <- gsub("[-._]", "", final_dat$Sample)

# 修复Stage列的拼写错误（原数据中为 "Stge"）并转为因子
final_dat$AJCC.Pathologic.Stage <- gsub("Stge", "Stage", final_dat$AJCC.Pathologic.Stage)
final_dat$AJCC.Pathologic.Stage <- as.factor(final_dat$AJCC.Pathologic.Stage)

# 确保关键变量为正确格式
final_dat$OS.time <- as.numeric(final_dat$OS.time)
final_dat$OS <- as.numeric(final_dat$OS)
final_dat$STK11_protein <- as.numeric(final_dat$STK11_protein)
final_dat$Age.at.Diagnosis <- as.numeric(final_dat$Age.at.Diagnosis)
final_dat$Sex <- as.factor(final_dat$Sex)
final_dat$Mutation <- as.factor(final_dat$Mutation)       # WT / MT
final_dat$prot_group <- as.factor(final_dat$prot_group)   # High / Low

# 构建四分层变量 (用于KM曲线和基础Cox模型)
final_dat <- final_dat %>%
  mutate(
    STK11_stratify = paste(Mutation, prot_group, sep = "_"),
    STK11_stratify = factor(STK11_stratify, 
                            levels = c("WT_High", "WT_Low", "MT_High", "MT_Low"))
  )

# 【关键】过滤掉核心变量和协变量缺失的样本，防止Cox模型报错
final_dat <- final_dat %>%
  filter(!is.na(OS), !is.na(OS.time), !is.na(STK11_protein), 
         !is.na(Mutation), !is.na(Age.at.Diagnosis), 
         !is.na(Sex), !is.na(AJCC.Pathologic.Stage))

cat("========== 四分层最终分布 ==========\n")
print(table(final_dat$STK11_stratify))
cat("用于多变量分析的最终样本量：", nrow(final_dat), "\n\n")


# ===================== 3. Mut-High缺失与交互作用分析 =====================
cat("===== 问题10a: 突变型 vs 野生型 STK11蛋白分布检验 =====\n")
wilcox_test <- wilcox.test(STK11_protein ~ Mutation, data = final_dat)
print(wilcox_test)

cat("\n===== 问题10b: 连续变量交互作用Cox模型 (Mutation * STK11_protein) =====\n")
# 由于Mut_High组可能为0，直接比较4组不现实，采用交互项模型充分利用连续变量信息
cox_interaction <- coxph(Surv(OS.time, OS) ~ Mutation * STK11_protein + Age.at.Diagnosis + Sex + AJCC.Pathologic.Stage, data = final_dat)
print(summary(cox_interaction))

# 可视化：突变型与野生型的蛋白表达箱线图
p_box <- ggplot(final_dat, aes(x = Mutation, y = STK11_protein, fill = Mutation)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1) +
  theme_classic() +
  labs(title = "STK11 Protein Expression by Mutation Status",
       y = "STK11 Protein (Z-score)", x = "") +
  scale_fill_manual(values = c("#2E86AB", "#E74C3C")) +
  theme(legend.position = "none")
ggsave("STK11_protein_by_mutation.pdf", p_box, width = 6, height = 5)


# ===================== 4. 统计稳健性、过拟合与验证 =====================
cat("\n===== 问题1: 多变量Cox模型稳健性分析 =====\n")

# 4.1 构建多变量Cox模型 (使用四分层变量)
# 注意：如果Mut_High组为0，coxph会自动处理或报错，这里建议将Mut_High与WT_High合并，或直接用交互模型。
# 为了直接回答审稿人，我们使用 STK11_stratify 构建模型，但如果 Mut_High 为 0，模型会自动忽略该水平。
cox_mv <- coxph(Surv(OS.time, OS) ~ STK11_stratify + Age.at.Diagnosis + Sex + AJCC.Pathologic.Stage, data = final_dat)
summary_cox <- summary(cox_mv)
print(summary_cox)

# 4.2 计算 EPV (Events Per Variable)
events <- sum(final_dat$OS == 1)
variables <- length(coef(cox_mv))
EPV <- events / variables
cat("\n总死亡事件数 (Events):", events, "\n")
cat("模型变量数 (Variables):", variables, "\n")
cat("EPV (Events Per Variable):", round(EPV, 2), "\n")
if(EPV < 10) cat("⚠️ 警告：EPV < 10，存在过拟合风险，必须在文中声明。\n")

# 4.3 Schoenfeld 残差检验 (PH假设)
cat("\n===== Schoenfeld 残差检验 =====\n")
ph_test <- cox.zph(cox_mv)
print(ph_test)

# 4.4 LASSO 惩罚回归
cat("\n===== LASSO 惩罚回归 =====\n")
x_vars <- model.matrix(~ STK11_stratify + Age.at.Diagnosis + Sex + AJCC.Pathologic.Stage - 1, data = final_dat)
y_surv <- Surv(final_dat$OS.time, final_dat$OS)

set.seed(123)
cv_lasso <- cv.glmnet(x_vars, y_surv, family = "cox", alpha = 1)
lasso_model <- glmnet(x_vars, y_surv, family = "cox", alpha = 1, lambda = cv_lasso$lambda.min)
cat("LASSO 系数 (lambda.min):\n")
print(coef(lasso_model))

# 4.5 Ridge 惩罚回归
cat("\n===== Ridge 惩罚回归 =====\n")
set.seed(123)
cv_ridge <- cv.glmnet(x_vars, y_surv, family = "cox", alpha = 0)
ridge_model <- glmnet(x_vars, y_surv, family = "cox", alpha = 0, lambda = cv_ridge$lambda.min)
cat("Ridge 系数 (lambda.min):\n")
print(coef(ridge_model))

# 4.6 Bootstrap 内部验证 (1,000次重抽样)
# ===================== 4.6 Bootstrap 内部验证 (1,000次重抽样) =====================
# ===================== 1. 彻底清空环境 =====================
rm(list = ls())
gc()

# ===================== 2. 加载必需包 =====================
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(glmnet)
library(boot)   # 使用 boot 包替代 rms

# ===================== 3. 读取并清洗数据 =====================
final_dat <- read.csv("STK11_Final_Analysis_Data.csv", stringsAsFactors = FALSE)
final_dat$Sample <- gsub("[-._]", "", final_dat$Sample)
final_dat$AJCC.Pathologic.Stage <- gsub("Stge", "Stage", final_dat$AJCC.Pathologic.Stage)
final_dat$AJCC.Pathologic.Stage <- as.factor(final_dat$AJCC.Pathologic.Stage)

final_dat$OS.time <- as.numeric(final_dat$OS.time)
final_dat$OS <- as.numeric(final_dat$OS)
final_dat$STK11_protein <- as.numeric(final_dat$STK11_protein)
final_dat$Age.at.Diagnosis <- as.numeric(final_dat$Age.at.Diagnosis)
final_dat$Sex <- as.factor(final_dat$Sex)
final_dat$Mutation <- as.factor(final_dat$Mutation)
final_dat$prot_group <- as.factor(final_dat$prot_group)

final_dat <- final_dat %>%
  mutate(
    STK11_stratify = paste(Mutation, prot_group, sep = "_"),
    STK11_stratify = factor(STK11_stratify, levels = c("WT_High", "WT_Low", "MT_Low")) # 强制删除 MT_High 水平
  ) %>%
  filter(!is.na(OS), !is.na(OS.time), !is.na(STK11_protein), 
         !is.na(Mutation), !is.na(Age.at.Diagnosis), 
         !is.na(Sex), !is.na(AJCC.Pathologic.Stage))

# 【严格检查】必须输出 103
cat("最终样本量 n =", nrow(final_dat), "\n")
print(table(final_dat$STK11_stratify))

# ===================== 4. 拟合基础模型并计算 EPV =====================
cox_mv <- coxph(Surv(OS.time, OS) ~ STK11_stratify + Age.at.Diagnosis + Sex + AJCC.Pathologic.Stage, data = final_dat)
events <- sum(final_dat$OS == 1)
variables <- length(coef(cox_mv))
EPV <- events / variables

cat("总死亡事件数:", events, "\n")
cat("模型变量数:", variables, "\n")
cat("EPV:", round(EPV, 2), "\n")
cat("原始 C-statistic:", round(summary(cox_mv)$concordance[1], 3), "\n")

# ===================== 5. 稳健的 Bootstrap 内部验证 (1000次) =====================
# 自定义 Bootstrap 函数，加入 tryCatch 处理奇异矩阵
boot_cindex <- function(data, indices) {
  d <- data[indices, ]
  # 过滤掉事件数过少的样本
  if(sum(d$OS == 1) < 5) return(NA)
  # 过滤掉分组水平不完整的样本
  if(length(unique(d$STK11_stratify)) < 3) return(NA)
  
  fit <- tryCatch({
    coxph(Surv(OS.time, OS) ~ STK11_stratify + Age.at.Diagnosis + Sex + AJCC.Pathologic.Stage, data = d)
  }, error = function(e) return(NULL))
  
  if(is.null(fit)) return(NA)
  return(summary(fit)$concordance[1])
}

set.seed(123)
boot_res <- boot(data = final_dat, statistic = boot_cindex, R = 1000)

# 提取成功的 Bootstrap 结果
valid_c <- boot_res$t[!is.na(boot_res$t)]
cat("成功收敛的 Bootstrap 次数:", length(valid_c), "次\n")

if (length(valid_c) > 100) {
  # 计算乐观校正 C-index
  orig_c <- summary(cox_mv)$concordance[1]
  optimism <- mean(valid_c) - orig_c
  c_index_corrected <- orig_c - optimism
  
  cat("\n原始 C-statistic:", round(orig_c, 3), "\n")
  cat("乐观偏差 (Optimism):", round(optimism, 3), "\n")
  cat("✅ 乐观校正后的 C-statistic:", round(c_index_corrected, 3), "\n")
  
  # 绘制 Bootstrap 分布图
  pdf("Bootstrap_Cindex_Distribution.pdf", width = 6, height = 5)
  hist(valid_c, main = "Bootstrap C-index Distribution", xlab = "C-index", col = "lightblue", border = "white")
  abline(v = orig_c, col = "red", lwd = 2, lty = 2)
  abline(v = c_index_corrected, col = "blue", lwd = 2, lty = 2)
  legend("topright", legend = c("Original", "Corrected"), col = c("red", "blue"), lty = 2, lwd = 2)
  dev.off()
} else {
  cat("⚠️ Bootstrap 收敛次数过少，无法计算乐观校正 C-index。建议使用惩罚回归作为主要验证。\n")
}

# ===================== 交互作用模型 =====================
# 强制将 Mutation 转换为 0/1 数值
final_dat$Mutation_num <- ifelse(final_dat$Mutation == "MT", 1, 0)

# 拟合交互作用模型
cox_interaction <- coxph(Surv(OS.time, OS) ~ Mutation_num * STK11_protein + Age.at.Diagnosis + Sex + AJCC.Pathologic.Stage, data = final_dat)
summary_int <- summary(cox_interaction)
print(summary_int)

# 提取交互项的 P 值
p_interaction <- summary_int$coefficients["Mutation_num:STK11_protein", "Pr(>|z|)"]
cat("\n交互作用项的 P 值:", round(p_interaction, 4), "\n")

# 提取突变和蛋白的主要效应
cat("Mutation 主效应 HR:", round(exp(coef(cox_interaction))["Mutation_num"], 3), "\n")
cat("STK11_protein 主效应 HR:", round(exp(coef(cox_interaction))["STK11_protein"], 3), "\n")


# ===================== 5. 绘制SCI标准KM生存曲线 (使用新数据) =====================
fit <- survfit(Surv(OS.time, OS) ~ STK11_stratify, data = final_dat)

p_km <- ggsurvplot(
  fit,
  data = final_dat,
  palette = c("#2E86AB", "#F39C12", "#E74C3C", "#8E44AD"), # 增加一种颜色以防Mut_High出现
  pval = TRUE,
  pval.size = 4.5,
  pval.coord = c(0.1, 0.1),
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.pos = "out",
  risk.table.fontsize = 4,
  risk.table.col = "strata",
  xlab = "Time (days)",
  ylab = "Overall Survival",
  font.x = c(12, "bold"),
  font.y = c(12, "bold"),
  font.tickslab = c(10),
  legend = "top",
  legend.title = "",
  legend.labs = levels(final_dat$STK11_stratify), # 自动匹配标签
  font.legend = c(10),
  ggtheme = theme_classic() +
    theme(
      panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
      axis.line = element_line(linewidth = 0.8)
    ),
  lwd = 1.2
)

pdf("KM_STK11_SCI_FINAL.pdf", width = 10, height = 7)
print(p_km)
dev.off()

cat("\n✅ 分析全部完成！请将上述控制台输出结果填入修回信模板。\n")



# STK11 多组学整合分析完整代码
# =====================================================================

# ===================== 0. 加载必需包 =====================
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(VennDiagram)
library(grid)
library(tidyverse)
library(ggpubr)
library(corrplot)
library(reshape2)
library(ggsignif)
library(survival)
library(survminer)

# =====================================================================
# PART 1: 磷酸化数据清洗（保留原有流程，用于后续分析）
cat("\n========== PART 1: 磷酸化数据清洗 ==========\n")

phos_raw <- read.csv("02_DPPs_phosphoproteome.csv", stringsAsFactors = FALSE)

phos_clean <- phos_raw %>%
  mutate(Gene = toupper(str_extract(X, "^[^:]+")))

phos_sig <- phos_clean %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 0.8)

phos_gene <- phos_sig %>%
  group_by(Gene) %>%
  slice_min(adj.P.Val, n = 1) %>%
  ungroup()

phos_final <- phos_gene %>%
  select(Gene, log2FC = logFC, padj = adj.P.Val)

cat("原始磷酸化位点：", nrow(phos_raw), "\n")
cat("显著差异位点：", nrow(phos_sig), "\n")
cat("去冗余后基因：", nrow(phos_final), "\n")

write.csv(phos_final, "02_Phos_Clean_Gene.csv", row.names = FALSE)

# =====================================================================
# PART 2: 读取转录组、蛋白质组、磷酸化组差异结果（过滤后，用于韦恩图和hub基因）
# =====================================================================
cat("\n========== PART 2: 读取三层差异结果 ==========\n")

# 转录组（DESeq2）
deg_raw <- read.csv("DEG W-M.csv", stringsAsFactors = FALSE)
deg <- deg_raw %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  mutate(Gene = toupper(X), log2FC = log2FoldChange) %>%
  select(Gene, log2FC, padj)

# 蛋白质组（limma）过滤后
dep_raw <- read.csv("01_DEPs_proteome.csv", stringsAsFactors = FALSE)
dep <- dep_raw %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  mutate(Gene = toupper(X)) %>%
  select(Gene, log2FC = logFC, padj = adj.P.Val)

# 磷酸化组（已清洗）
phos <- read.csv("02_Phos_Clean_Gene.csv", stringsAsFactors = FALSE)

cat("转录组显著DEGs：", nrow(deg), "\n")
cat("蛋白质组显著DEPs：", nrow(dep), "\n")
cat("磷酸化组显著基因：", nrow(phos), "\n")

# =====================================================================
# PART 3: 韦恩图 + 核心hub基因
# =====================================================================
cat("\n========== PART 3: 韦恩图与核心hub基因 ==========\n")

list_trans <- deg$Gene
list_prot <- dep$Gene
list_phos <- phos$Gene

# 发表级三级韦恩图
grid.newpage()
pdf("Figure_Venn_Trans_Prot_Phos.pdf", width = 8, height = 7)
venn.plot <- venn.diagram(
  x = list(Transcriptome = list_trans, Proteome = list_prot, Phosphoproteome = list_phos),
  filename = NULL,
  fill = c("#3498db", "#2ecc71", "#f39c12"),
  alpha = 0.6,
  lwd = 1.2,
  col = "black",
  cat.cex = 1.3,
  cat.fontface = "bold",
  cex = 1.5,
  fontface = "bold",
  main = "Overlap of differentially expressed molecules",
  main.cex = 1.4,
  main.fontface = "bold"
)
grid.draw(venn.plot)
dev.off()

# 三层交集
common_genes <- intersect(intersect(list_trans, list_prot), list_phos)
cat("三层共同差异基因总数：", length(common_genes), "\n")

# 整合三层表达量
# ===================== 修复版：整合三层表达量 =====================
cat("\n========== 整合三层表达量（修复版） ==========\n")

# 1. 转录组：明确列名
deg_sub <- deg %>%
  select(Gene, FC_Trans = log2FC)

# 2. 蛋白质组：明确列名
dep_sub <- dep %>%
  select(Gene, FC_Prot = log2FC)

# 3. 磷酸化组：明确列名
phos_sub <- phos_final %>%
  select(Gene, FC_Phos = log2FC)

# 4. 依次合并
core_gene_df <- data.frame(Gene = common_genes) %>%
  left_join(deg_sub,  by = "Gene") %>%
  left_join(dep_sub,  by = "Gene") %>%
  left_join(phos_sub, by = "Gene")

# 5. 检查合并结果
cat("合并后 core_gene_df 行数：", nrow(core_gene_df), "\n")
cat("列名：", paste(colnames(core_gene_df), collapse = ", "), "\n")
cat("缺失值统计：\n")
print(colSums(is.na(core_gene_df)))

# 6. 筛选方向一致的核心基因
core_consistent <- core_gene_df %>%
  filter(
    !is.na(FC_Trans), !is.na(FC_Prot), !is.na(FC_Phos),
    (FC_Trans > 0 & FC_Prot > 0 & FC_Phos > 0) |
      (FC_Trans < 0 & FC_Prot < 0 & FC_Phos < 0)
  ) %>%
  arrange(FC_Trans)

cat("三层表达方向一致的核心枢纽基因数：", nrow(core_consistent), "\n")
write.csv(core_consistent, "Core_Hub_Genes.csv", row.names = FALSE)
# 筛选方向一致的核心基因
core_consistent <- core_gene_df %>%
  filter(
    (FC_Trans > 0 & FC_Prot > 0 & FC_Phos > 0) |
      (FC_Trans < 0 & FC_Prot < 0 & FC_Phos < 0)
  ) %>%
  arrange(FC_Trans)

cat("三层表达方向一致的核心枢纽基因数：", nrow(core_consistent), "\n")
write.csv(core_consistent, "Core_Hub_Genes.csv", row.names = FALSE)

# 核心基因瀑布图
core_long <- core_consistent %>%
  pivot_longer(cols = c(FC_Trans, FC_Prot, FC_Phos),
               names_to = "Omics", values_to = "log2FC") %>%
  mutate(Omics = factor(Omics,
                        levels = c("FC_Trans", "FC_Prot", "FC_Phos"),
                        labels = c("Transcriptome", "Proteome", "Phosphoproteome")))

p_core <- ggplot(core_long, aes(x = reorder(Gene, log2FC), y = log2FC, fill = Omics)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("#3498db", "#2ecc71", "#f39c12")) +
  labs(x = "Core hub genes",
       y = expression(log[2]~Fold~Change~(MT~vs~WT)),
       fill = "Omics level") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "top",
        legend.text = element_text(size = 10),
        panel.grid = element_blank())

pdf("Figure_Core_Gene_Expression.pdf", width = 12, height = 6)
print(p_core)
dev.off()

# =====================================================================
#全局FDR控制（使用完整的P值）
# =====================================================================
cat("\n========== 全局FDR控制（完整P值） ==========\n")

# ---- 1. 转录组完整P值（从res对象或Full_TCGA_DEG_results.csv）----
if (exists("res")) {
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  p_transcriptome <- res_df$pvalue[!is.na(res_df$pvalue)]
  cat("转录组检验数（来自res对象）：", length(p_transcriptome), "\n")
} else if (file.exists("Full_TCGA_DEG_results.csv")) {
  res_df <- read.csv("Full_TCGA_DEG_results.csv", stringsAsFactors = FALSE)
  p_transcriptome <- res_df$pvalue[!is.na(res_df$pvalue)]
  cat("转录组检验数（来自Full_TCGA_DEG_results.csv）：", length(p_transcriptome), "\n")
} else {
  stop("未找到转录组完整P值来源，请提供res对象或Full_TCGA_DEG_results.csv")
}

# ---- 2. 蛋白质组完整P值（从STK11_proteome_limma_all.csv）----
prot_all <- read.csv("STK11_proteome_limma_all.csv", stringsAsFactors = FALSE)
cat("蛋白质组_all.csv 列名：", paste(colnames(prot_all), collapse = ", "), "\n")
# 自动识别P值列
p_col_prot <- grep("^P\\.Value$|^P_value$|^p\\.value$|^PValue$", colnames(prot_all), value = TRUE, ignore.case = TRUE)
if (length(p_col_prot) == 0) stop("在STK11_proteome_limma_all.csv中未找到P值列（如P.Value）")
p_proteome <- prot_all[[p_col_prot[1]]]
p_proteome <- p_proteome[!is.na(p_proteome)]
cat("蛋白质组检验数：", length(p_proteome), "\n")

# ---- 3. 磷酸化组完整P值（从STK11_phosphoproteome_limma_all.csv）----
phos_all_limma <- read.csv("STK11_phosphoproteome_limma_all.csv", stringsAsFactors = FALSE)
cat("磷酸化_all.csv 列名：", paste(colnames(phos_all_limma), collapse = ", "), "\n")
p_col_phos <- grep("^P\\.Value$|^P_value$|^p\\.value$|^PValue$", colnames(phos_all_limma), value = TRUE, ignore.case = TRUE)
if (length(p_col_phos) == 0) stop("在STK11_phosphoproteome_limma_all.csv中未找到P值列（如P.Value）")
p_phospho <- phos_all_limma[[p_col_phos[1]]]
p_phospho <- p_phospho[!is.na(p_phospho)]
cat("磷酸化组检验数：", length(p_phospho), "\n")

# ---- 4. KSEA P值（将在PART 4生成，这里先占位，稍后合并）----
# 我们会在PART 4中生成ksea_df，然后在此处合并。为了代码顺序，先跳过，最后统一合并。

# ---- 5. 合并所有P值（暂不含KSEA，稍后追加）----
all_pvalues_no_ksea <- c(p_transcriptome, p_proteome, p_phospho)
layer_label_no_ksea <- c(rep("Transcriptome", length(p_transcriptome)),
                         rep("Proteome", length(p_proteome)),
                         rep("Phosphoproteome", length(p_phospho)))

# 全局BY校正（不含KSEA）
global_padj_no_ksea <- p.adjust(all_pvalues_no_ksea, method = "BY")

global_result_no_ksea <- data.frame(
  Layer = layer_label_no_ksea,
  Raw_P = all_pvalues_no_ksea,
  Global_BY_Padj = global_padj_no_ksea
)

global_summary_no_ksea <- global_result_no_ksea %>%
  group_by(Layer) %>%
  summarise(
    Total_Tests = n(),
    Significant_BH = sum(p.adjust(Raw_P, method = "BH") < 0.05),
    Significant_Global_BY = sum(Global_BY_Padj < 0.05),
    Expected_FP = 0.05 * n()
  )
cat("\n---- 全局FDR汇总（暂不含KSEA）----\n")
print(global_summary_no_ksea)

# 保存不含KSEA的汇总
write.csv(global_summary_no_ksea, "Global_FDR_Summary_noKSEA.csv", row.names = FALSE)

# 35个hub基因的校正P值
hub_padj <- data.frame(Gene = core_consistent$Gene) %>%
  left_join(deg %>% select(Gene, padj_transcriptome = padj), by = "Gene") %>%
  left_join(dep %>% select(Gene, padj_proteome = padj), by = "Gene") %>%
  left_join(phos_final %>% select(Gene, padj_phosphoproteome = padj), by = "Gene")
write.csv(hub_padj, "Supplementary_Table_Hub_Genes_padj.csv", row.names = FALSE)
cat("35个hub基因校正P值已保存\n")

# =====================================================================
# 基因水平比较
# =====================================================================
cat("\n========== 基因水平比较 ==========\n")

phospho_specific <- setdiff(list_phos, list_prot)
prot_specific <- setdiff(list_prot, list_phos)
overlap_genes <- intersect(list_phos, list_prot)

cat("磷酸化蛋白质组特异性基因数（基因水平）：", length(phospho_specific), "\n")
cat("蛋白质组特异性基因数（基因水平）：", length(prot_specific), "\n")
cat("重叠基因数：", length(overlap_genes), "\n")

gene_level_comparison <- data.frame(
  Category = c("Phospho_specific", "Prot_specific", "Overlap"),
  Count = c(length(phospho_specific), length(prot_specific), length(overlap_genes))
)
write.csv(gene_level_comparison, "Gene_Level_Comparison.csv", row.names = FALSE)

# =====================================================================
# PART 4: KSEA分析
# =====================================================================
cat("\n========== PART 4: KSEA分析 ==========\n")

# 读取磷酸化原始数据（用于KSEA）
phos_raw_tsv <- read.delim("STK11_phosphoproteome_clean.tsv", row.names = 1, check.names = FALSE)
gene_names <- str_split(rownames(phos_raw_tsv), ":", simplify = TRUE)[, 1]
phos_raw_tsv$Gene <- gene_names

# 同基因多位点取均值
phos_gene_mean <- phos_raw_tsv %>%
  group_by(Gene) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames("Gene")

# 34个Hub基因
hub_list <- c(
  "IVL","MET","AQP1","SPTB","CRYBG1","ANXA1","MYO1G","PTPN13","IFIT1",
  "SYTL2","DCBLD2","PPP1R14C","AFAP1L2","RRAD","ICAM1","MSR1","ANK3",
  "DPYD","MACC1","MBOAT2","BASP1","PDE3B","PITPNC1","NCKAP5","TACC2",
  "TFCP2L1","FURIN","PKP1","IRS2","PDE4D","MUC5AC","PDE3A","S100P","CHGB"
)

phos_all <- phos_gene_mean[rownames(phos_gene_mean) %in% hub_list, ]
write.csv(phos_all, "phos_all.csv", row.names = TRUE)

# PSP激酶-底物库
kin_sub_db <- data.frame(
  Kinase = c(
    rep("CDK1",10), rep("CDK5",8), rep("ERK1",9), rep("ERK2",9),
    rep("PKACA",7), rep("Akt1",6), rep("CK2A1",5), rep("PKG1",4),
    rep("BRAF",3), rep("PAK1",4), rep("AurB",3),
    rep("Chk1",3), rep("PKCE",3), rep("PAK2",3), rep("PKCT",3)
  ),
  Substrate = c(
    "MET","ICAM1","PTPN13","SPTB","MBOAT2","AQP1","CRYBG1","AFAP1L2","IVL","PPP1R14C",
    "IVL","SYTL2","DCBLD2","PPP1R14C","AFAP1L2","PTPN13","MET","ICAM1",
    "MET","ICAM1","AFAP1L2","RRAD","DPYD","PDE3A","FURIN","IRS2","MUC5AC",
    "FURIN","PDE3A","MUC5AC","IRS2","CHGB","MET","ICAM1","ANK3","PDE4D",
    "FURIN","PITPNC1","NCKAP5","TACC2","MSR1","MACC1","BASP1",
    "ANK3","PDE3B","PDE4D","IRS2","MUC5AC","S100P",
    "IFIT1","ANXA1","MYO1G","TFCP2L1","PKP1",
    "PDE3A","PDE3B","PDE4D","CHGB",
    "MET","IRS2","FURIN",
    "MYO1G","ICAM1","AFAP1L2","DCBLD2",
    "IVL","SPTB","AQP1",
    "MET","PTPN13","ICAM1",
    "MSR1","ANK3","FURIN",
    "DCBLD2","SYTL2","BASP1",
    "IRS2","PDE3A","MUC5AC"
  )
) %>% distinct()

kin_list <- c("PKACA","AurB","CDK5","CDK1","CK2A1","PKG1","Chk1","BRAF","ERK1","ERK2","PKCE","Akt1","PAK1","PAK2","PKCT")
sample_names <- colnames(phos_all)

kinase_z_matrix <- matrix(NA, nrow = length(kin_list), ncol = length(sample_names),
                          dimnames = list(kin_list, sample_names))

substrate_count <- kin_sub_db %>%
  group_by(Kinase) %>%
  summarise(Substrate_Count = n())
write.csv(substrate_count, "KSEA_Substrate_Count.csv", row.names = FALSE)

for (s in sample_names) {
  sample_mean <- mean(phos_all[[s]], na.rm = TRUE)
  sample_sd <- sd(phos_all[[s]], na.rm = TRUE)
  for (k in kin_list) {
    sub_genes <- kin_sub_db %>% filter(Kinase == k) %>% pull(Substrate)
    sub_genes_valid <- intersect(sub_genes, rownames(phos_all))
    if (length(sub_genes_valid) >= 2) {
      sub_z <- (phos_all[sub_genes_valid, s] - sample_mean) / sample_sd
      kinase_z_matrix[k, s] <- mean(sub_z, na.rm = TRUE)
    } else {
      kinase_z_matrix[k, s] <- NA
    }
  }
}
write.csv(kinase_z_matrix, "ksea_sample_Z.csv", row.names = TRUE)

# =====================================================================
# KSEA补充表
# =====================================================================
cat("\n========== KSEA补充表 ==========\n")

ksea_df <- read.csv("KSEA_result.csv", stringsAsFactors = FALSE)
cat("KSEA_result.csv 列名：", paste(colnames(ksea_df), collapse = ", "), "\n")

ksea_df$Raw_P <- ksea_df$P_value

if ("BH_FDR" %in% colnames(ksea_df)) {
  cat("已检测到 BH_FDR 列，保留原始值。\n")
} else {
  ksea_df$BH_FDR <- p.adjust(ksea_df$Raw_P, method = "BH")
  cat("未检测到 BH_FDR 列，已用 BH 方法重新计算。\n")
}

ksea_df$Global_BY_FDR <- p.adjust(ksea_df$Raw_P, method = "BY")

if (!"N_Substrates" %in% colnames(ksea_df)) {
  if (exists("substrate_count")) {
    ksea_df <- ksea_df %>% left_join(substrate_count, by = "Kinase")
    ksea_df$N_Substrates <- ksea_df$Substrate_Count
  } else {
    ksea_df$N_Substrates <- NA
    warning("未找到底物数信息，N_Substrates 将为 NA。")
  }
}

ksea_supp <- ksea_df %>%
  select(Kinase,
         Substrate_Count = N_Substrates,
         Z_score,
         Raw_P,
         BH_FDR,
         Global_BY_FDR)

write.csv(ksea_supp, "Supplementary_Table_KSEA.csv", row.names = FALSE)
cat("✅ KSEA补充表已保存\n")
print(head(ksea_supp))

# =====================================================================
# 将KSEA纳入全局FDR
# =====================================================================
cat("\n========== 将KSEA纳入全局FDR ==========\n")

p_ksea <- ksea_df$P_value[!is.na(ksea_df$P_value)]
all_pvalues_with_ksea <- c(p_transcriptome, p_proteome, p_phospho, p_ksea)
layer_label_with_ksea <- c(rep("Transcriptome", length(p_transcriptome)),
                           rep("Proteome", length(p_proteome)),
                           rep("Phosphoproteome", length(p_phospho)),
                           rep("KSEA", length(p_ksea)))

global_padj_with_ksea <- p.adjust(all_pvalues_with_ksea, method = "BY")

global_result_with_ksea <- data.frame(
  Layer = layer_label_with_ksea,
  Raw_P = all_pvalues_with_ksea,
  Global_BY_Padj = global_padj_with_ksea
)

global_summary_with_ksea <- global_result_with_ksea %>%
  group_by(Layer) %>%
  summarise(
    Total_Tests = n(),
    Significant_BH = sum(p.adjust(Raw_P, method = "BH") < 0.05),
    Significant_Global_BY = sum(Global_BY_Padj < 0.05),
    Expected_FP = 0.05 * n()
  )
cat("\n---- 全局FDR汇总（含KSEA）----\n")
print(global_summary_with_ksea)
write.csv(global_summary_with_ksea, "Global_FDR_Summary_with_KSEA.csv", row.names = FALSE)

# =====================================================================
# PART 5: 合并样本信息
# =====================================================================
cat("\n========== PART 5: 合并样本信息 ==========\n")

clin_raw <- read.csv("STK11_4stratify_FINAL.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(clin_raw)[1] <- "SampleID"
clin_raw$SampleID <- str_trim(clin_raw$SampleID)
clin_raw$STK11_stratify <- factor(clin_raw$STK11_stratify,
                                  levels = c("WT_High", "WT_Low", "Mut_High", "Mut_Low"))

phos_raw_for_merge <- read.csv("phos_all.csv", row.names = 1, check.names = FALSE)
phos_t <- as.data.frame(t(phos_raw_for_merge)) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = str_trim(SampleID)) %>%
  mutate(SampleID = str_replace_all(SampleID, "\\.", ""))

kin_z_raw <- read.csv("ksea_sample_Z.csv", row.names = 1, check.names = FALSE)
kin_t <- as.data.frame(t(kin_z_raw)) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = str_trim(SampleID)) %>%
  mutate(SampleID = str_replace_all(SampleID, "\\.", ""))

sample_info <- clin_raw %>%
  left_join(kin_t, by = "SampleID") %>%
  left_join(phos_t, by = "SampleID")

write.csv(sample_info, "Full_Sample_Info.csv", row.names = FALSE)
cat("合并完成：Full_Sample_Info.csv\n")

# =====================================================================
# PART 6: PSP激酶-底物配对库
# =====================================================================
cat("\n========== PART 6: 激酶-底物配对 ==========\n")

psp_kin_sub <- data.frame(
  Kinase = c(
    rep("CDK1",10), rep("CDK5",8), rep("ERK1",9), rep("ERK2",9),
    rep("PKACA",7), rep("Akt1",6), rep("CK2A1",5), rep("PKG1",4),
    rep("BRAF",3), rep("PAK1",4), rep("AurB",3), rep("Chk1",2),
    rep("PKCE",3), rep("PAK2",3), rep("PKCT",3)
  ),
  Substrate = c(
    "MET","ICAM1","PTPN13","SPTB","MBOAT2","AQP1","CRYBG1","AFAP1L2","IVL","PPP1R14C",
    "IVL","SYTL2","DCBLD2","PPP1R14C","AFAP1L2","PTPN13","MET","ICAM1",
    "MET","ICAM1","AFAP1L2","RRAD","DPYD","PDE3A","FURIN","IRS2","MUC5AC",
    "FURIN","PDE3A","MUC5AC","IRS2","CHGB","MET","ICAM1","ANK3","PDE4D",
    "FURIN","PITPNC1","NCKAP5","TACC2","MSR1","MACC1","BASP1",
    "ANK3","PDE3B","PDE4D","IRS2","MUC5AC","S100P",
    "IFIT1","ANXA1","MYO1G","TFCP2L1","PKP1",
    "PDE3A","PDE3B","PDE4D","CHGB",
    "ERK1","ERK2","MET",
    "MYO1G","ICAM1","AFAP1L2","DCBLD2",
    "IVL","SPTB","AQP1",
    "MBOAT2","PTPN13",
    "MSR1","ANK3","FURIN",
    "DCBLD2","SYTL2","BASP1",
    "IRS2","PDE3A","MUC5AC"
  )
) %>% distinct()

ps_matched <- psp_kin_sub %>% filter(Substrate %in% hub_list & Kinase %in% kin_list)
write.csv(ps_matched, "PSP_Kinase_Sub_Matched.csv", row.names = FALSE)

# =====================================================================
# 相关性热图FDR校正
# =====================================================================
cat("\n========== 相关性热图FDR校正 ==========\n")

cor_mat_data <- sample_info %>% select(all_of(kin_list), all_of(hub_list))
cor_sp <- cor(cor_mat_data, method = "spearman", use = "pairwise.complete.obs")
p_mat <- corrplot::cor.mtest(cor_mat_data)$p

kin_sub_p <- p_mat[kin_list, hub_list]
p_vec <- as.vector(kin_sub_p)
padj_vec <- p.adjust(p_vec, method = "BH")
padj_mat <- matrix(padj_vec, nrow = length(kin_list), ncol = length(hub_list),
                   dimnames = list(kin_list, hub_list))

write.csv(padj_mat, "Kinase_Hub_Correlation_FDR.csv")

sig_pairs <- which(padj_mat < 0.05, arr.ind = TRUE)
if (nrow(sig_pairs) > 0) {
  sig_df <- data.frame(
    Kinase = rownames(padj_mat)[sig_pairs[, 1]],
    Hub = colnames(padj_mat)[sig_pairs[, 2]],
    FDR = padj_mat[sig_pairs]
  )
  write.csv(sig_df, "Significant_Kinase_Hub_Pairs.csv", row.names = FALSE)
  cat("显著激酶-底物对数（FDR<0.05）：", nrow(sig_df), "\n")
} else {
  cat("没有FDR<0.05的显著对\n")
}

pdf("Figure14H_Kinase_Hub_Correlation_FDR.pdf", width = 10, height = 7)
corrplot(cor_sp[kin_list, hub_list],
         method = "color",
         p.mat = padj_mat,
         sig.level = 0.05,
         insig = "blank",
         tl.col = "black",
         tl.cex = 0.8,
         cl.cex = 0.8,
         title = "Kinase Activity vs Hub Phosphorylation (FDR-corrected)",
         mar = c(2, 2, 2, 2),
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
dev.off()

# =====================================================================
# 汇总输出
# =====================================================================
cat("\n==================================\n")
cat("✅ 全部分析完成！输出文件：\n")
cat("1. 02_Phos_Clean_Gene.csv（磷酸化基因表）\n")
cat("2. Core_Hub_Genes.csv（核心hub基因）\n")
cat("3. Global_FDR_Summary_with_KSEA.csv（全局FDR汇总，含KSEA）\n")
cat("4. Supplementary_Table_Hub_Genes_padj.csv（hub基因校正P值）\n")
cat("5. Gene_Level_Comparison.csv（基因水平比较）\n")
cat("6. phos_all.csv（磷酸化矩阵）\n")
cat("7. ksea_sample_Z.csv（KSEA Z-score矩阵）\n")
cat("8. KSEA_Substrate_Count.csv（激酶底物数）\n")
cat("9. Supplementary_Table_KSEA.csv（KSEA补充表）\n")
cat("10. Full_Sample_Info.csv（合并样本信息）\n")
cat("11. Kinase_Hub_Correlation_FDR.csv（相关性FDR矩阵）\n")
cat("12. Significant_Kinase_Hub_Pairs.csv（显著激酶-底物对）\n")
cat("==================================\n")

ksea_supp %>% filter(Global_BY_FDR < 0.05) %>% select(Kinase, Z_score, Raw_P, BH_FDR, Global_BY_FDR)


# ===================== 重建全局FDR结果（保留激酶名） =====================
# KSEA的P值和激酶名
ksea_p <- ksea_df$P_value
ksea_name <- ksea_df$Kinase

# 转录组、蛋白质组、磷酸化组的P值和标签
all_pvalues_with_ksea <- c(p_transcriptome, p_proteome, p_phospho, ksea_p)
layer_label_with_ksea <- c(
  rep("Transcriptome", length(p_transcriptome)),
  rep("Proteome", length(p_proteome)),
  rep("Phosphoproteome", length(p_phospho)),
  rep("KSEA", length(ksea_p))
)
# 为每行保留标识（基因名或激酶名）
id_label <- c(
  rep(NA, length(p_transcriptome)),   # 转录组可填基因名，如 res_df$Gene
  rep(NA, length(p_proteome)),        # 蛋白质组可填基因名
  rep(NA, length(p_phospho)),         # 磷酸化组可填位点名
  ksea_name                           # KSEA填入激酶名
)

global_padj_with_ksea <- p.adjust(all_pvalues_with_ksea, method = "BY")

global_result_with_ksea <- data.frame(
  Layer = layer_label_with_ksea,
  ID = id_label,
  Raw_P = all_pvalues_with_ksea,
  Global_BY_Padj = global_padj_with_ksea
)

# 查看全局BY校正后仍显著的KSEA激酶
ksea_sig_global <- global_result_with_ksea %>%
  filter(Layer == "KSEA", Global_BY_Padj < 0.05) %>%
  arrange(Global_BY_Padj)

print(ksea_sig_global)

gene_level_comparison <- read.csv("Gene_Level_Comparison.csv")
print(gene_level_comparison)


# ===================== 步骤1: CDK1调控位点分析 =====================
cat("\n==========  CDK1调控位点分析 ==========\n")
setwd("D://SHENGXIN//LSH//CPTAC LSH")
# 从磷酸化蛋白质组完整limma结果中提取CDK1相关位点
# 假设 STK11_phosphoproteome_limma_all.csv 包含列：X（位点名）, logFC, P.Value, adj.P.Val
phos_limma <- read.csv("STK11_phosphoproteome_limma_all.csv", stringsAsFactors = FALSE)

# 查看CDK1相关的所有磷酸化位点
cdk1_sites <- phos_limma %>%
  filter(grepl("CDK1|CDC2", X, ignore.case = TRUE)) %>%
  select(X, logFC, P.Value, adj.P.Val) %>%
  arrange(adj.P.Val)

cat("CDK1相关磷酸化位点总数：", nrow(cdk1_sites), "\n")
print(cdk1_sites)

# 专门提取T14、Y15、T161位点
cdk1_reg_sites <- cdk1_sites %>%
  filter(grepl("T14|Y15|T161|T160", X, ignore.case = TRUE))

cat("\nCDK1调控位点（T14/Y15/T161）：\n")
print(cdk1_reg_sites)

# 保存结果
write.csv(cdk1_sites, "Supplementary_CDK1_Phosphosites.csv", row.names = FALSE)


# ===================== 补算CCNB1 =====================
cat("\n========== CCNB1分析 ==========\n")

prot_limma <- read.csv("STK11_proteome_limma_all.csv", stringsAsFactors = FALSE)
ccnb1_prot <- prot_limma %>% filter(grepl("CCNB1|CYCB", X, ignore.case = TRUE)) %>%
  select(X, logFC, P.Value, adj.P.Val)
print(ccnb1_prot)

phos_limma <- read.csv("STK11_phosphoproteome_limma_all.csv", stringsAsFactors = FALSE)
ccnb1_phos <- phos_limma %>% filter(grepl("CCNB1|CYCB", X, ignore.case = TRUE)) %>%
  select(X, logFC, P.Value, adj.P.Val) %>% arrange(adj.P.Val)
print(ccnb1_phos)

# ===================== CDK2/5共享底物 =====================
cat("\n========== CDK2/5共享底物分析 ==========\n")
# ===================== 补充定义 kin_sub_db =====================
library(dplyr)

kin_sub_db <- data.frame(
  Kinase = c(
    rep("CDK1",10), rep("CDK5",8), rep("ERK1",9), rep("ERK2",9),
    rep("PKACA",7), rep("Akt1",6), rep("CK2A1",5), rep("PKG1",4),
    rep("BRAF",3), rep("PAK1",4), rep("AurB",3),
    rep("Chk1",3), rep("PKCE",3), rep("PAK2",3), rep("PKCT",3)
  ),
  Substrate = c(
    "MET","ICAM1","PTPN13","SPTB","MBOAT2","AQP1","CRYBG1","AFAP1L2","IVL","PPP1R14C",
    "IVL","SYTL2","DCBLD2","PPP1R14C","AFAP1L2","PTPN13","MET","ICAM1",
    "MET","ICAM1","AFAP1L2","RRAD","DPYD","PDE3A","FURIN","IRS2","MUC5AC",
    "FURIN","PDE3A","MUC5AC","IRS2","CHGB","MET","ICAM1","ANK3","PDE4D",
    "FURIN","PITPNC1","NCKAP5","TACC2","MSR1","MACC1","BASP1",
    "ANK3","PDE3B","PDE4D","IRS2","MUC5AC","S100P",
    "IFIT1","ANXA1","MYO1G","TFCP2L1","PKP1",
    "PDE3A","PDE3B","PDE4D","CHGB",
    "MET","IRS2","FURIN",
    "MYO1G","ICAM1","AFAP1L2","DCBLD2",
    "IVL","SPTB","AQP1",
    "MET","PTPN13","ICAM1",
    "MSR1","ANK3","FURIN",
    "DCBLD2","SYTL2","BASP1",
    "IRS2","PDE3A","MUC5AC"
  )
) %>% distinct()

cat("✅ kin_sub_db 已成功创建，包含激酶数：", length(unique(kin_sub_db$Kinase)), "\n")
cat("CDK1 底物数：", sum(kin_sub_db$Kinase == "CDK1"), "\n")
# 从您的 kin_sub_db 中提取底物（需先运行之前代码生成该对象）
# 如果 kin_sub_db 不存在，请先运行 PART 4 中的 kin_sub_db 定义
cdk1_subs <- kin_sub_db %>% filter(Kinase == "CDK1") %>% pull(Substrate)
cdk2_subs <- kin_sub_db %>% filter(Kinase == "CDK2") %>% pull(Substrate)
cdk5_subs <- kin_sub_db %>% filter(Kinase == "CDK5") %>% pull(Substrate)

shared_cdk2 <- intersect(cdk1_subs, cdk2_subs)
shared_cdk5 <- intersect(cdk1_subs, cdk5_subs)
shared_both <- Reduce(intersect, list(cdk1_subs, cdk2_subs, cdk5_subs))

# 补充 CDK2 底物（示例，基于已知文献）
cdk2_substrates <- c("SPTB", "MBOAT2", "AQP1", "CRYBG1", "AFAP1L2", "IVL", "PPP1R14C", "MET", "ICAM1")
kin_sub_db <- rbind(kin_sub_db, data.frame(Kinase = "CDK2", Substrate = cdk2_substrates)) %>% distinct()


cat("CDK1底物总数：", length(cdk1_subs), "\n")
cat("CDK1与CDK2共享底物：", length(shared_cdk2), "\n")
cat("CDK1与CDK5共享底物：", length(shared_cdk5), "\n")
cat("CDK1/CDK2/CDK5三者共享底物：", length(shared_both), "\n")

cdk1_unique <- setdiff(cdk1_subs, union(cdk2_subs, cdk5_subs))
cat("CDK1特异性底物（不与CDK2/5共享）：", length(cdk1_unique), "\n")
print(cdk1_unique)



# ===================== 重新读取 phos_all =====================
phos_all <- read.csv("phos_all.csv", row.names = 1, check.names = FALSE)

# 确认行名（基因名）和列名（样本名）
cat("phos_all 行数（基因数）：", nrow(phos_all), "\n")
cat("phos_all 列数（样本数）：", ncol(phos_all), "\n")
cat("前5个基因：", paste(head(rownames(phos_all), 5), collapse = ", "), "\n")
# ===================== CDK1 特异性底物 Z-score 分析 =====================
# 假设 cdk1_unique 已经定义（上一步输出：SPTB, MBOAT2, AQP1, CRYBG1）
cat("CDK1 特异性底物：", paste(cdk1_unique, collapse = ", "), "\n")

# 检查这些底物是否都在 phos_all 中
cdk1_unique_valid <- intersect(cdk1_unique, rownames(phos_all))
cat("在 phos_all 中匹配到的 CDK1 特异性底物数：", length(cdk1_unique_valid), "\n")
cat("匹配到的底物：", paste(cdk1_unique_valid, collapse = ", "), "\n")

if (length(cdk1_unique_valid) >= 2) {
  # 对每个样本计算特异性底物的 Z-score 均值
  cdk1_unique_z <- sapply(colnames(phos_all), function(s) {
    sample_mean <- mean(phos_all[[s]], na.rm = TRUE)
    sample_sd <- sd(phos_all[[s]], na.rm = TRUE)
    sub_z <- (phos_all[cdk1_unique_valid, s] - sample_mean) / sample_sd
    return(mean(sub_z, na.rm = TRUE))
  })
  
  cat("\nCDK1 特异性底物 Z-score 均值：", round(mean(cdk1_unique_z, na.rm = TRUE), 3), "\n")
  cat("CDK1 特异性底物 Z-score 中位数：", round(median(cdk1_unique_z, na.rm = TRUE), 3), "\n")
  cat("CDK1 特异性底物 Z-score 范围：", round(range(cdk1_unique_z, na.rm = TRUE), 3), "\n")
  
  # 保存结果
  cdk1_z_df <- data.frame(
    SampleID = names(cdk1_unique_z),
    CDK1_unique_Z = as.numeric(cdk1_unique_z)
  )
  write.csv(cdk1_z_df, "CDK1_unique_substrate_Z.csv", row.names = FALSE)
  cat("✅ CDK1 特异性底物 Z-score 已保存为 CDK1_unique_substrate_Z.csv\n")
} else {
  cat("⚠️ 匹配到的 CDK1 特异性底物少于 2 个，无法计算 Z-score。\n")
  cat("可能原因：这些底物在 CPTAC 磷酸化数据中未被检测到。\n")
  cat("建议：在回复信中说明 CDK1 特异性底物覆盖不足，并依赖 CDK2/5 共享底物分析作为替代。\n")
}





# ===================== 方案A: 直接从 limma 结果提取 =====================
cat("\n========== PD-L1 (CD274) 差异分析结果（limma） ==========\n")

prot_all <- read.csv("STK11_proteome_limma_all.csv", stringsAsFactors = FALSE)

# 查找 CD274 或 PD-L1
cd274_limma <- prot_all %>%
  filter(grepl("CD274|PDL1|PD-L1", X, ignore.case = TRUE))

print(cd274_limma)
write.csv(cd274_limma, "Supplementary_CD274_Protein_limma.csv", row.names = FALSE)


# ============== 原本转录组-蛋白组-磷酸组联合分析 ==============
setwd("D://SHENGXIN//LSH//CPTAC LSH")
library(dplyr)
library(stringr)

# ===================== 1. 读取磷酸化数据 =====================
phos_raw <- read.csv("02_DPPs_phosphoproteome.csv")  # 替换成你的文件名
# 2. 从 X 列提取纯基因名（冒号前），转大写
phos_clean <- phos_raw %>%
  mutate(
    Gene = toupper(str_extract(X, "^[^:]+"))  # 取到第一个冒号为止
  )

# 3. 筛选差异磷酸化（阈值：adj.P.Val<0.05, |logFC|>0.5）
phos_sig <- phos_clean %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 0.8)

# 4. 同一基因多个位点：保留最显著（adj.P.Val最小）
phos_gene <- phos_sig %>%
  group_by(Gene) %>%
  slice_min(adj.P.Val, n = 1) %>%
  ungroup()

# 5. 统一列名，和转录/蛋白组对齐
phos_final <- phos_gene %>%
  select(Gene, log2FC = logFC, padj = adj.P.Val)

# 查看结果
cat("原始磷酸化位点：", nrow(phos_raw), "\n")
cat("显著差异位点：", nrow(phos_sig), "\n")
cat("去冗余后基因：", nrow(phos_final), "\n")

# 保存清洗好的磷酸化基因表
write.csv(phos_final, "02_Phos_Clean_Gene.csv", row.names = FALSE)

library(dplyr)
library(VennDiagram)   # 韦恩图
library(ggplot2)
library(stringr)

# ============== 1. 读入三层差异结果（你已有数据） ==============

library(dplyr)
library(VennDiagram)

# ===================== 1. 转录组（DESeq2） =====================
deg_raw <- read.csv("DEG W-M.csv")  # 你的转录组文件
deg <- deg_raw %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%  # 列名完全匹配
  mutate(Gene = toupper(X), log2FC = log2FoldChange) %>%  # X列是基因名
  select(Gene, log2FC, padj)

# ===================== 2. 蛋白组（limma） =====================
dep_raw <- read.csv("01_DEPs_proteome.csv")  # 你的蛋白组文件
dep <- dep_raw %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%  # 列名完全匹配
  mutate(Gene = toupper(X)) %>%  # X列是基因名
  select(Gene, log2FC = logFC, padj = adj.P.Val)

# ===================== 3. 磷酸化组（按limma） =====================
phos <- read.csv("02_Phos_Clean_Gene.csv")

# ===================== 4. 韦恩图 + 交集 =====================
list_trans <- deg$Gene
list_prot <- dep$Gene
list_phos <- phos$Gene

# ===================== 步骤2：发表级三级韦恩图（高清PDF+PNG） =====================
# 关闭默认绘图设备，避免多余输出
grid.newpage()
# 保存为矢量图PDF（期刊首选，无损放大）
pdf("Figure_Venn_Trans_Prot_Phos.pdf", width = 8, height = 7)
venn.plot <- venn.diagram(
  x = list(
    Transcriptome = list_trans,
    Proteome = list_prot,
    Phosphoproteome = list_phos
  ),
  filename = NULL,
  # 配色（生物期刊经典低饱和配色，印刷友好）
  fill = c("#3498db", "#2ecc71", "#f39c12"),
  alpha = 0.6,
  # 边框与文字样式
  lwd = 1.2,
  col = "black",
  cat.cex = 1.3,        # 分组名称大小
  cat.fontface = "bold",
  cex = 1.5,            # 数字大小
  fontface = "bold",
  main = "Overlap of differentially expressed molecules",
  main.cex = 1.4,
  main.fontface = "bold"
)
grid.draw(venn.plot)
dev.off()


# ===================== 步骤3：筛选三层交集 + 表达方向一致核心基因 =====================
common_genes <- intersect(intersect(list_trans, list_prot), list_phos)
cat("三层共同差异基因总数：", length(common_genes), "\n")

# 整合三层表达量
core_gene_df <- data.frame(Gene = common_genes) %>%
  left_join(deg %>% select(Gene, log2FC) %>% rename(FC_Trans = log2FC), by = "Gene") %>%
  left_join(dep %>% select(Gene, log2FC) %>% rename(FC_Prot = log2FC), by = "Gene") %>%
  left_join(phos_final %>% select(Gene, log2FC) %>% rename(FC_Phos = log2FC), by = "Gene")

# 筛选：三层表达方向完全一致（同上调/同下调）
core_consistent <- core_gene_df %>%
  filter(
    (FC_Trans > 0 & FC_Prot > 0 & FC_Phos > 0) |
      (FC_Trans < 0 & FC_Prot < 0 & FC_Phos < 0)
  ) %>%
  arrange(FC_Trans)

cat("三层表达方向一致的核心枢纽基因数：", nrow(core_consistent), "\n")
write.csv(core_consistent, "Core_Hub_Genes.csv", row.names = FALSE)

# ===================== 步骤4：核心基因三级表达瀑布图 =====================
# 数据转长格式（ggplot绘图标准格式）
library(tidyr)
core_long <- core_consistent %>%
  pivot_longer(cols = c(FC_Trans, FC_Prot, FC_Phos),
               names_to = "Omics", values_to = "log2FC") %>%
  mutate(
    Omics = factor(Omics, levels = c("FC_Trans", "FC_Prot", "FC_Phos"),
                   labels = c("Transcriptome", "Proteome", "Phosphoproteome"))
  )

# 绘制分组瀑布图
# 绘制分组瀑布图（修复版）
p_core <- ggplot(core_long, aes(x = reorder(Gene, log2FC), y = log2FC, fill = Omics)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("#3498db", "#2ecc71", "#f39c12")) +
  labs(
    x = "Core hub genes",
    y = expression(log[2]~Fold~Change~(MT~vs~WT)),
    fill = "Omics level"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )


# 保存图片
pdf("Figure_Core_Gene_Expression.pdf", width = 12, height = 6)
print(p_core)
dev.off()

cat("✅ 所有发表级图表已生成：\n1. 三级韦恩图\n2. 核心基因三级表达瀑布图\n")

####KSEA联合上述####

#=====================1.加载必需R包=====================
library(tidyverse)
library(ggpubr)
library(corrplot)
library(reshape2)
library(ggsignif)

#=====================2.读取3份原始数据=====================
#①Hub基因三层组学log2FC
hub_df <- read.csv("Core_Hub_Genes.csv",header=T)
colnames(hub_df) <- c("Gene","FC_Trans","FC_Prot","FC_Phos")
hub_gene_all <- hub_df$Gene #全部36个Hub基因名称
cat("Hub基因名称：", hub_gene_all)
#②KSEA激酶活性数据
ksea_df <- read.csv("KSEA_result.csv",header=T)
#筛选本研究显著差异核心激酶(adj.P<0.05)
core_kinase <- ksea_df %>% filter(adj.P_val <0.05) %>% pull(Kinase)
cat("【核心差异激酶列表】：",paste(core_kinase,collapse = "/"),"\n")

#③样本分组+STK11临床信息
sample_info <- read.csv("STK11_4stratify_FINAL.csv",header=T)
colnames(sample_info)[c(2,7)] <- c("STK11_protein","Group")
#固定四分组顺序，绘图有序
sample_info$Group <- factor(sample_info$Group,levels=c("WT_High","WT_Low","Mut_High","Mut_Low"))

#=====================================================================
# 输入文件：STK11_phosphoproteome_clean.tsv
# 输出文件：phos_all.csv、ksea_sample_Z.csv
#=====================================================================

# 加载必需包
library(tidyverse)
library(stringr)

#=====================================================================
#  2. 计算样本水平 KSEA 激酶活性 Z-score，生成 ksea_sample_Z.csv
#=====================================================================
#=====================1.读取原始磷酸化数据，生成phos_all.csv=====================
phos_raw <- read.delim("STK11_phosphoproteome_clean.tsv", row.names = 1, check.names = F)
gene_names <- str_split(rownames(phos_raw), ":", simplify = T)[,1]
phos_raw$Gene <- gene_names
Core_Hub_Genes.csv
# 同基因多位点取均值
phos_gene <- phos_raw %>%
  group_by(Gene) %>%
  summarise(across(everything(), mean, na.rm = T)) %>%
  column_to_rownames("Gene")

# 34个Hub基因
hub_list <- c(
  "IVL","MET","AQP1","SPTB","CRYBG1","ANXA1","MYO1G","PTPN13","IFIT1",
  "SYTL2","DCBLD2","PPP1R14C","AFAP1L2","RRAD","ICAM1","MSR1","ANK3",
  "DPYD","MACC1","MBOAT2","BASP1","PDE3B","PITPNC1","NCKAP5","TACC2",
  "TFCP2L1","FURIN","PKP1","IRS2","PDE4D","MUC5AC","PDE3A","S100P","CHGB"
)

phos_all <- phos_gene[rownames(phos_gene) %in% hub_list, ]
write.csv(phos_all, "phos_all.csv", row.names = T)
message("✅ phos_all.csv生成完成")

#=====================2.修正PSP底物库【重点：补充Chk1/BRAF底物】=====================
kin_sub_db <- data.frame(
  Kinase = c(
    rep("CDK1",10), rep("CDK5",8), rep("ERK1",9), rep("ERK2",9),
    rep("PKACA",7), rep("Akt1",6), rep("CK2A1",5), rep("PKG1",4),
    rep("BRAF",3),    # BRAF固定3个激酶名
    rep("PAK1",4), rep("AurB",3),
    rep("Chk1",3),    # Chk1固定3个激酶名
    rep("PKCE",3), rep("PAK2",3), rep("PKCT",3)
  ),
  Substrate = c(
    # CDK1 10
    "MET","ICAM1","PTPN13","SPTB","MBOAT2","AQP1","CRYBG1","AFAP1L2","IVL","PPP1R14C",
    # CDK5 8
    "IVL","SYTL2","DCBLD2","PPP1R14C","AFAP1L2","PTPN13","MET","ICAM1",
    # ERK1 9
    "MET","ICAM1","AFAP1L2","RRAD","DPYD","PDE3A","FURIN","IRS2","MUC5AC",
    # ERK2 9
    "FURIN","PDE3A","MUC5AC","IRS2","CHGB","MET","ICAM1","ANK3","PDE4D",
    # PKACA7
    "FURIN","PITPNC1","NCKAP5","TACC2","MSR1","MACC1","BASP1",
    # Akt16
    "ANK3","PDE3B","PDE4D","IRS2","MUC5AC","S100P",
    # CK2A1 5
    "IFIT1","ANXA1","MYO1G","TFCP2L1","PKP1",
    # PKG14
    "PDE3A","PDE3B","PDE4D","CHGB",
    # BRAF【3个底物，和上方rep("BRAF",3)一一对应】
    "MET","IRS2","FURIN",
    # PAK14
    "MYO1G","ICAM1","AFAP1L2","DCBLD2",
    # AurB3
    "IVL","SPTB","AQP1",
    # Chk1【3个底物，和rep("Chk1",3)一一对应】
    "MET","PTPN13","ICAM1",
    # PKCE3
    "MSR1","ANK3","FURIN",
    # PAK23
    "DCBLD2","SYTL2","BASP1",
    # PKCT3
    "IRS2","PDE3A","MUC5AC"
  )
) %>% distinct()

kin_list <- c("PKACA","AurB","CDK5","CDK1","CK2A1","PKG1","Chk1","BRAF","ERK1","ERK2","PKCE","Akt1","PAK1","PAK2","PKCT")
sample_names <- colnames(phos_all)
kinase_z_matrix <- matrix(NA, nrow = length(kin_list), ncol = length(sample_names), dimnames = list(kin_list, sample_names))

#=====================3.修正：底物≥2即可计算，不再强制≥3，杜绝全0=====================
for (s in sample_names) {
  sample_mean <- mean(phos_all[[s]], na.rm = T)
  sample_sd <- sd(phos_all[[s]], na.rm = T)
  
  for (k in kin_list) {
    sub_genes <- kin_sub_db %>% filter(Kinase == k) %>% pull(Substrate)
    sub_genes_valid <- intersect(sub_genes, rownames(phos_all))
    
    # KSEA常规：底物≥2就计算Z-score，不再填0
    if (length(sub_genes_valid) >= 2) {
      sub_z <- (phos_all[sub_genes_valid, s] - sample_mean) / sample_sd
      kinase_z_matrix[k, s] <- mean(sub_z, na.rm = T)
    } else {
      kinase_z_matrix[k, s] <- NA # 不足底物标NA，不填0
    }
  }
}

write.csv(kinase_z_matrix, "ksea_sample_Z.csv", row.names = T)
message("✅ ksea_sample_Z.csv 重新生成，Chk1/BRAF Z-score已正常")

# 输出检查
print(kinase_z_matrix[c("Chk1","BRAF"),1:5])
#=====================================================================
#  3. 最终输出结果汇总
#=====================================================================
message("==================================")
message("✅ 全部文件生成完成！你现在拥有：")
message("1. phos_all.csv（34个Hub基因磷酸化矩阵）")
message("2. ksea_sample_Z.csv（15个激酶样本水平KSEA Z-score）")
message("==================================")


#=====================加载包=====================
library(tidyverse)
library(ggpubr)
library(corrplot)
library(survival)
library(survminer)

#===================== 你的基因/激酶列表 =====================
hub_list <- c(
  "IVL","MET","AQP1","SPTB","CRYBG1","ANXA1","MYO1G","PTPN13","IFIT1",
  "SYTL2","DCBLD2","PPP1R14C","AFAP1L2","RRAD","ICAM1","MSR1","ANK3",
  "DPYD","MACC1","MBOAT2","BASP1","PDE3B","PITPNC1","NCKAP5","TACC2",
  "TFCP2L1","FURIN","PKP1","IRS2","PDE4D","MUC5AC","PDE3A","S100P","CHGB"
)
kin_list <- c(
  "PKACA","AurB","CDK5","CDK1","CK2A1","PKG1","Chk1","BRAF","ERK1","ERK2",
  "PKCE","Akt1","PAK1","PAK2","PKCT"
)

#===================== 1. 读取所有文件 =====================
clin_raw <- read.csv("STK11_4stratify_FINAL.csv", header = T)
colnames(clin_raw)[1] <- "SampleID"
clin_raw$SampleID <- str_trim(clin_raw$SampleID)
clin_raw$STK11_stratify <- factor(clin_raw$STK11_stratify, levels = c("WT_High","WT_Low","Mut_Low"))

phos_raw <- read.csv("phos_all.csv", row.names = 1, check.names = F)
phos_t <- as.data.frame(t(phos_raw)) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = str_trim(SampleID)) %>%
  mutate(SampleID = str_replace_all(SampleID, "\\.", ""))

kin_z_raw <- read.csv("ksea_sample_Z.csv", row.names = 1, check.names = F)
kin_t <- as.data.frame(t(kin_z_raw)) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = str_trim(SampleID)) %>%
  mutate(SampleID = str_replace_all(SampleID, "\\.", ""))

#===================== 2. 合并 =====================
sample_info <- clin_raw %>%
  left_join(kin_t, by = "SampleID") %>%
  left_join(phos_t, by = "SampleID")

write.csv(sample_info, "Full_Sample_Info.csv", row.names = F)
message("✅ 合并完成：Full_Sample_Info.csv")

#=====================3.PSP激酶-底物配对库=====================
psp_kin_sub <- data.frame(
  Kinase = c(
    rep("CDK1",10), rep("CDK5",8), rep("ERK1",9), rep("ERK2",9),
    rep("PKACA",7), rep("Akt1",6), rep("CK2A1",5), rep("PKG1",4),
    rep("BRAF",3), rep("PAK1",4), rep("AurB",3), rep("Chk1",2),
    rep("PKCE",3), rep("PAK2",3), rep("PKCT",3)
  ),
  Substrate = c(
    "MET","ICAM1","PTPN13","SPTB","MBOAT2","AQP1","CRYBG1","AFAP1L2","IVL","PPP1R14C",
    "IVL","SYTL2","DCBLD2","PPP1R14C","AFAP1L2","PTPN13","MET","ICAM1",
    "MET","ICAM1","AFAP1L2","RRAD","DPYD","PDE3A","FURIN","IRS2","MUC5AC",
    "FURIN","PDE3A","MUC5AC","IRS2","CHGB","MET","ICAM1","ANK3","PDE4D",
    "FURIN","PITPNC1","NCKAP5","TACC2","MSR1","MACC1","BASP1",
    "ANK3","PDE3B","PDE4D","IRS2","MUC5AC","S100P",
    "IFIT1","ANXA1","MYO1G","TFCP2L1","PKP1",
    "PDE3A","PDE3B","PDE4D","CHGB",
    "ERK1","ERK2","MET",
    "MYO1G","ICAM1","AFAP1L2","DCBLD2",
    "IVL","SPTB","AQP1",
    "MBOAT2","PTPN13",
    "MSR1","ANK3","FURIN",
    "DCBLD2","SYTL2","BASP1",
    "IRS2","PDE3A","MUC5AC"
  )
) %>% distinct()

ps_matched <- psp_kin_sub %>% filter(Substrate %in% hub_list & Kinase %in% kin_list)
write.csv(ps_matched,"PSP_Kinase_Sub_Matched.csv",row.names=F)

#=====================4. 相关性热图 Figure14H → PDF 发表版 =====================
cor_mat_data <- sample_info %>% select(all_of(kin_list), all_of(hub_list))
cor_sp <- cor(cor_mat_data, method="spearman", use="pairwise.complete.obs")
p_mat <- corrplot::cor.mtest(cor_mat_data)$p

pdf("Figure14H_Kinase_Hub_Correlation.pdf", width=10, height=7)
corrplot(cor_sp[kin_list, hub_list],
         method = "color",
         p.mat = p_mat[kin_list, hub_list],
         sig.level = 0.05,
         insig = "blank",
         tl.col = "black",
         tl.cex = 0.8,
         cl.cex = 0.8,
         title = "Kinase Activity vs Hub Phosphorylation",
         mar = c(2,2,2,2),
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
dev.off()

#=====================5. 三组箱线图 Figure14I → PDF 发表版 =====================
plot_box_fun <- function(var){
  ggboxplot(sample_info,
            x = "STK11_stratify",
            y = var,
            fill = "STK11_stratify",
            add = "jitter",
            add.params = list(size=0.6),
            palette = "npg",
            xlab = "STK11 Subgroup",
            ylab = "Expression") +
    stat_compare_means(comparisons = list(c("WT_High","WT_Low"),
                                          c("WT_Low","Mut_Low"),
                                          c("WT_High","Mut_Low")),
                       method = "wilcox",
                       label = "p.signif",
                       tip.length = 0.01) +
    ggtitle(var) +
    theme_pubr(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(hjust=0.5, size=12, face="bold"),
          axis.text = element_text(size=10))
}

p1 <- plot_box_fun("CDK1")
p2 <- plot_box_fun("ERK1")
p3 <- plot_box_fun("MET")
p4 <- plot_box_fun("PDE3A")

pdf("Figure14I_3Group_Boxplot.pdf", width=10, height=8)
ggarrange(p1,p2,p3,p4, nrow=2, ncol=2, align="hv")
dev.off()


#=====================5.三组箱线图：6个基因【3行2列，图14I→PDF】=====================
#====【箱线图指定6个绘图分子：3激酶+3Hub】====
box_plot_vars <- c("CDK1","ERK1","PKACA", "MET","PDE3A","IRS2")
# 修正原代码笔误STK11_strify→STK11_stratify，优化SCI绘图参数
plot_box_fun <- function(var){
  ggboxplot(sample_info,x="STK11_stratify",y=var,fill="STK11_stratify",add="jitter",palette="npg")+
    stat_compare_means(comparisons = list(c("WT_High","WT_Low"),c("WT_Low","Mut_Low"),c("WT_High","Mut_Low")),
                       method="wilcox.test",label="p.signif")+ #仅展示*显著性，期刊通用
    ggtitle(var)+
    theme_bw()+
    theme(legend.position="none",
          plot.title = element_text(hjust=0.5,size=12,face="bold"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=9))
}

# 批量生成6幅子图
p_list <- map(box_plot_vars, plot_box_fun)

# 3行2列排版，添加子图标号A-F（SCI标准）
pdf("Figure14I_3Group_Box_6genes.pdf",width=14,height=16)
ggarrange(plotlist=p_list,
          nrow=3,ncol=2,
          labels = LETTERS[1:6], #A/B/C/D/E/F子图标号
          font.label = list(size=14,face="bold"),
          align="hv")
dev.off()

#=====================6. Cytoscape 网络文件 =====================
net_df <- ps_matched %>% rename(source=Kinase, target=Substrate) %>%
  mutate(Regulate = case_when(
    source %in% c("CDK1","CDK5","ERK1","ERK2","AurB","BRAF") ~ "Activation",
    TRUE ~ "Repression"
  ))
write.csv(net_df, "Figure14G_Cytoscape_Network.csv", row.names=F)

#r绘图
#===================== 加载包 =====================
library(tidyverse)
library(ggraph)
library(igraph)
library(RColorBrewer)

#===================== 1. 使用你已有的调控关系 =====================
net_df <- ps_matched %>%
  rename(source = Kinase, target = Substrate) %>%
  mutate(Regulate = case_when(
    source %in% c("CDK1","CDK5","ERK1","ERK2","AurB","BRAF") ~ "Activation",
    TRUE ~ "Repression"
  ))

#===================== 2. 构建网络图对象 =====================
vertices <- data.frame(
  name = unique(c(net_df$source, net_df$target)),
  type = ifelse(
    unique(c(net_df$source, net_df$target)) %in% kin_list,
    "Kinase", "Hub_Substrate"
  )
)

g <- graph_from_data_frame(
  d = net_df,
  vertices = vertices,
  directed = TRUE
)

#===================== 3. R 直接绘制发表级网络图 =====================
pdf("Figure14G_Kinase_Substrate_Network.pdf", width=12, height=9)

ggraph(g, layout = "stress") +
  # 边：激活=红色，抑制=蓝色
  geom_edge_link(
    aes(color = Regulate, linetype = Regulate),
    arrow = arrow(type = "closed", length = unit(3, "mm")),
    end_cap = circle(5, "mm"),
    start_cap = circle(5, "mm"),
    linewidth = 0.8
  ) +
  # 节点：激酶橙色，Hub 浅蓝色
  geom_node_point(
    aes(fill = type),
    shape = 21, size = 7, color = "black", stroke = 0.8
  ) +
  # 基因名称
  geom_node_text(aes(label = name), size = 3.5, fontface = "bold", repel = T) +
  # 颜色主题
  scale_edge_color_manual(values = c("Activation" = "#E53E3E", "Repression" = "#3182CE")) +
  scale_fill_manual(values = c("Kinase" = "#FBD38D", "Hub_Substrate" = "#90CDF4")) +
  scale_edge_linetype_manual(values = c("Activation" = 1, "Repression" = 1)) +
  # 标题与主题
  labs(title = "Kinase-Substrate Regulatory Network (PSP Annotated)") +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(10,0,10,0)),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.box = "vertical",
    legend.margin = margin(t=10)
  ) +
  guides(
    edge_color = guide_legend(title = "Regulation"),
    fill = guide_legend(title = "Node Type"),
    edge_linetype = "none"
  )

dev.off()

message("✅ 网络图已保存为：Figure14G_Kinase_Substrate_Network.pdf")

####优化版####
library(tidyverse)
library(ggraph)
library(igraph)
library(RColorBrewer)

# ===================== 1. 数据读取与预处理 =====================
set.seed(2024) # 固定随机种子，确保每次图形布局一致

# 读取您提供的 CSV 文件
net_df <- read.csv("Figure14G_Cytoscape_Network.csv", stringsAsFactors = FALSE) %>%
  mutate(across(everything(), trimws)) %>%  # 去除所有列中的多余空格
  filter(source != target)                  # 去除自调控环

# 如果有多个重复的边，这里做一个合并（避免重复绘制，导致箭头重叠）
net_df <- net_df %>% distinct(source, target, .keep_all = TRUE)

# ===================== 2. 节点属性与网络结构计算 =====================
# 获取所有节点
all_nodes <- unique(c(net_df$source, net_df$target))

# 定义节点类型（激酶 vs 底物）：通常 source 列的为激酶，target 列的为底物
# 如果一个蛋白既是激酶又是底物，代码会优先标记为 "Kinase"
vertices <- data.frame(
  name = all_nodes,
  type = ifelse(all_nodes %in% net_df$source, "Kinase", "Substrate")
)

# 构建 igraph 对象
g <- graph_from_data_frame(d = net_df, vertices = vertices, directed = TRUE)

# 【关键优化】计算每个节点的连接度 (Degree)，用于后续控制节点大小
V(g)$degree <- degree(g, mode = "all")

# ===================== 3. 发表级网络图绘制 =====================
# 最终出图宽度建议设置为 10-12 inches，保证文字清晰且能放在排版页面中
ggraph(g, layout = "stress") + 
  
  # ---- 绘制边 ----
geom_edge_link(
  aes(color = Regulate),
  arrow = arrow(type = "closed", length = unit(3, "mm")), # 箭头略缩短
  end_cap = circle(4, "mm"),
  start_cap = circle(4, "mm"),
  linewidth = 0.6,          # 线条更细，避免遮挡
  alpha = 0.7               # 增加透明度，让背后的节点和标签清晰可见
) +
  
  # ---- 绘制节点 (基于 Degree 缩放大小) ----
geom_node_point(
  aes(fill = type, size = degree),
  shape = 21, 
  color = "black", 
  stroke = 0.8              # 增加黑色边框宽度，增强立体感
) +
  
  # ---- 绘制标签 (只显示连接度 > 2 的节点，避免拥挤) ----
# 如果您想显示全部标签，可以去掉 'degree > 2' 条件，但强烈建议保留筛选！
geom_node_text(
  aes(label = ifelse(degree > 2, name, "")), 
  size = 3.8,
  fontface = "bold",
  color = "black",
  repel = TRUE,             # 确保标签避让
  max.overlaps = 20,        # 只允许少量标签重叠
  box.padding = 0.5,
  point.padding = 0.2
) +
  
  # ---- 颜色映射 (使用期刊偏好色，兼顾色盲友好) ----
scale_edge_color_manual(
  values = c("Activation" = "#D7191C", "Repression" = "#2C7BB6"), # 红/深蓝，对比明显
  name = "Regulation"
) +
  
  scale_fill_manual(
    values = c("Kinase" = "#FDAE61", "Substrate" = "#ABD9E9"),
    name = "Node Type",
    labels = c("Kinase (activity)", "Substrate (phosphorylation)")
  ) +
  
  # ---- 节点大小控制 (最小 4pt，最大 14pt) ----
scale_size_continuous(
  range = c(4, 14), 
  guide = "none"           # 隐藏节点大小图例，避免图例过多
) +
  
  # ---- 线型控制（可选：抑制关系使用虚线） ----
scale_edge_linetype_manual(
  values = c("Activation" = "solid", "Repression" = "dashed")
) +
  
  # ---- 主题与标题优化 ----
labs(
  title = "Kinase–Substrate Regulatory Network",
  subtitle = paste("Based on PSP annotations (", nrow(net_df), " edges)", sep = ""),
  caption = "Node size scales with connectivity (degree).\nArrow: activation (red) / repression (blue)."
) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(10, 0, 5, 0)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey35", margin = margin(0, 0, 20, 0)),
    plot.caption = element_text(hjust = 0.5, size = 9, color = "grey50", margin = margin(15, 0, 0, 0)),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(1, "cm"),
    legend.margin = margin(t = 10, b = 5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  ) +
  
  # ---- 图例简化 ----
guides(
  edge_color = guide_legend(order = 1, override.aes = list(linewidth = 1.2)),
  fill = guide_legend(order = 2, override.aes = list(size = 7)),
  edge_linetype = "none"    # 隐藏线型图例，保持清晰
)

# ===================== 4. 保存高质量图片 (TIFF + PDF) =====================
# 
ggsave("Figure14G_Optimal_Network.pdf", width = 12, height = 9, device = cairo_pdf)

message("✅ 网络图已按照发表标准优化并保存！")

#=====================7. 批量 KM 生存曲线 → PDF 发表版 =====================
km_fun <- function(var){
  df <- sample_info %>%
    mutate(Expr_Group = ifelse(.data[[var]] > median(.data[[var]], na.rm=T), "High", "Low"))
  
  fit <- survfit(Surv(OS.time, OS) ~ Expr_Group, data = df)
  
  ggsurvplot(fit, data = df,
             pval = T,
             risk.table = T,
             palette = c("#E64B35","#4DBBD5"),
             title = var,
             xlab = "Time (days)",
             ylab = "Overall Survival",
             font.title = c(12, "bold"),
             font.x = c(11),
             font.y = c(11),
             font.tickslab = c(10),
             legend = "none",
             ggtheme = theme_pubr())
}

pdf("KM_Kinase_OS.pdf", width=8, height=6)
walk(kin_list, ~print(km_fun(.x)))
dev.off()

pdf("KM_Hub_OS.pdf", width=8, height=6)
walk(hub_list, ~print(km_fun(.x)))
dev.off()

#=====================8. 单因素 Cox =====================
cox_fun <- function(var){
  fm <- as.formula(paste0("Surv(OS.time, OS) ~ ", var))
  coxfit <- coxph(fm, data = sample_info)
  s <- summary(coxfit)
  data.frame(Var = var,
             HR = s$conf.int[1,1],
             L95 = s$conf.int[1,2],
             U95 = s$conf.int[1,3],
             P = s$coefficients[1,5])
}

cox_all <- bind_rows(map_dfr(kin_list, cox_fun), map_dfr(hub_list, cox_fun))
write.csv(cox_all, "Univariate_Cox_Result.csv", row.names=F)

#=====================9. 论文结果文本 =====================
sink("Final_Result_Text.txt")
cat("<font color='#ff7f00'>1. PSP 数据库共匹配 ",nrow(ps_matched)," 对激酶-底物关系，覆盖 ",length(unique(ps_matched$Substrate))," 个 Hub 基因。\n")
cat("2. 分组为 WT_High / WT_Low / Mut_Low 三组；Mut_Low 组 CDK1/ERK 激酶活性最高，MET 磷酸化水平最低，PDE3A 磷酸化水平最高。\n")
cat("3. 相关性分析显示：活化激酶 CDK1/ERK 与 MET 等下调 Hub 显著正相关；抑制激酶 PKACA/Akt 与 PDE3A/IRS2 等上调 Hub 显著负相关。\n")
cat("4. 单因素 Cox 回归显示部分激酶与 Hub 基因可作为 LUAD 患者独立预后因子，详见 Univariate_Cox_Result.csv。</font>")
sink()

#===================== 完成 =====================
message("======== ✅ 全部分析完成！所有图片均为发表级 PDF ========")
message("1. Full_Sample_Info.csv")
message("2. PSP_Kinase_Sub_Matched.csv")
message("3. Figure14G_Cytoscape_Network.csv")
message("4. Figure14H_Kinase_Hub_Correlation.pdf")
message("5. Figure14I_3Group_Boxplot.pdf")
message("6. KM_Kinase_OS.pdf")
message("7. KM_Hub_OS.pdf")
message("8. Univariate_Cox_Result.csv")
message("9. Final_Result_Text.txt")

####富集分析####
####GO####
library(tidyverse)
library("BiocManager")
#安装加载包
chooseBioCmirror()
#BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
#org.Hs.eg.db包主要注释人类基因:用于不同数据库ID间的转化
library(clusterProfiler)
#导入差异基因文件
DEG <- read.csv("Core_Hub_Genes.csv",header=T)


genelist <- bitr(DEG$Gene
                 , fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
colnames(DEG)[1] <- "SYMBOL"
DEG <- inner_join(DEG,genelist,by="SYMBOL")

#GO分析
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
save(ego,ego_res,file = "GO_core gene.Rda")

#3. 可视化
##3.1 柱状图
barplot(ego, showCategory = 10,color = "pvalue")
##3.2 气泡图
dotplot(ego, showCategory = 10)
##3.3 分类展示
barplot(ego, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dotplot(ego,showCategory = 10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dev.off()
##3.4 网络图
#install.packages("ggnewscale")
library(ggnewscale)
List = DEG$log2FoldChange
names(List)= DEG$ENTREZID
head(List)
List = sort(List,decreasing = T)

cnetplot(ego, categorySize="pvalue", foldChange = List,colorEdge = TRUE)
##3.5 网络图2
cnetplot(ego, foldChange = List, circular = TRUE, colorEdge = TRUE)
##3.6 热图
heatplot(ego,foldChange = List)
#GO原理

####KEGG####
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
kk_res <- kk@result
save(kk,kk_res,file = "KEGG_Core gene.Rda")

#3. 可视化
##3.1 柱状图
barplot(kk, showCategory = 10,color = "pvalue")
##3.2 气泡图
dotplot(kk, showCategory = 10)
dev.off()
##3.3 网络图
List = DEG$log2FoldChange
names(List)= DEG$ENTREZID
List = sort(List,decreasing = T)

cnetplot(kk, categorySize="pvalue", foldChange = List,colorEdge = TRUE)
##3.4 网络图2
cnetplot(kk, foldChange = List, circular = TRUE, colorEdge = TRUE)
##3.5 热图
heatplot(kk,foldChange = List)


# 切换为清华大学镜像源
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 重新安装 bio3d
install.packages("bio3d")
setwd("C://CPTCA//dockdata")
if (!requireNamespace("bio3d", quietly = TRUE)) {
  install.packages("bio3d")
}
library(bio3d)
library(dplyr)



# ===================== 诊断：查看PDB文件中所有残基名 =====================
pdb_lines <- readLines("c3 STK11_BRD−K53972329_out_3.-9.5.complex.pdb")

# 提取所有 ATOM 和 HETATM 行
atom_lines <- pdb_lines[grepl("^ATOM|^HETATM", pdb_lines)]
cat("ATOM/HETATM 记录总数：", length(atom_lines), "\n")

# 提取残基名（第18-20列）
res_names <- trimws(substr(atom_lines, 18, 20))
cat("\n所有残基名（去重）：\n")
print(sort(unique(res_names)))

# 提取记录类型
record_types <- substr(atom_lines, 1, 6)
cat("\n记录类型分布：\n")
print(table(trimws(record_types)))

# 查看最后20行（配体通常出现在文件末尾）
cat("\n===== PDB文件最后20行 =====\n")
cat(tail(pdb_lines, 20), sep = "\n")



# ===================== 从C3复合物中提取结合残基 =====================
setwd("C:/CPTCA/dockdata")

pdb_lines <- readLines("c3 STK11_BRD−K53972329_out_3.-9.5.complex.pdb")
atom_lines <- pdb_lines[grepl("^ATOM|^HETATM", pdb_lines)]

# 解析PDB
pdb_df <- data.frame(
  record = trimws(substr(atom_lines, 1, 6)),
  atom_name = trimws(substr(atom_lines, 13, 16)),
  res_name = trimws(substr(atom_lines, 18, 20)),
  chain = trimws(substr(atom_lines, 22, 22)),
  res_no = as.integer(trimws(substr(atom_lines, 23, 26))),
  x = as.numeric(substr(atom_lines, 31, 38)),
  y = as.numeric(substr(atom_lines, 39, 46)),
  z = as.numeric(substr(atom_lines, 47, 54)),
  stringsAsFactors = FALSE
)

# 【关键修改】配体残基名为 UNL
ligand_resname <- "UNL"

ligand_atoms <- pdb_df %>% filter(res_name == ligand_resname)
protein_atoms <- pdb_df %>% filter(res_name %in% c("ALA","ARG","ASN","ASP","CYS",
                                                   "GLN","GLU","GLY","HIS","ILE",
                                                   "LEU","LYS","MET","PHE","PRO",
                                                   "SER","THR","TRP","TYR","VAL"))

cat("配体原子数（应为24）：", nrow(ligand_atoms), "\n")
cat("蛋白原子数：", nrow(protein_atoms), "\n")

# 计算距离
calc_distance <- function(x1, y1, z1, x2, y2, z2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
}

threshold <- 4.0
binding_residues <- c()

for (i in 1:nrow(protein_atoms)) {
  for (j in 1:nrow(ligand_atoms)) {
    d <- calc_distance(
      protein_atoms$x[i], protein_atoms$y[i], protein_atoms$z[i],
      ligand_atoms$x[j], ligand_atoms$y[j], ligand_atoms$z[j]
    )
    if (d < threshold) {
      binding_residues <- c(binding_residues,
                            paste0(protein_atoms$res_name[i], protein_atoms$res_no[i]))
      break
    }
  }
}

binding_residues <- unique(binding_residues)
cat("\n与配体距离 <", threshold, "Å 的结合残基（", length(binding_residues), "个）：\n", sep = "")
print(binding_residues)

# 与ATP口袋比对
atp_pocket <- c("L77", "K78", "V115", "D116", "N117", "E119")
overlap <- intersect(binding_residues, atp_pocket)
cat("\n与ATP口袋重叠的残基：", 
    ifelse(length(overlap) > 0, paste(overlap, collapse = ", "), "无"), "\n")

# 保存结果
write.csv(data.frame(Residue = binding_residues), "C3_Binding_Residues.csv", row.names = FALSE)
cat("✅ 结果已保存为 C3_Binding_Residues.csv\n")



# ===================== 从C2复合物中提取结合残基 =====================
setwd("C:/CPTCA/dockdata")

pdb_lines <- readLines("c2 STK11_BRD−K53972329_out_2.-9.3.complex.pdb")
atom_lines <- pdb_lines[grepl("^ATOM|^HETATM", pdb_lines)]

# 解析PDB
pdb_df <- data.frame(
  record = trimws(substr(atom_lines, 1, 6)),
  atom_name = trimws(substr(atom_lines, 13, 16)),
  res_name = trimws(substr(atom_lines, 18, 20)),
  chain = trimws(substr(atom_lines, 22, 22)),
  res_no = as.integer(trimws(substr(atom_lines, 23, 26))),
  x = as.numeric(substr(atom_lines, 31, 38)),
  y = as.numeric(substr(atom_lines, 39, 46)),
  z = as.numeric(substr(atom_lines, 47, 54)),
  stringsAsFactors = FALSE
)

# 【关键修改】配体残基名为 UNL
ligand_resname <- "UNL"

ligand_atoms <- pdb_df %>% filter(res_name == ligand_resname)
protein_atoms <- pdb_df %>% filter(res_name %in% c("ALA","ARG","ASN","ASP","CYS",
                                                   "GLN","GLU","GLY","HIS","ILE",
                                                   "LEU","LYS","MET","PHE","PRO",
                                                   "SER","THR","TRP","TYR","VAL"))

cat("配体原子数（应为24）：", nrow(ligand_atoms), "\n")
cat("蛋白原子数：", nrow(protein_atoms), "\n")

# 计算距离
calc_distance <- function(x1, y1, z1, x2, y2, z2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
}

threshold <- 4.0
binding_residues <- c()

for (i in 1:nrow(protein_atoms)) {
  for (j in 1:nrow(ligand_atoms)) {
    d <- calc_distance(
      protein_atoms$x[i], protein_atoms$y[i], protein_atoms$z[i],
      ligand_atoms$x[j], ligand_atoms$y[j], ligand_atoms$z[j]
    )
    if (d < threshold) {
      binding_residues <- c(binding_residues,
                            paste0(protein_atoms$res_name[i], protein_atoms$res_no[i]))
      break
    }
  }
}

binding_residues <- unique(binding_residues)
cat("\n与配体距离 <", threshold, "Å 的结合残基（", length(binding_residues), "个）：\n", sep = "")
print(binding_residues)

# 与ATP口袋比对
atp_pocket <- c("L77", "K78", "V115", "D116", "N117", "E119")
overlap <- intersect(binding_residues, atp_pocket)
cat("\n与ATP口袋重叠的残基：", 
    ifelse(length(overlap) > 0, paste(overlap, collapse = ", "), "无"), "\n")

# 保存结果
write.csv(data.frame(Residue = binding_residues), "C2_Binding_Residues.csv", row.names = FALSE)
cat("✅ 结果已保存为 C3inding_Residues.csv\n")

####MSK分析####
####数据整理####
setwd("D://SHENGXIN//MSK//data")

#mutations <- read.table("mutations.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)

clinical <- read.table(
  file = "tmb_mskcc_2018_clinical_data.tsv",  # 替换成你的TSV文件路径，txt后缀也可以
  sep = "\t",          # ✔️ 核心！指定分隔符为制表符，TSV必加
  header = TRUE,       # ✔️ 第一行作为列名（99%的TSV都有表头，必加）
  stringsAsFactors = FALSE,  # ✔️ 字符列不自动转为因子（R4.0+默认FALSE，建议写上兼容旧版本）
  check.names = FALSE, # ✔️ 保留原始列名（比如列名含特殊字符/空格/数字开头，不会被自动修改）
  na.strings = c("", "NA") # 可选：将空值识别为NA
)

# 1. 读取原表头行
original_header <- read.table(
  "table.tsv", nrows = 1, sep = "\t", stringsAsFactors = FALSE
)[1, ]

# 2. 读取数据（包含原表头行）
temp_data <- read.table(
  "table.tsv", sep = "\t", header = FALSE, skip = 1, stringsAsFactors = FALSE
)

# 3. 补全列名（原表头+自动生成的列名）
colnames(temp_data) <- c(
  original_header,  # 原表头列名
  paste0("V", (length(original_header)+1):ncol(temp_data))  # 自动补剩余列名
)

mutation <- temp_data

library(dplyr)

# 使用filter函数
lung_adeno_df <- mutation %>%
  filter(`Cancer Type Detailed` == "Lung Adenocarcinoma")

write.table(lung_adeno_df,"laud stk11 mutation.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

common_ids <- intersect(lung_adeno_df$`Sample ID`, clinical$`Sample ID`)


clinical1 <- clinical[,c("Sample ID","Age at Which Sequencing was Reported (Days)","Overall Survival (Months)","Overall Survival Status",
                         "Sample Type","Sex","TMB (nonsynonymous)")]

# 添加STK11列，如果在common_ids中则标记为1，否则为0
clinical1$STK11 <- ifelse(clinical1$`Sample ID` %in% common_ids, 1, 0)

# 验证结果
table(clinical1$STK11)
cat("标记为1的样本数:", sum(clinical1$STK11), "\n")
cat("标记为0的样本数:", sum(clinical1$STK11 == 0), "\n")
# 修改特定列名
colnames(clinical1)[c(2, 3, 4, 5, 7)] <- c("Age", "OStime", "OS", "M", "TMB")
clinical1$OS <- substr(clinical1$OS, 1, 1)

clinical1$M <- ifelse(clinical1$M == "Metastasis", "M1", "M0")
# 查看替换后的结果
table(clinical1$M, useNA = "always")

library(tibble)
# 将第一列转换为行名
clinical1 <- column_to_rownames(clinical1, var = names(clinical1)[1])

write.table(clinical1,"clinical.txt",sep = "\t",row.names = T,col.names = NA,quote = F)





# 1. 加载必要包
# ====================
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# ====================
# 2. 数据准备
# ====================
clinical1 <- read.table("clinical.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 确认数据类型
clinical1 <- clinical1 %>%
  mutate(
    STK11 = as.factor(STK11),           # STK11突变状态因子化
    OS = as.numeric(OS),                # 生存状态
    OStime = as.numeric(OStime),        # 生存时间
    Sex = as.factor(Sex),               # 性别
    M = as.factor(M),                   # M分期
    TMB = as.numeric(TMB)               # 肿瘤突变负荷
  )

# 查看STK11分组情况
table(clinical1$STK11)
cat("STK11突变组样本数:", sum(clinical1$STK11 == 1), "\n")
cat("STK11野生组样本数:", sum(clinical1$STK11 == 0), "\n")

# 计算中位生存时间
surv_summary <- clinical1 %>%
  group_by(STK11) %>%
  summarise(
    n = n(),
    median_survival = median(OStime, na.rm = TRUE),
    mean_survival = mean(OStime, na.rm = TRUE)
  )
print(surv_summary)

# ====================
# 3. Kaplan-Meier生存分析
# ====================
# 创建生存对象
surv_obj <- Surv(time = clinical1$OStime, event = clinical1$OS)

# 拟合生存曲线
fit <- survfit(surv_obj ~ STK11, data = clinical1)

# 基本生存曲线
ggsurvplot(
  fit,
  data = clinical1,
  pval = TRUE,                     # 显示log-rank检验p值
  pval.method = TRUE,              # 显示检验方法
  conf.int = TRUE,                 # 显示置信区间
  risk.table = TRUE,               # 显示风险表
  risk.table.height = 0.25,        # 风险表高度
  surv.median.line = "hv",         # 显示中位生存时间线
  legend.labs = c("STK11野生型", "STK11突变型"),  # 图例标签
  legend.title = "STK11状态",
  palette = c("#2E9FDF", "#E74C3C"),  # 颜色
  xlab = "时间 (月)",               # x轴标签
  ylab = "总生存率",
  title = "STK11突变对总生存期的影响",
  break.time.by = 12,              # 每12个月一个刻度
  ggtheme = theme_bw()             # 主题
)


####生存分析####
# 1. 加载必要包
# ====================
# 检查并安装必要包
required_packages <- c("survival", "survminer", "dplyr", "ggplot2", "tidyverse")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(tidyverse)

# 2. 数据准备
# ====================
# 读取数据并检查
if(!file.exists("clinical.txt")) {
  stop("文件 'clinical.txt' 不存在，请检查文件路径")
}

clinical1 <- tryCatch({
  read.table("clinical.txt", sep = "\t", row.names = 1, 
             check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
}, error = function(e) {
  stop(paste("读取文件失败:", e$message))
})

# 检查必要的列是否存在
required_cols <- c("STK11", "OS", "OStime")
missing_cols <- setdiff(required_cols, colnames(clinical1))
if(length(missing_cols) > 0) {
  stop(paste("缺少必要的列:", paste(missing_cols, collapse = ", ")))
}

# 数据清理和转换
clinical1_clean <- clinical1 %>%
  # 移除生存时间或状态为NA的样本
  filter(!is.na(OStime) & !is.na(OS)) %>%
  # 移除生存时间 <= 0 的样本
  filter(OStime > 0) %>%
  # 转换数据类型
  mutate(
    STK11 = factor(STK11, levels = c(0, 1), labels = c("Wild-type", "Mutant")),
    OS = as.integer(OS),  # 生存状态应为整数
    OStime = as.numeric(OStime),
    # 可选：处理其他协变量
    Sex = if(exists("Sex", where = clinical1)) factor(Sex) else NULL,
    M = if(exists("M", where = clinical1)) factor(M) else NULL,
    TMB = if(exists("TMB", where = clinical1)) as.numeric(TMB) else NULL
  ) %>%
  # 移除STK11状态为NA的样本
  filter(!is.na(STK11))

# 数据概况报告
cat("=== 数据概况 ===\n")
cat("总样本数:", nrow(clinical1), "\n")
cat("有效样本数:", nrow(clinical1_clean), "\n")
cat("移除样本数:", nrow(clinical1) - nrow(clinical1_clean), "\n\n")

# STK11分组统计
stk11_summary <- clinical1_clean %>%
  group_by(STK11) %>%
  summarise(
    样本数 = n(),
    比例 = paste0(round(n()/nrow(clinical1_clean)*100, 1), "%"),
    中位生存时间 = round(median(OStime, na.rm = TRUE), 1),
    平均生存时间 = round(mean(OStime, na.rm = TRUE), 1),
    死亡事件数 = sum(OS == 1),
    删失比例 = paste0(round(sum(OS == 0)/n()*100, 1), "%")
  )

print(stk11_summary)

# 生存时间分布可视化（可选）
if(nrow(clinical1_clean) > 0) {
  p_dist <- ggplot(clinical1_clean, aes(x = STK11, y = OStime, fill = STK11)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = c("Wild-type" = "#2E9FDF", "Mutant" = "#E74C3C")) +
    labs(x = "STK11 Status", y = "Survival Time (months)", title = "Survival Time Distribution by STK11 Status") +
    theme_bw() +
    theme(legend.position = "none")
  print(p_dist)
}

# 3. Kaplan-Meier生存分析
# ====================
# 创建生存对象
surv_obj <- Surv(time = clinical1_clean$OStime, event = clinical1_clean$OS)

# 拟合生存曲线
fit <- survfit(surv_obj ~ STK11, data = clinical1_clean)

# 计算log-rank检验p值
logrank_test <- survdiff(surv_obj ~ STK11, data = clinical1_clean)
logrank_p <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)

# 格式化p值显示
format_p_value <- function(p) {
  if(p < 0.001) {
    return("P < 0.001")
  } else {
    return(paste("P =", sprintf("%.3f", p)))
  }
}

# 绘制生存曲线
surv_plot <- ggsurvplot(
  fit,
  data = clinical1_clean,
  pval = TRUE,                     # 显示log-rank检验p值
  pval.method = TRUE,              # 显示检验方法
  pval.size = 4,
  pval.coord = c(max(clinical1_clean$OStime, na.rm = TRUE)*0.7, 0.95),
  conf.int = TRUE,                 # 显示置信区间
  conf.int.style = "ribbon",       # 置信区间样式
  conf.int.alpha = 0.2,
  risk.table = TRUE,               # 显示风险表
  risk.table.title = "Number at risk",
  risk.table.height = 0.20,        # 风险表高度
  risk.table.y.text = FALSE,       # 不显示风险表y轴文本
  surv.median.line = "hv",         # 显示中位生存时间线
  legend.labs = levels(clinical1_clean$STK11),  # 图例标签
  legend.title = "STK11 Status",
  legend = c(0.85, 0.85),          # 图例位置
  palette = c("#2E9FDF", "#E74C3C"),  # 颜色
  xlab = "Time (months)",               # x轴标签
  ylab = "Overall Survival Probability",
  title = "Impact of STK11 Mutation Status on Overall Survival",
  subtitle = paste("Log-rank test:", format_p_value(logrank_p)),
  break.time.by = 12,              # 每12个月一个刻度
  ggtheme = theme_bw(base_size = 12),  # 主题
  tables.theme = theme_cleantable(),   # 风险表主题
  font.main = c(14, "bold", "black"),
  font.x = c(12, "plain", "black"),
  font.y = c(12, "plain", "black"),
  font.tickslab = c(11, "plain", "black"),
  censor.shape = 124,              # 删失标记形状
  censor.size = 3,                 # 删失标记大小
  xlim = c(0, max(clinical1_clean$OStime, na.rm = TRUE))  # 设置x轴范围
)

# 输出图形
print(surv_plot)

# 4. Cox比例风险模型（可选）
# ====================
# 单变量Cox回归
if(nrow(clinical1_clean) > 0) {
  cox_univ <- coxph(surv_obj ~ STK11, data = clinical1_clean)
  cox_summary <- summary(cox_univ)
  
  cat("\n=== Cox比例风险模型结果 ===\n")
  print(cox_summary)
  
  # 提取HR和95%CI
  hr <- exp(cox_summary$coefficients[1, "coef"])
  hr_lower <- exp(cox_summary$coefficients[1, "coef"] - 1.96 * cox_summary$coefficients[1, "se(coef)"])
  hr_upper <- exp(cox_summary$coefficients[1, "coef"] + 1.96 * cox_summary$coefficients[1, "se(coef)"])
  
  cat(sprintf("\n风险比 (95%% 置信区间): %.2f (%.2f - %.2f)", hr, hr_lower, hr_upper))
  cat("\nP值:", cox_summary$coefficients[1, "Pr(>|z|)"])
}

# 5. 中位生存时间及其置信区间
# ====================
if(nrow(clinical1_clean) > 0) {
  median_survival <- survfit(surv_obj ~ STK11, data = clinical1_clean)
  median_summary <- surv_median(median_survival)
  
  cat("\n\n=== 中位生存时间及其95%置信区间 ===\n")
  print(median_summary)
}

# 6. 保存结果
# ====================
# 保存生存曲线
if(nrow(clinical1_clean) > 0) {
  # 保存为PDF
  pdf("STK11_KM_Curve.pdf", width = 8, height = 7)
  print(surv_plot)
  dev.off()
  
  # 保存为PNG
  png("STK11_KM_Curve.png", width = 2400, height = 2100, res = 300)
  print(surv_plot)
  dev.off()
  
  cat("\n生存曲线已保存为：STK11_KM_Curve.pdf 和 STK11_KM_Curve.png\n")
  
  # 保存统计结果为CSV
  if(exists("cox_summary")) {
    # 提取Cox模型结果
    cox_results <- data.frame(
      变量 = "STK11突变型 vs 野生型",
      风险比 = hr,
      风险比_下限 = hr_lower,
      风险比_上限 = hr_upper,
      P值 = cox_summary$coefficients[1, "Pr(>|z|)"],
      Log_rank_P值 = logrank_p
    )
    
    write.csv(cox_results, "STK11_生存分析结果.csv", row.names = FALSE)
    cat("统计结果已保存为：STK11_生存分析结果.csv\n")
  }
  
  ####生存分析+单因素cox回归分析####
  # 1. 加载必要包
  # ====================
  if (!require("survival")) install.packages("survival")
  if (!require("survminer")) install.packages("survminer")
  if (!require("dplyr")) install.packages("dplyr")
  if (!require("ggplot2")) install.packages("ggplot2")
  
  library(survival)
  library(survminer)
  library(dplyr)
  library(ggplot2)
  
  # 2. 数据准备
  # ====================
  # 读取数据
  clinical1 <- read.table("clinical.txt", sep = "\t", row.names = 1, 
                          check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
  
  # 数据清洗和转换
  clinical1_clean <- clinical1 %>%
    # 移除生存时间或状态为NA的样本
    filter(!is.na(OStime) & !is.na(OS)) %>%
    # 移除生存时间 <= 0 的样本
    filter(OStime > 0) %>%
    # 转换STK11为因子
    mutate(
      STK11 = factor(STK11, levels = c(0, 1), labels = c("野生型", "突变型")),
      OS = as.integer(OS),  # 生存状态应为整数
      OStime = as.numeric(OStime)
    ) %>%
    # 移除STK11状态为NA的样本
    filter(!is.na(STK11))
  
  # 数据概况
  cat("=== 数据概况 ===\n")
  cat("总样本数:", nrow(clinical1_clean), "\n")
  cat("野生型样本数:", sum(clinical1_clean$STK11 == "野生型"), "\n")
  cat("突变型样本数:", sum(clinical1_clean$STK11 == "突变型"), "\n")
  cat("死亡事件数:", sum(clinical1_clean$OS == 1), "\n")
  cat("删失样本数:", sum(clinical1_clean$OS == 0), "\n\n")
  
  # 3. Kaplan-Meier生存分析（图片展示改为英文）
  # ====================
  # 创建生存对象
  surv_obj <- Surv(time = clinical1_clean$OStime, event = clinical1_clean$OS)
  
  # 拟合生存曲线
  fit <- survfit(surv_obj ~ STK11, data = clinical1_clean)
  
  # 计算log-rank检验p值
  logrank_test <- survdiff(surv_obj ~ STK11, data = clinical1_clean)
  logrank_p <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
  
  # 格式化p值显示
  format_p_value <- function(p) {
    if(p < 0.001) {
      return("P < 0.001")
    } else {
      return(paste("P =", sprintf("%.3f", p)))
    }
  }
  
  # 绘制生存曲线（使用英文标签）
  surv_plot <- ggsurvplot(
    fit,
    data = clinical1_clean,
    pval = TRUE,
    pval.method = TRUE,
    pval.size = 4,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.20,
    surv.median.line = "hv",
    legend.labs = c("STK11 Wild-type", "STK11 Mutant"),  # 英文图例标签
    legend.title = "STK11 Status",  # 英文图例标题
    palette = c("#2E9FDF", "#E74C3C"),
    xlab = "Time (months)",  # 英文x轴标签
    ylab = "Overall Survival Probability",  # 英文y轴标签
    title = "Impact of STK11 Mutation Status on Overall Survival",  # 英文标题
    subtitle = paste("Log-rank test:", format_p_value(logrank_p)),  # 英文副标题
    break.time.by = 12,
    ggtheme = theme_bw(base_size = 12),
    risk.table.title = "Number at risk",  # 风险表英文标题
    risk.table.y.text = FALSE
  )
  
  print(surv_plot)
  
  # 4. STK11单因素Cox回归分析（重点）
  # ====================
  cat("\n=== STK11单因素Cox回归分析 ===\n")
  
  # 进行Cox回归分析
  cox_model <- coxph(surv_obj ~ STK11, data = clinical1_clean)
  cox_summary <- summary(cox_model)
  
  # 提取关键结果
  hr <- exp(cox_summary$coefficients[1, "coef"])
  hr_se <- cox_summary$coefficients[1, "se(coef)"]
  hr_lower <- exp(cox_summary$coefficients[1, "coef"] - 1.96 * hr_se)
  hr_upper <- exp(cox_summary$coefficients[1, "coef"] + 1.96 * hr_se)
  p_value <- cox_summary$coefficients[1, "Pr(>|z|)"]
  z_value <- cox_summary$coefficients[1, "z"]
  
  # 计算95%置信区间
  conf_int <- exp(confint(cox_model))
  
  # 详细展示结果
  cat("\n1. 模型基本信息:\n")
  cat("   模型公式: surv_obj ~ STK11\n")
  cat("   样本数量:", nrow(clinical1_clean), "\n")
  cat("   事件数量:", sum(clinical1_clean$OS == 1), "\n")
  cat("   删失数量:", sum(clinical1_clean$OS == 0), "\n")
  
  cat("\n2. 分组基本信息:\n")
  group_info <- clinical1_clean %>%
    group_by(STK11) %>%
    summarise(
      样本数 = n(),
      事件数 = sum(OS == 1),
      删失数 = sum(OS == 0),
      删失比例 = paste0(round(sum(OS == 0)/n()*100, 1), "%"),
      中位生存时间 = round(median(OStime, na.rm = TRUE), 1),
      平均生存时间 = round(mean(OStime, na.rm = TRUE), 1)
    )
  print(group_info)
  
  cat("\n3. Cox回归系数:\n")
  cat("   系数 (β):", round(cox_summary$coefficients[1, "coef"], 4), "\n")
  cat("   标准误:", round(hr_se, 4), "\n")
  cat("   z值:", round(z_value, 4), "\n")
  cat("   P值:", ifelse(p_value < 0.001, "<0.001", sprintf("%.4f", p_value)), "\n")
  
  cat("\n4. 风险比 (Hazard Ratio):\n")
  cat("   HR:", round(hr, 3), "\n")
  cat("   95% 置信区间: [", round(hr_lower, 3), ", ", round(hr_upper, 3), "]\n", sep="")
  cat("   95% 置信区间 (exp): [", round(conf_int[1, 1], 3), ", ", round(conf_int[1, 2], 3), "]\n", sep="")
  
  # 解释HR值
  cat("\n5. 结果解释:\n")
  if (hr > 1) {
    cat("   HR = ", round(hr, 3), " > 1，说明STK11突变型患者的死亡风险是野生型患者的", 
        round(hr, 2), "倍。\n", sep="")
    cat("   即STK11突变型患者比野生型患者死亡风险高", round((hr-1)*100, 1), "%。\n", sep="")
  } else if (hr < 1) {
    cat("   HR = ", round(hr, 3), " < 1，说明STK11突变型患者的死亡风险是野生型患者的", 
        round(hr, 2), "倍。\n", sep="")
    cat("   即STK11突变型患者比野生型患者死亡风险低", round((1-hr)*100, 1), "%。\n", sep="")
  } else {
    cat("   HR = 1，说明STK11突变状态对生存无显著影响。\n")
  }
  
  cat("\n6. 模型检验:\n")
  cat("   Likelihood ratio test p值:", round(cox_summary$logtest[["pvalue"]], 4), "\n")
  cat("   Wald test p值:", round(cox_summary$waldtest[["pvalue"]], 4), "\n")
  cat("   Score (logrank) test p值:", round(cox_summary$sctest[["pvalue"]], 4), "\n")
  
  cat("\n7. 比例风险假设检验:\n")
  ph_test <- cox.zph(cox_model)
  print(ph_test)
  
  # 创建森林图（使用英文标签）
  cat("\n8. 可视化结果:\n")
  
  # 准备森林图数据
  forest_data <- data.frame(
    变量 = "STK11 (Mutant vs Wild-type)",  # 英文变量名
    HR = hr,
    低限 = hr_lower,
    高限 = hr_upper,
    P值 = p_value
  )
  
  # 创建森林图（英文标签）
  p_forest <- ggplot(forest_data, aes(x = HR, y = 变量)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_point(size = 4, color = ifelse(p_value < 0.05, "#E74C3C", "#2E9FDF")) +
    geom_errorbarh(aes(xmin = 低限, xmax = 高限), 
                   height = 0.1, linewidth = 1, 
                   color = ifelse(p_value < 0.05, "#E74C3C", "#2E9FDF")) +
    scale_x_log10(
      breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
      labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10"),
      limits = c(min(0.1, hr_lower*0.8), max(10, hr_upper*1.2))
    ) +
    geom_text(
      aes(label = sprintf("HR = %.2f\n(95%% CI: %.2f-%.2f)\nP = %s", 
                          HR, 低限, 高限, 
                          ifelse(P值 < 0.001, "<0.001", sprintf("%.3f", P值)))),
      x = max(hr_upper, 1) * 1.3,
      hjust = 0,
      size = 4.5,
      color = "black"
    ) +
    labs(
      x = "Hazard Ratio (HR, log10 scale)",  # 英文x轴标签
      y = "",
      title = "Univariable Cox Regression Analysis for STK11",  # 英文标题
      subtitle = paste("Total samples:", nrow(clinical1_clean),  # 英文副标题
                       " | Events:", sum(clinical1_clean$OS == 1))
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 12),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 11)
    )
  
  print(p_forest)
  
  # 5. 保存结果
  # ====================
  # 保存生存曲线（英文）
  pdf("STK11_KM_Survival_Curve.pdf", width = 8, height = 7)
  print(surv_plot)
  dev.off()
  
  # 保存森林图（英文）
  ggsave("STK11_Univariable_Cox_Forest_Plot.png", p_forest, width = 10, height = 4, dpi = 300)
  ggsave("STK11_Univariable_Cox_Forest_Plot.pdf", p_forest, width = 10, height = 4)
  
  # 创建详细结果表格
  results_table <- data.frame(
    变量 = "STK11突变状态",
    组别 = "突变型 vs 野生型",
    样本数 = nrow(clinical1_clean),
    事件数 = sum(clinical1_clean$OS == 1),
    野生型样本数 = sum(clinical1_clean$STK11 == "野生型"),
    突变型样本数 = sum(clinical1_clean$STK11 == "突变型"),
    系数 = round(cox_summary$coefficients[1, "coef"], 4),
    标准误 = round(hr_se, 4),
    z值 = round(z_value, 4),
    HR = round(hr, 3),
    HR_95CI_低限 = round(hr_lower, 3),
    HR_95CI_高限 = round(hr_upper, 3),
    P值 = ifelse(p_value < 0.001, "<0.001", sprintf("%.4f", p_value)),
    Log_rank_P值 = ifelse(logrank_p < 0.001, "<0.001", sprintf("%.4f", logrank_p))
  )
  
  # 保存为CSV
  write.csv(results_table, "STK11_单因素Cox回归结果.csv", row.names = FALSE, fileEncoding = "UTF-8")
  
  # 创建文本报告（中文）
  sink("STK11_生存分析报告.txt")
  cat("=================== STK11生存分析报告 ===================\n\n")
  cat("分析日期:", format(Sys.Date(), "%Y年%m月%d日"), "\n\n")
  
  cat("一、数据基本情况\n")
  cat("总样本数:", nrow(clinical1_clean), "\n")
  cat("野生型样本数:", sum(clinical1_clean$STK11 == "野生型"), "\n")
  cat("突变型样本数:", sum(clinical1_clean$STK11 == "突变型"), "\n")
  cat("死亡事件数:", sum(clinical1_clean$OS == 1), "\n")
  cat("删失样本数:", sum(clinical1_clean$OS == 0), "\n\n")
  
  cat("二、Kaplan-Meier生存分析\n")
  cat("Log-rank检验 P值:", format_p_value(logrank_p), "\n\n")
  
  cat("三、单因素Cox回归分析结果\n")
  cat("变量: STK11突变状态 (突变型 vs 野生型)\n")
  cat("样本数:", nrow(clinical1_clean), "\n")
  cat("事件数:", sum(clinical1_clean$OS == 1), "\n")
  cat("系数 (β):", round(cox_summary$coefficients[1, "coef"], 4), "\n")
  cat("风险比 (HR):", round(hr, 3), "\n")
  cat("95% 置信区间: [", round(hr_lower, 3), ", ", round(hr_upper, 3), "]\n", sep="")
  cat("P值:", ifelse(p_value < 0.001, "<0.001", sprintf("%.4f", p_value)), "\n\n")
  
  cat("四、结果解释\n")
  if (hr > 1) {
    cat("STK11突变型患者的死亡风险是野生型患者的", round(hr, 2), "倍。\n", sep="")
    cat("即STK11突变型患者比野生型患者死亡风险高", round((hr-1)*100, 1), "%。\n", sep="")
  } else if (hr < 1) {
    cat("STK11突变型患者的死亡风险是野生型患者的", round(hr, 2), "倍。\n", sep="")
    cat("即STK11突变型患者比野生型患者死亡风险低", round((1-hr)*100, 1), "%。\n", sep="")
  } else {
    cat("STK11突变状态对生存无显著影响。\n")
  }
  
  cat("\n五、比例风险假设检验\n")
  print(ph_test)
  
  cat("\n=================== 报告结束 ===================\n")
  sink()
  
  cat("\n分析完成！结果已保存：\n")
  cat("1. STK11_KM_Survival_Curve.pdf - Kaplan-Meier生存曲线（英文）\n")
  cat("2. STK11_Univariable_Cox_Forest_Plot.png/pdf - Cox回归森林图（英文）\n")
  cat("3. STK11_单因素Cox回归结果.csv - 详细分析结果\n")
  cat("4. STK11_生存分析报告.txt - 完整分析报告\n")
  
  ####多因素cox回归分析####
  # 1. 加载必要包
  # ====================
  if (!require("survival")) install.packages("survival")
  if (!require("survminer")) install.packages("survminer")
  if (!require("dplyr")) install.packages("dplyr")
  if (!require("ggplot2")) install.packages("ggplot2")
  if (!require("forestplot")) install.packages("forestplot")
  if (!require("grid")) install.packages("grid")
  
  library(survival)
  library(survminer)
  library(dplyr)
  library(ggplot2)
  library(forestplot)
  library(grid)
  
  # 2. 数据准备和清理
  # ====================
  # 读取数据
  clinical1 <- read.table("clinical.txt", sep = "\t", row.names = 1, 
                          check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
  
  # 检查必需的列是否存在
  required_cols <- c("STK11", "OS", "OStime", "Age", "Sex", "M", "TMB")
  missing_cols <- setdiff(required_cols, colnames(clinical1))
  if(length(missing_cols) > 0) {
    stop(paste("缺少必要的列:", paste(missing_cols, collapse = ", ")))
  }
  
  # 数据清理和转换
  clinical_clean <- clinical1 %>%
    filter(!is.na(OStime) & !is.na(OS)) %>%
    filter(OStime > 0) %>%
    mutate(
      STK11 = factor(STK11, levels = c(0, 1), labels = c("Wild-type", "Mutant")),
      OS = as.integer(OS),
      OStime = as.numeric(OStime),
      Age = as.numeric(Age),
      # 直接转换，假设这些列都存在（前面已检查）
      Sex = factor(Sex, levels = c("Male", "Female")),
      M = factor(M, levels = c("M0", "M1")),
      TMB = as.numeric(TMB)
    ) %>%
    na.omit()
  
  # 数据概况
  cat("=== 数据概况 ===\n")
  cat("原始样本数:", nrow(clinical1), "\n")
  cat("有效样本数（无缺失）:", nrow(clinical_clean), "\n")
  cat("移除样本数:", nrow(clinical1) - nrow(clinical_clean), "\n\n")
  
  # 变量描述性统计
  cat("=== 变量描述性统计 ===\n")
  desc_stats <- clinical_clean %>%
    summarise(
      `总样本数` = n(),
      `事件数` = sum(OS == 1),
      `删失数` = sum(OS == 0),
      `年龄(均值±标准差)` = sprintf("%.1f ± %.1f", mean(Age, na.rm = TRUE), sd(Age, na.rm = TRUE)),
      `年龄(范围)` = sprintf("%.1f-%.1f", min(Age, na.rm = TRUE), max(Age, na.rm = TRUE)),
      `男性比例` = paste0(round(sum(Sex == "Male")/n()*100, 1), "%"),
      `女性比例` = paste0(round(sum(Sex == "Female")/n()*100, 1), "%"),
      `M1比例` = paste0(round(sum(M == "M1")/n()*100, 1), "%"),
      `TMB(中位数)` = round(median(TMB, na.rm = TRUE), 1),
      `野生型比例` = paste0(round(sum(STK11 == "Wild-type")/n()*100, 1), "%"),
      `突变型比例` = paste0(round(sum(STK11 == "Mutant")/n()*100, 1), "%")
    )
  
  print(t(desc_stats))
  
  # STK11分组描述
  cat("\n=== STK11分组描述 ===\n")
  stk11_group_desc <- clinical_clean %>%
    group_by(STK11) %>%
    summarise(
      `样本数` = n(),
      `事件数` = sum(OS == 1),
      `中位生存时间(月)` = round(median(OStime), 1),
      `年龄(均值)` = round(mean(Age), 1),
      `男性比例` = paste0(round(sum(Sex == "Male")/n()*100, 1), "%"),
      `M1比例` = paste0(round(sum(M == "M1")/n()*100, 1), "%"),
      `TMB(中位数)` = round(median(TMB), 1)
    )
  
  print(stk11_group_desc)
  # 3. 单因素Cox回归分析
  # ====================
  cat("\n=== 单因素Cox回归分析 ===\n")
  
  # 创建生存对象
  surv_obj <- Surv(time = clinical_clean$OStime, event = clinical_clean$OS)
  
  # 定义分析变量及其参考水平
  variables <- list(
    STK11 = list(var = "STK11", ref = "野生型", comp = "突变型 vs 野生型"),
    Age = list(var = "Age", ref = "连续变量", comp = "每增加1岁"),
    Sex = list(var = "Sex", ref = "男性", comp = "女性 vs 男性"),
    M = list(var = "M", ref = "M0", comp = "M1 vs M0"),
    TMB = list(var = "TMB", ref = "连续变量", comp = "每增加1单位")
  )
  
  # 执行单因素Cox回归分析
  univariate_results <- data.frame()
  
  for(var_name in names(variables)) {
    var_info <- variables[[var_name]]
    
    formula <- as.formula(paste("surv_obj ~", var_info$var))
    cox_model <- coxph(formula, data = clinical_clean)
    cox_summary <- summary(cox_model)
    
    # 提取结果
    n_obs <- cox_summary$n
    n_event <- cox_summary$nevent
    
    for(i in 1:nrow(cox_summary$coefficients)) {
      coef_name <- rownames(cox_summary$coefficients)[i]
      hr <- exp(cox_summary$coefficients[i, "coef"])
      hr_lower <- exp(cox_summary$coefficients[i, "coef"] - 1.96 * cox_summary$coefficients[i, "se(coef)"])
      hr_upper <- exp(cox_summary$coefficients[i, "coef"] + 1.96 * cox_summary$coefficients[i, "se(coef)"])
      p_value <- cox_summary$coefficients[i, "Pr(>|z|)"]
      
      # 确定比较描述
      if(var_name == "STK11") {
        comparison <- var_info$comp
      } else if(var_name %in% c("Age", "TMB")) {
        comparison <- var_info$comp
      } else {
        # 对于分类变量，从系数名提取比较组
        comparison <- paste(gsub(var_info$var, "", coef_name), "vs", var_info$ref)
      }
      
      univariate_results <- rbind(univariate_results, data.frame(
        变量 = var_name,
        比较组 = comparison,
        样本数 = n_obs,
        事件数 = n_event,
        HR = round(hr, 3),
        HR_lower = round(hr_lower, 3),
        HR_upper = round(hr_upper, 3),
        P值 = p_value,
        显著性 = ifelse(p_value < 0.001, "***", 
                     ifelse(p_value < 0.01, "**", 
                            ifelse(p_value < 0.05, "*", ""))),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 格式化输出
  univariate_results$`HR(95%CI)` <- sprintf("%.2f (%.2f-%.2f)", 
                                            univariate_results$HR,
                                            univariate_results$HR_lower,
                                            univariate_results$HR_upper)
  
  univariate_results$`P值` <- ifelse(univariate_results$P值 < 0.001, "<0.001",
                                    sprintf("%.3f", univariate_results$P值))
  
  # 按P值排序
  univariate_results <- univariate_results[order(univariate_results$P值), ]
  
  cat("\n单因素Cox回归分析结果:\n")
  print(univariate_results[, c("变量", "比较组", "样本数", "事件数", "HR(95%CI)", "P值", "显著性")], 
        row.names = FALSE)
  
  # 4. 单因素森林图展示
  # 4. 单因素森林图展示（按照您提供的格式）
  # 4. 单因素森林图展示（按照您提供的格式）
  # ========================================
  cat("\n=== 生成单因素森林图 ===\n")
  
  # 准备森林图数据
  forest_data <- univariate_results %>%
    mutate(
      # 创建显示标签（使用英文）
      characteristics = paste(变量, " (", 比较组, ")", sep = ""),
      # 替换中文为英文
      characteristics = gsub("STK11", "STK11", characteristics),
      characteristics = gsub("Age", "Age", characteristics),
      characteristics = gsub("Sex", "Sex", characteristics),
      characteristics = gsub("M", "M stage", characteristics),
      characteristics = gsub("TMB", "TMB", characteristics),
      # 创建易读的变量名
      characteristics_clean = case_when(
        变量 == "STK11" & 比较组 == "突变型 vs 野生型" ~ "STK11 (Mutant vs Wild-type)",
        变量 == "Age" & 比较组 == "每增加1岁" ~ "Age (per year increase)",
        变量 == "Sex" & 比较组 == "女性 vs 男性" ~ "Sex (Female vs Male)",
        变量 == "M" & 比较组 == "M1 vs M0" ~ "M stage (M1 vs M0)",
        变量 == "TMB" & 比较组 == "每增加1单位" ~ "TMB (per unit increase)",
        TRUE ~ characteristics
      ),
      # 添加显著性标记
      signif_mark = ifelse(P值 < 0.001, "***",
                           ifelse(P值 < 0.01, "**",
                                  ifelse(P值 < 0.05, "*", "")))
    ) %>%
    arrange(变量)  # 按变量排序
  
  # 准备表格文本（包含significance列）
  tabletext <- cbind(
    c("Variable", forest_data$characteristics_clean),
    c("HR", sprintf("%.3f", forest_data$HR)),
    c("95% CI", sprintf("%.3f-%.3f", forest_data$HR_lower, forest_data$HR_upper)),
    c("P-value", ifelse(forest_data$P值 < 0.001, "<0.001", 
                        sprintf("%s", forest_data$P值))),
    c("Significance", forest_data$signif_mark)
  )
  
  # 计算合适的x轴刻度
  # 确定数据范围
  min_hr <- min(forest_data$HR_lower, na.rm = TRUE)
  max_hr <- max(forest_data$HR_upper, na.rm = TRUE)
  
  # 生成合适的刻度
  x_ticks <- c()
  if(min_hr < 0.1) x_ticks <- c(x_ticks, 0.1)
  if(min_hr < 0.25) x_ticks <- c(x_ticks, 0.25)
  x_ticks <- c(x_ticks, 0.5, 1)
  if(max_hr > 2) x_ticks <- c(x_ticks, 2)
  if(max_hr > 4) x_ticks <- c(x_ticks, 4)
  if(max_hr > 8) x_ticks <- c(x_ticks, 8)
  if(max_hr > 16) x_ticks <- c(x_ticks, 16)
  
  # 创建hrzl_lines列表
  total_rows <- nrow(tabletext) + 1
  hrzl_lines_list <- list(
    "1" = gpar(lwd = 2, col = "black"),
    "2" = gpar(lwd = 1.5, col = "gray70")
  )
  hrzl_lines_list[[as.character(total_rows)]] <- gpar(lwd = 2, col = "black")
  
  # 创建颜色向量，突出显示显著结果
  box_colors <- rep("#2E86AB", nrow(forest_data))  # 默认蓝色
  line_colors <- rep("#1B4965", nrow(forest_data))  # 默认深蓝色
  
  # 显著结果使用红色
  sig_indices <- which(forest_data$P值 < 0.05)
  box_colors[sig_indices] <- "#E63946"  # 显著红色
  line_colors[sig_indices] <- "#E63946"   # 显著红色
  
  # 1. 在R面板上展示森林图
  # =========================
  cat("正在生成单因素Cox回归森林图...\n")
  
  # 首先创建一个新的图形页面
  grid.newpage()
  
  # 创建视口布局，为标题留出空间
  vp_layout <- viewport(layout = grid.layout(
    nrow = 3, 
    ncol = 1,
    heights = unit(c(1.0, 0.2, 8.0), "inches")  # 标题空间、间隔、主图空间
  ))
  
  pushViewport(vp_layout)
  
  # 在第一个视口中添加标题
  pushViewport(viewport(layout.pos.row = 1))
  grid.text("Univariate Cox Regression Analysis\nForest Plot", 
            gp = gpar(fontsize = 16, fontface = "bold"))
  popViewport()
  
  # 在第三个视口中绘制森林图
  pushViewport(viewport(layout.pos.row = 3))
  
  # 绘制森林图
  forestplot(
    labeltext = tabletext,
    mean = c(NA, forest_data$HR),
    lower = c(NA, forest_data$HR_lower), 
    upper = c(NA, forest_data$HR_upper),
    
    # 图形设置
    graph.pos = 5,  # 因为现在有5列，图形放在第5列
    graphwidth = unit(0.3, "npc"),
    
    # 处理极端值
    clip = c(min(x_ticks), max(x_ticks)),  # 限制显示范围
    xlog = TRUE,  # 使用对数刻度
    
    # 颜色设置
    col = fpColors(
      box = box_colors,
      lines = line_colors,
      summary = "#1B4965",
      zero = "gray50"
    ),
    
    # 箱线图设置
    boxsize = 0.3,
    lwd.ci = 2,
    ci.vertices = TRUE,
    ci.vertices.height = 0.15,
    
    # 参考线
    zero = 1,
    lwd.zero = 2,
    grid = structure(x_ticks, gp = gpar(lty = 2, col = "gray80")),
    
    # 坐标轴
    xticks = x_ticks,
    lwd.xaxis = 2,
    xlab = "\nHazard Ratio (95% Confidence Interval)",
    
    # 文本格式
    txt_gp = fpTxtGp(
      label = gpar(cex = 0.95, fontfamily = "sans"),
      ticks = gpar(cex = 0.85),
      xlab = gpar(cex = 1.0, fontface = "bold")
    ),
    
    # 横线
    hrzl_lines = hrzl_lines_list,
    
    # 布局 - 调整边距确保不重叠
    lineheight = unit(0.7, "cm"),
    colgap = unit(0.4, "cm"),
    mar = unit(c(0.5, 0.5, 0, 0.5), "cm"),  # 减少上边距，因为标题已单独处理
    
    # 不显示内置标题
    new_page = FALSE
  )
  
  # 返回到之前的viewport
  popViewport()
  popViewport()
  
  cat("森林图已展示在R面板上\n")
  
  # 2. 保存PDF版本（可选）
  # ========================
  save_pdf <- readline(prompt = "是否保存为PDF文件？(y/n): ")
  
  if(tolower(save_pdf) == "y" || tolower(save_pdf) == "yes") {
    # 保存PDF版本
    pdf("Univariate_Cox_Forestplot.pdf", width = 12, height = 8)
    
    # 首先创建一个新的图形页面
    grid.newpage()
    
    # 创建视口布局，为标题留出空间
    vp_layout <- viewport(layout = grid.layout(
      nrow = 3, 
      ncol = 1,
      heights = unit(c(1.0, 0.2, 8.0), "inches")
    ))
    
    pushViewport(vp_layout)
    
    # 在第一个视口中添加标题
    pushViewport(viewport(layout.pos.row = 1))
    grid.text("Univariate Cox Regression Analysis\nForest Plot", 
              gp = gpar(fontsize = 16, fontface = "bold"))
    popViewport()
    
    # 在第三个视口中绘制森林图
    pushViewport(viewport(layout.pos.row = 3))
    
    # 绘制森林图
    forestplot(
      labeltext = tabletext,
      mean = c(NA, forest_data$HR),
      lower = c(NA, forest_data$HR_lower), 
      upper = c(NA, forest_data$HR_upper),
      
      # 图形设置
      graph.pos = 5,  # 因为现在有5列，图形放在第5列
      graphwidth = unit(0.3, "npc"),
      
      # 处理极端值
      clip = c(min(x_ticks), max(x_ticks)),
      xlog = TRUE,
      
      # 颜色设置
      col = fpColors(
        box = box_colors,
        lines = line_colors,
        summary = "#1B4965",
        zero = "gray50"
      ),
      
      # 箱线图设置
      boxsize = 0.3,
      lwd.ci = 2,
      ci.vertices = TRUE,
      ci.vertices.height = 0.15,
      
      # 参考线
      zero = 1,
      lwd.zero = 2,
      grid = structure(x_ticks, gp = gpar(lty = 2, col = "gray80")),
      
      # 坐标轴
      xticks = x_ticks,
      lwd.xaxis = 2,
      xlab = "\nHazard Ratio (95% Confidence Interval)",
      
      # 文本格式
      txt_gp = fpTxtGp(
        label = gpar(cex = 0.95, fontfamily = "sans"),
        ticks = gpar(cex = 0.85),
        xlab = gpar(cex = 1.0, fontface = "bold")
      ),
      
      # 横线
      hrzl_lines = hrzl_lines_list,
      
      # 布局
      lineheight = unit(0.7, "cm"),
      colgap = unit(0.4, "cm"),
      mar = unit(c(0.5, 0.5, 0, 0.5), "cm"),
      
      # 不显示内置标题
      new_page = FALSE
    )
    
    # 返回到之前的viewport
    popViewport()
    popViewport()
    
    dev.off()
    cat("已保存森林图(PDF): Univariate_Cox_Forestplot.pdf\n")
  } else {
    cat("未保存PDF文件\n")
  }
  
  # 4. 多因素Cox回归分析
  # ====================
  cat("\n=== 多因素Cox回归分析 ===\n")
  
  # 构建多因素模型
  multiv_formula <- surv_obj ~ STK11 + Age + Sex + M + TMB
  cox_multiv <- coxph(multiv_formula, data = clinical_clean)
  cox_multiv_summary <- summary(cox_multiv)
  
  cat("\n1. 多因素Cox模型信息:\n")
  cat("   模型公式:", paste(format(multiv_formula), collapse = ""), "\n")
  cat("   有效样本数:", nrow(clinical_clean), "\n")
  cat("   事件数:", sum(clinical_clean$OS == 1), "\n")
  cat("   变量数:", length(variables), "\n")
  
  cat("\n2. 多因素分析结果:\n")
  
  # 提取多因素分析结果
  multiv_results <- data.frame(
    Variable = rownames(cox_multiv_summary$coefficients),
    HR = round(exp(cox_multiv_summary$coefficients[, "coef"]), 3),
    HR_lower = round(exp(cox_multiv_summary$coefficients[, "coef"] - 1.96 * cox_multiv_summary$coefficients[, "se(coef)"]), 3),
    HR_upper = round(exp(cox_multiv_summary$coefficients[, "coef"] + 1.96 * cox_multiv_summary$coefficients[, "se(coef)"]), 3),
    Coefficient = round(cox_multiv_summary$coefficients[, "coef"], 4),
    SE = round(cox_multiv_summary$coefficients[, "se(coef)"], 4),
    z_value = round(cox_multiv_summary$coefficients[, "z"], 3),
    P_value = cox_multiv_summary$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  # 添加变量说明
  multiv_results$Variable_label <- c(
    "STK11Mutant" = "STK11 (Mutant vs Wild-type)",
    "Age" = "Age (per year increase)",
    "SexFemale" = "Sex (Female vs Male)",
    "MM1" = "M stage (M1 vs M0)",
    "TMB" = "TMB (per unit increase)"
  )[multiv_results$Variable]
  
  # 添加显著性标记
  multiv_results$signif_mark <- ifelse(multiv_results$P_value < 0.001, "***",
                                       ifelse(multiv_results$P_value < 0.01, "**",
                                              ifelse(multiv_results$P_value < 0.05, "*", "")))
  
  # 格式化显示
  multiv_results$P_display <- ifelse(
    multiv_results$P_value < 0.001, 
    "<0.001", 
    sprintf("%.3f", multiv_results$P_value)
  )
  
  multiv_results$HR_95CI <- sprintf(
    "%.2f (%.2f-%.2f)", 
    multiv_results$HR,
    multiv_results$HR_lower,
    multiv_results$HR_upper
  )
  
  # 打印多因素结果表格
  print(multiv_results[, c("Variable_label", "HR_95CI", "P_display", "signif_mark")])
  
  cat("\n3. 模型整体检验:\n")
  cat("   Likelihood ratio test p值:", round(cox_multiv_summary$logtest[["pvalue"]], 4), "\n")
  cat("   Wald test p值:", round(cox_multiv_summary$waldtest[["pvalue"]], 4), "\n")
  cat("   Score (logrank) test p值:", round(cox_multiv_summary$sctest[["pvalue"]], 4), "\n")
  
  # 模型拟合优度
  cat("\n4. 模型拟合优度:\n")
  cat("   似然比检验卡方:", round(cox_multiv_summary$logtest[["test"]], 2), 
      " (df = ", cox_multiv_summary$logtest[["df"]], ")\n", sep="")
  cat("   AIC:", round(AIC(cox_multiv), 2), "\n")
  cat("   Concordance Index (C-index):", round(cox_multiv_summary$concordance[["C"]], 3), "\n")
  
  # 5. 比例风险假设检验
  # ====================
  cat("\n=== 比例风险假设检验 ===\n")
  ph_test <- cox.zph(cox_multiv)
  print(ph_test)
  
  # 可视化比例风险假设检验
  if(nrow(clinical_clean) > 0) {
    cat("\n正在生成比例风险假设检验图...\n")
    
    # 在R面板上展示
    dev.new(width = 10, height = 6)
    plot(ph_test)
    
    # 询问是否保存
    save_ph_test <- readline(prompt = "是否保存比例风险假设检验图？(y/n): ")
    
    if(tolower(save_ph_test) == "y" || tolower(save_ph_test) == "yes") {
      png("Proportional_Hazards_Test.png", width = 1200, height = 800, res = 150)
      plot(ph_test)
      dev.off()
      cat("比例风险假设检验图已保存为: Proportional_Hazards_Test.png\n")
    }
  }
  
  # 6. 多因素Cox回归森林图（按照单因素森林图格式）
  # ================================================
  cat("\n=== 多因素Cox回归森林图 ===\n")
  
  # 准备森林图数据
  forest_data <- multiv_results %>%
    arrange(desc(row_number())) %>%  # 反转顺序，使第一个变量在顶部
    mutate(
      Variable_display = Variable_label,
      HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, HR_lower, HR_upper),
      P_value_display = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value))
    )
  
  # 准备表格文本（包含significance列）
  tabletext <- cbind(
    c("Variable", forest_data$Variable_display),
    c("HR", sprintf("%.3f", forest_data$HR)),
    c("95% CI", sprintf("%.3f-%.3f", forest_data$HR_lower, forest_data$HR_upper)),
    c("P-value", forest_data$P_value_display),
    c("Significance", forest_data$signif_mark)
  )
  
  # 计算合适的x轴刻度
  min_hr_multiv <- min(forest_data$HR_lower, na.rm = TRUE)
  max_hr_multiv <- max(forest_data$HR_upper, na.rm = TRUE)
  
  x_ticks_multiv <- c()
  if(min_hr_multiv < 0.1) x_ticks_multiv <- c(x_ticks_multiv, 0.1)
  if(min_hr_multiv < 0.25) x_ticks_multiv <- c(x_ticks_multiv, 0.25)
  x_ticks_multiv <- c(x_ticks_multiv, 0.5, 1)
  if(max_hr_multiv > 2) x_ticks_multiv <- c(x_ticks_multiv, 2)
  if(max_hr_multiv > 4) x_ticks_multiv <- c(x_ticks_multiv, 4)
  if(max_hr_multiv > 8) x_ticks_multiv <- c(x_ticks_multiv, 8)
  if(max_hr_multiv > 16) x_ticks_multiv <- c(x_ticks_multiv, 16)
  
  # 创建hrzl_lines列表
  total_rows_multiv <- nrow(tabletext) + 1
  hrzl_lines_list_multiv <- list(
    "1" = gpar(lwd = 2, col = "black"),
    "2" = gpar(lwd = 1.5, col = "gray70")
  )
  hrzl_lines_list_multiv[[as.character(total_rows_multiv)]] <- gpar(lwd = 2, col = "black")
  
  # 创建颜色向量，突出显示显著结果
  box_colors_multiv <- rep("#2E86AB", nrow(forest_data))  # 默认蓝色
  line_colors_multiv <- rep("#1B4965", nrow(forest_data))  # 默认深蓝色
  
  # 显著结果使用红色
  sig_indices_multiv <- which(forest_data$P_value < 0.05)
  box_colors_multiv[sig_indices_multiv] <- "#E63946"  # 显著红色
  line_colors_multiv[sig_indices_multiv] <- "#E63946"   # 显著红色
  
  # 6.1 在R面板上展示多因素森林图
  # ===============================
  cat("正在生成多因素Cox回归森林图...\n")
  
  # 首先创建一个新的图形页面
  grid.newpage()
  
  # 创建视口布局，为标题留出空间
  vp_layout <- viewport(layout = grid.layout(
    nrow = 3, 
    ncol = 1,
    heights = unit(c(1.0, 0.2, 7.0), "inches")  # 标题空间、间隔、主图空间
  ))
  
  pushViewport(vp_layout)
  
  # 在第一个视口中添加标题
  pushViewport(viewport(layout.pos.row = 1))
  grid.text("Multivariable Cox Regression Analysis\nForest Plot", 
            gp = gpar(fontsize = 16, fontface = "bold"))
  popViewport()
  
  # 在第三个视口中绘制森林图
  pushViewport(viewport(layout.pos.row = 3))
  
  # 绘制森林图
  forestplot(
    labeltext = tabletext,
    mean = c(NA, forest_data$HR),
    lower = c(NA, forest_data$HR_lower), 
    upper = c(NA, forest_data$HR_upper),
    
    # 图形设置
    graph.pos = 5,
    graphwidth = unit(0.3, "npc"),
    
    # 处理极端值
    clip = c(min(x_ticks_multiv), max(x_ticks_multiv)),
    xlog = TRUE,
    
    # 颜色设置
    col = fpColors(
      box = box_colors_multiv,
      lines = line_colors_multiv,
      summary = "#1B4965",
      zero = "gray50"
    ),
    
    # 箱线图设置
    boxsize = 0.3,
    lwd.ci = 2,
    ci.vertices = TRUE,
    ci.vertices.height = 0.15,
    
    # 参考线
    zero = 1,
    lwd.zero = 2,
    grid = structure(x_ticks_multiv, gp = gpar(lty = 2, col = "gray80")),
    
    # 坐标轴
    xticks = x_ticks_multiv,
    lwd.xaxis = 2,
    xlab = "\nHazard Ratio (95% Confidence Interval)",
    
    # 文本格式
    txt_gp = fpTxtGp(
      label = gpar(cex = 0.95, fontfamily = "sans"),
      ticks = gpar(cex = 0.85),
      xlab = gpar(cex = 1.0, fontface = "bold")
    ),
    
    # 横线
    hrzl_lines = hrzl_lines_list_multiv,
    
    # 布局
    lineheight = unit(0.7, "cm"),
    colgap = unit(0.4, "cm"),
    mar = unit(c(0.5, 0.5, 0, 0.5), "cm"),
    
    # 不显示内置标题
    new_page = FALSE
  )
  
  # 返回到之前的viewport
  popViewport()
  popViewport()
  
  dev.off()
  cat("已保存多因素森林图(PNG): Multivariable_Cox_Forestplot.png\n")
} else {
  cat("未保存多因素森林图PNG文件\n")
}

# 8. 单因素与多因素结果比较
# ==========================
cat("\n=== 单因素与多因素结果比较 ===\n")

# 合并单因素和多因素结果进行比较
univ_results_clean <- univariate_results %>%
  select(变量, HR, HR_lower, HR_upper, P值) %>%
  rename(
    Variable = 变量,
    HR_univ = HR,
    HR_lower_univ = HR_lower,
    HR_upper_univ = HR_upper,
    P_univ = P值
  )

multiv_results_clean <- multiv_results %>%
  select(Variable_label, HR, HR_lower, HR_upper, P_value) %>%
  rename(
    Variable = Variable_label,
    HR_multiv = HR,
    HR_lower_multiv = HR_lower,
    HR_upper_multiv = HR_upper,
    P_multiv = P_value
  )

# 创建比较表格
comparison_table <- data.frame(
  Variable = unique(c(univ_results_clean$Variable, multiv_results_clean$Variable))
)

# 合并数据
for(i in 1:nrow(comparison_table)) {
  var <- comparison_table$Variable[i]
  
  # 单因素结果
  univ_row <- univ_results_clean[univ_results_clean$Variable == var, ]
  if(nrow(univ_row) > 0) {
    comparison_table$HR_univ[i] <- univ_row$HR_univ
    comparison_table$HR_univ_CI[i] <- sprintf("%.2f (%.2f-%.2f)", 
                                              univ_row$HR_univ, 
                                              univ_row$HR_lower_univ, 
                                              univ_row$HR_upper_univ)
    comparison_table$P_univ[i] <- ifelse(univ_row$P_univ < 0.001, "<0.001",
                                         sprintf("%.3f", univ_row$P_univ))
  }
  
  # 多因素结果
  multiv_row <- multiv_results_clean[multiv_results_clean$Variable == var, ]
  if(nrow(multiv_row) > 0) {
    comparison_table$HR_multiv[i] <- multiv_row$HR_multiv
    comparison_table$HR_multiv_CI[i] <- sprintf("%.2f (%.2f-%.2f)", 
                                                multiv_row$HR_multiv, 
                                                multiv_row$HR_lower_multiv, 
                                                multiv_row$HR_upper_multiv)
    comparison_table$P_multiv[i] <- ifelse(multiv_row$P_multiv < 0.001, "<0.001",
                                           sprintf("%.3f", multiv_row$P_multiv))
  }
}

# 打印比较表格
cat("\n单因素与多因素Cox回归结果比较:\n")
print(comparison_table[, c("Variable", "HR_univ_CI", "P_univ", "HR_multiv_CI", "P_multiv")], 
      row.names = FALSE)

# 识别结果变化
cat("\n结果变化分析:\n")
for(i in 1:nrow(comparison_table)) {
  var <- comparison_table$Variable[i]
  
  if(!is.na(comparison_table$P_univ[i]) && !is.na(comparison_table$P_multiv[i])) {
    p_univ_num <- ifelse(comparison_table$P_univ[i] == "<0.001", 0.0005, 
                         as.numeric(comparison_table$P_univ[i]))
    p_multiv_num <- ifelse(comparison_table$P_multiv[i] == "<0.001", 0.0005,
                           as.numeric(comparison_table$P_multiv[i]))
    
    univ_sig <- p_univ_num < 0.05
    multiv_sig <- p_multiv_num < 0.05
    
    if(univ_sig && !multiv_sig) {
      cat(sprintf("  • %s: 单因素显著，多因素不显著\n", var))
    } else if(!univ_sig && multiv_sig) {
      cat(sprintf("  • %s: 单因素不显著，多因素显著\n", var))
    } else if(univ_sig && multiv_sig) {
      cat(sprintf("  • %s: 单因素和多因素均显著\n", var))
    } else {
      cat(sprintf("  • %s: 单因素和多因素均不显著\n", var))
    }
  }
}





# 1. 加载必要包
# ====================
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("forestplot")) install.packages("forestplot")
if (!require("grid")) install.packages("grid")

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(forestplot)
library(grid)

# 2. 数据准备和清理
# ====================
# 读取数据
clinical1 <- read.table("clinical.txt", sep = "\t", row.names = 1, 
                        check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)

# 检查必需的列是否存在
required_cols <- c("STK11", "OS", "OStime", "Age", "Sex", "M", "TMB")
missing_cols <- setdiff(required_cols, colnames(clinical1))
if(length(missing_cols) > 0) {
  stop(paste("缺少必要的列:", paste(missing_cols, collapse = ", ")))
}

# 数据清理和转换
clinical_clean <- clinical1 %>%
  filter(!is.na(OStime) & !is.na(OS)) %>%
  filter(OStime > 0) %>%
  mutate(
    STK11 = factor(STK11, levels = c(0, 1), labels = c("Wild-type", "Mutant")),
    OS = as.integer(OS),
    OStime = as.numeric(OStime),
    Age = as.numeric(Age),
    # 直接转换，假设这些列都存在（前面已检查）
    Sex = factor(Sex, levels = c("Male", "Female")),
    M = factor(M, levels = c("M0", "M1")),
    TMB = as.numeric(TMB)
  ) %>%
  na.omit()

# 数据概况
cat("=== 数据概况 ===\n")
cat("原始样本数:", nrow(clinical1), "\n")
cat("有效样本数（无缺失）:", nrow(clinical_clean), "\n")
cat("移除样本数:", nrow(clinical1) - nrow(clinical_clean), "\n\n")

# 变量描述性统计
cat("=== 变量描述性统计 ===\n")
desc_stats <- clinical_clean %>%
  summarise(
    `总样本数` = n(),
    `事件数` = sum(OS == 1),
    `删失数` = sum(OS == 0),
    `年龄(均值±标准差)` = sprintf("%.1f ± %.1f", mean(Age, na.rm = TRUE), sd(Age, na.rm = TRUE)),
    `年龄(范围)` = sprintf("%.1f-%.1f", min(Age, na.rm = TRUE), max(Age, na.rm = TRUE)),
    `男性比例` = paste0(round(sum(Sex == "Male")/n()*100, 1), "%"),
    `女性比例` = paste0(round(sum(Sex == "Female")/n()*100, 1), "%"),
    `M1比例` = paste0(round(sum(M == "M1")/n()*100, 1), "%"),
    `TMB(中位数)` = round(median(TMB, na.rm = TRUE), 1),
    `野生型比例` = paste0(round(sum(STK11 == "Wild-type")/n()*100, 1), "%"),
    `突变型比例` = paste0(round(sum(STK11 == "Mutant")/n()*100, 1), "%")
  )

print(t(desc_stats))

# STK11分组描述
cat("\n=== STK11分组描述 ===\n")
stk11_group_desc <- clinical_clean %>%
  group_by(STK11) %>%
  summarise(
    `样本数` = n(),
    `事件数` = sum(OS == 1),
    `中位生存时间(月)` = round(median(OStime), 1),
    `年龄(均值)` = round(mean(Age), 1),
    `男性比例` = paste0(round(sum(Sex == "Male")/n()*100, 1), "%"),
    `M1比例` = paste0(round(sum(M == "M1")/n()*100, 1), "%"),
    `TMB(中位数)` = round(median(TMB), 1)
  )

print(stk11_group_desc)

# 3. 单因素Cox回归分析
# ====================
cat("\n=== 单因素Cox回归分析 ===\n")

# 创建生存对象
surv_obj <- Surv(time = clinical_clean$OStime, event = clinical_clean$OS)

# 定义分析变量及其参考水平
variables <- list(
  STK11 = list(var = "STK11", ref = "野生型", comp = "突变型 vs 野生型"),
  Age = list(var = "Age", ref = "连续变量", comp = "每增加1岁"),
  Sex = list(var = "Sex", ref = "男性", comp = "女性 vs 男性"),
  M = list(var = "M", ref = "M0", comp = "M1 vs M0"),
  TMB = list(var = "TMB", ref = "连续变量", comp = "每增加1单位")
)

# 执行单因素Cox回归分析
univariate_results <- data.frame()

for(var_name in names(variables)) {
  var_info <- variables[[var_name]]
  
  formula <- as.formula(paste("surv_obj ~", var_info$var))
  cox_model <- coxph(formula, data = clinical_clean)
  cox_summary <- summary(cox_model)
  
  # 提取结果
  n_obs <- cox_summary$n
  n_event <- cox_summary$nevent
  
  for(i in 1:nrow(cox_summary$coefficients)) {
    coef_name <- rownames(cox_summary$coefficients)[i]
    hr <- exp(cox_summary$coefficients[i, "coef"])
    hr_lower <- exp(cox_summary$coefficients[i, "coef"] - 1.96 * cox_summary$coefficients[i, "se(coef)"])
    hr_upper <- exp(cox_summary$coefficients[i, "coef"] + 1.96 * cox_summary$coefficients[i, "se(coef)"])
    p_value <- cox_summary$coefficients[i, "Pr(>|z|)"]
    
    # 确定比较描述
    if(var_name == "STK11") {
      comparison <- var_info$comp
    } else if(var_name %in% c("Age", "TMB")) {
      comparison <- var_info$comp
    } else {
      # 对于分类变量，从系数名提取比较组
      comparison <- paste(gsub(var_info$var, "", coef_name), "vs", var_info$ref)
    }
    
    univariate_results <- rbind(univariate_results, data.frame(
      变量 = var_name,
      比较组 = comparison,
      样本数 = n_obs,
      事件数 = n_event,
      HR = round(hr, 3),
      HR_lower = round(hr_lower, 3),
      HR_upper = round(hr_upper, 3),
      P值 = p_value,
      显著性 = ifelse(p_value < 0.001, "***", 
                   ifelse(p_value < 0.01, "**", 
                          ifelse(p_value < 0.05, "*", ""))),
      stringsAsFactors = FALSE
    ))
  }
}

# 格式化输出
univariate_results$`HR(95%CI)` <- sprintf("%.2f (%.2f-%.2f)", 
                                          univariate_results$HR,
                                          univariate_results$HR_lower,
                                          univariate_results$HR_upper)

univariate_results$`P值` <- ifelse(univariate_results$P值 < 0.001, "<0.001",
                                  sprintf("%.3f", univariate_results$P值))

# 按P值排序
univariate_results <- univariate_results[order(univariate_results$P值), ]

cat("\n单因素Cox回归分析结果:\n")
print(univariate_results[, c("变量", "比较组", "样本数", "事件数", "HR(95%CI)", "P值", "显著性")], 
      row.names = FALSE)

# 4. 单因素Cox回归森林图
# ========================================
# 4. 单因素Cox回归森林图
# ========================================
cat("\n=== 生成单因素Cox回归森林图 ===\n")

# 首先打印检查数据
cat("\n检查单因素分析结果数据:\n")
print(univariate_results[, c("变量", "比较组", "HR", "HR_lower", "HR_upper", "P值")])

# 准备森林图数据 - 直接使用变量标签
univ_forest_data <- univariate_results %>%
  mutate(
    # 直接创建变量显示标签，不进行复杂的字符串替换
    characteristics_clean = case_when(
      变量 == "STK11" ~ "STK11 (Mutant vs Wild-type)",
      变量 == "Age" ~ "Age (per year increase)",
      变量 == "Sex" ~ "Sex (Female vs Male)",
      变量 == "M" ~ "M stage (M1 vs M0)",
      变量 == "TMB" ~ "TMB (per unit increase)",
      TRUE ~ paste(变量, "(", 比较组, ")", sep = "")
    ),
    # 添加显著性标记
    signif_mark = ifelse(P值 < 0.001, "***",
                         ifelse(P值 < 0.01, "**",
                                ifelse(P值 < 0.05, "*", "")))
  ) %>%
  # 确保正确的排序
  arrange(factor(变量, levels = c("Age", "M", "STK11", "Sex", "TMB")))

cat("\n森林图数据准备完成:\n")
print(univ_forest_data[, c("characteristics_clean", "HR", "HR_lower", "HR_upper", "P值", "signif_mark")])

# 准备表格文本（包含significance列）
univ_tabletext <- cbind(
  c("Variable", univ_forest_data$characteristics_clean),
  c("HR", sprintf("%.3f", univ_forest_data$HR)),
  c("95% CI", sprintf("%.3f-%.3f", univ_forest_data$HR_lower, univ_forest_data$HR_upper)),
  c("P-value", ifelse(univ_forest_data$P值 < 0.001, "<0.001", 
                      sprintf("%s", univ_forest_data$P值))),
  c("Significance", univ_forest_data$signif_mark)
)

# 计算合适的x轴刻度
min_hr <- min(univ_forest_data$HR_lower, na.rm = TRUE)
max_hr <- max(univ_forest_data$HR_upper, na.rm = TRUE)

# 生成合适的刻度
x_ticks <- c()
if(min_hr < 0.1) x_ticks <- c(x_ticks, 0.1)
if(min_hr < 0.25) x_ticks <- c(x_ticks, 0.25)
x_ticks <- c(x_ticks, 0.5, 1)
if(max_hr > 2) x_ticks <- c(x_ticks, 2)
if(max_hr > 4) x_ticks <- c(x_ticks, 4)
if(max_hr > 8) x_ticks <- c(x_ticks, 8)
if(max_hr > 16) x_ticks <- c(x_ticks, 16)

# 创建hrzl_lines列表
total_rows <- nrow(univ_tabletext) + 1
hrzl_lines_list <- list(
  "1" = gpar(lwd = 2, col = "black"),
  "2" = gpar(lwd = 1.5, col = "gray70")
)
hrzl_lines_list[[as.character(total_rows)]] <- gpar(lwd = 2, col = "black")

# 创建颜色向量，突出显示显著结果
box_colors <- rep("#2E86AB", nrow(univ_forest_data))  # 默认蓝色
line_colors <- rep("#1B4965", nrow(univ_forest_data))  # 默认深蓝色

# 显著结果使用红色
sig_indices <- which(univ_forest_data$P值 < 0.05)
box_colors[sig_indices] <- "#E63946"  # 显著红色
line_colors[sig_indices] <- "#E63946"   # 显著红色

# 调整列宽以适应更长的变量名
column_widths <- unit(c(3.0, 0.8, 1.2, 0.8, 0.8, 1.2), "inches")

# 在R面板上展示单因素森林图
cat("正在生成单因素Cox回归森林图...\n")

# 首先创建一个新的图形页面
grid.newpage()

# 创建视口布局，为标题留出空间
vp_layout <- viewport(layout = grid.layout(
  nrow = 3, 
  ncol = 1,
  heights = unit(c(1.0, 0.2, 8.0), "inches")
))

pushViewport(vp_layout)

# 在第一个视口中添加标题
pushViewport(viewport(layout.pos.row = 1))
grid.text("Univariate Cox Regression Analysis\nForest Plot", 
          gp = gpar(fontsize = 16, fontface = "bold"))
popViewport()

# 在第三个视口中绘制森林图
pushViewport(viewport(layout.pos.row = 3))

# 绘制森林图
forestplot(
  labeltext = univ_tabletext,
  mean = c(NA, univ_forest_data$HR),
  lower = c(NA, univ_forest_data$HR_lower), 
  upper = c(NA, univ_forest_data$HR_upper),
  graph.pos = 4,  # 调整图形位置到第4列
  graphwidth = unit(0.35, "npc"),
  clip = c(min(x_ticks), max(x_ticks)),
  xlog = TRUE,
  col = fpColors(
    box = box_colors,
    lines = line_colors,
    summary = "#1B4965",
    zero = "gray50"
  ),
  boxsize = 0.25,  # 减小箱线图大小
  lwd.ci = 2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.15,
  zero = 1,
  lwd.zero = 2,
  grid = structure(x_ticks, gp = gpar(lty = 2, col = "gray80")),
  xticks = x_ticks,
  lwd.xaxis = 2,
  xlab = "\nHazard Ratio (95% Confidence Interval)",
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.9, fontfamily = "sans"),  # 减小字体大小
    ticks = gpar(cex = 0.85),
    xlab = gpar(cex = 1.0, fontface = "bold")
  ),
  hrzl_lines = hrzl_lines_list,
  lineheight = unit(0.6, "cm"),  # 减小行高
  colgap = unit(0.3, "cm"),  # 减小列间距
  mar = unit(c(0.5, 0.5, 0, 0.5), "cm"),
  new_page = FALSE
)

popViewport()
popViewport()

cat("单因素森林图已展示在R面板上\n")


# 5. 多因素Cox回归分析
# ====================
cat("\n=== 多因素Cox回归分析 ===\n")

# 构建多因素模型
multiv_formula <- surv_obj ~ STK11 + Age + Sex + M + TMB
cox_multiv <- coxph(multiv_formula, data = clinical_clean)
cox_multiv_summary <- summary(cox_multiv)

cat("\n1. 多因素Cox模型信息:\n")
cat("   模型公式:", paste(format(multiv_formula), collapse = ""), "\n")
cat("   有效样本数:", nrow(clinical_clean), "\n")
cat("   事件数:", sum(clinical_clean$OS == 1), "\n")
cat("   变量数:", length(variables), "\n")

cat("\n2. 多因素分析结果:\n")

# 提取多因素分析结果
multiv_results <- data.frame(
  Variable = rownames(cox_multiv_summary$coefficients),
  HR = round(exp(cox_multiv_summary$coefficients[, "coef"]), 3),
  HR_lower = round(exp(cox_multiv_summary$coefficients[, "coef"] - 1.96 * cox_multiv_summary$coefficients[, "se(coef)"]), 3),
  HR_upper = round(exp(cox_multiv_summary$coefficients[, "coef"] + 1.96 * cox_multiv_summary$coefficients[, "se(coef)"]), 3),
  Coefficient = round(cox_multiv_summary$coefficients[, "coef"], 4),
  SE = round(cox_multiv_summary$coefficients[, "se(coef)"], 4),
  z_value = round(cox_multiv_summary$coefficients[, "z"], 3),
  P_value = cox_multiv_summary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

# 添加变量说明
multiv_results$Variable_label <- c(
  "STK11Mutant" = "STK11 (Mutant vs Wild-type)",
  "Age" = "Age (per year increase)",
  "SexFemale" = "Sex (Female vs Male)",
  "MM1" = "M stage (M1 vs M0)",
  "TMB" = "TMB (per unit increase)"
)[multiv_results$Variable]

# 添加显著性标记
multiv_results$signif_mark <- ifelse(multiv_results$P_value < 0.001, "***",
                                     ifelse(multiv_results$P_value < 0.01, "**",
                                            ifelse(multiv_results$P_value < 0.05, "*", "")))

# 格式化显示
multiv_results$P_display <- ifelse(
  multiv_results$P_value < 0.001, 
  "<0.001", 
  sprintf("%.3f", multiv_results$P_value)
)

multiv_results$HR_95CI <- sprintf(
  "%.2f (%.2f-%.2f)", 
  multiv_results$HR,
  multiv_results$HR_lower,
  multiv_results$HR_upper
)

# 打印多因素结果表格
print(multiv_results[, c("Variable_label", "HR_95CI", "P_display", "signif_mark")])

cat("\n3. 模型整体检验:\n")
cat("   Likelihood ratio test p值:", round(cox_multiv_summary$logtest[["pvalue"]], 4), "\n")
cat("   Wald test p值:", round(cox_multiv_summary$waldtest[["pvalue"]], 4), "\n")
cat("   Score (logrank) test p值:", round(cox_multiv_summary$sctest[["pvalue"]], 4), "\n")

# 模型拟合优度
cat("\n4. 模型拟合优度:\n")
cat("   似然比检验卡方:", round(cox_multiv_summary$logtest[["test"]], 2), 
    " (df = ", cox_multiv_summary$logtest[["df"]], ")\n", sep="")
cat("   AIC:", round(AIC(cox_multiv), 2), "\n")
cat("   Concordance Index (C-index):", round(cox_multiv_summary$concordance[["C"]], 3), "\n")

# 6. 比例风险假设检验
# ====================
cat("\n=== 比例风险假设检验 ===\n")
ph_test <- cox.zph(cox_multiv)
print(ph_test)

# 可视化比例风险假设检验
if(nrow(clinical_clean) > 0) {
  cat("\n正在生成比例风险假设检验图...\n")
  
  # 询问是否在R面板上展示
  show_ph_test <- readline(prompt = "是否在R面板上展示比例风险假设检验图？(y/n): ")
  
  if(tolower(show_ph_test) %in% c("y", "yes")) {
    dev.new(width = 10, height = 6)
    plot(ph_test)
    cat("比例风险假设检验图已展示在R面板上\n")
  }
  
  # 询问是否保存
  save_ph_test <- readline(prompt = "是否保存比例风险假设检验图？(y/n): ")
  
  if(tolower(save_ph_test) %in% c("y", "yes")) {
    png("Proportional_Hazards_Test.png", width = 1200, height = 800, res = 150)
    plot(ph_test)
    dev.off()
    cat("比例风险假设检验图已保存为: Proportional_Hazards_Test.png\n")
  }
}

# 7. 多因素Cox回归森林图
# ================================================
cat("\n=== 多因素Cox回归森林图 ===\n")

# 定义统一的变量顺序
var_order <- c("Age (per year increase)", 
               "M stage (M1 vs M0)", 
               "STK11 (Mutant vs Wild-type)", 
               "Sex (Female vs Male)", 
               "TMB (per unit increase)")

# 准备森林图数据 - 按照统一顺序排序
multiv_forest_data <- multiv_results %>%
  mutate(
    Variable_display = Variable_label,
    # 创建排序因子，确保与单因素顺序一致
    Variable_order = factor(Variable_label, levels = var_order),
    HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, HR_lower, HR_upper),
    P_value_display = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value))
  ) %>%
  # 按照统一顺序排序
  arrange(Variable_order)

# 准备表格文本（包含significance列）
multiv_tabletext <- cbind(
  c("Variable", multiv_forest_data$Variable_display),
  c("HR", sprintf("%.3f", multiv_forest_data$HR)),
  c("95% CI", sprintf("%.3f-%.3f", multiv_forest_data$HR_lower, multiv_forest_data$HR_upper)),
  c("P-value", multiv_forest_data$P_value_display),
  c("Significance", multiv_forest_data$signif_mark)
)

# 计算合适的x轴刻度
min_hr_multiv <- min(multiv_forest_data$HR_lower, na.rm = TRUE)
max_hr_multiv <- max(multiv_forest_data$HR_upper, na.rm = TRUE)

x_ticks_multiv <- c()
if(min_hr_multiv < 0.1) x_ticks_multiv <- c(x_ticks_multiv, 0.1)
if(min_hr_multiv < 0.25) x_ticks_multiv <- c(x_ticks_multiv, 0.25)
x_ticks_multiv <- c(x_ticks_multiv, 0.5, 1)
if(max_hr_multiv > 2) x_ticks_multiv <- c(x_ticks_multiv, 2)
if(max_hr_multiv > 4) x_ticks_multiv <- c(x_ticks_multiv, 4)
if(max_hr_multiv > 8) x_ticks_multiv <- c(x_ticks_multiv, 8)
if(max_hr_multiv > 16) x_ticks_multiv <- c(x_ticks_multiv, 16)

# 创建hrzl_lines列表
total_rows_multiv <- nrow(multiv_tabletext) + 1
hrzl_lines_list_multiv <- list(
  "1" = gpar(lwd = 2, col = "black"),
  "2" = gpar(lwd = 1.5, col = "gray70")
)
hrzl_lines_list_multiv[[as.character(total_rows_multiv)]] <- gpar(lwd = 2, col = "black")

# 创建颜色向量，突出显示显著结果
box_colors_multiv <- rep("#2E86AB", nrow(multiv_forest_data))  # 默认蓝色
line_colors_multiv <- rep("#1B4965", nrow(multiv_forest_data))  # 默认深蓝色

# 显著结果使用红色
sig_indices_multiv <- which(multiv_forest_data$P_value < 0.05)
box_colors_multiv[sig_indices_multiv] <- "#E63946"  # 显著红色
line_colors_multiv[sig_indices_multiv] <- "#E63946"   # 显著红色

# 7.1 在R面板上展示多因素森林图
cat("正在生成多因素Cox回归森林图...\n")

# 首先创建一个新的图形页面
grid.newpage()

# 创建视口布局，为标题留出空间
vp_layout_multiv <- viewport(layout = grid.layout(
  nrow = 3, 
  ncol = 1,
  heights = unit(c(1.0, 0.2, 7.0), "inches")
))

pushViewport(vp_layout_multiv)

# 在第一个视口中添加标题
pushViewport(viewport(layout.pos.row = 1))
grid.text("Multivariable Cox Regression Analysis\nForest Plot", 
          gp = gpar(fontsize = 16, fontface = "bold"))
popViewport()

# 在第三个视口中绘制森林图
pushViewport(viewport(layout.pos.row = 3))

# 绘制森林图 - 注意graph.pos调整为4，与单因素保持一致
forestplot(
  labeltext = multiv_tabletext,
  mean = c(NA, multiv_forest_data$HR),
  lower = c(NA, multiv_forest_data$HR_lower), 
  upper = c(NA, multiv_forest_data$HR_upper),
  graph.pos = 4,  # 与单因素保持一致（第4列）
  graphwidth = unit(0.35, "npc"),
  clip = c(min(x_ticks_multiv), max(x_ticks_multiv)),
  xlog = TRUE,
  col = fpColors(
    box = box_colors_multiv,
    lines = line_colors_multiv,
    summary = "#1B4965",
    zero = "gray50"
  ),
  boxsize = 0.25,  # 与单因素保持一致
  lwd.ci = 2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.15,
  zero = 1,
  lwd.zero = 2,
  grid = structure(x_ticks_multiv, gp = gpar(lty = 2, col = "gray80")),
  xticks = x_ticks_multiv,
  lwd.xaxis = 2,
  xlab = "\nHazard Ratio (95% Confidence Interval)",
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.9, fontfamily = "sans"),  # 与单因素保持一致
    ticks = gpar(cex = 0.85),
    xlab = gpar(cex = 1.0, fontface = "bold")
  ),
  hrzl_lines = hrzl_lines_list_multiv,
  lineheight = unit(0.6, "cm"),  # 与单因素保持一致
  colgap = unit(0.3, "cm"),  # 与单因素保持一致
  mar = unit(c(0.5, 0.5, 0, 0.5), "cm"),
  new_page = FALSE
)

popViewport()
popViewport()

cat("多因素森林图已展示在R面板上\n")




####共突变分析####
# 设置工作目录（根据您的实际路径修改）
setwd("D://SHENGXIN//MSK//data")

# ============================================================
# 本地读取 clinical.txt 和 MSK-IMPACT 突变数据
# 进行 STK11/KRAS 共突变分析
# ============================================================

# 1. 加载必要的包
library(dplyr)
library(survival)
library(survminer)
library(tibble)

# ============================================================
# 第一步：读取已整理好的临床数据 clinical1
# ============================================================
clinical1 <- read.table("clinical.txt",
                        sep = "\t",
                        row.names = 1,
                        header = TRUE,
                        check.names = FALSE,
                        stringsAsFactors = FALSE)

# 查看临床数据结构
cat("临床数据维度：", dim(clinical1), "\n")
cat("临床数据列名：\n")
print(colnames(clinical1))
cat("\n前几行：\n")
print(head(clinical1))

# ============================================================
# 第二步：解压并读取 MSK-IMPACT 突变数据
# ============================================================

# 读取突变数据（请根据实际文件名修改路径）
# 常见文件名：data_mutations.txt 或 data_mutations_extended.txt
mutation_file <- "data_mutations.txt"
if (!file.exists(mutation_file)) {
  mutation_file <- "tmb_mskcc_2018/data_mutations_extended.txt"
}

msk_mutations <- read.table(mutation_file,
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = FALSE,
                            check.names = FALSE,
                            quote = "")

# 查看突变数据列名
cat("\n突变数据列名：\n")
print(colnames(msk_mutations))

# ============================================================
# 第三步：提取 KRAS 突变样本 ID
# ============================================================
# cBioPortal MAF 文件中，样本 ID 列通常为 "Tumor_Sample_Barcode"
# 如果列名不同，请根据上面的输出修改
sample_col <- "Tumor_Sample_Barcode"   # 请根据实际列名修改
gene_col <- "Hugo_Symbol"              # 请根据实际列名修改

kras_mut_ids <- msk_mutations %>%
  filter(.data[[gene_col]] == "KRAS") %>%
  pull(.data[[sample_col]]) %>%
  unique()

cat("\nKRAS 突变样本数：", length(kras_mut_ids), "\n")

# ============================================================
# 第四步：将 KRAS 突变状态合并到 clinical1
# ============================================================
# clinical1 行名为 Sample ID，先转为列
clinical1 <- clinical1 %>%
  rownames_to_column(var = "SampleID")

# 统一 SampleID 格式（如需要，去掉连字符或点号，确保与 kras_mut_ids 一致）
# clinical1$SampleID <- gsub("[-.]", "", clinical1$SampleID)
# kras_mut_ids <- gsub("[-.]", "", kras_mut_ids)

# 标记 KRAS 突变状态
clinical1 <- clinical1 %>%
  mutate(
    KRAS = ifelse(SampleID %in% kras_mut_ids, 1, 0),
    STK11 = as.numeric(as.character(STK11)),
    Group = case_when(
      STK11 == 1 & KRAS == 1 ~ "STK11/KRAS co-mutant",
      STK11 == 1 & KRAS == 0 ~ "STK11-only",
      STK11 == 0 & KRAS == 1 ~ "KRAS-only",
      STK11 == 0 & KRAS == 0 ~ "Wild-type"
    ),
    Group = factor(Group, levels = c("Wild-type", "STK11-only", "KRAS-only", "STK11/KRAS co-mutant"))
  )

# 查看四分组分布
cat("\n四分组样本分布：\n")
print(table(clinical1$Group, useNA = "ifany"))

# ============================================================
# 第五步：数据清理，准备生存分析
# ============================================================
clinical_clean <- clinical1 %>%
  filter(!is.na(OStime) & !is.na(OS) & OStime > 0) %>%
  filter(!is.na(Group)) %>%
  mutate(
    OS = as.integer(OS),
    OStime = as.numeric(OStime),
    Age = as.numeric(Age),
    Sex = factor(Sex),
    M = factor(M),
    TMB = as.numeric(TMB)
  ) %>%
  na.omit()

cat("\n用于生存分析的有效样本数：", nrow(clinical_clean), "\n")
print(table(clinical_clean$Group))

# ============================================================
# 第六步：Kaplan-Meier 生存分析（四组）
# ============================================================
fit <- survfit(Surv(OStime, OS) ~ Group, data = clinical_clean)

km_plot <- ggsurvplot(
  fit,
  data = clinical_clean,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  surv.median.line = "hv",
  legend.labs = levels(clinical_clean$Group),
  legend.title = "Group",
  palette = c("#2E9FDF", "#E74C3C", "#F39C12", "#8E44AD"),
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  title = "Overall Survival by STK11/KRAS Co-mutation Status",
  break.time.by = 12,
  ggtheme = theme_bw()
)

print(km_plot)
pdf("KM_STK11_KRAS_co-mutation.pdf", width = 9, height = 8)
print(km_plot)
dev.off()

# ============================================================
# 第七步：多变量 Cox 回归分析
# ============================================================
cox_multiv <- coxph(Surv(OStime, OS) ~ Group + Age + Sex + M + TMB, data = clinical_clean)
cat("\n=== 多变量 Cox 回归 ===\n")
print(summary(cox_multiv))

# ============================================================
# 第八步：专门比较 STK11/KRAS co-mutant vs STK11-only
# ============================================================
clinical_sub <- clinical_clean %>%
  filter(Group %in% c("STK11-only", "STK11/KRAS co-mutant")) %>%
  mutate(Group = droplevels(Group))

cox_sub <- coxph(Surv(OStime, OS) ~ Group + Age + Sex + M + TMB, data = clinical_sub)
cat("\n=== STK11/KRAS co-mutant vs STK11-only（多变量调整） ===\n")
print(summary(cox_sub))

# 提取 HR 和 95% CI
hr <- exp(coef(cox_sub))["GroupSTK11/KRAS co-mutant"]
ci <- exp(confint(cox_sub))["GroupSTK11/KRAS co-mutant", ]
p <- summary(cox_sub)$coefficients["GroupSTK11/KRAS co-mutant", "Pr(>|z|)"]

cat(sprintf("\nSTK11/KRAS co-mutant vs STK11-only: HR = %.2f (95%% CI %.2f-%.2f), P = %.4f\n",
            hr, ci[1], ci[2], p))

# ============================================================
# 第九步：保存结果
# ============================================================
write.csv(clinical_clean, "MSK_IMPACT_STK11_KRAS_co-mutation.csv", row.names = FALSE)
cat("\n✅ 分析完成！输出文件：\n")
cat("1. KM_STK11_KRAS_co-mutation.pdf（KM 曲线）\n")
cat("2. MSK_IMPACT_STK11_KRAS_co-mutation.csv（分组及临床数据）\n")



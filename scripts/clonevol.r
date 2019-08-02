library(clonevol)
library(viridis)

args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  return(res)
}

argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$help)) {help=FALSE} else {help=TRUE}

if(is.null(args$pyclone_res_folder) | is.null(args$tumor1_reformated) | is.null(args$tumor1_id) | is.null(args$tumor2_id) | help) {
  cat("

      reformat4pyclone.r: build a tsv file to input to pyclone from one multisample VCF from mutect2 and one CNV.txt file from facets 

      Mandatory arguments:
      --pyclone_res_folder=folder_path             - path to the pyclone output folder
      --tumor1_reformated=file                     - Reformated file for tumor 1 passed to pyclone
      --tumor1_id=id                               - ID of tumor 1
      --tumor2_id=id                               - ID of tumor 2

      Optional arguments:
      --kept_genes=list                            - list of kept genes separated by commas
      --help                                       - print this text \n\n")
  q(save="no")
}
tumor1_id = args$tumor1_id
tumor2_id = args$tumor2_id
pyclone_res = read.table(paste(args$pyclone_res_folder,"/tables/loci.tsv",sep=""), quote="\"", 
                         stringsAsFactors=F, sep="\t", header=T, comment.char = "")

# add the gene information
T1_reformated4pyclone = read.table(args$tumor1_reformated, quote="\"", stringsAsFactors=F, sep="\t", header=T, comment.char = "")

pyclone_res$gene = unlist(lapply(pyclone_res$mutation_id, function(m){
  T1_reformated4pyclone[which(T1_reformated4pyclone$mutation_id==m),"gene"]
}))
  
# recurrent genes
if(is.null(args$kept_genes)) {kept_genes = c()} else {kept_genes = unlist(strsplit(args$kept_genes,","))}

for(i in (1:nrow(pyclone_res))[c(T,F)]){
  cluster = pyclone_res[i,"cluster_id"]
  gene = pyclone_res[i,"gene"]
  VAF_T1 = pyclone_res[i,"cellular_prevalence"]
  VAF_T2 = pyclone_res[i+1,"cellular_prevalence"]
  is.recurrent = ifelse(gene %in% kept_genes, TRUE, FALSE)
  if(i==1) res = data.frame(cluster, VAF_T1, VAF_T2, gene, is.recurrent)
  if(i!=1) res = rbind(res, data.frame(cluster, VAF_T1, VAF_T2, gene, is.recurrent))
}

colnames(res)[which(colnames(res)=="VAF_T1")] = paste("VAF",tumor1_id,sep="_")
colnames(res)[which(colnames(res)=="VAF_T2")] = paste("VAF",tumor2_id,sep="_")
  
vaf.col.names <- colnames(res)[grepl("VAF",colnames(res))]
sample.groups <- vaf.col.names
names(sample.groups) <- vaf.col.names

# clonevol requires vaf in [0:100] scale
for(v in vaf.col.names){
  res[,v] = res[,v] * 100
}

# setup the order of clusters to display in various plots (later)
x <- res[order(res$cluster),]
x$cluster = x$cluster + 1
kept_clusters = names(table(x$cluster))[which(table(x$cluster)>=10)]
x = x[which(x$cluster %in% kept_clusters),]

# choosing color for the clones
clone.colors <- viridis(length(table(x$cluster))+1)[(length(table(x$cluster))+1):2] #c('#999793', '#8d4891', '#f8e356', '#fe9536', '#d7352e')

# infer clonal evolution trees
cs = 1:length(table(x$cluster))
names(cs) = sort(unique(x$cluster))
x$cluster = as.numeric(cs[as.character(x$cluster)])

founding.cluster = as.numeric(which.max(table(x$cluster)))
if(founding.cluster!=1){
  fc = x$cluster == founding.cluster
  c1 = x$cluster == 1
  x$cluster[c1] = founding.cluster
  x$cluster[fc] = 1
  x = x[order(x$cluster),]
}

# plot the trees 

y = infer.clonal.models(variants = x,
                        cluster.col.name = 'cluster',
                        vaf.col.names = vaf.col.names,
                        sample.groups = sample.groups,
                        cancer.initiation.model='monoclonal',
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 1000,
                        founding.cluster = 1,
                        cluster.center = 'mean',
                        ignore.clusters = NULL,
                        clone.colors = clone.colors,
                        min.cluster.vaf = 0.05,
                        # min probability that CCF(clone) is non-negative
                        sum.p = 0.05,
                        # alpha level in confidence interval estimate for CCF(clone)
                        alpha = 0.05)

if(nrow(x[x$is.recurrent,])>0){
  y <- transfer.events.to.consensus.trees(y,
                                          x[x$is.recurrent,],
                                          cluster.col.name = 'cluster',
                                          event.col.name = 'gene')
}

y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')

plot.clonal.models(y, out.dir = paste("output_pyclone",tumor1_id,sep="_"))


# plot the clusters
pdf(paste("output_pyclone_",tumor1_id,"/clusters.pdf",sep=""), width = 5, height = 5, useDingbats = FALSE, title='')
pp <- plot.variant.clusters(x,
                            cluster.col.name = 'cluster',
                            show.cluster.size = FALSE,
                            cluster.size.text.color = 'blue',
                            vaf.col.names = vaf.col.names,
                            vaf.limits = 100,
                            sample.title.size = 20,
                            violin = FALSE,
                            box = FALSE,
                            jitter = TRUE,
                            jitter.shape = 1,
                            jitter.color = clone.colors,
                            jitter.size = 3,
                            jitter.alpha = 1,
                            jitter.center.method = 'median',
                            jitter.center.size = 1,
                            jitter.center.color = 'darkgray',
                            jitter.center.display.value = 'none',
                            highlight = 'is.recurrent',
                            highlight.shape = 21,
                            highlight.color = 'blue',
                            highlight.fill.color = 'green',
                            highlight.note.col.name = 'gene',
                            highlight.note.size = 2,
                            order.by.total.vaf = FALSE)
dev.off()

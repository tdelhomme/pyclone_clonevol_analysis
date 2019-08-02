args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  return(res)
}

argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$help)) {help=FALSE} else {help=TRUE}

if(is.null(args$VCF) | is.null(args$CNV_tumor1) | is.null(args$CNV_tumor2) | is.null(args$normal_id) | is.null(args$tumor1_id)
   | is.null(args$tumor2_id)| help) {
  cat("

      reformat4pyclone.r: build a tsv file to input to pyclone from one multisample VCF from mutect2 and one CNV.txt file from facets 

      Mandatory arguments:
      --VCF=file_name             - multisample VCF file from mutect2
      --CNV_tumor1=file_name      - CNV file from facets tool for tumor 1
      --CNV_tumor2=file_name      - CNV file from facets tool for tumor 2
      --normal_id=id              - ID of the normal sample in VCF
      --tumor1_id=id              - ID of the tumor 1 sample in VCF
      --tumor2_id=id              - ID of the tumor 2 sample in VCF

      Optional arguments:
      --CNV_normal=file_name      - CNV file from facets tool for normal
      --min_dp=value              - minimum dp to consider a mutation (default: 20)
      --min_cnv_length=value      - minimum length in pb to consider a copy number (default:1)
      --output_folder=path        - folder to save the models (default: .)
      --chunk_size=value          - size of reading (default: 1000)
      --kept_genes=list           - list of kept genes separated by commas
      --help                      - print this text \n\n")
  q(save="no")
}

library(VariantAnnotation)

if(is.null(args$output_folder)) {output_folder="."} else {output_folder=args$output_folder}
if(is.null(args$chunk_size)) {chunk_size = 100000} else {chunk_size = as.numeric(args$chunk_size)}
if(is.null(args$min_dp)) {min_dp = 20} else {min_dp = as.numeric(args$min_dp)}
if(is.null(args$min_cnv_length)) {min_cnv_length = 1} else {min_cnv_length = as.numeric(args$min_cnv_length)}
if(is.null(args$kept_genes)) {kept_genes = c()} else {kept_genes = unlist(strsplit(args$kept_genes,","))}

normal_id = args$normal_id
tumor1_id = args$tumor1_id
tumor2_id = args$tumor2_id

CNV_tumor1 = read.table(args$CNV_tumor1, quote="\"", stringsAsFactors=F, sep="\t", header=T)
CNV_tumor2 = read.table(args$CNV_tumor2, quote="\"", stringsAsFactors=F, sep="\t", header=T)

# only consider CNV with a sufficient length to de-noise the data
CNV_tumor1$seglength = CNV_tumor1$end - CNV_tumor1$start + 1
CNV_tumor2$seglength = CNV_tumor2$end - CNV_tumor2$start + 1
CNV_tumor1 = CNV_tumor1[which(CNV_tumor1$seglength > min_cnv_length),]
CNV_tumor2 = CNV_tumor2[which(CNV_tumor2$seglength > min_cnv_length),]

VCF = read.table(args$VCF, quote="\"", stringsAsFactors=F, sep="\t", header=T)
if(!file.exists(paste(args$VCF,".tbi",sep=""))){
  system(paste("tabix -p vcf ", args$VCF, sep=""))
}

vcf <- open(VcfFile(args$VCF,  yieldSize=chunk_size))
vcf_chunk = readVcf(vcf, "hg19")

while(dim(vcf_chunk)[1] != 0) {
  muts = unlist(lapply(rownames(vcf_chunk), function(s) unlist(strsplit(s,"_"))[[2]]))
  vcf_chunk = vcf_chunk[which(nchar(muts)==3),]
  AD_matrix = geno(vcf_chunk,"AD")
  GT_matrix = geno(vcf_chunk,"GT")
  DP_matrix = geno(vcf_chunk,"DP")
  gene = unlist(info(vcf_chunk)$Gene.refGene)
  tumor1_mutid = tumor2_mutid = as.numeric(which( (GT_matrix[,normal_id]=="0/0" | GT_matrix[,normal_id]=="0|0") &
                                                ( (DP_matrix[,tumor1_id]>min_dp & DP_matrix[,tumor2_id]>min_dp)
                                                  | gene %in% kept_genes)))

  mutid=rownames(AD_matrix)
  # for the tumor 1 sample
  rc = unlist(lapply(AD_matrix[tumor1_mutid,tumor1_id], function(x) x[1]))
  vc = unlist(lapply(AD_matrix[tumor1_mutid,tumor1_id], function(x) x[2]))
  gene = unlist(info(vcf_chunk)$Gene.refGene)[tumor1_mutid]
  gene = unlist(lapply(gene, function(g) unlist(strsplit(g,"[\\]"))[1]))
  if(!exists("tumor1_reformated4pyclone")){
    tumor1_reformated4pyclone = data.frame(mutation_id=mutid[tumor1_mutid], ref_counts=rc, var_counts=vc, normal_cn=NA, minor_cn=NA, major_cn=NA, gene=gene)
  } else {
    tumor1_reformated4pyclone = rbind(tumor1_reformated4pyclone,
           data.frame(mutation_id=mutid[tumor1_mutid], ref_counts=rc, var_counts=vc, normal_cn=NA, minor_cn=NA, major_cn=NA, gene=gene))
  }
  # for the tumor 2 sample
  rc = unlist(lapply(AD_matrix[tumor2_mutid,tumor2_id], function(x) x[1]))
  vc = unlist(lapply(AD_matrix[tumor2_mutid,tumor2_id], function(x) x[2]))
  gene = unlist(info(vcf_chunk)$Gene.refGene)[tumor2_mutid]
  if(!exists("tumor2_reformated4pyclone")){
    tumor2_reformated4pyclone = data.frame(mutation_id=mutid[tumor2_mutid], ref_counts=rc, var_counts=vc, normal_cn=NA, minor_cn=NA, major_cn=NA, gene=gene)
  } else {
    tumor2_reformated4pyclone = rbind(tumor2_reformated4pyclone,
                                      data.frame(mutation_id=mutid[tumor2_mutid], ref_counts=rc, var_counts=vc, normal_cn=NA, minor_cn=NA, major_cn=NA, gene=gene))
  }
  vcf_chunk = readVcf(vcf, "hg19")
}

get_CN <-function(mutid, CNV_data, type="tumor"){
  mut_start = as.numeric(strsplit(strsplit(mutid, ":")[[1]][2],"_")[[1]][1])
  mut_chr = as.numeric(gsub("chr","",strsplit(mutid, ":")[[1]][1]))
  id = which(CNV_data$chrom == mut_chr & CNV_data$start <= mut_start & CNV_data$end>= mut_start)
  if(length(id)>0){
    if(type=="normal") {return(list(normal_cn = CNV_data[id,"tcn.em"]))}
    if(type=="tumor") {return(list(minor_cn = ifelse(is.na(CNV_data[id,"lcn.em"]),0,CNV_data[id,"lcn.em"]), 
                                   major_cn = ifelse(is.na(CNV_data[id,"lcn.em"]),CNV_data[id,"tcn.em"],
                                                     (CNV_data[id,"tcn.em"] - CNV_data[id,"lcn.em"]))))}
  } else {
    if(type=="normal") {return(list(normal_cn = 2))}
    if(type=="tumor") {return(list(minor_cn = 2, major_cn = 2))}
  }
}

if(is.null(args$CNV_normal)){tumor1_reformated4pyclone$normal_cn = tumor2_reformated4pyclone$normal_cn = 2} else{
  CNV_normal = read.table(args$CNV_normal, quote="\"", stringsAsFactors=F, sep="\t", header=T)
  res_n = get_CN(as.character(tumor1_reformated4pyclone[i,"mutation_id"]), CNV_normal, type="normal")
}

for(i in 1:nrow(tumor1_reformated4pyclone)){
  res = get_CN(as.character(tumor1_reformated4pyclone[i,"mutation_id"]), CNV_tumor1)
  tumor1_reformated4pyclone[i,"minor_cn"] = res$minor_cn
  tumor1_reformated4pyclone[i,"major_cn"] = res$major_cn
  if(!is.null(args$CNV_normal)) { tumor1_reformated4pyclone[i,"normal_cn"] = res_n$major_cn }
}

for(i in 1:nrow(tumor2_reformated4pyclone)){
  res = get_CN(as.character(tumor2_reformated4pyclone[i,"mutation_id"]), CNV_tumor2)
  tumor2_reformated4pyclone[i,"minor_cn"] = res$minor_cn
  tumor2_reformated4pyclone[i,"major_cn"] = res$major_cn
  if(!is.null(args$CNV_normal)) { tumor2_reformated4pyclone[i,"normal_cn"] = res_n$major_cn }
}

# for the moment remove mutations in region where CN was estimated as 0
rm_mut = unique(c(as.character(tumor1_reformated4pyclone[which(tumor1_reformated4pyclone$major_cn == 0),"mutation_id"]),
                  as.character(tumor2_reformated4pyclone[which(tumor2_reformated4pyclone$major_cn == 0),"mutation_id"])))
tumor1_reformated4pyclone = tumor1_reformated4pyclone[which(!tumor1_reformated4pyclone$mutation_id %in% rm_mut),]
tumor2_reformated4pyclone = tumor2_reformated4pyclone[which(!tumor2_reformated4pyclone$mutation_id %in% rm_mut),]
# keep only in the known chromosomes
tumor1_reformated4pyclone = tumor1_reformated4pyclone[which(grepl("chr",tumor1_reformated4pyclone$mutation_id)),]
tumor2_reformated4pyclone = tumor2_reformated4pyclone[which(grepl("chr",tumor2_reformated4pyclone$mutation_id)),]

write.table(tumor1_reformated4pyclone, file=paste(output_folder,"/",tumor1_id,"_reformated4pyclone.csv",sep=""), quote=F, sep="\t", row.names = F)
write.table(tumor2_reformated4pyclone, file=paste(output_folder,"/",tumor2_id,"_reformated4pyclone.csv",sep=""), quote=F, sep="\t", row.names = F)

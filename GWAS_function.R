library(rrBLUP)
library(BGLR)
library(DT)
library(SNPRelate)#BiocManager::install("SNPRelate")
library(dplyr)
library(qqman)
library(grid)
library(gridGraphics)
library(ggplot2)
library(readr)

plink_dir = '/home/starjjbang/an_proj/4_GWAS/plink'
bcftools_dir = '/home/starjjbang/an_proj/4_GWAS/bcftools/bcftools'
#/home/starjjbang/an_proj/4_GWAS/bcftools/bcftools isec -p /home/starjjbang/an_proj/4_GWAS/result/linear -n=2 G1-1.vcf.gz G1-2.vcf.gz
#/home/starjjbang/an_proj/4_GWAS/bcftools/bcftools stats G1-1.vcf.gz > cot.txt


qq_fun<-function(in_dat_f, re_fig_lab,cut_p){
  #re_dat = read_table(in_dat_f)
  re_dat = data.table::fread(in_dat_f)
  m(x = re_dat, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
    col = c("#BBC814",  "#0A4A61"), 
    suggestiveline = cut_p, 
    logp = TRUE,annotatePval = 5e-05,
    cex.lab=1.5,snp_font_size=0.8,snp_font_color  = "#FEA609")
  p <- recordPlot()
  g <- grid.grabExpr(grid.echo(p))
  ggsave(re_fig_lab, g, bg = "transparent", width = 10, height = 8,dpi=300)
}


runplink <- function(opts,re_f_name){
  q = paste0(plink_dir, opts, re_f_name)
  system(q)
}

runplink_q <- function(opts,re_f_name){
  q = paste0(plink_dir, opts, re_f_name)
  return(q)
}


vcf_merge_fun <-function(f_lst,vcf_dir,ref_name){
  in_f_dir = "/home/starjjbang/an_proj/4_GWAS/HN00159051/"
  for(i in 1:length(f_lst)){
    fn = gsub("-","-",f_lst[i])
    f_dir = paste0(in_f_dir,f_lst[i],'/',fn,'.final.vcf')
    q = paste0('cp ',f_dir, ' ',vcf_dir,fn,'.vcf')
    system(q)
    q1 = paste0('bgzip ',vcf_dir,fn,'.vcf')
    system(q1)
    q2 = paste0('tabix -p vcf ',vcf_dir,fn,'.vcf.gz')
    system(q2)
  }
  q3 = paste0('vcf-merge ',vcf_dir,'*.vcf.gz > ',vcf_dir,ref_name) 
  system(q3)  
  #q4 = paste0('bgzip ', vcf_dir,cot_lab,'_',case_lab,'.vcf.gz')
  #system(q4)
}
Geno_dat <-function(ped_f,fam_f,map_f,rds_dir){
  Geno <- read_ped(ped_f) 
  FAM <- read.csv(fam_f,header=FALSE,sep=' ')
  MAP <- read.table(map_f)
  p = Geno$p
  n = Geno$n
  Geno = Geno$x
  #.ped (conatins the marker allele data in 0, 1, 2, and 3 format; 2 represents missing data). 
  # The files are plink converted files.
  #1 represents heterozygous, 0 and 3 represents homozygous for major and minor allele.
  Geno[Geno == 2] <- NA  # Converting missing data to NA
  Geno[Geno == 0] <- 0  # Converting 0 data to 0
  Geno[Geno == 1] <- 1  # 1 represents heterozygous
  Geno[Geno == 3] <- 2  # Converting 3 to 2
  Geno <- matrix(Geno, nrow = p, ncol = n, byrow = TRUE)
  Geno <- t(Geno)
  rownames(Geno) <- FAM$V2
  dim(Geno); head(t(Geno))  
  saveRDS(Geno,paste0(rds_dir,'1_Geno.rds'))
  
  #y <- matrix(FAM$V6)  # # use the first trait 
  #rownames(y) <- FAM$V2
  #index <- !is.na(y)
  #y <- y[index, 1, drop = FALSE]  
  #Geno <- Geno[index, ] 
  return(Geno)
}

imputation_fun <-function(Geno){
  re_dat = Geno
  for (j in 1:ncol(re_dat)) {
    missing_data_cal = (length(which(is.na(re_dat[,j])))/length(re_dat[,j])) * 100
    if(missing_data_cal <= 10){
      re_dat[, j] <- ifelse(is.na(re_dat[, j]), mean(re_dat[, j], na.rm = TRUE), re_dat[, j]) 
    }
  }
  return(re_dat)
}

pca_run <-function(in_Geno,FAM,MAP,p_dir,gds_name){
  gds_name = "pca.gds"
  #remove duplication
  MAP_1 = MAP[-which(duplicated(MAP$V2)),]
  Geno_1 = in_Geno[,-which(duplicated(MAP$V2))]
  #Population structure
  # Create geno matrix file and assign the row and column names from fam and
  # map files
  Geno1 <- as.matrix(Geno_1); rownames(Geno_1) <- FAM$V2
  sample <- row.names(Geno_1); colnames(Geno_1) <- MAP_1$V2
  snp.id <- colnames(Geno_1)
  snpgdsCreateGeno(paste0(p_dir,gds_name), genmat = Geno1, sample.id = sample, snp.id = snp.id, 
                   snp.chromosome = MAP_1$V1, snp.position = MAP_1$V4, snpfirstdim = FALSE)
  geno_total <- snpgdsOpen(paste0(p_dir,gds_name))
  snpgdsSummary(paste0(p_dir,gds_name))
  pca1 <- snpgdsPCA(geno_total, snp.id = colnames(Geno1), num.thread=2)
  pca <- data.frame(sample.id = row.names(Geno1), 
                    EV1 = pca1$eigenvect[, 1], 
                    EV2 = pca1$eigenvect[, 2], 
                    EV3 = pca1$eigenvect[, 3], 
                    EV3 = pca1$eigenvect[, 4], 
                    stringsAsFactors = FALSE)
  snpgdsClose(geno_total)  
  return(pca)
}
# Plot the PCA
#par(mfrow = c(1,1))
#plot(pca_po$EV1, pca_po$EV2, xlab = "PC1", ylab = "PC2", col = c(1:4)[factor(pca_po$group)])
#legend(x = "topright", legend = levels(factor(pca_po$group)), col = c(1:4),pch = 1, cex = 0.6)

#pc.percent <- pca$varprop*100
#lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
#pairs(pca1$eigenvect[,1:4], col=c(1:4)[factor(pca_po$group)], labels=lbls)

#fig <- plot_ly(pca_po, x = ~EV1, y = ~EV2, z = ~EV3, color = ~pca_po$group, colors = c('#FBADAD','#FDFEB6','#E4F1EE',"#D9EDF8") )%>%
#  add_markers(size = 12)
#fig <- fig %>% layout(title = "test",scene = list(bgcolor = "white"))
#fig


num_lab = "2";cot_lab = "G1";case_lab = "G3"



GWAS_pair<-function(num_lab,cot_lab,case_lab){
  re_dir = paste0('/home/starjjbang/an_proj/4_GWAS/result/',num_lab,'_',cot_lab,'_',case_lab,'/')
  p_dir = paste0(re_dir, 'plink_result/')
  vcf_dir = paste0(re_dir,'vcf/')
  rds_dir = paste0(re_dir,'rds/')
  fig_dir = paste0(re_dir,'fig/')
  system(paste0('mkdir ',re_dir))
  system(paste0('mkdir ',p_dir))
  system(paste0('mkdir ',vcf_dir))
  system(paste0('mkdir ',rds_dir))
  system(paste0('mkdir ',fig_dir))
  total_vcf_dir = paste0(vcf_dir,cot_lab,'_',case_lab,'.vcf') 
  re_f_name = "total"
  ######################################
  #. vcf merge
  ######################################
  source('/home/starjjbang/an_proj/4_GWAS/code/GWAS_function.R')
  f_lst = list.files('/home/starjjbang/an_proj/4_GWAS/HN00159051/',pattern = paste0(cot_lab,"|",case_lab))
  ref_name = paste0(cot_lab,'_',case_lab,'.vcf')
  vcf_merge_fun(f_lst,vcf_dir,ref_name)
  ######################################
  #. make plink file format
  # - remove indel 
  ######################################
  #['#Family_ID', 'Individual_ID', 'Father_ID', 'Mother_ID', 'Sex','Phenotype','Genotype’]
  q = paste0(bcftools_dir," view --max-alleles 2 --exclude-types indels ",total_vcf_dir, ' > ', p_dir,'filter_total.vcf')
  system(q)
  runplink(paste0(" --vcf ",paste0(p_dir,'filter_total.vcf') ," --recode ped --out "),paste0(p_dir,"2_filter_indel"))
  runplink(paste0(" --vcf ",paste0(p_dir,'filter_total.vcf') ," --make-bed --out "),paste0(p_dir,"2_filter_indel"))
 
  ######################################
  #pedfile_fix --> insert phenotype
  ######################################
  runplink(paste0(" --file ",paste0(p_dir,"2_filter_indel"),
                  " --pheno /home/starjjbang/an_proj/4_GWAS/1_input_result/phe-1.txt --recode --out "),
           paste0(p_dir,"3_phe_update"))
  runplink(paste0(" --file ",paste0(p_dir,"2_filter_indel"),
                  " --pheno /home/starjjbang/an_proj/4_GWAS/1_input_result/phe-1.txt --make-bed --out "),
           paste0(p_dir,"3_phe_update"))
  
  tmp_fam = read.csv(paste0(p_dir,"3_phe_update.fam"),sep=' ',header = FALSE)
  phe_dat = read.csv('/home/starjjbang/an_proj/4_GWAS/1_input_result/phe-2.txt',sep='\t',header=TRUE)
  tmp_fam$V5 =  sapply(tmp_fam$V1, function(x) phe_dat[which( phe_dat$FID  == x),5])
  tmp_fam$V6[which(tmp_fam$V6 != 1)]=2
  write.table(tmp_fam,paste0(p_dir,"3_phe_update.fam"),sep=' ',quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  ######################################
  #Quality Control
  #Missingness per SNP: 0.1 --geno
  #Missingness per individual: 0.1 --mind
  #Minor allele frequency: 0.05 --maf
  #Hardy-Weinberg threshold: 0.0000001 --hwe
  ######################################
  runplink(paste0(" --bfile ",p_dir,"3_phe_update",
                  " --autosome --output-missing-genotype N --geno 0.8 --maf 0.05  --hwe 1e-6 --nonfounders --recode ped --out "),
           paste0(p_dir,"4_filter"))
  
  runplink(paste0(" --bfile ",p_dir,"3_phe_update",
                  " --autosome --output-missing-genotype N --geno 0.8 --maf 0.05  --hwe 1e-6 --nonfounders --make-bed --out "),
           paste0(p_dir,"4_filter"))
  #runplink(paste0(" --bfile ",re_dir,"2_filter_indel --filter-cases --make-bed -out "),"test01_case")
  ######################################
  #PCA
  ######################################
  runplink(paste0(" --bfile ", paste0(p_dir,"4_filter"),
                  " --allow-extra-chr --pca var-wts --recode --out "),
           paste0(p_dir,"pca"))
  vec = read.table(paste0(p_dir,"pca.eigenvec"),header=F,stringsAsFactors=F)[,c(1:6)]
  val = read.table(paste0(p_dir,"pca.eigenval"),header=F,stringsAsFactors=F)
  val = round(val / sum(val) * 100, 2)
  vec$V2 =  sapply( vec$V2,function(x) unlist(strsplit(x,'-'))[1])
  lab =  sapply( vec$V1,function(x) paste0(unlist(strsplit(x,'-'))[1:2],collapse = '-'))
  p = ggplot(data = vec, aes(x = -V3, y = -V4, fill = V2)) +
    geom_point(size=6, pch=21,stroke = 0) + 
    labs(x = paste0('PC1(',val$V1[1],'%)'),
         y = paste0('PC2(',val$V1[2],'%)')) + 
    theme(axis.title = element_text(size=7),
          axis.text = element_text(size=7),
          legend.title = element_text(size=7,face="bold"),
          panel.background = element_rect(fill = "white",color = "black",size=0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0, linetype = "dashed",color = "Dim gray"),
          panel.grid.minor = element_line(size = 0, linetype = "dashed",color = "gray")) + 
    scale_fill_manual(values = c("G1" = "#FDFEB9", "G2" = "#5C8984","G3" = "#F2D8D8", "G4" = "#374259")) +
    geom_text(label=lab,color="#1A2537",size = 2)
  
  
  png(filename = paste0(fig_dir,"1_PCA.png"), width = 2000, height = 1000,units = "px", bg = "white", res = 300)
  print(p)
  dev.off()
  ######################################
  #GWAS 
  ######################################
  runplink(paste0(" --file ",p_dir,'4_filter'," --assoc fisher --ci 0.95 --adjust --out "),
           paste0(p_dir,"5_fisher_",re_f_name))
  
  runplink(paste0(" --file ",p_dir,'4_filter'," --logistic --ci 0.95 --adjust --out "),
           paste0(p_dir,"5_logistic_",re_f_name))
  #qq_fun(paste0(p_dir,"5_fisher_total.assoc.fisher"), paste0(fig_dir,"2_fisher_MH.png" ))
  #qq_fun(paste0(p_dir,"5_logistic_total.assoc.logistic"), paste0(fig_dir,"3_logistic_MH.png" ))
}

m <- function(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", "gray60"), 
              chrlabs = NULL, suggestiveline = -log10(0.00001), 
              genomewideline = -log10(0.00000005), highlight = NULL, logp = TRUE, 
              annotatePval = NULL, annotateTop = TRUE,snp_font_size,snp_font_color, ...) {
  library(calibrate )
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  if (!is.null(x[[snp]])) 
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                   pos = NA, index = NA, SNP = x[[snp]], stringsAsFactors = FALSE)
  else d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                      pos = NA, index = NA)
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, 
                                                             d$CHR, length))
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    xlabel = paste("Chromosome", unique(d$CHR), "position")
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + max(d[d$index == (i - 1), 
                                    "BP"])
        d[d$index == i, "BP"] = d[d$index == i, "BP"] - 
          min(d[d$index == i, "BP"]) + 1
        d[d$index == i, "pos"] = d[d$index == i, "BP"] + 
          lastbase
      }
    }
    ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp[!is.na(d$logp)]) )), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep_len(col, max(d$index))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      points(d[d$index == i, "pos"], d[d$index == i, "logp"], 
             col = col[icol], pch = 20, ...)
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = "#0F7878",lwd=2)
  if (genomewideline) 
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = "green3", pch = 20, 
                             ...))
  }
  if (!is.null(annotatePval)) {
    if (logp) {
      topHits = subset(d, P <= annotatePval)
    }
    else topHits = subset(d, P >= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      if (logp) {
        with(subset(d, P <= annotatePval), textxy(pos, 
                                                  -log10(P), offset = 0.625, labs = topHits$SNP, 
                                                  cex = 0.45), ...)
      }
      else with(subset(d, P >= annotatePval), textxy(pos, 
                                                     P, offset = 0.625, labs = topHits$SNP, cex =0.45), 
                ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      if (logp) {
        textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
               labs = topSNPs$SNP, cex = snp_font_size,col = snp_font_color, ...)
      }
      else textxy(topSNPs$pos, topSNPs$P, offset = 0.625, 
                  labs = topSNPs$SNP, cex = snp_font_size,col = snp_font_color, ...)
    }
  }
  par(xpd = FALSE)
}



# library(qqman)
# re_dat = read_table('/home/starjjbang/an_proj/4_GWAS/result/1_G1_G2/plink_result/5_fisher_total.assoc.fisher')
# manhattan(x = re_dat, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
#           col = c("#A9B8CE",  "#091251"), suggestiveline = -log10(1e-03), logp = TRUE,annotatePval = 0.001)
# 
# #https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
# library(tidyverse)
# library(ggtext)
# library(normentR)
# sig_data <- re_dat %>% subset(P < 0.05)
# notsig_data <- re_dat %>% subset(P >= 0.05) %>% group_by(CHR) %>%  sample_frac(0.1)
# gwas_data <- bind_rows(sig_data, notsig_data)
# data_cum <- gwas_data %>% group_by(CHR) %>% summarise(max_bp = max(BP)) %>% 
#   mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
#   select(CHR, bp_add)
# gwas_data <- gwas_data %>% inner_join(data_cum, by = "CHR") %>% mutate(bp_cum = BP + bp_add)
# axis_set <- gwas_data %>% group_by(CHR) %>% summarize(center = mean(bp_cum))
# ylim <- gwas_data %>% filter(P == min(P)) %>% mutate(ylim = abs(floor(log10(P))) + 2) %>% pull(ylim)
# sig <- 1e-04
# manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P),  color = as_factor(CHR), size = -log10(P))) +
#   geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
#   geom_point(alpha = 0.75) +
#   scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
#   scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$CHR)))) +
#   scale_size_continuous(range = c(0.5,3)) +
#   labs(x = NULL, y = "-log<sub>10</sub>(p)") + 
#   theme_minimal() +
#   theme( 
#     legend.position = "none",
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     axis.title.y = element_markdown(),
#     axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
#   )
# manhplot
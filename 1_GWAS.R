#sample information 
#G1 : the first medical examination
#G2 : improve
#G3 : exacerbation
#G4 : Stop progression

#https://varamos.tistory.com/11
#https://plink.readthedocs.io/en/latest/analysis_bash/

#phenotype data make
phe_dat = read.csv('/home/starjjbang/an_proj/4_GWAS/1_input_result/phe.txt',sep='\t',header=TRUE)
phe_tmp = data.frame("FID" = paste0(phe_dat$Data.identifire ,'-2'), "IID" = paste0(phe_dat$Data.identifire ,'-2'),"CLS" = rep(c(1,2,3,4),each=5))
write.table(phe_tmp, '/home/starjjbang/an_proj/4_GWAS/1_input_result/phe-1.txt',sep='\t',row.names=FALSE,quote = FALSE)
data_dir = '/home/starjjbang/an_proj/3_hurb/1_clinical_data/1_data/'
phe_dat2 = read.csv(paste0(data_dir,'6_0331_merge_dat.txt'),sep='\t',header=TRUE)
head(phe_dat2)
phe_dat2 = phe_dat2 %>% select(-row.names)
head(phe_dat2)
add_data = cbind(c('yj191','yj192','yj193','yj194','yj195'),
      c(55,59,57,67,54),
      c(0,0,0,0,1),
      c(1,1,2,2,1),
      matrix(rep(0,5*28), nrow = 5, ncol = 28))
colnames(add_data) = colnames(phe_dat2)
phe_dat2 = rbind(phe_dat2, add_data)
write.table(phe_dat2,paste0(data_dir,'7_0807_add_merge_dat.txt'),sep='\t',quote = FALSE,row.names = FALSE)  

phe_dat = read.csv('/home/starjjbang/an_proj/4_GWAS/1_input_result/phe-1.txt',sep='\t')
phe_dat1 = read.csv('/home/starjjbang/an_proj/4_GWAS/1_input_result/phe.txt',sep='\t',header=TRUE)
phe_total = read.csv(paste0(data_dir,'7_0807_add_merge_dat.txt'),sep='\t',header=TRUE)

phe_dat$key = sapply(phe_dat$FID, function(x) paste0(unlist(strsplit(x,'-'))[1:2],collapse = '-'))
m_dat = merge(phe_dat, phe_dat1, by.x = "key",by.y = "Data.identifire",all.x=TRUE)
m_dat1 = merge(m_dat, phe_total, by.x = "Clinical.Patient.Number",by.y = "P_id",all.x=TRUE)
m_dat1  = m_dat1 %>% select( FID,IID,CLS,age, sex)
m_dat1$sex[which(m_dat1$sex=="1")] = 1
m_dat1$sex[which(m_dat1$sex=="0")] = 2
m_dat1$sex[is.na(m_dat1$sex)] = 0
write.table(m_dat1,'/home/starjjbang/an_proj/4_GWAS/1_input_result/phe-2.txt',sep='\t',quote = FALSE,row.names = FALSE) 
##################################################################################################################

rm(list=ls())
source('/home/starjjbang/an_proj/4_GWAS/code/GWAS_function.R')

GWAS_run <- function(exp_no,lab_1,lab_2){
  GWAS_pair(exp_no,lab_1,lab_2)
  fig_lab = paste0(exp_no,'_',lab_1,'_',lab_2)
  qq_fun(paste0(pt_dir, fig_lab,'/plink_result/',"5_fisher_total.assoc.fisher"), 
         paste0(pt_dir, fig_lab,"/fig/2_",fig_lab,"_fisher_MH.png" ))  
}
#GWAS_run("1","G1","G2")
#GWAS_run("3","G1","G3")
#GWAS_run("3","G1","G4")
#GWAS_run("4","G2","G3")
#GWAS_run("4","G2","G4")
#GWAS_run("4","G3","G4")

exp_no = "4"
lab_1 = "G2"
lab_2 = "G3"
#ref_na = "test01rest.assoc"
ref_na = "5_fisher_total.assoc.fisher"

fig_lab = paste0(exp_no,'_',lab_1,'_',lab_2)
re_dir = paste0('/home/starjjbang/an_proj/4_GWAS/result/',fig_lab,'/')
in_dat_f = paste0(re_dir,'plink_result/',ref_na)
re_dat = data.table::fread(in_dat_f)
m(x = re_dat, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
  col = c("#BBC814",  "#0A4A61"), 
  suggestiveline = -log10(5e-05), 
  logp = TRUE,annotatePval = 5e-05,
  cex.lab=1.5,snp_font_size=0.8,snp_font_color  = "#FEA609")
p <- recordPlot()
g <- grid.grabExpr(grid.echo(p))
ggsave(paste0(re_dir,'fig/MH.png'), g, bg = "transparent", width = 10, height = 8,dpi=300)



in_dat_f = paste0(pt_dir, fig_lab,'/plink_result/',"5_fisher_total.assoc.fisher")
re_dat = data.table::fread(in_dat_f)

m(x = re_dat, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
  col = c("#BBC814",  "#0A4A61"), 
  suggestiveline = -log10(5e-04),
  logp = TRUE,
  annotatePval = 5e-04,
  cex.lab=1.5,
  snp_font_size=0.8,
  snp_font_color  = "#FEA609",ylim = c(0,5))

p <- recordPlot()
g <- grid.grabExpr(grid.echo(p))




GWAS_pair("2","G1","G3")
GWAS_pair("3","G1","G4")



p_dir = '/home/starjjbang/an_proj/4_GWAS/result/1_G1_G2/plink_result/'
qq_fun(paste0(p_dir,"5_fisher_total.assoc.fisher"), paste0(fig_dir,"2_fisher_MH.png" ))




#plink --noweb --bfile [mydata] --extract [Gene_SNP_list] --recodeHV --out [Gene_haploview]
#runplink(paste0(" --noweb --bfile ",p_dir,'4_phe_update'," --extract  --out "),paste0(p_dir,"5_logistic_",re_f_name))

##################################################################################################################
#linear
##################################################################################################################
rm(list=ls())
source('/home/starjjbang/an_proj/4_GWAS/code/GWAS_function.R')
re_dir = paste0('/home/starjjbang/an_proj/4_GWAS/result/linear/')
in_f_dir = "/home/starjjbang/an_proj/4_GWAS/HN00159051/"
ref_name ='total.vcf'

GWAS_total <-function(re_dir,in_f_dir,ref_name){
  p_dir = paste0(re_dir, 'plink_result/')
  vcf_dir = paste0(re_dir,'vcf/')
  rds_dir = paste0(re_dir,'rds/')
  system(paste0('mkdir ',re_dir))
  system(paste0('mkdir ',p_dir))
  system(paste0('mkdir ',vcf_dir))
  system(paste0('mkdir ',rds_dir))
  ######################################
  #. vcf merge
  ######################################
  f_lst = list.files(in_f_dir)
  vcf_merge_fun(f_lst,vcf_dir,ref_name)
  total_vcf_dir = paste0(vcf_dir,ref_name) 
  ######################################
  #. make plink file format
  # - remove indel 
  ######################################
  #['#Family_ID', 'Individual_ID', 'Father_ID', 'Mother_ID', 'Sex','Phenotype','Genotype’]
  q = paste0(bcftools_dir," view --max-alleles 2 --exclude-types indels ",total_vcf_dir, ' > ', p_dir,'1_filter_total.vcf')
  system(q)
  runplink(paste0(" --vcf ",paste0(p_dir,'1_filter_total.vcf') ," --recode ped --out "),paste0(p_dir,"2_filter_indel"))
  runplink(paste0(" --vcf ",paste0(p_dir,'1_filter_total.vcf') ," --make-bed --out "),paste0(p_dir,"2_filter_indel"))
 
  ######################################
  #pedfile_fix --> insert phenotype
  ######################################
  runplink(paste0(" --file ",paste0(p_dir,"2_filter_indel"),
                  " --pheno /home/starjjbang/an_proj/4_GWAS/1_input_result/phe-2.txt --recode --out "),
           paste0(p_dir,"3_phe_update_",ref_name))
  runplink(paste0(" --file ",paste0(p_dir,"2_filter_indel"),
                  " --pheno /home/starjjbang/an_proj/4_GWAS/1_input_result/phe-2.txt --make-bed --out "),
           paste0(p_dir,"3_phe_update_",ref_name))
  
  tmp_fam = read.csv(paste0(p_dir,"3_phe_update_total.vcf.fam"),sep=' ',header = FALSE)
  phe_dat = read.csv('/home/starjjbang/an_proj/4_GWAS/1_input_result/phe-2.txt',sep='\t',header=TRUE)
  tmp_fam$V5 =  sapply(tmp_fam$V1, function(x) phe_dat[which( phe_dat$FID  == x),5])
  write.table(tmp_fam,paste0(p_dir,"3_phe_update_total.vcf.fam"),sep=' ',quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  ######################################
  #Quality Control
  #Missingness per SNP: 0.1 --geno
  #Missingness per individual: 0.1 --mind
  #Minor allele frequency: 0.05 --maf
  #Hardy-Weinberg threshold: 0.0000001 --hwe
  ######################################
  runplink(paste0(" --bfile ",p_dir,"3_phe_update_total.vcf",
                  " --autosome --output-missing-genotype N --geno 0.6 --maf 0.05  --hwe 1e-6 --nonfounders --recode ped --out "),
           paste0(p_dir,"4_filter_",ref_name))
  
  runplink(paste0(" --bfile ",p_dir,"3_phe_update_total.vcf",
                  " --autosome --output-missing-genotype N --geno 0.6 --maf 0.05  --hwe 1e-6 --nonfounders --make-bed --out "),
           paste0(p_dir,"4_filter_",ref_name))
  

  ######################################
  #PCA
  ######################################
  runplink(paste0(" --bfile ", paste0(p_dir,"4_filter_",ref_name),
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
  
}





#fig <- plot_ly(vec, x = ~V3, y = ~V4, z = ~V5, color = ~V2, colors = c('#FDFEB9','#5C8984','#F2D8D8',"#374259") )%>%
#  add_markers(marker=list(sizeref=0.02))
#fig <- fig %>% layout(title = "test",scene = list(bgcolor = "white"))
#fig
######################################
#GWAS 
######################################
runplink(paste0(" --file ",p_dir,"4_filter_",ref_name," --linear --ci 0.95 --adjust --out "),
         paste0(p_dir,"5_linear_",ref_name))
library(qqman)
re_dat = read_table(paste0(p_dir,'5_linear_total.vcf.assoc.linear'))
f_re_dat = re_dat[which(re_dat$P < 0.05),]
dim(f_re_dat)


p_dir = paste0(re_dir, 'plink_result/')
runplink_q(paste0(" --ped ",p_dir,'4_filter'," --assoc --adjust --out "),
           paste0(p_dir,"test01",re_f_name))
in_dat_f = paste0(pt_dir, fig_lab,'/plink_result/',"test_logistic_rest.assoc")#chi-square test,



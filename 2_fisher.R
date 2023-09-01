#http://zhanxw.github.io/rvtests/
rm(list=ls())
source('/home/starjjbang/an_proj/4_GWAS/code/GWAS_function.R')
f_lst = list.files('/home/starjjbang/an_proj/4_GWAS/HN00159051/')
re_vcf_dir = '/home/starjjbang/an_proj/4_GWAS/result/5_freq/vcf/'
bcftools_dir = '/home/starjjbang/an_proj/4_GWAS/bcftools/bcftools'
plink_dir = '/home/starjjbang/an_proj/4_GWAS/plink'
re_dir = '/home/starjjbang/an_proj/4_GWAS/result/5_freq/2_rare_test/'
bcftools_dir = '/home/starjjbang/an_proj/4_GWAS/bcftools/bcftools'


#total VCF merge
cot_lab = "G1"
case_lab = "G2"
re_f_dir = '/home/starjjbang/an_proj/4_GWAS/result/5_rare_test/'
vcf_dir = '/home/starjjbang/an_proj/4_GWAS/result/5_freq/vcf/'
filter_vcf_dir = paste0(re_dir,'filter_total.vcf')
total_vcf_dir = paste0(re_dir,cot_lab,'-',case_lab,'.vcf')

#VCF merge
f_lst = list.files(vcf_dir,pattern = paste0(cot_lab,"|",case_lab))
lab_1 = paste(f_lst[grep("*.vcf.gz$",f_lst)],collapse = " ")
re_dir = paste0(re_f_dir,cot_lab,'_',case_lab,'/')
system(paste0('mkdir ',re_dir))
q3 = paste0('vcf-merge ',lab_1,' > ',re_dir,cot_lab,'-',case_lab,'.vcf') 
print(q3)

#remove indel 
q = paste0(bcftools_dir," view --max-alleles 2 --exclude-types indels ",total_vcf_dir, 
           ' > ',re_dir,'filter_total.vcf')
system(q)

#make bed, fam, bim
q1 = paste0(plink_dir, ' --vcf ',filter_vcf_dir,' --out ',re_dir,'total --make-bed' )
system(q1)

# make ped, map
#['#Family_ID', 'Individual_ID', 'Father_ID', 'Mother_ID', 'Sex','Phenotype','Genotype’]
q2 = paste0(plink_dir, ' --vcf ',filter_vcf_dir,' --out ',re_dir,'total --recode ped' )
system(q2) 

runplink(paste0(" --file ",paste0(re_dir,"total"),
                " --pheno /home/starjjbang/an_proj/4_GWAS/1_input_result/phe-2.txt",
                " --recode ped --out "),
         paste0(re_dir,"3_phe_update"))


#Rare GWAS test
rvt_dir = '/home/starjjbang/an_proj/4_GWAS/rvtests/executable/rvtest'
input.vcf = paste0(re_dir, 'filter_total.vcf')
ped_dir = paste0(re_dir, '3_phe_update.ped')
q4 = paste0(rvt_dir,' --inVcf ',input.vcf,' --pheno ', ped_dir, 
            ' --out ',re_dir,'GWAS_test --single wald,score,exact' )
system(q4)

dat = read.csv(paste0(re_dir, 'GWAS_test.FisherExact.assoc'),sep='\t')
head(dat)
unique(dat$Pvalue)
which(dat$Pvalue < 0.05)


in_lab = "G4"
f_lst = list.files(re_vcf_dir,pattern = in_lab)
lab_1 = paste(f_lst[grep("*.vcf.gz$",f_lst)],collapse = " ")
vcf_dir = '/home/starjjbang/an_proj/4_GWAS/result/5_freq/vcf/'
q3 = paste0('vcf-merge ',vcf_dir,lab_1,' > ',
            '/home/starjjbang/an_proj/4_GWAS/result/5_freq/1_merge_VCF/',in_lab,'.vcf') 

q = paste0(bcftools_dir," view --max-alleles 2 --exclude-types indels ",
           paste0(p_dir,in_lab,'.vcf'), ' > ', p_dir,'1_filter_',in_lab,'.vcf')
system(q)
runplink(paste0(" --vcf ",p_dir,'1_filter_',in_lab,'.vcf',
                " --make-bed --out "), paste0(p_dir,'2_filter_',in_lab))
runplink(paste0(" --bfile ",p_dir,'2_filter_',in_lab," --freq --out "),
         paste0(p_dir,'3_filter_',in_lab,'.vcf'))


p_dir = '/home/starjjbang/an_proj/4_GWAS/result/5_freq/1_merge_VCF/'

g1_dat = read_table(paste0(p_dir,'3_filter_G1.vcf.frqx'))
g2_dat = read_table(paste0(p_dir,'3_filter_G2.vcf.frqx'))

for_tmp = c()
for(i in 1:nrow(g2_dat)){
  if(length(which(g1_dat$SNP == g2_dat$SNP[i])) > 0){
    for_tmp = rbind(for_tmp, cal_fisher(g2_dat$SNP[i],g2_dat,g1_dat)) 
  }
  print(i)
}

for_tmp[which(as.numeric(for_tmp[,2]) < 0.05),]

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


/home/starjjbang/an_proj/4_GWAS/rvtests/executable/rvtest --inVcf input.vcf --pheno phenotype.ped --out output --single wald,score




cal_fisher <-function(v_id,case_dat, control_dat){
  tmp = rbind(case_dat[which(case_dat$SNP == v_id),5:7],control_dat[which(control_dat$SNP == v_id),5:7])
  t_dat = rbind(c( as.numeric((tmp[1,1]*2) + tmp[1,2]),as.numeric((tmp[1,3]*2) + tmp[1,2]  ) ),
                c( as.numeric((tmp[2,1]*2) + tmp[2,2]),as.numeric((tmp[2,3]*2) + tmp[2,2]  ) ) )
  t_re = fisher.test(t_dat)
  re_dat = c(v_id,t_re$p.value, t_re$estimate)
  return(re_dat)
}



# q4 = paste0(bcftools_dir," query -f '%CHROM\t%POS0\t%END\n' ",
#             "/home/starjjbang/an_proj/4_GWAS/result/5_freq/vcf/G1.vcf > G1.bed")
# bedtools fisher -a a.bed -b b.bed
# sort -k1,1 hg38.chrom.sizes > sorted_hg38.chrom.sizes
# sort -k1,1 -k2,2n G1.bed > G1_sorted.bed
# sort -k1,1 -k2,2n G2.bed > G2_sorted.bed
# bedtools fisher -a G1_sorted.bed -b G2_sorted.bed -g sorted_hg38.chrom.sizes


vcf_merge_fun <-function(f_lst,vcf_dir){
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
  #q3 = paste0('vcf-merge ',vcf_dir,'*.vcf.gz > ',vcf_dir,ref_name) 
  #system(q3)  
  #q4 = paste0('bgzip ', vcf_dir,cot_lab,'_',case_lab,'.vcf.gz')
  #system(q4)
}

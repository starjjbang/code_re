source('/home/starjjbang/an_proj/4_GWAS/code/GWAS_function.R')
in_lab = "G1"
vcf_dir = '/home/starjjbang/an_proj/4_GWAS/result/5_freq/vcf/'
re_dir = '/home/starjjbang/an_proj/4_GWAS/result/7_MAF_test/'
bcftools_dir = '/home/starjjbang/an_proj/4_GWAS/bcftools/bcftools'
p_dir = paste0(re_dir,'1_merge_VCF/')
  
f_lst = list.files(vcf_dir,pattern = in_lab)
lab_1 = paste(f_lst[grep("*.vcf.gz$",f_lst)],collapse = " ")
q3 = paste0('vcf-merge ',vcf_dir,lab_1,' > ',
            re_dir,'1_merge_VCF/',in_lab,'.vcf') 
print(q3)
q = paste0(bcftools_dir," view --max-alleles 2 --exclude-types indels ",
           paste0(p_dir,in_lab,'.vcf'), ' > ', p_dir,'1_filter_',in_lab,'.vcf')
system(q)
runplink(paste0(" --vcf ",p_dir,'1_filter_',in_lab,'.vcf',
                " --make-bed --out "), paste0(p_dir,'2_filter_',in_lab))
runplink(paste0(" --bfile ",p_dir,'2_filter_',in_lab," --freq --out "),
         paste0(p_dir,'3_filter_',in_lab,'.vcf'))
#https://blog.naver.com/sw4r/221541700667
#MAF = 0.5 : 4 명 8 allele 중 T가 4명이 다 있었음
#A1에 대한 MAF https://mopipe.tistory.com/13


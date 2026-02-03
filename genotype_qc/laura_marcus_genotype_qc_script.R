## Scripts for Tobi Genotyping QC ##
working.dir <- "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/data/genotypes/dna_genotypes/2023_Nov2/"
# setwd("/lustre/scratch126/humgen/projects/sc-eqtl-ibd/data/genotypes/genotypes_qc_imputed_june_2023/")
setwd(working.dir)

# 23rd Nov 2023
# Merge plates across the runs

# NOTE Have to first edit VCF to remove space from sample name: Pos Con -> PosCon2

# Then bcf tools filter out the positive controls
system()

system("/software/team152/plink_linux_x86_64_20181202/plink --vcf UKX68_All.vcf.gz --allow-extra-chr --output-chr MT --make-bed --out merged_plates")

#####################
# 1.- REMOVE INDELS #
#####################

#Read in .bim file to get list of SNPs
bim<-read.table("merged_plates.bim")

#Get list of indels
indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]

#Get list of SNPs on non 1-22 & sex chr's (usually mitochondrial)
otherchr <- bim[ !(bim$V1 %in% c(1:22,"X","Y")),]

#Get joint list of SNPs to remove
all_remove<-rbind(indels,otherchr)
#Remove duplicate SNPs in list
all_remove<-all_remove[!duplicated(all_remove),]

#Generate list of SNPs to remove
write.table(all_remove[,"V2"],"list_indel_var_exclude", col.names=F,row.names=F,quote=F,sep="\t")

#Remove SNPs in this list from merged file:
system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates --allow-extra-chr --output-chr MT --exclude list_indel_var_exclude --make-bed --out merged_plates_ATCG")


#########################
# 4.- ALIGN TO + STRAND #
#########################
system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG --allow-extra-chr --output-chr MT --allow-no-sex --recode vcf --out merged_plates_ATCG")

system("export BCFTOOLS_PLUGINS=/software/hgi/installs/anaconda3/envs/hgi_base/libexec/bcftools/; /software/team152/bcftools-1.9/./bcftools +fixref merged_plates_ATCG.vcf -Oz -o merged_plates_hg19_posstrandaligned.vcf.gz -- -f /lustre/scratch125/humgen/resources/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta  -m top")

# Remove mismatch alleles, they cause imputation server to crash:
system("/software/team152/bcftools-1.9/./bcftools norm --check-ref x -f /lustre/scratch125/humgen/resources/ref/Homo_sapiens/1000Genomes/human_g1k_v37.fasta merged_plates_hg19_posstrandaligned.vcf.gz -Oz -o merged_plates_hg19_posstrandaligned_2.vcf.gz")

# After alignment, format vcf to bed file
system("/software/team152/plink_linux_x86_64_20181202/./plink --vcf merged_plates_hg19_posstrandaligned_2.vcf.gz --double-id --make-bed --out merged_plates_hg19_noind_posstr")
system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr --missing --out merged_plates_hg19_noind_posstr")

################################################################################
# 5.- UPDATE NAME VARIANTS TO CHR:POSITION_REF_ALT; REMOVE DUPLICATED VARIANTS #
################################################################################

# Update the IDs to chr:pos:ref:alt
system("zcat merged_plates_hg19_posstrandaligned_2.vcf.gz | cut -f '1-5' | awk '{print $3,$1\":\"$2\"_\"$4\"_\"$5}' > list_variants_merged_plates_hg19_posstrandaligned")

###### R
bim<-read.table("merged_plates_hg19_noind_posstr.bim",sep="\t",head=F)
ids<-read.table("list_variants_merged_plates_hg19_posstrandaligned",sep=" ",head=F,skip=35)

varmiss<-read.table("merged_plates_hg19_noind_posstr.lmiss",head=T)

colnames(ids)[2]<-"ids"


#the major allele is set to A2 by default by Plink, keep ids with real ref/alt as in vcf using the ids file
bim.1<-cbind(bim,ids[,"ids",drop=F])

bim.1<-merge(bim.1,varmiss[,c("SNP","F_MISS")],by.x="V2",by.y="SNP",sort=F)

bim.1$ids<-as.character(bim.1$ids)

# identify duplicated variants (same chr position ref and alt)
dups<-bim.1[which(duplicated(bim.1$ids)),"ids"]

### Not done because no duplicates ###

# for (i in 1:length(dups)){
#
#   tmp<-bim.1[which(bim.1$ids %in% dups[i]),]
#
#   keep<-tmp[which(tmp$F_MISS==min(tmp$F_MISS)),]
#
#   if (nrow(keep)>1){
#     keep<-keep[1,]
#   }
#
#   exclude<-tmp[which(!tmp$V2 %in% keep$V2),]
#
#   bim.1$ids[which(bim.1$V2 %in% exclude$V2)]<-paste(bim.1$ids[which(bim.1$V2 %in% exclude$V2)],"_rm",sep="")
#
# }

write.table(duplicated_variants,"list_duplicated_var_exclude",col.names=F,row.names=F,quote=F,sep="\t")
write.table(bim.1[,c(2,7,3:6)],"merged_plates_hg19_noind_posstr_edited.bim",col.names=F,row.names=F,quote=F,sep="\t")

system("/software/team152/plink_linux_x86_64_20181202/./plink --bed merged_plates_hg19_noind_posstr.bed --fam merged_plates_hg19_noind_posstr.fam --bim merged_plates_hg19_noind_posstr_edited.bim --split-x 'b37' 'no-fail' --make-bed --out merged_plates_hg19_noind_posstr_nodup")

##############################################
# 6.- COMPARE ALLELE FREQUENCIES WITH 1000GP #
##############################################

system("cat list_variants_merged_plates_hg19_posstrandaligned | cut -f 2 > list_variants_merged_plates_hg19_noind_posstr_nodup")

for (i in 1:23) {
  system(sprintf("/software/team152/plink_linux_x86_64_20181202/./plink --bfile /lustre/scratch123/hgi/projects/ibdgwas_ukibdgc/resources/1000GP/1000GP_EUR_chr%s_b37 --extract list_variants_merged_plates_hg19_noind_posstr_nodup --make-bed --out ./1KG_data/1000GP_EUR_chr%s_b37_merged_plates_variants",i,i))
}

dat<-matrix(ncol=3,nrow=22)
dat<-as.data.frame(dat)
for (i in 1:nrow(dat)){
  dat[i,1]<-paste(working.dir,"/1KG_data/1000GP_EUR_chr",i+1,"_b37_merged_plates_variants.bed",sep="")
  dat[i,2]<-paste(working.dir,"/1KG_data/1000GP_EUR_chr",i+1,"_b37_merged_plates_variants.bim",sep="")
  dat[i,3]<-paste(working.dir,"/1KG_data/1000GP_EUR_chr",i+1,"_b37_merged_plates_variants.fam",sep="")
}

write.table(dat,"./1KG_data/list_1000GP_files_merged_plates_variants_tomerge.txt",col.names=F,row.names=F,quote=F,sep="\t")

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile ./1KG_data/1000GP_EUR_chr1_b37_merged_plates_variants --merge-list ./1KG_data/list_1000GP_files_merged_plates_variants_tomerge.txt --allow-no-sex --make-bed --out ./1000GP_EUR_b37_merged_plates_variants")

# 1000GP
system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile 1000GP_EUR_b37_merged_plates_variants --freq --out 1000GP_EUR_b37_merged_plates_variants")

# Genotypes
system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup --freq --out merged_plates_hg19_noind_posstr_nodup")

#### R
gp<-read.table("1000GP_EUR_b37_merged_plates_variants.frq",head=T)
g1<-read.table("merged_plates_hg19_noind_posstr_nodup.frq",head=T)

colnames(gp)[3:6]<-paste(colnames(gp)[3:6],"_gp",sep="")
colnames(g1)[3:6]<-paste(colnames(g1)[3:6],"_g1",sep="")

all<-merge(g1,gp[,2:6],by="SNP",all=T)

check<-all[which(all$A1_g1!=all$A1_gp),]

# keep only A/T C/G
check<-check[which( (check$A1_g1=="G" & check$A2_g1=="C") | (check$A1_g1=="C" & check$A2_g1=="G") | (check$A1_g1=="A" & check$A2_g1=="T") | (check$A1_g1=="T" & check$A2_g1=="A")),]

# LIST OF VARIANTS TO REMOVE, WE CANNOT REALLY BE SURE WHETHER THERE IS STRAND ISSUE OR NOT
remove<-check[which(check$MAF_g1>=0.45),]

# LIST OF VARIANTS TO FLIP:
flip<-check[which(check$MAF_g1<0.45),]

flip<-flip[order(flip$MAF_g1,decreasing=T),]

remove_2<-flip[which(flip$MAF_g1>0.2 & flip$MAF_gp<0.1),]
remove_3<-flip[which(flip$MAF_g1<0.1 & flip$MAF_gp>0.2),]

remove<-rbind(remove,remove_2,remove_3)

write.table(remove[,"SNP"],"list_variants_to_remove_AT_CG",col.names=F,row.names=F,quote=F,sep="\t")

flip<-flip[which(!flip$SNP %in% remove$SNP),]

write.table(flip[,"SNP"],"list_variants_to_flip_AT_CG",col.names=F,row.names=F,quote=F,sep="\t")

###############
# 6.5 REMOVE VARIANTS WE CANNOT BE SURE ARE IN THE RIGHT STRAND, AND FLIP THE OTHERS
################

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup --exclude list_variants_to_remove_AT_CG --flip list_variants_to_flip_AT_CG --make-bed --out merged_plates_hg19_noind_posstr_nodup_flip")

###################
# 8.- REMOVE ChrY #
###################

system("awk '{if($1==24)print $2}' < merged_plates_hg19_noind_posstr_nodup_flip.bim > list_chry_variants_toremove")
system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_hg19_noind_posstr_nodup_flip --allow-no-sex --exclude list_chry_variants_toremove --make-bed --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry")

##########
# 9 - QC #
##########

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry --allow-no-sex --missing --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry")

### R

sample_miss<-read.table("merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.imiss",head=T)
var_miss<-read.table("merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.lmiss",head=T)

####################################
# 9.1 - REMOVE SAMPLES CallPP <0.80%

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry --mind 0.20  --allow-no-sex --make-bed --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8")

####################################
# 9.2- REMOVE VARIANTS CallPP <0.80%

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8 --geno 0.20 --allow-no-sex --make-bed --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8")

####################################
# 9.3 - REMOVE SAMPLES CallPP <0.95%

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8 --mind 0.05 --allow-no-sex --make-bed --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95")

#####################################
# 9.4 - REMOVE VARIANTS CallPP <0.95%

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95 --geno 0.05 --allow-no-sex --make-bed --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95")

#######################################################
# 9.5 - REMOVE VARIANTS FREQ <0.05 AND CallPP <0.98%

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 --allow-no-sex --missing --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95")

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 --allow-no-sex --freq --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95")

###### R

sample_miss<-read.table("merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.imiss",head=T)
var_miss<-read.table("merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.lmiss",head=T)
frq<-read.table("merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.frq",head=T)

var<-merge(frq[,c(2:6)],var_miss,by="SNP")
var.1<-var[which(var$MAF<0.05 & var$F_MISS>0.02),]

write.table(var.1[,"SNP",drop=F],"list_monomorphic_vcr098_var_exclude",col.names=F,row.names=F,quote=F,sep="\t")

##############

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 --exclude list_monomorphic_vcr098_var_exclude --allow-no-sex --make-bed --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05")


###############################################
# 12.- REMOVE INTRA-COHORT DUPLICATED SAMPLES #
###############################################

system("/software/team152/./king -b merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05.bed --kinship --prefix kinship")
system("/software/team152/./king -b merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05.bed --related --prefix test")
kin<-read.table("test.kin0",head=T)

#Extract the duplicates
kin_duplicates = kin[which(kin$Kinship > 0.354),]
# NOTE: Decided to drop 

#Remove the duplicate with lower callrate
system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05 --missing")
missing = read.table("plink.imiss", header = T)

#Loop through each pair and extract the sample with low call rate to exclude:
samples_to_exclude = c()
for (i in 1:nrow(kin_duplicates)) {
  sample_A = kin_duplicates[i,]$FID1
  sample_B = kin_duplicates[i,]$FID2
  #Find these samples in the missing table
  sample_A_missing = (missing[(which(sample_A == missing$FID)),]$F_MISS)
  sample_B_missing = (missing[(which(sample_B == missing$FID)),]$F_MISS)

  #Find sample with lower call rate to exclude
  sample_remove_index = which.min(c(sample_A_missing,sample_B_missing))
  if (sample_remove_index == 1) {
    samples_to_exclude[i] = sample_A
  }
  else{samples_to_exclude[i] = sample_B}
}

write.table(cbind(samples_to_exclude,samples_to_exclude),file = "samples_to_exclude_duplicate", quote = F, col.names = F, row.names = F)

#Exclude with PLINK
system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05 --remove samples_to_exclude_duplicate --make-bed --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05_nodupsample")


#### PCA Check ####
system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05_nodupsample --indep-pairwise 50 10 0.1 --out pruned_study")
system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05_nodupsample --extract pruned_study.prune.in --make-bed --pca --out pca_pruned")
#Read in eigenvecs
pca_eigenvecs = read.table("pca_pruned.eigenvec", header = F)

# Why are these samples outliers? 
#3.9815767767e+012.CEL_3.9815767767e+012.CEL 1
#3.9815768208e+012.CEL_3.9815768208e+012.CEL 2
#3.9815768007e+012.CEL_3.9815768007e+012.CEL 3
#3.9815767947e+012.CEL_3.9815767947e+012.CEL 4
# 1 and 2 are non-EUR (see below analyses), 3 and 4 are parent-child (see King relatedness)

#### Project into 1000G PC space and predict ancestry #####
system("/software/team152/./king -b king/KGref_varid2.bed,merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05.bed --pca --projection --rgplot --prefix 1000g_projected")
projected_data <- read.table(file="1000g_projectedpc.txt", header=T)
# ggplot(projected_data, aes(x=PC1, y=PC2, colour=as.factor(AFF), alpha=AFF/2)) +
#   geom_point() + theme_classic() + scale_colour_manual(values = c("grey", "red")) +
#   theme(legend.position = "none")

# Drop based on predicted ancestry
inferred_ancestry <- read.table(file='1000g_projected_InferredAncestry.txt', header=T)
inferred_non_eurs <- inferred_ancestry[inferred_ancestry$Ancestry != 'EUR',]$FID

write.table(cbind(inferred_non_eurs,inferred_non_eurs),"samples_non_eur_remove", quote = F, col.names = F, row.names = F)
#Exclude with PLINK
system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05_nodupsample --remove samples_non_eur_remove --make-bed --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05_nodupsample_eur")

#######################
# 15 - HETEROZYGOSITY #
#######################

system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05_nodupsample_eur --hwe 1e-6 --make-bed --out merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05_nodupsample_eur_hwe")


# #Filter heterozygosity
# filter_heterozygosity = function(){
#
#   system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.05_nodup_sample_eur --allow-no-sex --chr 1-22 --het --out autosomal_het_nodup")
#   het<-read.table("autosomal_het_nodup.het",head=T)
#   het$het<-(het$N.NM.-het$O.HOM)/het$N.NM.
#
#   #No samples were removed
#
# }

######################
# 16.- UPDATE TO B38 #
######################

# NOTE: data in build 37, needs to be liftover first:
system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05_nodupsample_eur_hwe --recode tab --out merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.05_nodup_sample_eur")

#### R convert map into bed:

map<-read.table("merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.05_nodup_sample_eur.map",head=F)

map$chr<-paste("chr",map$V1,sep="")
map$chr[which(map$V1=="23")]<-"chrX"
map$chr[which(map$V1=="25")]<-"chrX"
map$chr[which(map$V1=="24")]<-"chrY"
map$chromStart<-map$V4
map$chromEnd<-map$V4+1

map$chromStart<-format(map$chromStart, scientific=F)
map$chromEnd<-format(map$chromEnd, scientific=F)


write.table(map[,c(5:7,2)],"merged_plates_hg19_postqc.bed",col.names=F,row.names=F,quote=F,sep="\t")

#### lift positions

system("/software/team152/liftover/liftOver merged_plates_hg19_postqc.bed /software/team152/liftover/hg19ToHg38.over.chain.gz merged_plates_hg19_postqc_lifted_hg38 merged_plates_hg19_postqc_no_lifted_hg38")

system("cut -f 4 merged_plates_hg19_postqc_no_lifted_hg38 | sed \"/^#/d\" > nonlifted_hg38_variants_to_exclude_tmp.dat")

system("grep \"alt\" merged_plates_hg19_postqc_lifted_hg38 | cut -f 4 | cat - nonlifted_hg38_variants_to_exclude_tmp.dat > nonlifted_hg38_variants_to_exclude.dat")

### exclude non lifted:

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.05_nodupsample_eur_hwe --exclude nonlifted_hg38_variants_to_exclude.dat --recode tab --out merged_plates_hg19_postqc_lifted_hg38_liftedvariants")

## R

bed_lifted<-read.table("merged_plates_hg19_postqc_lifted_hg38",head=F)
excluded<-read.table("nonlifted_hg38_variants_to_exclude.dat",head=F)
bed_lifted<-bed_lifted[which(!bed_lifted$V4 %in% excluded$V1),]
map<-read.table("merged_plates_hg19_postqc_lifted_hg38_liftedvariants.map",head=F)

map$pos<-bed_lifted$V2

write.table(map[,c(1:3,5)],"merged_plates_hg19_postqc_lifted_hg38_liftedvariants_edited.map",col.names=F,row.names=F,quote=F,sep="\t")

########

# put together updated .ped plus lifted.map
system("/software/team152/plink_linux_x86_64_20181202/./plink --ped merged_plates_hg19_postqc_lifted_hg38_liftedvariants.ped --map merged_plates_hg19_postqc_lifted_hg38_liftedvariants_edited.map --merge-x 'no-fail' --make-bed --out merged_plates_postqc_lifted_hg38")

######################################################
# 16.1 - REMOVE MONOMORPHIC VARIANTS

### R

bim<-read.table("merged_plates_postqc_lifted_hg38.bim",head=F)

monom<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]

write.table(monom[,"V2",drop=F],"list_indel_var_exclude_2",col.names=F,row.names=F,quote=F,sep="\t")

######################################################
# 16.2 - CREATE FILE FORCING A1 ALLELE TO BE REFERENCE

# force A1 and A2 to be ref and alt alleles:

system("zcat merged_plates_hg19_posstrandaligned.vcf.gz | cut -f '1-5' | awk '{print $1\":\"$2\"_\"$4\"_\"$5,$4}' | sed '/^##/d' > list_variants_merged_plates_hg19_posstrandaligned_with_A1")

system("sort list_variants_merged_plates_hg19_posstrandaligned_with_A1 | uniq > list_variants_merged_plates_hg19_posstrandaligned_with_A1_ed")

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_postqc_lifted_hg38 --exclude list_indel_var_exclude_2 --make-bed --out merged_plates_postqc_lifted_hg38_nomonom")

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_postqc_lifted_hg38_nomonom --allow-no-sex --a2-allele list_variants_merged_plates_hg19_posstrandaligned_with_A1_ed --make-bed --out merged_plates_postqc_lifted_hg38_nomonom_RefAlt")

##############################################################
# 16.1 DOUBLE CHECK REF ALT ALLELE WITH HG38 REF SEQUENCE:

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_postqc_lifted_hg38_nomonom_RefAlt --allow-no-sex --keep-allele-order --output-chr M --recode vcf-iid --out merged_plates_postqc_lifted_hg38_nomonom_RefAlt")

system("export BCFTOOLS_PLUGINS=/software/hgi/installs/anaconda3/envs/hgi_base/libexec/bcftools/; /software/team152/bcftools-1.9/./bcftools +fixref merged_plates_postqc_lifted_hg38_nomonom_RefAlt.vcf -Oz -o merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.vcf.gz -- -f /lustre/scratch123/hgi/projects/ibdgwas_ukibdgc/resources/hg38/hg38_edited.fa -m top")

# REMOVE MISSMATCH ALLELES, THEY  IMPUTATION SERVER TO CRASH:

#### VCF to BED

system("/software/team152/plink_linux_x86_64_20181202/./plink --vcf merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.vcf.gz --keep-allele-order --allow-no-sex --double-id --make-bed --out merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned")

##############################################################
# 16.2 UPDATE NAME VARIANTS TO CHR:POSITION_REF_ALT in b38

system("zcat merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.vcf.gz | cut -f '1-5' | awk '{print $3,$1\":\"$2\"_\"$4\"_\"$5}' > list_variants_merged_plates_hg38_posstrandaligned")

### R
bim<-read.table("merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.bim",sep="\t",head=F)
ids<-read.table("list_variants_merged_plates_hg38_posstrandaligned",sep=" ",head=F,skip=32)

colnames(ids)[2]<-"ids"

#the major allele is set to A2 by default by Plink, keep ids with real ref/alt as in vcf using the ids file
bim.1<-cbind(bim,ids[,"ids",drop=F])

bim.1$ids<-gsub("X:","23:",bim.1$ids)
write.table(bim.1[,c(1,7,3:6)],"merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_edited.bim",col.names=F,row.names=F,quote=F,sep="\t")

##############################################################
# 16.3 RENAME VARIANTS USING B38

system("/software/team152/plink_linux_x86_64_20181202/./plink --bed merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.bed --bim merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_edited.bim --fam merged_plates_postqc_lifted_hg38_nomonom_RefAlt.fam --keep-allele-order --allow-no-sex --freq counts --make-bed --out merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated")

###################################################################################################################################################################################

##################################
# 17.- KEEP ONLY TOPMed VARIANTS #
##################################

MEM=40000;bsub -J"TOPMed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 -q normal -e stderr_TOPMed_subset -o stdout_TOPMed_subset "awk 'NR==FNR{vals[\$2];next} (\$9) in vals' merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated.bim <(zcat /lustre/scratch123/hgi/projects/ibdgwas_ukibdgc/resources/TOPMed/PASS.Variantsbravo-dbsnp-all_edited_3.vcf.gz) > merged_plates_TOPMed_variants"

system("cat merged_plates_TOPMed_variants | cut -f 9 > merged_plates_TOPMed_variants_ids")

### R

bim<-read.table("merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated.bim",head=F)
tm<-read.table("merged_plates_TOPMed_variants",head=F)
frq<-read.table("merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated.frq.counts",head=T)


frq$MAF<-frq$C1/(frq$C1+frq$C2)
frq$MAF[which(frq$MAF>0.5)]<-1-frq$MAF[which(frq$MAF>0.5)]

#######

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated --keep-allele-order --allow-no-sex --extract merged_plates_TOPMed_variants_ids --make-bed --out merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed")


###################################################################################################################################################################################


####################################
# 18.- DOUBLE CHEK A/T C/G ALLELES #
####################################


system("awk -v OFS='\t' '{print $2,$6}' merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim > list_variants_merged_plates_hg38_posstrandaligned_with_A1")

### EUR ctr

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed --keep-allele-order --allow-no-sex --freq counts --out merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur_ctr")


# GNOMAD:
MEM=40000; bsub -J "gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 -q normal -e stderr_gnomad_subset -o stdout_gnomad_subset "awk 'NR==FNR{vals[\$2];next} (\$1) in vals' merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim <( zcat /lustre/scratch123/hgi/projects/ibdgwas_ukibdgc/resources/gnomad/gnomad_freq_edited.gz) > merged_plates_gnomad_variants"

### R

gnomad<-read.table("merged_plates_gnomad_variants",head=F)
colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))

vfreq<-read.table("merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur_ctr.frq.counts",head=T)
vfreq$freq_alt<-vfreq$C1/(vfreq$C1+vfreq$C2)
vfreq$freq_alt<-as.numeric(vfreq$freq_alt)

all<-merge(vfreq,gnomad,by="SNP",sort=F)

remove1<-all[which( (all$freq_alt>=0.45 & all$freq_alt<=0.55) &
                      ( (all$A1=="G" & all$A2=="C") | (all$A1=="C" & all$A2=="G") |
                          (all$A1=="A" & all$A2=="T") | (all$A1=="T" & all$A2=="A") ) ),]

# remove A/T G/C that cannot be evaluated:
remove2<-vfreq[which(!(vfreq$SNP %in% gnomad$SNP) & ( (vfreq$A1=="G" & vfreq$A2=="C") | (vfreq$A1=="C" & vfreq$A2=="G") |
                                                        (vfreq$A1=="A" & vfreq$A2=="T") | (vfreq$A1=="T" & vfreq$A2=="A") ) ),]

remove<-rbind(remove1[,"SNP",drop=F],remove2[,"SNP",drop=F])

flip<-all[which( (all$freq_alt<0.45 | all$freq_alt>0.55) &
                   ( (all$A1=="G" & all$A2=="C") | (all$A1=="C" & all$A2=="G") |
                       (all$A1=="A" & all$A2=="T") | (all$A1=="T" & all$A2=="A") ) ),]
flip<-flip[which( (flip$freq_alt>0.5 & flip$AF_nfe<0.5) | (flip$freq_alt<0.5 & flip$AF_nfe>0.5) ),]

write.table(remove[,"SNP",drop=F],"list_variants_to_remove_AT_CG_merged_plates",col.names=F,row.names=F,quote=F,sep="\t")
write.table(flip[,"SNP",drop=F],"list_variants_to_flip_AT_CG_merged_plates",col.names=F,row.names=F,quote=F,sep="\t")


##########################################################################################
# 18.2 REMOVE VARIANTS WE CANNOT BE SURE ARE IN THE RIGHT STRAND, AND FLIP THE OTHERS

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed --allow-no-sex --exclude list_variants_to_remove_AT_CG_merged_plates --flip list_variants_to_flip_AT_CG_merged_plates --make-bed --out merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip")

system("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip --allow-no-sex --a2-allele list_variants_merged_plates_hg38_posstrandaligned_with_A1 --make-bed --out merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2")

#########
## Prepare TOPMED server
#########

for (i in seq(1,23)) {
  system(sprintf("/software/team152/plink_linux_x86_64_20181202/./plink --bfile merged_plates_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2 --allow-no-sex --keep-allele-order --chr %s --output-chr chr26 --recode vcf-iid --out ./submit_for_imputation/merged_plates_chr%s",i,i))
}

for (i in seq(1,23)) {
  system(sprintf("/software/team152/bcftools-1.9/./bcftools sort ./submit_for_imputation/merged_plates_chr%s.vcf -Oz -o ./submit_for_imputation/merged_plates_chr%s.vcf.gz",i,i))
}

############
# SUBMIT TO TOPMED
#######



########
# FINAL STEPS
########

# #Re-name the samples with mapping as before:
# #Read in the sample names from fam file
# sample_names = read.table("merged_plates_hg19_noind_posstr.fam",header = F)[1]
# sample_names[2] = gsub('[.]CEL.*', '', sample_names[1])
# gut_metadata = read.csv("GUT_scRNAseq_metadata - GUT_scRNAseq-cleaned.csv")
# prev_mapping = read.table("plate_1-3_mapping.txt")
# mapped_names = merge(sample_names, prev_mapping, by = 'V1',all.x = T)
# names(mapped_names) = paste0('V',1:4)
# 
# mapped_names$V4 <- ifelse(is.na(mapped_names$V4),
#                           mapped_names$V2, mapped_names$V4)
# 
# write.table(mapped_names[c(1,4)], 'plate_1-4_mapping.txt',
#             row.names = FALSE,col.names = FALSE, quote = FALSE)

#Merge the chr's together after imputation
system("/software/team152/bcftools-1.9/bcftools concat $(for i in {1..22} X; do echo \"chr$i.dose.vcf.gz\"; done) -Oz -o merged_chr.vcf.gz")
sytem("/software/team152/bcftools-1.9/bcftools reheader -s ../plate_1-4_mapping.txt merged_chr.vcf.gz | /software/team152/bcftools-1.9/bcftools view -s ^hb15_D45,hb15_D46,hb15_D47 -Oz -o CCF_OTAR-N_individuals-imputed-all_chr.vcf.gz")
system("mv CCF_OTAR-N_individuals-imputed-all_chr.vcf.gz CCF_OTAR-$(/software/team152/bcftools-1.9/bcftools query -l CCF_OTAR-N_individuals-imputed-all_chr.vcf.gz | wc -l)_individuals-imputed-all_chr.vcf.gz")




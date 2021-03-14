#! /bin/bash
#$ -V
#$ -cwd
#$ -o "log_master_amend.txt"
#$ -j y
#$ -l mem_free=20G,h_vmem=20G
#$ -pe threaded 8

module load eagle
module unload plink
module load plink2
module load minimac3
module load vcftools
module load R
module load tabix

#### Code modified for specific imputation desired

##############################
#### Creating Directories ####
##############################
mkdir Norwegian/chr${1}
mkdir NonNorwegian/chr${1}
mkdir SICCA/chr${1}
mkdir Phase1/chr${1}
mkdir Reference/Full/chr${1}
mkdir Reference/North/chr${1}
mkdir Reference/South/chr${1}

###############################
#### Reference Panel Setup ####
###############################
### convert reference panel files to plink files (ped/map)
vcftools --gzvcf /.../1000g/phase3_20130502/ALL.chr${1}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --plink --chr ${1} --out Reference/Full/chr${1}/ref_chr${1}

### read plink reference panel and make binary files
plink --file Reference/Full/chr${1}/ref_chr${1} --make-bed --out Reference/Full/chr${1}/ref_chr${1}
### zip files
plink --bfile Reference/Full/chr${1}/ref_chr${1} --recode vcf-iid --out  Reference/Full/chr${1}/ref_chr${1}
bgzip Reference/Full/chr${1}/ref_chr${1}.vcf
tabix Reference/Full/chr${1}/ref_chr${1}.vcf.gz -f

### filter reference panel for North European and South European samples
### North
plink --bfile Reference/Full/chr${1}/ref_chr${1} --keep 1000g_north_eur_samps.txt --make-bed --out Reference/North/chr${1}/north_eur_ref_chr${1}
### zip files
plink --bfile Reference/North/chr${1}/north_eur_ref_chr${1} --recode vcf-iid --out  Reference/North/chr${1}/north_eur_ref_chr${1}
bgzip Reference/North/chr${1}/north_eur_ref_chr${1}.vcf
tabix Reference/North/chr${1}/north_eur_ref_chr${1}.vcf.gz -f

### South
plink --bfile Reference/Full/chr${1}/ref_chr${1} --keep 1000g_south_eur_samps.txt --make-bed --out Reference/South/chr${1}/south_eur_ref_chr${1}
### zip files
plink --bfile Reference/South/chr${1}/south_eur_ref_chr${1} --recode vcf-iid --out  Reference/South/chr${1}/south_eur_ref_chr${1}
bgzip Reference/South/chr${1}/south_eur_ref_chr${1}.vcf
tabix Reference/South/chr${1}/south_eur_ref_chr${1}.vcf.gz -f

#########################
#### GWAS Data Setup ####
#########################

### Read in datasets to be imputed ###
### Norwegian GWAS
plink --bfile ../01_GWAS_Data/02_GWAS_Norwegian/SS_GWAS_Norwegian_Merged   --chr ${1} --make-bed --out Norwegian/chr${1}/chr${1}_nor
### Remove SNPs with mendelian errors
plink --bfile Norwegian/chr${1}/chr${1}_nor --mendel --out Norwegian/nor${1}
awk '$3>0' Norwegian/nor${1}.lmendel > Norwegian/mend${1}_snps_nor.txt
plink --bfile Norwegian/chr${1}/chr${1}_nor --exclude Norwegian/mend${1}_snps_nor.txt --make-bed --out Norwegian/chr${1}/chr${1}_no_mend_nor

### Non-Norwegian GWAS 
plink --bfile ../01_GWAS_Data/03_GWAS_Non_Norwegian/SS_GWAS_Non_Norwegian_Merged   --chr ${1} --make-bed --out NonNorwegian/chr${1}/chr${1}_non_nor
### Remove SNPs with mendelian errors
plink --bfile NonNorwegian/chr${1}/chr${1}_non_nor --mendel --out NonNorwegian/non_nor${1}
awk '$3>0' NonNorwegian/non_nor${1}.lmendel > NonNorwegian/mend${1}_snps_non_nor.txt
plink --bfile NonNorwegian/chr${1}/chr${1}_non_nor --exclude NonNorwegian/mend${1}_snps_non_nor.txt --make-bed --out NonNorwegian/chr${1}/chr${1}_no_mend_non_nor

### SICCA GWAS 
plink --bfile ../01_GWAS_Data/04_GWAS_SICCA/SS_GWAS_SICCA   --chr ${1} --make-bed --out SICCA/chr${1}/chr${1}_sicca
### Remove SNPs with mendelian errors
plink --bfile SICCA/chr${1}/chr${1}_sicca --mendel --out SICCA/sicca${1}
awk '$3>0' SICCA/sicca${1}.lmendel > SICCA/mend${1}_snps_sicca.txt
plink --bfile SICCA/chr${1}/chr${1}_sicca --exclude SICCA/mend${1}_snps_sicca.txt --make-bed --out SICCA/chr${1}/chr${1}_no_mend_sicca

### Phase1 GWAS
plink --bfile ../01_GWAS_Data/01_GWAS_Phase1/SS_GWAS_Phase1   --chr ${1} --make-bed --out Phase1/chr${1}/chr${1}_p1
### Remove SNPs with mendelian errors
plink --bfile Phase1/chr${1}/chr${1}_p1 --mendel --out Phase1/p1${1}
awk '$3>0' Phase1/p1${1}.lmendel > Phase1/mend${1}_snps_p1.txt
plink --bfile Phase1/chr${1}/chr${1}_p1 --exclude Phase1/mend${1}_snps_p1.txt --make-bed --out Phase1/chr${1}/chr${1}_no_mend_p1

### Execute R script to identify snps to have their strand flipped so they match reference panel
Rscript --vanilla pipeline3_amend.R ${1}

### Flip snps identified in previous line and remove snps mismatched to reference panel ###
### Norwegian
plink --bfile Norwegian/chr${1}/chr${1}_no_mend_nor --flip Norwegian/chr${1}/nor_snps_to_flip${1}.txt --extract Norwegian/chr${1}/nor_extract${1}.txt --make-bed --out Norwegian/chr${1}/chr${1}_update_nor

### Non-Norwegian
plink --bfile NonNorwegian/chr${1}/chr${1}_no_mend_non_nor --flip NonNorwegian/chr${1}/non_nor_snps_to_flip${1}.txt --extract NonNorwegian/chr${1}/non_nor_extract${1}.txt --make-bed --out NonNorwegian/chr${1}/chr${1}_update_non_nor

###SICCA
plink --bfile SICCA/chr${1}/chr${1}_no_mend_sicca --flip SICCA/chr${1}/sicca_snps_to_flip${1}.txt --extract SICCA/chr${1}/sicca_extract${1}.txt --make-bed --out SICCA/chr${1}/chr${1}_update_sicca

###Phase1
plink --bfile Phase1/chr${1}/chr${1}_no_mend_p1 --flip Phase1/chr${1}/p1_snps_to_flip${1}.txt --extract Phase1/chr${1}/p1_extract${1}.txt --make-bed --out Phase1/chr${1}/chr${1}_update_p1

### Update SNP names in all data sets to match reference panel ###
### Norwegian
plink --bfile Norwegian/chr${1}/chr${1}_update_nor --update-name Norwegian/chr${1}/nor_upSNP${1}.txt --make-bed --out Norwegian/chr${1}/chr${1}_nor2
### Zip files
plink --bfile Norwegian/chr${1}/chr${1}_nor2 --recode vcf-iid --out  Norwegian/chr${1}/chr${1}_nor2
bgzip Norwegian/chr${1}/chr${1}_nor2.vcf
tabix Norwegian/chr${1}/chr${1}_nor2.vcf.gz -f

### Non-Norwegian
plink --bfile NonNorwegian/chr${1}/chr${1}_update_non_nor --update-name NonNorwegian/chr${1}/non_nor_upSNP${1}.txt --make-bed --out NonNorwegian/chr${1}/chr${1}_non_nor2
### Zip files
plink --bfile NonNorwegian/chr${1}/chr${1}_non_nor2 --recode vcf-iid --out  NonNorwegian/chr${1}/chr${1}_non_nor2
bgzip NonNorwegian/chr${1}/chr${1}_non_nor2.vcf
tabix NonNorwegian/chr${1}/chr${1}_non_nor2.vcf.gz -f

### SICCA
plink --bfile SICCA/chr${1}/chr${1}_update_sicca --update-name SICCA/chr${1}/sicca_upSNP${1}.txt --make-bed --out SICCA/chr${1}/chr${1}_sicca2
### Zip files
plink --bfile SICCA/chr${1}/chr${1}_sicca2 --recode vcf-iid --out  SICCA/chr${1}/chr${1}_sicca2
bgzip SICCA/chr${1}/chr${1}_sicca2.vcf
tabix SICCA/chr${1}/chr${1}_sicca2.vcf.gz -f

### Phase1
plink --bfile Phase1/chr${1}/chr${1}_update_p1 --update-name Phase1/chr${1}/p1_upSNP${1}.txt --make-bed --out Phase1/chr${1}/chr${1}_p12
### Zip files
plink --bfile Phase1/chr${1}/chr${1}_p12 --recode vcf-iid --out  Phase1/chr${1}/chr${1}_p12
bgzip Phase1/chr${1}/chr${1}_p12.vcf
tabix Phase1/chr${1}/chr${1}_p12.vcf.gz -f

### Create txt file to change snp name from chr:pos to rsID ###
### Norwegian
awk '{print $1":"$4"\t"$2}' Norwegian/chr${1}/chr${1}_no_mend_nor.bim > Norwegian/chr${1}/rename_snps_chr${1}.txt

### Non-Norwegian
awk '{print $1":"$4"\t"$2}' NonNorwegian/chr${1}/chr${1}_no_mend_non_nor.bim > NonNorwegian/chr${1}/rename_snps_chr${1}.txt

### SICCA
awk '{print $1":"$4"\t"$2}' SICCA/chr${1}/chr${1}_no_mend_sicca.bim > SICCA/chr${1}/rename_snps_chr${1}.txt

### Phase1
awk '{print $1":"$4"\t"$2}' Phase1/chr${1}/chr${1}_no_mend_p1.bim > Phase1/chr${1}/rename_snps_chr${1}.txt

######################
#### Phasing Data ####
######################
### Reference panel phasing ###
### Full
eagle --vcfTarget=Reference/Full/chr${1}/ref_chr${1}.vcf.gz --vcfRef=Reference/Full/chr${1}/ref_chr${1}.vcf.gz --geneticMapFile=/Volumes/hts_core/Shared/eagle-tables/genetic_map_hg19_withX.txt.gz --outPrefix=Reference/Full/chr${1}/ref_chr${1}_phased --numThreads=8 --vcfOutFormat=v --noImpMissing
### Convert the phased North reference panel data to VCF format (needed for minimac3)
bgzip Reference/Full/chr${1}/ref_chr${1}_phased.vcf
tabix Reference/Full/chr${1}/ref_chr${1}_phased.vcf.gz -f

### North
eagle --vcfTarget=Reference/North/chr${1}/north_eur_ref_chr${1}.vcf.gz --vcfRef=Reference/North/chr${1}/north_eur_ref_chr${1}.vcf.gz --geneticMapFile=/Volumes/hts_core/Shared/eagle-tables/genetic_map_hg19_withX.txt.gz --outPrefix=Reference/North/chr${1}/north_eur_ref_chr${1}_phased --numThreads=8 --vcfOutFormat=v --noImpMissing
### Convert the phased North reference panel data to VCF format (needed for minimac3)
bgzip Reference/North/chr${1}/north_eur_ref_chr${1}_phased.vcf
tabix Reference/North/chr${1}/north_eur_ref_chr${1}_phased.vcf.gz -f

###South
eagle --vcfTarget=Reference/South/chr${1}/south_eur_ref_chr${1}.vcf.gz --vcfRef=Reference/South/chr${1}/south_eur_ref_chr${1}.vcf.gz --geneticMapFile=/Volumes/hts_core/Shared/eagle-tables/genetic_map_hg19_withX.txt.gz --outPrefix=Reference/South/chr${1}/south_eur_ref_chr${1}_phased --numThreads=8 --vcfOutFormat=v --noImpMissing
### Convert the phased South reference panel data to VCF format (needed for minimac3)
bgzip Reference/South/chr${1}/south_eur_ref_chr${1}_phased.vcf
tabix Reference/South/chr${1}/south_eur_ref_chr${1}_phased.vcf.gz -f

### GWAS phasing ###
### Norwegian
eagle --vcfTarget=Norwegian/chr${1}/chr${1}_nor2.vcf.gz --vcfRef=Reference/North/chr${1}/north_eur_ref_chr${1}_phased.vcf.gz --geneticMapFile=/Volumes/hts_core/Shared/eagle-tables/genetic_map_hg19_withX.txt.gz --outPrefix=Norwegian/chr${1}/chr${1}_nor_phased --numThreads=8 --vcfOutFormat=v --noImpMissing
### Convert the phased GWAS data to VCF format (needed for minimac3)
bgzip Norwegian/chr${1}/chr${1}_nor_phased.vcf
tabix Norwegian/chr${1}/chr${1}_nor_phased.vcf.gz -f

### Non-Norwegian
eagle --vcfTarget=NonNorwegian/chr${1}/chr${1}_non_nor2.vcf.gz --vcfRef=Reference/South/chr${1}/south_eur_ref_chr${1}_phased.vcf.gz --geneticMapFile=/Volumes/hts_core/Shared/eagle-tables/genetic_map_hg19_withX.txt.gz --outPrefix=NonNorwegian/chr${1}/chr${1}_non_nor_phased --numThreads=8 --vcfOutFormat=v --noImpMissing
### Convert the phased GWAS data to VCF format (needed for minimac3)
bgzip NonNorwegian/chr${1}/chr${1}_non_nor_phased.vcf
tabix NonNorwegian/chr${1}/chr${1}_non_nor_phased.vcf.gz -f

### SICCA
eagle --vcfTarget=SICCA/chr${1}/chr${1}_sicca2.vcf.gz --vcfRef=Reference/Full/chr${1}/ref_chr${1}_phased.vcf.gz --geneticMapFile=/Volumes/hts_core/Shared/eagle-tables/genetic_map_hg19_withX.txt.gz --outPrefix=SICCA/chr${1}/chr${1}_sicca_phased --numThreads=8 --vcfOutFormat=v --noImpMissing
### Convert the phased GWAS data to VCF format (needed for minimac3)
bgzip SICCA/chr${1}/chr${1}_sicca_phased.vcf
tabix SICCA/chr${1}/chr${1}_sicca_phased.vcf.gz -f

### Phase1
eagle --vcfTarget=Phase1/chr${1}/chr${1}_p12.vcf.gz --vcfRef=Reference/Full/chr${1}/ref_chr${1}_phased.vcf.gz --geneticMapFile=/Volumes/hts_core/Shared/eagle-tables/genetic_map_hg19_withX.txt.gz --outPrefix=Phase1/chr${1}/chr${1}_p1_phased --numThreads=8 --vcfOutFormat=v --noImpMissing
### Convert the phased GWAS data to VCF format (needed for minimac3)
bgzip Phase1/chr${1}/chr${1}_p1_phased.vcf
tabix Phase1/chr${1}/chr${1}_p1_phased.vcf.gz -f

####################
#### Imputation ####
####################

### Norwegian data using north reference panel
minimac3-omp --refHaps Reference/North/chr${1}/north_eur_ref_chr${1}_phased.vcf.gz --haps Norwegian/chr${1}/chr${1}_nor_phased.vcf.gz --prefix Norwegian/chr${1}/chr${1}_imputed --cpus 8 --noPhoneHome

### Non-Norwegian data using south reference panel
minimac3-omp --refHaps Reference/South/chr${1}/south_eur_ref_chr${1}_phased.vcf.gz --haps NonNorwegian/chr${1}/chr${1}_non_nor_phased.vcf.gz --prefix NonNorwegian/chr${1}/chr${1}_imputed --cpus 8 --noPhoneHome

### SICCA data using full reference panel
minimac3-omp --refHaps Reference/Full/chr${1}/ref_chr${1}_phased.vcf.gz --haps SICCA/chr${1}/chr${1}_sicca_phased.vcf.gz --prefix SICCA/chr${1}/chr${1}_imputed --cpus 8 --noPhoneHome

### Phase1 data using full referene panel
minimac3-omp --refHaps Reference/Full/chr${1}/ref_chr${1}_phased.vcf.gz --haps Phase1/chr${1}/chr${1}_p1_phased.vcf.gz --prefix Phase1/chr${1}/chr${1}_imputed --cpus 8 --noPhoneHome

### Convert imputed data to transposed plink format ###
###Norwegian
vcftools --gzvcf Norwegian/chr${1}/chr${1}_imputed.dose.vcf.gz --plink-tped --out Norwegian/chr${1}/chr${1}_imputed
plink --tfile Norwegian/chr${1}/chr${1}_imputed --make-bed --out Norwegian/chr${1}/chr${1}_imputed

### Non-Norwegian
vcftools --gzvcf NonNorwegian/chr${1}/chr${1}_imputed.dose.vcf.gz --plink-tped --out NonNorwegian/chr${1}/chr${1}_imputed
plink --tfile NonNorwegian/chr${1}/chr${1}_imputed --make-bed --out NonNorwegian/chr${1}/chr${1}_imputed

### SICCA
vcftools --gzvcf SICCA/chr${1}/chr${1}_imputed.dose.vcf.gz --plink-tped --out SICCA/chr${1}/chr${1}_imputed
plink --tfile SICCA/chr${1}/chr${1}_imputed --make-bed --out SICCA/chr${1}/chr${1}_imputed

### Phase1
vcftools --gzvcf Phase1/chr${1}/chr${1}_imputed.dose.vcf.gz --plink-tped --out Phase1/chr${1}/chr${1}_imputed
plink --tfile Phase1/chr${1}/chr${1}_imputed --make-bed --out Phase1/chr${1}/chr${1}_imputed

############################
#### Creating txt files ####
############################
### create txt files to correct FID, PID, MID, SEX, and Phenotype to match what originally appears in the data ###
### Norwegian
awk '{print $2, $2, $1, $2}' Norwegian/chr${1}/chr${1}_nor2.fam > Norwegian/chr${1}/upIDs.txt
awk '{print $1, $2, $5}' Norwegian/chr${1}/chr${1}_nor2.fam > Norwegian/chr${1}/upSex.txt
awk '{print $1, $2, $3, $4}' Norwegian/chr${1}/chr${1}_nor2.fam > Norwegian/chr${1}/upPar.txt
awk '{print $1, $2, $6}' Norwegian/chr${1}/chr${1}_nor2.fam > Norwegian/chr${1}/upPheno.txt

### Non-Norwegian
awk '{print $2, $2, $1, $2}' NonNorwegian/chr${1}/chr${1}_non_nor2.fam > NonNorwegian/chr${1}/upIDs.txt
awk '{print $1, $2, $5}' NonNorwegian/chr${1}/chr${1}_non_nor2.fam > NonNorwegian/chr${1}/upSex.txt
awk '{print $1, $2, $3, $4}' NonNorwegian/chr${1}/chr${1}_non_nor2.fam > NonNorwegian/chr${1}/upPar.txt
awk '{print $1, $2, $6}' NonNorwegian/chr${1}/chr${1}_non_nor2.fam > NonNorwegian/chr${1}/upPheno.txt

### SICCA
awk '{print $2, $2, $1, $2}' SICCA/chr${1}/chr${1}_sicca2.fam > SICCA/chr${1}/upIDs.txt
awk '{print $1, $2, $5}' SICCA/chr${1}/chr${1}_sicca2.fam > SICCA/chr${1}/upSex.txt
awk '{print $1, $2, $3, $4}' SICCA/chr${1}/chr${1}_sicca2.fam > SICCA/chr${1}/upPar.txt
awk '{print $1, $2, $6}' SICCA/chr${1}/chr${1}_sicca2.fam > SICCA/chr${1}/upPheno.txt

### Phase1
awk '{print $2, $2, $1, $2}' Phase1/chr${1}/chr${1}_p12.fam > Phase1/chr${1}/upIDs.txt
awk '{print $1, $2, $5}' Phase1/chr${1}/chr${1}_p12.fam > Phase1/chr${1}/upSex.txt
awk '{print $1, $2, $3, $4}' Phase1/chr${1}/chr${1}_p12.fam > Phase1/chr${1}/upPar.txt
awk '{print $1, $2, $6}' Phase1/chr${1}/chr${1}_p12.fam > Phase1/chr${1}/upPheno.txt

### update info on imputed data ###
###Remove duplicates from imputed and rename_snps for all datasets
Rscript --vanilla Duplicate_IDs.R ${1}

### Norwegian
plink --bfile Norwegian/chr${1}/chr${1}_imputed --exclude Norwegian/chr${1}/duplicates.txt --make-bed --out Norwegian/chr${1}/chr${1}_imputed_temp
# Update IDS
plink --bfile Norwegian/chr${1}/chr${1}_imputed_temp --update-ids Norwegian/chr${1}/upIDs.txt --make-bed --out Norwegian/chr${1}/chr${1}_imputed_temp
# Update sex, parents, pheno, snp names
plink --bfile Norwegian/chr${1}/chr${1}_imputed_temp --update-sex Norwegian/chr${1}/upSex.txt --update-parents Norwegian/chr${1}/upPar.txt --pheno Norwegian/chr${1}/upPheno.txt --update-name Norwegian/chr${1}/rename_snps_chr${1}_2.txt --make-bed --out Norwegian/chr${1}/chr${1}_imputed_updated

### Non-Norwegian
plink --bfile NonNorwegian/chr${1}/chr${1}_imputed --exclude NonNorwegian/chr${1}/duplicates.txt --make-bed --out NonNorwegian/chr${1}/chr${1}_imputed_temp
#update IDS
plink --bfile NonNorwegian/chr${1}/chr${1}_imputed_temp --update-ids NonNorwegian/chr${1}/upIDs.txt --make-bed --out NonNorwegian/chr${1}/chr${1}_imputed_temp
#update sex, parents, pheno, snp names
plink --bfile NonNorwegian/chr${1}/chr${1}_imputed_temp --update-sex NonNorwegian/chr${1}/upSex.txt --update-parents NonNorwegian/chr${1}/upPar.txt --pheno NonNorwegian/chr${1}/upPheno.txt --update-name NonNorwegian/chr${1}/rename_snps_chr${1}_2.txt --make-bed --out NonNorwegian/chr${1}/chr${1}_imputed_updated

### SICCA
plink --bfile SICCA/chr${1}/chr${1}_imputed --exclude SICCA/chr${1}/duplicates.txt --make-bed --out SICCA/chr${1}/chr${1}_imputed_temp
# update IDS
plink --bfile SICCA/chr${1}/chr${1}_imputed_temp --update-ids SICCA/chr${1}/upIDs.txt --make-bed --out SICCA/chr${1}/chr${1}_imputed_temp
#update sex, parents, pheno, snp names
plink --bfile SICCA/chr${1}/chr${1}_imputed_temp --update-sex SICCA/chr${1}/upSex.txt --update-parents SICCA/chr${1}/upPar.txt --pheno SICCA/chr${1}/upPheno.txt --update-name SICCA/chr${1}/rename_snps_chr${1}_2.txt --make-bed --out SICCA/chr${1}/chr${1}_imputed_updated

### Phase1
plink --bfile Phase1/chr${1}/chr${1}_imputed --exclude Phase1/chr${1}/duplicates.txt --make-bed --out Phase1/chr${1}/chr${1}_imputed_temp
#update IDS
plink --bfile Phase1/chr${1}/chr${1}_imputed_temp --update-ids Phase1/chr${1}/upIDs.txt --make-bed --out Phase1/chr${1}/chr${1}_imputed_temp
#update sex, parents, pheno, snp names
plink --bfile Phase1/chr${1}/chr${1}_imputed_temp --update-sex Phase1/chr${1}/upSex.txt --update-parents Phase1/chr${1}/upPar.txt --pheno Phase1/chr${1}/upPheno.txt --update-name Phase1/chr${1}/rename_snps_chr${1}_2.txt --make-bed --out Phase1/chr${1}/chr${1}_imputed_updated

#################
#### Storage ####
#################
### Norwegian
cd Completed_Imputation/Norwegian

cp ../../Norwegian/chr${1}/chr${1}_imputed_updated.bed chr${1}_imputed_updated.bed
cp ../../Norwegian/chr${1}/chr${1}_imputed_updated.bim chr${1}_imputed_updated.bim
cp ../../Norwegian/chr${1}/chr${1}_imputed_updated.fam chr${1}_imputed_updated.fam

### Non-Norwegian
cd ../NonNorwegian

cp ../../NonNorwegian/chr${1}/chr${1}_imputed_updated.bed chr${1}_imputed_updated.bed
cp ../../NonNorwegian/chr${1}/chr${1}_imputed_updated.bim chr${1}_imputed_updated.bim
cp ../../NonNorwegian/chr${1}/chr${1}_imputed_updated.fam chr${1}_imputed_updated.fam

### SICCA
cd ../SICCA

cp ../../SICCA/chr${1}/chr${1}_imputed_updated.bed chr${1}_imputed_updated.bed
cp ../../SICCA/chr${1}/chr${1}_imputed_updated.bim chr${1}_imputed_updated.bim
cp ../../SICCA/chr${1}/chr${1}_imputed_updated.fam chr${1}_imputed_updated.fam

### Phase1
cd ../Phase1

cp ../../Phase1/chr${1}/chr${1}_imputed_updated.bed chr${1}_imputed_updated.bed
cp ../../Phase1/chr${1}/chr${1}_imputed_updated.bim chr${1}_imputed_updated.bim
cp ../../Phase1/chr${1}/chr${1}_imputed_updated.fam chr${1}_imputed_updated.fam

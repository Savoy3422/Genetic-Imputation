#### Extract command line argument of chromosome number ####
args=commandArgs(trailingOnly=T)        #Only arguments after 'args' be returned (${1} value)
chr=as.numeric(args[1])                 #Create chr number as numeric using the previous command and the first argument

######################################################
#### Reading in datasets with no mendelian errors ####
######################################################
#### Reference Panels ####
ref=read.table(paste0("/Volumes/qac/prj/sjogrens/Imputation/Reference/Full/chr",chr,"/ref_chr",chr,".bim"), header=F)               
dupes=which(duplicated(ref$V4))              #Locating duplicates at base-pair position on chr
if(length(dupes)>0) ref=ref[-dupes,]   #Removing duplicates from dataset

ref_north=read.table(paste0("/Volumes/qac/prj/sjogrens/Imputation/Reference/North/chr",chr,"/north_eur_ref_chr",chr,".bim"), header=F)               
dupes=which(duplicated(ref_north$V4))              #Locating duplicates at base-pair position on chr
if(length(dupes)>0) ref_north=ref_north[-dupes,]   #Removing duplicates from dataset

ref_south=read.table(paste0("/Volumes/qac/prj/sjogrens/Imputation/Reference/South/chr",chr,"/south_eur_ref_chr",chr,".bim"), header=F)               
dupes=which(duplicated(ref_south$V4))              #Locating duplicates at base-pair position on chr
if(length(dupes)>0) ref_south=ref_south[-dupes,]   #Removing duplicates from dataset

#### GWAS Datasets ####
nor=read.table(paste0("/Volumes/qac/prj/sjogrens/Imputation/Norwegian/chr",chr,"/chr",chr,"_no_mend_nor.bim"), header=F)
dupes=which(duplicated(nor$V4))             #Locating duplicates at base-pair position on chr
if(length(dupes)>0) nor=nor[-dupes,]        #Removing duplicates from dataset

non_nor=read.table(paste0("/Volumes/qac/prj/sjogrens/Imputation/NonNorwegian/chr",chr,"/chr",chr,"_no_mend_non_nor.bim"), header=F)
dupes=which(duplicated(non_nor$V4))             #Locating duplicates at base-pair position on chr
if(length(dupes)>0) non_nor=non_nor[-dupes,]        #Removing duplicates from dataset

sicca=read.table(paste0("/Volumes/qac/prj/sjogrens/Imputation/SICCA/chr",chr,"/chr",chr,"_no_mend_sicca.bim"), header=F)
dupes=which(duplicated(sicca$V4))             #Locating duplicates at base-pair position on chr
if(length(dupes)>0) sicca=sicca[-dupes,]        #Removing duplicates from dataset

p1=read.table(paste0("/Volumes/qac/prj/sjogrens/Imputation/Phase1/chr",chr,"/chr",chr,"_no_mend_p1.bim"), header=F)
dupes=which(duplicated(p1$V4))             #Locating duplicates at base-pair position on chr
if(length(dupes)>0) p1=p1[-dupes,]        #Removing duplicates from dataset

################################################
#### Comparing position SNP reference names ####
################################################
#### Norwegian ####
ref_north_nor1=ref_north[which(ref_north$V4 %in% nor$V4),]          #Keep positions in reference that are in norwegian dataset
ref_north_nor2=nor[which(nor$V4 %in% ref_north$V4),]                #Keep positions in norwegian dataset that are in reference
ref_north_nor=cbind(ref_north_nor2[,c(2,4)],ref_north_nor1[,c(2)],ref_north_nor2[,c(5,6)],ref_north_nor1[,c(5,6)]) #Combine columns into new dataset
names(ref_north_nor)=c("SNP_nor","Pos","SNP_ref_north","nor_ref","nor_alt","ref_ref","ref_alt")    #Name columns
ref_north_nor$ref_ref=as.character(ref_north_nor$ref_ref)                #Making all columns charaters
ref_north_nor$ref_alt=as.character(ref_north_nor$ref_alt)
ref_north_nor$nor_ref=as.character(ref_north_nor$nor_ref)
ref_north_nor$nor_alt=as.character(ref_north_nor$nor_alt)
ref_north_nor$SNP_nor=as.character(ref_north_nor$SNP_nor)
ref_north_nor$SNP_ref_north=as.character(ref_north_nor$SNP_ref_north)

#### Non-Norwegian ####
ref_south_non1=ref_south[which(ref_south$V4 %in% non_nor$V4),]          #Keep positions in reference that are in norwegian dataset
ref_south_non2=non_nor[which(non_nor$V4 %in% ref_south$V4),]                #Keep positions in norwegian dataset that are in reference
ref_south_non=cbind(ref_south_non2[,c(2,4)],ref_south_non1[,c(2)],ref_south_non2[,c(5,6)],ref_south_non1[,c(5,6)]) #Combine columns into new dataset
names(ref_south_non)=c("SNP_non","Pos","SNP_ref_south","non_ref","non_alt","ref_ref","ref_alt")    #Name columns
ref_south_non$ref_ref=as.character(ref_south_non$ref_ref)                #Making all columns charaters
ref_south_non$ref_alt=as.character(ref_south_non$ref_alt)
ref_south_non$non_ref=as.character(ref_south_non$non_ref)
ref_south_non$non_alt=as.character(ref_south_non$non_alt)
ref_south_non$SNP_non=as.character(ref_south_non$SNP_non)
ref_south_non$SNP_ref_south=as.character(ref_south_non$SNP_ref_south)

#### SICCA ####
ref_full_sicca1=ref[which(ref$V4 %in% sicca$V4),]          #Keep positions in reference that are in norwegian dataset
ref_full_sicca2=sicca[which(sicca$V4 %in% ref$V4),]                #Keep positions in norwegian dataset that are in reference
ref_full_sicca=cbind(ref_full_sicca2[,c(2,4)],ref_full_sicca1[,c(2)],ref_full_sicca2[,c(5,6)],ref_full_sicca1[,c(5,6)]) #Combine columns into new dataset
names(ref_full_sicca)=c("SNP_sicca","Pos","SNP_ref","sicca_ref","sicca_alt","ref_ref","ref_alt")    #Name columns
ref_full_sicca$ref_ref=as.character(ref_full_sicca$ref_ref)                #Making all columns charaters
ref_full_sicca$ref_alt=as.character(ref_full_sicca$ref_alt)
ref_full_sicca$sicca_ref=as.character(ref_full_sicca$sicca_ref)
ref_full_sicca$sicca_alt=as.character(ref_full_sicca$sicca_alt)
ref_full_sicca$SNP_sicca=as.character(ref_full_sicca$SNP_sicca)
ref_full_sicca$SNP_ref=as.character(ref_full_sicca$SNP_ref)

#### Phase1 ####
ref_full_p11=ref[which(ref$V4 %in% p1$V4),]          #Keep positions in reference that are in norwegian dataset
ref_full_p12=sicca[which(p1$V4 %in% ref$V4),]                #Keep positions in norwegian dataset that are in reference
ref_full_p1=cbind(ref_full_p12[,c(2,4)],ref_full_p11[,c(2)],ref_full_p12[,c(5,6)],ref_full_p11[,c(5,6)]) #Combine columns into new dataset
names(ref_full_p1)=c("SNP_p1","Pos","SNP_ref","p1_ref","p1_alt","ref_ref","ref_alt")    #Name columns
ref_full_p1$ref_ref=as.character(ref_full_p1$ref_ref)                #Making all columns charaters
ref_full_p1$ref_alt=as.character(ref_full_p1$ref_alt)
ref_full_p1$p1_ref=as.character(ref_full_p1$p1_ref)
ref_full_p1$p1_alt=as.character(ref_full_p1$p1_alt)
ref_full_p1$SNP_p1=as.character(ref_full_p1$SNP_p1)
ref_full_p1$SNP_ref=as.character(ref_full_p1$SNP_ref)

opp=function(x){                                    #Defining function to correct allele changes
  if(x=="A") return("T") else
    if(x=="T") return("A") else
      if(x=="C") return("G") else
        if(x=="G") return("C") else
          return("")
}

'%!in%'=function(x,y) !('%in%'(x,y))                 #Defining function for 'not in this set'

#################################
#### Correcting allele flips ####
#################################
#### Norwegian ####
at=which((ref_north_nor$nor_ref=="A" & ref_north_nor$nor_alt=="T") | (ref_north_nor$ref_ref=="T" & ref_north_nor$ref_alt=="A"))   #
ta=which((ref_north_nor$nor_ref=="T" & ref_north_nor$nor_alt=="A") | (ref_north_nor$ref_ref=="A" & ref_north_nor$ref_alt=="T"))
cg=which((ref_north_nor$nor_ref=="C" & ref_north_nor$nor_alt=="G") | (ref_north_nor$ref_ref=="G" & ref_north_nor$ref_alt=="C"))
gc=which((ref_north_nor$nor_ref=="G" & ref_north_nor$nor_alt=="C") | (ref_north_nor$ref_ref=="C" & ref_north_nor$ref_alt=="G"))
remove_snps=unique(c(at,ta,cg,gc))                       #List SNPs where solution isn't apparent
remove_snps=remove_snps[order(remove_snps)]              #Ordering SNPs to remove
if(length(remove_snps)>0) ref_north_nor=ref_north_nor[-remove_snps,]  #Removing problem SNPs #Format note:No SNPs to remove, empty vector causes problems
match=which(ref_north_nor$ref_ref==ref_north_nor$nor_ref & ref_north_nor$ref_alt==ref_north_nor$nor_alt) #SNPs where alleles match
minor=which(ref_north_nor$ref_ref==ref_north_nor$nor_alt & ref_north_nor$ref_alt==ref_north_nor$nor_ref) #SNPs where major and minor alleles are switched
flip=vector()                                            #Creating empty vector for SNPs that should be flipped
one=rep(0,nrow(ref_north_nor))
two=rep(0,nrow(ref_north_nor))
for(i in 1:nrow(ref_north_nor)){
  if(ref_north_nor$ref_ref[i]==opp(ref_north_nor$nor_ref[i]) & ref_north_nor$ref_alt[i]==opp(ref_north_nor$nor_alt[i])) one[i]=T else one[i]=F
  if(ref_north_nor$ref_ref[i]==opp(ref_north_nor$nor_alt[i]) & ref_north_nor$ref_alt[i]==opp(ref_north_nor$nor_ref[i])) two[i]=T else two[i]=F
  if(one[i] | two[i]) flip=c(flip,i)                     #If either flip condition is true, place in flip vector
}

flip2=ref_north_nor[flip,]                     #Selecting rows to flip
keep=unique(c(match,minor,flip))               #Identifying rows to keep. Rows that match alleles, major and minor are flipped, and those that need flippping
keep=keep[order(keep)]                         #Ordering keep vector
ref_north_nor=ref_north_nor[keep,]             #Pruning dataset for all that are kept

#############################################################################
#### Writing out SNPs to change name of, SNPs to flip, and SNPs to remove####
#############################################################################
write.table(ref_north_nor[which(ref_north_nor$SNP_nor!=ref_north_nor$SNP_ref_north),c(1,3)],paste0("/Volumes/qac/prj/sjogrens/Imputation/Norwegian/chr",chr,"/nor_upSNP",chr,".txt"),col.names=F,row.names=F,quote=F)
write.table(flip2$SNP_nor,paste0("/Volumes/qac/prj/sjogrens/Imputation/Norwegian/chr",chr,"/nor_snps_to_flip",chr,".txt"),col.names=F,row.names=F,quote=F)
write.table(ref_north_nor$SNP_nor,paste0("/Volumes/qac/prj/sjogrens/Imputation/Norwegian/chr",chr,"/nor_extract",chr,".txt"),col.names=F,row.names=F,quote=F)

#### Non-Norwegian ####
at=which((ref_south_non$non_ref=="A" & ref_south_non$non_alt=="T") | (ref_south_non$ref_ref=="T" & ref_south_non$ref_alt=="A"))   #
ta=which((ref_south_non$non_ref=="T" & ref_south_non$non_alt=="A") | (ref_south_non$ref_ref=="A" & ref_south_non$ref_alt=="T"))
cg=which((ref_south_non$non_ref=="C" & ref_south_non$non_alt=="G") | (ref_south_non$ref_ref=="G" & ref_south_non$ref_alt=="C"))
gc=which((ref_south_non$non_ref=="G" & ref_south_non$non_alt=="C") | (ref_south_non$ref_ref=="C" & ref_south_non$ref_alt=="G"))
remove_snps=unique(c(at,ta,cg,gc))                       #List SNPs where solution isn't apparent
remove_snps=remove_snps[order(remove_snps)]              #Ordering SNPs to remove
if(length(remove_snps)>0) ref_south_non=ref_south_non[-remove_snps,]  #Removing problem SNPs #Format note:No SNPs to remove, empty vector causes problems
match=which(ref_south_non$ref_ref==ref_south_non$non_ref & ref_south_non$ref_alt==ref_south_non$non_alt) #SNPs where alleles match
minor=which(ref_south_non$ref_ref==ref_south_non$non_alt & ref_south_non$ref_alt==ref_south_non$non_ref) #SNPs where major and minor alleles are switched
flip=vector()                                            #Creating empty vector for SNPs that should be flipped
one=rep(0,nrow(ref_south_non))
two=rep(0,nrow(ref_south_non))
for(i in 1:nrow(ref_south_non)){
  if(ref_south_non$ref_ref[i]==opp(ref_south_non$non_ref[i]) & ref_south_non$ref_alt[i]==opp(ref_south_non$non_alt[i])) one[i]=T else one[i]=F
  if(ref_south_non$ref_ref[i]==opp(ref_south_non$non_alt[i]) & ref_south_non$ref_alt[i]==opp(ref_south_non$non_ref[i])) two[i]=T else two[i]=F
  if(one[i] | two[i]) flip=c(flip,i)                     #If either flip condition is true, place in flip vector
}

flip2=ref_south_non[flip,]                     #Selecting rows to flip
keep=unique(c(match,minor,flip))               #Identifying rows to keep. Rows that match alleles, major and minor are flipped, and those that need flippping
keep=keep[order(keep)]                         #Ordering keep vector
ref_south_non=ref_south_non[keep,]             #Pruning dataset for all that are kept

#############################################################################
#### Writing out SNPs to change name of, SNPs to flip, and SNPs to remove####
#############################################################################
write.table(ref_south_non[which(ref_south_non$SNP_non!=ref_south_non$SNP_ref_south),c(1,3)],paste0("/Volumes/qac/prj/sjogrens/Imputation/NonNorwegian/chr",chr,"/non_nor_upSNP",chr,".txt"),col.names=F,row.names=F,quote=F)
write.table(flip2$SNP_non,paste0("/Volumes/qac/prj/sjogrens/Imputation/NonNorwegian/chr",chr,"/non_nor_snps_to_flip",chr,".txt"),col.names=F,row.names=F,quote=F)
write.table(ref_south_non$SNP_non,paste0("/Volumes/qac/prj/sjogrens/Imputation/NonNorwegian/chr",chr,"/non_nor_extract",chr,".txt"),col.names=F,row.names=F,quote=F)

#### SICCA ####
at=which((ref_full_sicca$sicca_ref=="A" & ref_full_sicca$sicca_alt=="T") | (ref_full_sicca$ref_ref=="T" & ref_full_sicca$ref_alt=="A"))   #
ta=which((ref_full_sicca$sicca_ref=="T" & ref_full_sicca$sicca_alt=="A") | (ref_full_sicca$ref_ref=="A" & ref_full_sicca$ref_alt=="T"))
cg=which((ref_full_sicca$sicca_ref=="C" & ref_full_sicca$sicca_alt=="G") | (ref_full_sicca$ref_ref=="G" & ref_full_sicca$ref_alt=="C"))
gc=which((ref_full_sicca$sicca_ref=="G" & ref_full_sicca$sicca_alt=="C") | (ref_full_sicca$ref_ref=="C" & ref_full_sicca$ref_alt=="G"))
remove_snps=unique(c(at,ta,cg,gc))                       #List SNPs where solution isn't apparent
remove_snps=remove_snps[order(remove_snps)]              #Ordering SNPs to remove
if(length(remove_snps)>0) ref_full_sicca=ref_full_sicca[-remove_snps,]               #Removing problem SNPs #Format note:No SNPs to remove, empty vector causes problems
match=which(ref_full_sicca$ref_ref==ref_full_sicca$sicca_ref & ref_full_sicca$ref_alt==ref_full_sicca$sicca_alt) #SNPs where alleles match
minor=which(ref_full_sicca$ref_ref==ref_full_sicca$sicca_alt & ref_full_sicca$ref_alt==ref_full_sicca$sicca_ref) #SNPs where major and minor alleles are switched
flip=vector()                                            #Creating empty vector for SNPs that should be flipped
one=rep(0,nrow(ref_full_sicca))
two=rep(0,nrow(ref_full_sicca))
for(i in 1:nrow(ref_full_sicca)){
  if(ref_full_sicca$ref_ref[i]==opp(ref_full_sicca$sicca_ref[i]) & ref_full_sicca$ref_alt[i]==opp(ref_full_sicca$sicca_alt[i])) one[i]=T else one[i]=F
  if(ref_full_sicca$ref_ref[i]==opp(ref_full_sicca$sicca_alt[i]) & ref_full_sicca$ref_alt[i]==opp(ref_full_sicca$sicca_ref[i])) two[i]=T else two[i]=F
  if(one[i] | two[i]) flip=c(flip,i)                     #If either flip condition is true, place in flip vector
}

flip2=ref_full_sicca[flip,]                     #Selecting rows to flip
keep=unique(c(match,minor,flip))               #Identifying rows to keep. Rows that match alleles, major and minor are flipped, and those that need flippping
keep=keep[order(keep)]                         #Ordering keep vector
ref_full_sicca=ref_full_sicca[keep,]             #Pruning dataset for all that are kept

#############################################################################
#### Writing out SNPs to change name of, SNPs to flip, and SNPs to remove####
#############################################################################
write.table(ref_full_sicca[which(ref_full_sicca$SNP_sicca!=ref_full_sicca$SNP_ref),c(1,3)],paste0("/Volumes/qac/prj/sjogrens/Imputation/SICCA/chr",chr,"/sicca_upSNP",chr,".txt"),col.names=F,row.names=F,quote=F)
write.table(flip2$SNP_sicca,paste0("/Volumes/qac/prj/sjogrens/Imputation/SICCA/chr",chr,"/sicca_snps_to_flip",chr,".txt"),col.names=F,row.names=F,quote=F)
write.table(ref_full_sicca$SNP_sicca,paste0("/Volumes/qac/prj/sjogrens/Imputation/SICCA/chr",chr,"/sicca_extract",chr,".txt"),col.names=F,row.names=F,quote=F)

#### Phase1 ####
at=which((ref_full_p1$p1_ref=="A" & ref_full_p1$p1_alt=="T") | (ref_full_p1$ref_ref=="T" & ref_full_p1$ref_alt=="A"))   #
ta=which((ref_full_p1$p1_ref=="T" & ref_full_p1$p1_alt=="A") | (ref_full_p1$ref_ref=="A" & ref_full_p1$ref_alt=="T"))
cg=which((ref_full_p1$p1_ref=="C" & ref_full_p1$p1_alt=="G") | (ref_full_p1$ref_ref=="G" & ref_full_p1$ref_alt=="C"))
gc=which((ref_full_p1$p1_ref=="G" & ref_full_p1$p1_alt=="C") | (ref_full_p1$ref_ref=="C" & ref_full_p1$ref_alt=="G"))
remove_snps=unique(c(at,ta,cg,gc))                       #List SNPs where solution isn't apparent
remove_snps=remove_snps[order(remove_snps)]              #Ordering SNPs to remove
if(length(remove_snps)>0) ref_full_p1 = ref_full_p1[-remove_snps,]  #Removing problem SNPs #Format note:No SNPs to remove, empty vector causes problems
match=which(ref_full_p1$ref_ref==ref_full_p1$p1_ref & ref_full_p1$ref_alt==ref_full_p1$p1_alt) #SNPs where alleles match
minor=which(ref_full_p1$ref_ref==ref_full_p1$p1_alt & ref_full_p1$ref_alt==ref_full_p1$p1_ref) #SNPs where major and minor alleles are switched
flip=vector()                                            #Creating empty vector for SNPs that should be flipped
one=rep(0,nrow(ref_full_p1))
two=rep(0,nrow(ref_full_p1))
for(i in 1:nrow(ref_full_p1)){
  if(ref_full_p1$ref_ref[i]==opp(ref_full_p1$p1_ref[i]) & ref_full_p1$ref_alt[i]==opp(ref_full_p1$p1_alt[i])) one[i]=T else one[i]=F
  if(ref_full_p1$ref_ref[i]==opp(ref_full_p1$p1_alt[i]) & ref_full_p1$ref_alt[i]==opp(ref_full_p1$p1_ref[i])) two[i]=T else two[i]=F
  if(one[i] | two[i]) flip=c(flip,i)                     #If either flip condition is true, place in flip vector
}

flip2=ref_full_p1[flip,]                     #Selecting rows to flip
keep=unique(c(match,minor,flip))               #Identifying rows to keep. Rows that match alleles, major and minor are flipped, and those that need flippping
keep=keep[order(keep)]                         #Ordering keep vector
ref_full_p1=ref_full_p1[keep,]             #Pruning dataset for all that are kept

#############################################################################
#### Writing out SNPs to change name of, SNPs to flip, and SNPs to remove####
#############################################################################
write.table(ref_full_p1[which(ref_full_p1$SNP_p1!=ref_full_p1$SNP_ref),c(1,3)],paste0("/Volumes/qac/prj/sjogrens/Imputation/Phase1/chr",chr,"/p1_upSNP",chr,".txt"),col.names=F,row.names=F,quote=F)
write.table(flip2$SNP_p1,paste0("/Volumes/qac/prj/sjogrens/Imputation/Phase1/chr",chr,"/p1_snps_to_flip",chr,".txt"),col.names=F,row.names=F,quote=F)
write.table(ref_full_p1$SNP_p1,paste0("/Volumes/qac/prj/sjogrens/Imputation/Phase1/chr",chr,"/p1_extract",chr,".txt"),col.names=F,row.names=F,quote=F)



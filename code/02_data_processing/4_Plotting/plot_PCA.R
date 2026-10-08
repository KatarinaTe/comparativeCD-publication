####load eigenvec - PCA - new 2024-12-02####
#Load file created in 3_Merging_filtering/01_gencove_axiom_merging_filtering.sh
pca <- read.csv("prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.eigenvec", sep= " ", 
                header=F,col.names =c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7",
                                      "PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16",
                                      "PC17","PC18","PC19","PC20") )
head(pca)
nrow(pca)
pca$IID <- as.character(pca$IID)

data6 <- read.csv("data6_2024-10-14.txt", sep= "\t", header=T)
data6$IID <- as.character(data6$IID)

fam_pc <- full_join(pca, data6,  by="IID")
nrow(fam_pc)
head(fam_pc)

eigenval <- read.csv("prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.eigenval", sep= " ", 
                header=F,col.names =c("eigenvals"))
head(eigenval)
                                      
#####Get breeds for ALL FAM dogs#####
dog <- read.csv("/path/to/DarwinsArk_20220715_dogs.csv", sep= ",")
head(dog)
nrow(dog)
dog$IID <- as.character(dog$dog)
length(unique(dog$dog))

#subset(dog, dog$IID=158)


fam_breed <- left_join(fam_pc, dog,  by="IID")
head(fam_breed)
table(fam_breed$breed1_inputted, exclude=NULL)
unique(fam_breed$breed1_inputted, exclude=NULL)
table(fam_breed$purebred,fam_breed$sex.y, exclude=NULL)
#208/1219
#221/1219

table(fam_breed$purebred, fam_breed$breed1_inputted, exclude=NULL) 
nrow(fam_breed)

not_purebred <- subset(fam_breed, fam_breed$purebred == "no")
table(not_purebred$breed1_inputted, exclude=NULL) #NA=525
nrow(not_purebred)

purebred <- subset(fam_breed, fam_breed$purebred == "yes")
unique(purebred$breed1_inputted, exclude=NULL) #105 different breeds
nrow(purebred)
table(purebred$sex.y, exclude=NULL)


pch.vector = rep(1,nrow(fam_breed))
pch.vector[is.na(fam_breed$purebred) ] = 21
pch.vector[fam_breed$breed1_inputted == "Wolf"] = 8

plot(fam_breed$PC1, fam_breed$PC2,  pch=pch.vector, #col=col.vector,
     xlab="PCA1 (Variance = 37.3%)", ylab = "PCA2 (Variance = 24.7%)")
legend("topright", c( "other breed", "Labrador Retriever" ), 
       col=c("red", "blue"), pch=17, cex=0.8 )
legend("bottomright", c("purebred","not purebred", "purebred=NA", "wolf"), pch=c(17,2,21, 8), cex=0.8 )


####plot pca and breeds####
table(fam_breed$breed1_inputted, fam_breed$purebred, exclude=NULL)
#write.table(fam_breed, file= "fam_breed.txt", row.names =F, sep = ",", col.names= T, quote = F )

##plot breeds in PCAplot
col.vector = rep(1,nrow(fam_breed))
col.vector[fam_breed$purebred=="yes"] = "black"
col.vector[is.na(fam_breed$purebred) ] = "lightgrey"
col.vector[fam_breed$purebred=="no" ] = "lightgrey"
col.vector[fam_breed$breed1_inputted == "Labrador Retriever" & fam_breed$purebred=="yes"] = "red"
col.vector[fam_breed$breed1_inputted == "Golden Retriever"  & fam_breed$purebred=="yes" ] = "blue"
col.vector[fam_breed$breed1_inputted == "German Shepherd Dog"  & fam_breed$purebred=="yes"] = "purple"
col.vector[fam_breed$breed1_inputted == "Border Collie"  & fam_breed$purebred=="yes"] = "green"
col.vector[fam_breed$breed1_inputted == "Australian Shepherd"  & fam_breed$purebred=="yes"] = "yellow"
col.vector[fam_breed$breed1_inputted == "American Pit Bull Terrier"  & fam_breed$purebred=="yes"] = "orange"
col.vector[fam_breed$breed1_inputted == "Boxer"  & fam_breed$purebred=="yes"] = "lightblue3"
col.vector[fam_breed$breed1_inputted == "Siberian Husky"  & fam_breed$purebred=="yes"] = "darkred"
col.vector[fam_breed$breed1_inputted == "Chow Chow"  & fam_breed$purebred=="yes"] = "red3"
col.vector[fam_breed$breed1_inputted == "Shiba Inu"  & fam_breed$purebred=="yes"] = "green4"
col.vector[fam_breed$breed1_inputted == "Wolf"] = "black"
col.vector[fam_breed$breed1_inputted == "New Guinea Singing Dog"] = "yellow3"


pch.vector = rep(1,nrow(fam_breed))
pch.vector[fam_breed$purebred=="no"] = 2
pch.vector[fam_breed$purebred=="yes"] = 17
pch.vector[is.na(fam_breed$purebred) ] = 21
pch.vector[fam_breed$breed1_inputted == "Wolf"] = 8

plot(fam_breed$PC1, fam_breed$PC2, col=col.vector, pch=pch.vector, 
     xlab="PCA1 (Variance = 37.3%)", ylab = "PCA2 (Variance = 24.7%)")
legend("topright", c( "other 'purebreed' ", "Labrador Retriever","Golden Retriever","German Shepherd Dog", "Border Collie",
                      "Australian Shepherd","American Pit Bull Terrier", "Boxer", "Siberian Husky", "Chow Chow","Shiba Inu",
                      "New Guinea Singing Dog" ), 
       col=c("black", "red", "blue", "purple","green", "yellow","orange", "lightblue3", "darkred","red3","green4","yellow3"), 
       pch=17, cex=0.8 )
legend("bottomright", c("purebred","not purebred", "purebred=NA", "wolf"), pch=c(17,2,21, 8), cex=0.8 )


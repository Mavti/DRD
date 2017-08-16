setwd("~/Desktop/DRD_2017_LAB/Report_input")
baseDir <- ("~/Desktop/DRD_2017_LAB/Report_input")
library(minfi)
list.files(baseDir)
targets <- read.metharray.sheet(baseDir) #contains phenotypic info
targets

#Create an object of class RGChannelSet
RGset <- read.metharray.exp(targets = targets,extended=T) 
save(RGset,file="RGset.RData")

Red <- data.frame(getRed(RGset))
dim(Red)

Green <- data.frame(getGreen(RGset))
dim(Green)


load('~/Desktop/DRD_2017_LAB/R_files/Illumina450Manifest.RData')

# Get the manifest object
?getManifest
getManifest(RGset)
getManifestInfo(RGset)

df <- data.frame(getProbeInfo(RGset))
dim(df) 
df_II <- data.frame(getProbeInfo(RGset, type = "II"))
dim(df_II)


### Extract methylated and unmethylated signals
MSet.raw <- preprocessRaw(RGset)

### QC plot
qc <- getQC(MSet.raw)
plotQC(qc)

### Control probes
df_TypeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
str(df_TypeControl)
df_TypeControl$Type <- factor(df_TypeControl$Type)
str(df_TypeControl)
controlStripPlot(RGset, controls="NEGATIVE")

library(shinyMethyl)
vignette("shinyMethyl")
summary <- shinySummarize(RGset)


### Detection p-value
?detectionP
detP <- detectionP(RGset)
str(detP)
dim(detP)
head(detP) #row names are the probes. if p is value is low high prob that signal is different from background.
?detectionP #works on a RGchannelset object

failed <- detP>0.01
head(failed) #dataframe with trues and flases.
dim(failed)
table(failed)
summary(failed) #i have number of true and flase per sample


### Samples cleaning: filter the ones that have too many failed probes
# Fraction of failed postions per sample
head(colMeans(failed)) #each value corresponds to the ratio (number of TRUE)/(number of FALSE) for each column (that is, for each sample)
means_of_columns <- colMeans(failed) 
means_of_columns

### Probes cleaning
# Fraction of samples with failed postions per probe
head(rowMeans(failed)) #each value corresponds to the ratio (number of TRUE)/(number of FALSE) for each row (that is, for each probe)
# cfr with head(failed)
head(failed)
# 6/8=0.75
# 7/8=0.875
#5/8= 0.625 | 4/8=0.5


#to remove the probes having > 1 % of samples with a detection p-value greater than 0.01
means_of_rows <- rowMeans(failed)
head(means_of_rows)
str(means_of_rows)
probes_to_retain <- means_of_rows<0.01
table(probes_to_retain)
names_probes_to_remove <- names(probes_to_retain)[probes_to_retain==FALSE]
names_probes_to_remove
MSet.raw_filtered <- MSet.raw[probes_to_retain,]
MSet.raw
MSet.raw_filtered

# probes we have removed are in the MSet.raw_filtered object???
Unmeth_filtered <- getUnmeth(MSet.raw_filtered)
Meth_filtered <- getMeth(MSet.raw_filtered)
intersect(rownames(Unmeth_filtered), names_probes_to_remove)
intersect(rownames(Meth_filtered), names_probes_to_remove)


library(wateRmelon)
wateRmelon_filtered=pfilter(RGset)
wateRmelon_filtered



### beta-values and M-values
head(Unmeth_filtered)
head(Meth_filtered)

beta <- getBeta(MSet.raw_filtered)
str(beta)
summary(beta)
head(beta)

M <- getM(MSet.raw_filtered)
str(M)
summary(M)
head(M)



### Density plots of Beta and M values
mean_of_beta <- apply(beta,1,mean)
head(beta)
head(mean_of_beta)
?density
d_mean_of_beta <- density(mean_of_beta,na.rm=T)
plot(d_mean_of_beta,main="Density of Beta Values",col="blue")

mean_of_M <- apply(M,1,mean)
head(M)
head(mean_of_M)
?density
d_mean_of_M <- density(mean_of_M,na.rm=T)
plot(d_mean_of_M,main="Density of M Values", col="red")


### Check the distributions of Type I and Type II probes
load('../R_files/Illumina450Manifest_clean.RData')
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)

beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
dim(beta_I)
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]
dim(beta_II)

mean_of_beta_I <- apply(beta_I,1,mean)
mean_of_beta_II <- apply(beta_II,1,mean)

d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)

plot(d_mean_of_beta_I,main="Density of Beta Values, Type I and Type II", col="blue")
lines(d_mean_of_beta_II,col="red")

M_I <- M[rownames(M) %in% dfI$IlmnID,]
dim(M_I)
M_II <- M[rownames(M) %in% dfII$IlmnID,]
dim(M_II)

mean_of_M_I <- apply(M_I,1,mean)
mean_of_M_II <- apply(M_II,1,mean)

d_mean_of_M_I <- density(mean_of_M_I,na.rm=T)
d_mean_of_M_II <- density(mean_of_M_II,na.rm=T)

plot(d_mean_of_M_II,main="Density of M Values, Type I and Type II", col="red",xlim=c(-6, 7))
lines(d_mean_of_M_I,col="blue" )

#### standard deviations
sd_of_beta_I <- apply(beta_I,1,sd)
sd_of_beta_II <- apply(beta_II,1,sd)
d_sd_of_beta_I <- density(sd_of_beta_I,na.rm=T)
d_sd_of_beta_II <- density(sd_of_beta_II,na.rm=T)

plot(d_sd_of_beta_I,col="blue",main="raw sd")
lines(d_sd_of_beta_II,col="red")
boxplot(beta, main="Boxplots of beta values")

library(wateRmelon)

# betaqn
?betaqn
beta_betaqn <- betaqn(MSet.raw)
beta_betaqn_I <- beta_betaqn[rownames(beta_betaqn) %in% dfI$IlmnID,]
beta_betaqn_II <- beta_betaqn[rownames(beta_betaqn) %in% dfII$IlmnID,]
mean_of_beta_betaqn_I <- apply(beta_betaqn_I,1,mean)
mean_of_beta_betaqn_II <- apply(beta_betaqn_II,1,mean)
d_mean_of_beta_betaqn_I <- density(mean_of_beta_betaqn_I,na.rm=T)
d_mean_of_beta_betaqn_II <- density(mean_of_beta_betaqn_II,na.rm=T)
sd_of_beta_betaqn_I <- apply(beta_betaqn_I,1,sd)
sd_of_beta_betaqn_II <- apply(beta_betaqn_II,1,sd)
d_sd_of_beta_betaqn_I <- density(sd_of_beta_betaqn_I,na.rm=T)
d_sd_of_beta_betaqn_II <- density(sd_of_beta_betaqn_II,na.rm=T)
pdf("betaqn.pdf",width=10,height=5)
par(mfrow=c(1,3))
plot(d_mean_of_beta_betaqn_I,col="blue",main="betaqn beta")
lines(d_mean_of_beta_betaqn_II,col="red")
plot(d_sd_of_beta_betaqn_I,col="blue",main="betaqn sd")
lines(d_sd_of_beta_betaqn_II,col="red")
boxplot(beta_betaqn, main=" betqn boxplot")
dev.off() #shuts down the graphical device


# fuks
?fuks
beta_fuks <- fuks(MSet.raw)
beta_fuks_I <- beta_fuks[rownames(beta_fuks) %in% dfI$IlmnID,]
beta_fuks_II <- beta_fuks[rownames(beta_fuks) %in% dfII$IlmnID,]
mean_of_beta_fuks_I <- apply(beta_fuks_I,1,mean)
mean_of_beta_fuks_II <- apply(beta_fuks_II,1,mean)
d_mean_of_beta_fuks_I <- density(mean_of_beta_fuks_I,na.rm=T)
d_mean_of_beta_fuks_II <- density(mean_of_beta_fuks_II,na.rm=T)
sd_of_beta_fuks_I <- apply(beta_fuks_I,1,sd)
sd_of_beta_fuks_II <- apply(beta_fuks_II,1,sd)
d_sd_of_beta_fuks_I <- density(sd_of_beta_fuks_I,na.rm=T)
d_sd_of_beta_fuks_II <- density(sd_of_beta_fuks_II,na.rm=T)
pdf("fuks.pdf",width=10,height=5)
par(mfrow=c(1,3))
plot(d_mean_of_beta_fuks_I,col="blue",main="fuks beta")
lines(d_mean_of_beta_fuks_II,col="red")
plot(d_sd_of_beta_fuks_I,col="blue",main="fuks sd")
lines(d_sd_of_beta_fuks_II,col="red")
boxplot(beta_fuks, main="boxplot")
dev.off() #shuts down the graphical device

# dasen
?dasen
beta_dasen <- dasen(MSet.raw)
beta_dasen <- getBeta(beta_dasen)
beta_dasen <- na.omit(beta_dasen)
beta_dasen_I <- beta_dasen[rownames(beta_dasen) %in% dfI$IlmnID,]
beta_dasen_II <- beta_dasen[rownames(beta_dasen) %in% dfII$IlmnID,]
mean_of_beta_dasen_I <- apply(beta_dasen_I,1,mean)
mean_of_beta_dasen_II <- apply(beta_dasen_II,1,mean)
d_mean_of_beta_dasen_I <- density(mean_of_beta_dasen_I,na.rm=T)
d_mean_of_beta_dasen_II <- density(mean_of_beta_dasen_II,na.rm=T)
sd_of_beta_dasen_I <- apply(beta_dasen_I,1,sd)
sd_of_beta_dasen_II <- apply(beta_dasen_II,1,sd)
d_sd_of_beta_dasen_I <- density(sd_of_beta_dasen_I,na.rm=T)
d_sd_of_beta_dasen_II <- density(sd_of_beta_dasen_II,na.rm=T)
pdf("dasen.pdf",width=10,height=5)
par(mfrow=c(1,3))
plot(d_mean_of_beta_dasen_I,col="blue",main="dasen beta")
lines(d_mean_of_beta_dasen_II,col="red")
plot(d_sd_of_beta_dasen_I,col="blue",main="dasen sd")
lines(d_sd_of_beta_dasen_II,col="red")
boxplot(beta_dasen, main="boxplot")
dev.off() #shuts down the graphical device

# preprocessQuantile
?preprocessQuantile
preprocessQuantile_results <- preprocessQuantile(RGset)
str(preprocessQuantile_results)
class(preprocessQuantile_results)
preprocessQuantile_results
beta_preprocessQuantile <- getBeta(preprocessQuantile_results)
head(beta_preprocessQuantile)
preprocessQuantile_results <- preprocessQuantile(MSet.raw)
str(preprocessQuantile_results)
class(preprocessQuantile_results)
preprocessQuantile_results
beta_preprocessQuantile <- getBeta(preprocessQuantile_results)
head(beta_preprocessQuantile)
beta_preprocessQuantile_I <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfI$IlmnID,]
beta_preprocessQuantile_II <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfII$IlmnID,]
mean_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I,1,mean)
mean_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II,1,mean)
d_mean_of_beta_preprocessQuantile_I <- density(mean_of_beta_preprocessQuantile_I,na.rm=T)
d_mean_of_beta_preprocessQuantile_II <- density(mean_of_beta_preprocessQuantile_II,na.rm=T)
sd_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I,1,sd)
sd_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II,1,sd)
d_sd_of_beta_preprocessQuantile_I <- density(sd_of_beta_preprocessQuantile_I,na.rm=T)
d_sd_of_beta_preprocessQuantile_II <- density(sd_of_beta_preprocessQuantile_II,na.rm=T)
pdf("preprocessQuantile.pdf",width=10,height=5)
par(mfrow=c(1,3))
plot(d_mean_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile beta")
lines(d_mean_of_beta_preprocessQuantile_II,col="red")
plot(d_sd_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile sd")
lines(d_sd_of_beta_preprocessQuantile_II,col="red")
boxplot(beta_preprocessQuantile, main="boxplot")
dev.off() #shuts down the graphical device


### HOMOSCEDASTICITY and HETEROSCEDASTICITY

summary(beta_betaqn) #i have no 0 nor -inf, but i have missing values
beta_betaqn_noNA <- na.omit(beta_betaqn)

beta_betaqn_mean <- apply(beta_betaqn_noNA,1,mean)
beta_betaqn_stdev <- apply(beta_betaqn_noNA,1, sd)

M_betaqn <- Beta2M(beta_betaqn_noNA)
summary(M_betaqn) 
M_betaqn_mean <- apply(M_betaqn,1, mean)
M_betaqn_stdev <- apply(M_betaqn,1, sd)

par(mfrow=c(1,2))
smoothScatter(beta_betaqn_mean, beta_betaqn_stdev)
lines(lowess(beta_betaqn_mean, beta_betaqn_stdev), col="red") # lowess carries out a locally weighted regression of y on x
smoothScatter(M_betaqn_mean, M_betaqn_stdev)
lines(lowess(M_betaqn_mean, M_betaqn_stdev), col="red") # lowess line (x,y)


###IDENTIFICATION OF DMPs =Differentially Methilated Positions
#==> ALL PROBES ####

### ANOVA ###
#function to perform Anova on all probes according to groups and BMI
MYlmFunction <- function(x) {
  anova_test <- aov(x~ pheno$Group+pheno$BMI)
  return(summary(anova_test)[[1]][[5]][1])
}

pValuesAnova <- apply(M_betaqn,1, MYlmFunction)
# We can create a data.frame with all the M values and the pValue column
final_Anova <- data.frame(M_betaqn, pValuesAnova)
head(final_Anova)
str(final_Anova)
# We can order the probes on the basis of pValuesAnova (from the smallest to the larget value)
final_Anova <- final_Anova[order(final_Anova$pValuesAnova),] #the comma means that i want to order the rows, it's like a subset
head(final_Anova)
dim(final_Anova)
# How many probes have a pValue<=0.05?
final_Anova_0.05 <- final_Anova[final_Anova$pValuesAnova<=0.05,] #remember that when you add a vector to a dataframe, it automatically add as col name, the name of the vector! in this case pValuesAnova
dim(final_Anova_0.05)
head(final_Anova_0.05)
#49598

### WILCOXON ###
str(wilcox)
MYWilcoxFunction <- function(x) {
  wilcox <- wilcox.test(x~ pheno$Group)
  return(wilcox[[3]])
} 

pValuesWilcox <- apply(M_betaqn,1, MYWilcoxFunction)
final_Wilcox <- data.frame(M_betaqn, pValuesWilcox)
head(final_Wilcox)
final_Wilcox <- final_Wilcox[order(pValuesWilcox),]
final_Wilcox_0.05 <- final_Wilcox[final_Wilcox$pValuesWilcox<=0.05,]
dim(final_Wilcox_0.05)
#34367

### how many probes are differentially methylated according to both anova and wilcoxon?
intersection <- intersect(rownames(final_Anova_0.05),rownames(final_Wilcox_0.05))
length(intersection)
#30044

############################################
# MULTIPLE TEST CORRECTION
############################################
final_Anova <- final_Anova[order(final_Anova$pValuesAnova),]
head(final_Anova)
raw_pValues <- final_Anova[,9]
?p.adjust
corrected_pValues_BH <- p.adjust(raw_pValues,"BH")
corrected_pValues_Bonf <- p.adjust(raw_pValues,"bonferroni")
final_Anova_corrected <- data.frame(final_Anova, corrected_pValues_BH, corrected_pValues_Bonf) #attatch corrections to original dataframe
head(final_Anova_corrected)

pval_bonf<- final_Anova_corrected[final_Anova_corrected$corrected_pValues_Bonf<=0.05,]
pval_bh<- final_Anova_corrected[final_Anova_corrected$corrected_pValues_BH<=0.05,]
dim(pval_bh)
dim(pval_bonf)

#bonf made huge correction! values equal to 1
#pdf('p_value_corrected_anova_cov')
boxplot(final_Anova_corrected[,9:11], main="Multiple Test Correction")
#dev.off()

### ANNOTATION ###
##################

final_Anova_corrected <- data.frame(rownames(final_Anova_corrected), final_Anova_corrected)
head(final_Anova_corrected)
colnames(final_Anova_corrected)[1] <- "IlmnID"
head(final_Anova_corrected)

# merge the two dataframes
final_TOT_Anova_corrected_annotated <- merge(final_Anova_corrected, Illumina450Manifest_clean,by="IlmnID")
dim(final_TOT_Anova_corrected_annotated)
head(final_TOT_Anova_corrected_annotated)
colnames(final_TOT_Anova_corrected_annotated)
#discard the design and keep only the annotation part
final_TOT_Anova_corrected_annotated <- final_TOT_Anova_corrected_annotated[,c(1:12,23:43)]
str(final_TOT_Anova_corrected_annotated)
final_TOT_Anova_corrected_annotated <- droplevels(final_TOT_Anova_corrected_annotated)
str(final_TOT_Anova_corrected_annotated)

############################################
# VOLCANO PLOTS
############################################
# VISUALIZE DNA METHYLATION
beta_values <- M2Beta(final_TOT_Anova_corrected_annotated[,2:9])
head(beta_values)
final_TOT_Anova_corrected_annotated <- data.frame(final_TOT_Anova_corrected_annotated[,1], beta_values, final_TOT_Anova_corrected_annotated[,10:33])
head(final_TOT_Anova_corrected_annotated)
colnames(final_TOT_Anova_corrected_annotated) [1] <- "IlmnID"

#Volcano plot
beta_groupA <- beta_values[,pheno$Group=="A"]
mean_beta_groupA <- apply(beta_groupA,1,mean)
beta_groupB <- beta_values[,pheno$Group=="B"]
mean_beta_groupB <- apply(beta_groupB,1,mean)
delta <- mean_beta_groupB-mean_beta_groupA
head(delta)

toVolcPlot <- data.frame(delta, -log10(final_TOT_Anova_corrected_annotated$pValuesAnova))
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16, main="Volcano Plot")
# Add the threshold of significan pvalue
abline(a=-log10(0.01),b=0,col="red")
toHighlight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.01)),]
head(toHighlight)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col='yellow')

# MANHATTAN PLOTS
############################################
library(gap)
#we need a 3 col table:
db=data.frame(final_TOT_Anova_corrected_annotated$CHR, final_TOT_Anova_corrected_annotated$MAPINFO, final_TOT_Anova_corrected_annotated$pValuesAnova)

dim(db)
db <- droplevels(db)
str(db)
palette <- rainbow(24)
palette
mhtplot(db,control=mht.control(colors=palette))

palette <- c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green")
mhtplot(db,control=mht.control(colors=palette))
axis(2,cex=0.5)
abline(a=-log10(0.05),b=0)

############################################
# PCA and MDS
############################################

?prcomp
#transpose col with be val, scale=T means that val should be scaled.. by default is false
pca <- prcomp(t(final_TOT_Anova_corrected_annotated[,2:9]),scale=T)
print(summary(pca)) # print variance accounted for each component
str(pca)
pca$x
plot(pca$x[,1], pca$x[,2])
plot(pca$x[,1], pca$x[,2],cex=1) #cex increase dimension of dots
text(pca$x[,1], pca$x[,2],labels=rownames(pca$x),pos = 3, cex = 0.5, col = "blue") #text function adds labels

### MDS 
final_TOT_Anova_corrected_annotated_sign <- final_TOT_Anova_corrected_annotated[final_TOT_Anova_corrected_annotated$pValuesAnova <0.01,]
dim(final_TOT_Anova_corrected_annotated_sign)
final_TOT_Anova_corrected_annotated_sign <- droplevels(final_TOT_Anova_corrected_annotated_sign)

#euclidean distaces
dist.euc=dist(t(final_TOT_Anova_corrected_annotated_sign[,2:9])) # euclidean distances between the rows
mds=cmdscale(dist.euc)
mds
plot(mds,type="n", main="MDS PLOT")
points(mds[pheno$Group=="A",],col="red",cex=2)
points(mds[pheno$Group=="B",],col="orange",cex=2)
legend("topright",c("A","B"),pch=1,col=c("red","orange"))


############################################
# HEATMAP
############################################
library(gplots)

# Heatmap on top 100 most significant CpG probes
final_TOT_Anova_corrected_annotated_sign <- final_TOT_Anova_corrected_annotated_sign[order(final_TOT_Anova_corrected_annotated_sign$pValuesAnova),]

matrix=as.matrix(final_TOT_Anova_corrected_annotated_sign[1:100,2:9])
pheno$Group
colorbar <- c("green","green","orange","orange","green","green","orange","orange")

# Let's compare the 3 linkage methods (complete, single, average)

?rainbow

# Complete (default options)
heatmap.2(matrix,col=terrain.colors(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F, main="Complete linkage")


# Single
heatmap.2(main="Single Linkage", matrix,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'single'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)


# Average
heatmap.2(main="Average linkage", matrix,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)






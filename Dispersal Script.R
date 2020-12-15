library(dartR)
library(SNPRelate)
library(sm)

###read in DArT file to Genlight object
Lepdatanests <- gl.read.dart(filename="Report_DLepo19-4747_1_moreOrders_SNP_1.csv", ind.metafile="dartmetadatanests.csv")
AridOnly <- gl.read.dart(filename="ARID ONLY SNPS.csv", ind.metafile="aridmetadatanestsSNPSEX.csv")

###read in Genlight
Lepdatanests <- readRDS("Lepdatanests.Rdata")
AridUnfiltered <- readRDS("AridOnlyUnfiltered.Rdata")

###convert genlight to genind/bayescan
ARgi <- gl2gi(AR, v = 1)
genomic_converter(lepfilter, output = "bayescan", filename = "lepfilteredbayescan.txt")

###filter out monomorphs
aridfilter <- gl.filter.monomorphs(AridOnly)

###filter out SNPs with call rate <0.70
aridfilter <- gl.filter.callrate(aridfilter, threshold = 0.70)

#merge pops to calculate heterozygosit
aridfilterallpops <- gl.merge.pop(aridfilter, old = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","19","NA"), new="ALL")
gl.report.heterozygosity(aridfilterallpops)

###Below not run
###filter out MAFs <0.40
#aridfilter <- gl.filter.maf(aridfilter, threshold = 0.40)

###filter out secondaries
#aridfilter <-gl.filter.secondaries(aridfilter)

###filter based on repeatability
#lepfilter <- gl.filter.repavg(lepfilter, threshold = 0.9)

###SNPRelate Identity-By-Descent Analysis(KING method of moment)
ibdking <- snpgdsIBDKING(ARgds, sample.id=NULL, snp.id=NULL, autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, type=c("KING-robust", "KING-homo"), family.id=NULL, num.thread=1L, useMatrix=FALSE, verbose=TRUE)
dat <- snpgdsIBDSelection(ibdking , 1/32)
plot(dat$IBS0, dat$kinship, xlab="Proportion of Zero IBS",
     +      ylab="Estimated Kinship Coefficient (KING-robust)")
write.csv(dat, "AR Kinship Analysis.csv")

###Convert to COLONY file
library(radiator)
artidy <- tidy_genomic_data(ARfilter)
colony <- write_colony(artidy, sample.markers = 500)

###Construct network
links <- read.table("AR Kinship Analysis CutSNPSEXAdultCaptures.csv", sep=',', dec='.', header=T)
nodes <- read.table("aridmetadatanestsSNPSEXAdultCaptures.csv", sep=',', dec='.', header=T)
network <- graph_from_data_frame(d=links,vertices=nodes,directed = F)
c20 <- c("dodgerblue2", "#E31A1C","green4","#6A3D9A", "#FF7F00", "orchid1", "gold1", "skyblue2", "#FB9A99",  "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "deeppink1", "blue1", "darkturquoise", "green1", "brown")
shapes <- c("circle","square")
size <- c(7,7)
my_color <- c20[as.numeric(as.factor(V(network)$pop))]
my_shapes <- shapes[as.numeric(as.factor(V(network)$sex))]
my_size <- size[as.numeric(as.factor(V(network)$sex))]
layout <- layout_with_graphopt(network)
plot(network, layout=layout.fruchterman.reingold,vertex.color=my_color, vertex.shape=my_shapes, vertex.size=my_size ,vertex.label.color="black",
     vertex.label.cex=0.6, vertex.label.dist=0,edge.width=E(network)$kinship*20, main='Stick Nest Rat Relatedness',cex=1.5)

###Statistical test for relatedness within/between nests
kinship <- read.csv("AR Kinship Analysis CutSNPSEXAdultCaptures.csv", header = T)
Females <- kinship[ which(kinship$Sex== 'FF'),]
Femaleswithin <- Females[ which(Females$WB== 'W'),]
Femalesbetween <- Females[ which(Females$WB== 'B'),]
Males <- kinship[ which(kinship$Sex== 'MM'),]
Maleswithin <- Males[ which(Males$WB== 'W'),]
Malesbetween <- Males[ which(Males$WB== 'B'),]
MaleFemale <- kinship[ which(kinship$Sex== 'MF'),]
MFwithin <- MaleFemale[ which(MaleFemale$WB== 'W'),]
MFbetween <- MaleFemale[ which(MaleFemale$WB== 'B'),]
#test for normal distribution
shapiro.test(Females$kinship)
shapiro.test(Males$kinship)
#distribution is not normal (p<0.05) so we will use Wilcoxon test (non parametric)
wilcoxtestF <- wilcox.test(Femaleswithin$kinship, Femalesbetween$kinship)
wilcoxtestM <- wilcox.test(Maleswithin$kinship, Malesbetween$kinship)
wilcoxtestMF <- wilcox.test(Maleswithin$kinship, Femaleswithin$kinship)
#visualise with plots
bpF <- boxplot(Femalesbetween$kinship, Femaleswithin$kinship, main = "Females", ylab = "Kinship Coefficient", names = c("Between Nests", "Within Nests"), cex.main=2, cex.lab=1.5, cex.axis=1.5)
bpM <- boxplot(Malesbetween$kinship, Maleswithin$kinship, main = "Males", ylab = "Kinship Coefficient", names = c("Between Nests", "Within Nests"), cex.main=2, cex.lab=1.5, cex.axis=1.5)
bpMF <- boxplot(Maleswithin$kinship, Femaleswithin$kinship, main = "Male vs Female Within Nest Relatedness", ylab = "Kinship Coefficient", names = c("Male", "Female"), cex.main=2, cex.lab=1.5, cex.axis=1.5)
bpMF2 <- boxplot(MFwithin$kinship, MFbetween$kinship, main = "Male-Female Nest Relatedness", ylab = "Kinship Coefficient", names = c("Within Nests", "Between Nests"), cex.main=2, cex.lab=1.5, cex.axis=1.5)
d <- density(Femalesbetween$kinship)
plot(d)
Femkinship <- Females$kinship[-c(24,33,34)] #NAs Removed
WB <- Females$WB
#remove NAs
WB <- WB[-c(24,33,34)]
WB <- as.factor(WB)
compare <- sm.density.compare(Femkinship, WB, xlab = "Kinship Coefficient")
title(main="Female Kinship Within and Between Nest Sites", cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
#compare weights at first capture - high kinship & large mass difference probably means parent/offspring
plot(kinship$kinship, kinship$Body.Mass.Diff)
text(labels=rownames(kinship), data=kinship)


##Violin plots
FemNAsRem <- Females[-c(24,33,34),]
p <- ggplot(FemNAsRem, aes(x=WB, y=kinship, fill=WB)) + geom_violin() + labs(title="Females", x="", y="Kinship Coefficient", names = c("Between Nests", "Within Nests"))
p + geom_boxplot(width=0.1, fill="white") + labs(fill="") + theme_classic() + theme(plot.title = element_text(size = 30), axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 20), axis.text.x = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size = 20)) + scale_fill_discrete(labels=c("Between Nests", "Within Nests"))

MaleNAsRem <- Males[-c(14,15,19,20),]
p <- ggplot(MaleNAsRem, aes(x=WB, y=kinship, fill=WB)) + geom_violin() + labs(title="Males", x="", y="Kinship Coefficient", names = c("Between Nests", "Within Nests"))
p + geom_boxplot(width=0.1, fill="white") + labs(fill="") + theme_classic() + theme(plot.title = element_text(size = 30), axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 20), axis.text.x = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size = 20)) + scale_fill_discrete(labels=c("Between Nests", "Within Nests"))

malevsfemale <- read.csv("malevsfemale.csv")
p <- ggplot(malevsfemale, aes(x=Sex, y=Kinship.Coefficient, fill=Sex)) + geom_violin() + labs(title="Female vs Male Cohabiting Pairings", x="", y="Kinship Coefficient", names = c("Female", "Male"))
p + geom_boxplot(width=0.1, fill="white") + labs(fill="") + theme_classic() + theme(plot.title = element_text(size = 30), axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), legend.text = element_text(size = 20)) + scale_fill_discrete(labels=c("Female", "Male")) + scale_x_discrete(labels=c("Cohabiting Females", "Cohabiting Males"))

#For poster (larger text)
p <- ggplot(malevsfemale, aes(x=Sex, y=Kinship.Coefficient, fill=Sex)) + geom_violin() + labs(title="", x="", y="Kinship Coefficient", names = c("Female", "Male"))
p + geom_boxplot(width=0.1, fill="white") + labs(fill="") + theme_classic() + theme(plot.title = element_text(size = 1), axis.title.y = element_text(size = 35), axis.text.y = element_text(size = 25), axis.text.x = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size = 35)) + scale_fill_discrete(labels=c("Female", "Male"))

##Descriptive stats
Summarize(Maleswithin$kinship)
Summarize(Malesbetween$kinship)
Summarize(Femaleswithin$kinship)
Summarize(Femalesbetween$kinship)
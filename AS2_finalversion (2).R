########################################################################################
##Notes from KH:

##This code encapsulates all steps post sample PAO1 alignment, including options for DEseq2 normalizing and feature counts.
## that is, aside from a couple steps which are finished in a data editor like Excel, which I wrote notes about.
##From what I have discovered, to get the most similar counts to both me and to the Lewin et. al paper-
##use the referenced annotation and tigrfam files from supplementary. 

##You may not be processing all the data at once.
##I have included 2 "ladder" options.

#The first starts at 221 and allows you to read in a raw count, combined feature counts for all relevant samples (including Clinical), which you can then normalize.
#The second starts at 235 and allows you to read in individually completed raw feature counts for each group of samples, which you can then normalize.

#These steps do not streamline, but they should make it easier to work through the code in chunks.
#That being said, anytime you download and read something back in, there is room for error.
#ALWAYS CHECK YOUR ENVIRONMENT
#####Look for too few columns or rows, or names that don't match or aren't formatted correctly
#####Check that there are column names (not just the first row as column names)
#It is tedious, but keeping track of the environment at each step makes troubleshooting much easier

###If you still complete everything, but your values don't match exactly, make sure you used the provided alignment files in Bowtie prior to this.
###Below is a list of package versions

#R version 4.2.2 (2022-10-31 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 22631)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
#[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
#[5] LC_TIME=English_United States.utf8    
#
#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] scico_1.5.0                 zeallot_0.1.0               cowplot_1.1.1              
#[4] lubridate_1.9.2             forcats_1.0.0               stringr_1.5.0              
#[7] dplyr_1.1.1                 purrr_1.0.1                 readr_2.1.4                
#[10] tidyr_1.3.0                 tibble_3.2.1                tidyverse_2.0.0            
#[13] DESeq2_1.38.3               SummarizedExperiment_1.28.0 Biobase_2.58.0             
#[16] MatrixGenerics_1.10.0       matrixStats_0.63.0          GenomicRanges_1.50.2       
#[19] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2           
#[22] BiocGenerics_0.44.0         ggsunburst_0.3.0            ggplot2_3.4.1              
#[25] reticulate_1.28             devtools_2.4.5              usethis_2.1.6              

        
#####
###Script is updated to include more code and comments
#script calculates AS2 to compare a model to human transcriptome and then graphs output using sunburst plots
#   1. Calculates AS2 or zscores for all samples in a model and plots sunburst graph
#   2. Calculates AS2 and zscores after subsampling to keep input numbers consistent across different conditions
#   3. Calculates zscores (penalties) for each individual replicate of a model

##########################################################################################
#set parameters for all analyses
#also need to run functions at the end of the script


#set working directory
wd = "C:/Users/Gina/Documents/PA/"
setwd(wd)

library(devtools)
install_github("didacs/ggsunburst")
library(ggsunburst)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(zeallot)
library(ggsunburst)
library(scico)


################################################################################
#run these functions first
#functions used above, run each
score_target_vs_model_PA <- function(counts_DF, Target_samplenames, Model_samplenames) # given a counts DF, lists of target samples, and lists of model samples, will return a score for each gene
{
  
  all_tested_samples_list <- names(counts_DF)[-1]
  names(counts_DF)
  
  conditions_list <- all_tested_samples_list
  dim(counts_DF %>% select(locus_tag, all_tested_samples_list) %>% column_to_rownames("locus_tag"))
  
  df1 <- counts_DF
  setdiff(Target_samplenames, names(df1))
  DF_target <- df1  %>% select(locus_tag, Target_samplenames) %>% gather(key=sample_name, value=expression, -locus_tag) %>% group_by(locus_tag) %>% dplyr::summarize(target_mean = mean(expression), target_SD = sd(expression))
  
  df_allzscores <- DF_target %>% inner_join(df1 %>% select(locus_tag, Model_samplenames)) %>% mutate_at(.vars=vars(-locus_tag, -target_mean, -target_SD), .funs=funs((.-target_mean)/target_SD))
  mean_modelZscoreDF <- df_allzscores %>% select(locus_tag, Model_samplenames)  %>% transmute(locus_tag, penalty_temp = pmap_dbl(.[c(-1)], function(...)  mean(c(...)))) %>% mutate(penalty = round(penalty_temp, digits = 4)) %>% select(locus_tag, penalty) 
  
  # below should get rid of NAs
  
  if(length(which(is.na(mean_modelZscoreDF$penalty)))>0)
  {
    mean_modelZscoreDF[which(is.na(mean_modelZscoreDF$penalty)),]$penalty <- 0
  }  
  final_DF <- mean_modelZscoreDF  
  
  return(final_DF)
}

draw_sunburst_STANDARD_PA <- function(PA_df_counts, TIGR_PA, sb, target_list, model_list, maxlim, node_labs)
{
  
  TIGR_PA <- TIGR_PA  %>% mutate(sub_role = if_else(sub_role == "Other", paste0(sub_role, "@", main_role), sub_role))
  PA_DE <- score_target_vs_model_PA(PA_df_counts, target_list, model_list)
  
  subrole_penalty <- TIGR_PA %>% left_join(PA_DE)  %>% group_by(sub_role) %>%  select(locus_tag, sub_role, penalty) %>% distinct %>% dplyr::summarize(penalty = mean(abs(penalty)<maxlim, na.rm=TRUE))
  main_role_penalty <- TIGR_PA %>% left_join(PA_DE) %>% group_by(main_role) %>% select(locus_tag, main_role, penalty) %>% distinct %>% dplyr::summarize(penalty = mean(abs(penalty)<maxlim, na.rm=TRUE))
  Meta_penalty <- TIGR_PA %>% left_join(PA_DE) %>% group_by(Meta) %>% select(locus_tag, Meta, penalty) %>% distinct %>% dplyr::summarize(penalty = mean(abs(penalty)<maxlim, na.rm=TRUE))
  subrole_penalty %>% arrange(penalty)
  sum(abs(PA_DE$penalty) < 2)
  TIGR_PA %>% left_join(PA_DE)  %>% arrange(abs(penalty)) #%>% filter(sub_role == "Sulfur metabolism")
  
  subrole_penalty %>% filter()
  
  names(Meta_penalty)[1] <- "role"
  names(subrole_penalty)[1] <- "role"
  names(main_role_penalty)[1] <- "role"
  
  all_penalty_levels <- rbind(Meta_penalty, main_role_penalty, subrole_penalty)
  all_penalty_levels$role <- all_penalty_levels$role %>% as.character
  
  sb$rects$name <- sb$rects$name %>% as.character
  temp1 <- sb$rects %>% left_join(all_penalty_levels, by = c("name"="role"))
  
  sb$rects$color <- temp1$penalty
  sb$rects$leaf <- TRUE
  
  
  tempinnerplot <- sunburst(sb, node_labels = F, node_labels.min = .2, node_labels.size = .8, leaf_labels = FALSE, rects.fill.aes = "color", rects.size  = 0.1) + scico::scale_fill_scico(palette = "romaO", limits=c(-.1,1))
  d=data.frame(x1=.5, x2=1, y1=.5, y2=1)
  backplot_1 <- ggplot() + geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill= sum(PA_DE$penalty %>% abs < maxlim, na.rm = TRUE)/length(PA_DE$penalty))) + scico::scale_fill_scico(palette = "romaO", limits=c(-.1,1)) + theme(legend.position="none")
  prow <- plot_grid(tempinnerplot,
                    backplot_1,
                    align = 'v',
                    labels = c("inner", "middle"),
                    ncol = 1
  )
  
  plot_grid(prow)
  p <- plot_grid(prow)
  
  return(list(p, PA_DE))
}
################################################################################

#the following is the original feature counts code that would be done in cmd line, however I find
#it easier to work in R, included here if you prefer to work in cmd.

#featureCounts(T 24 -a ~/LewinFiles/Pseudomonas_aeruginosa_PAO1_107_chr_ed_Dan.gff 
#-O -s 1 -g locus -t CDS -o /storage/home/hcoda1/4/glewin3/scratch/Model_Improvement/MW79_NB_2.fastq.mapped_toPAO1_107.05_10_2022.forward.O.count.txt /storage/home/hcoda1/4/glewin3/scratch/Model_Improvement/MW79_NB_2.fastq.mapped_toPAO1_107.05_10_2022.sam.gz

#if you are use bam files just alter the next line of code to reflect that.
files<-list.files("~/LewinFiles/singlestrand/samfiles", pattern="\\.gz$")
files

####featurecounts finished sam or bam files
##both the clinical data and Bomberger data is single strand, so I process both at the same time.
Clinical_Bomberger<-featureCounts(files=files, annot.ext ="~/LewinFiles/Pseudomonas_aeruginosa_PAO1_107_chr_ed_Dan.gff", isGTFAnnotationFile = T, allowMultiOverlap = T, GTF.featureType="CDS", GTF.attrType = "locus", strandSpecific = 1)
Clinical_Bomberger_counts<-as.data.frame(Clinical_Bomberger$counts)
Clinical_Bomberger_counts$LocusID<-row.names(Clinical_Bomberger_counts)


files<-list.files("~/LewinFiles/Pairedend?samfiles", pattern="\\.gz$")
files

HunterPAO1<-featureCounts(files=files, annot.ext ="~/LewinFiles/Pseudomonas_aeruginosa_PAO1_107_chr_ed_Dan.gff", isGTFAnnotationFile = T, allowMultiOverlap = T, GTF.featureType="CDS", GTF.attrType = "locus", strandSpecific = 0, isPairedEnd = T)
HunterPAO1counts<-as.data.frame(HunterPAO1$counts)
HunterPAO1counts$LocusID<-row.names(HunterPAO1counts)


CHBcombinedcore<-as.data.frame(CHBcombinedcore)
CHBcombinedcore<-subset(CHBcombinedcore,select=-7)
CHBcombinedcore$locus_tag<-row.names(CHBcombinedcore)
CHBcombinedcore2<-inner_join(HSP_13rawcounts,CHBcombinedcore, by="locus_tag")
CHBcombinedcore2<-as.data.frame(CHBcombinedcore2)
row.names(CHBcombinedcore2)<-CHBcombinedcore2$locus_tag
CHBcombinedcore2<-subset(CHBcombinedcore2, select=-c(1,50))

#Combine paired end and single strand reads
CHBcombined<-inner_join(Clinical_Bomberger_counts,HunterPAO1counts, by = "LocusID")
#relocate LocusID to front
CHBcombined<-relocate(CHBcombined, LocusID)
#Remove all but core genes
Lewin95genes<-read.csv("Lewin95genes.csv")
Lewin95genes<-Lewin95_genes %>% filter(Lewin95_genes$`CORE (in 95% of genomes or 95% of meta-transcriptomes; 5147 genes)` =="CORE")
CHBcombinedcore<-inner_join(CHBcombined,Lewin95genes, by="LocusID")
#check column names
colnames(CHBcombinedcore)
#drop core gene marker
CHBcombinedcore<-subset(CHBcombinedcore,select= -(51))
#check it is gone
colnames(CHBcombinedcore)
#change LocusID back to rownames
row.names(CHBcombinedcore)<-CHBcombinedcore$LocusID
#drop LocusID column
CHBcombinedcore<-subset(CHBcombinedcore,select= -(1))
write.csv(CHBcombinedcore,"~/LewinFiles/CHBcombinedcore.csv")

#I am doing this part on my computer RStudio because installing new packages in MSI portal is a nightmare
#Download DESeq
#Normalizing data
CHBcombinedcore<-read.csv("CHBcombinedcore.csv")
CHBcombinedcore<-as.data.frame(CHBcombinedcore)
#may or not need this line depending on your csv and version of R
row.names(CHBcombinedcore)<-CHBcombinedcore$...1
CHBcombinedcore<-subset(CHBcombinedcore, select=-(1))
CHBcombinedcore<-as.matrix(CHBcombinedcore)
#double check that all samples are present prior to normalizing
View(CHBcombinedcore)

CHBcombinedcorenorm<-vst(CHBcombinedcore2, blind=T)
write.csv(CHBcombinedcorenorm, "CHBcombinedcorenorm.csv")    
#
#Use if pulling entire combined count csv back in at a later date, otherwise see next comment
CHBcombinedcorenorm<-subset(CHBcombinedcorenorm, select=-1)
CHBcombinedcorenorm<-as.data.frame(CHBcombinedcorenorm)

#Use if continuing at same date
CHBcombinedcorenorm<-as.data.frame(CHBcombinedcorenorm)
CHBcombinedcorenorm$locus_tag<-row.names(CHBcombinedcorenorm)
row.names(CHBcombinedcorenorm)<-NULL
CHBcombinedcorenorm<-as.numeric(CHBcombinedcorenorm)

#input data. 
inputversion <- "input1"


#Use these lines if you previously used featurecounts and have saved individual csv files with count info
#skip if you just did feature counts, or read in a combined raw/normalized count file.
###organize csv files

ClinicalPAO1<-as.data.frame(ClinicalPAO1counts$counts)

ClinicalPAO1<-read.csv("ClinicalPAO1counts.csv")
HunterPAO1<-read.csv("HunterPAO1counts.csv")
BombergerPAO1<-read.csv("BombergerPAO1.csv", sep="")


#PAO1
BombergerPAO1<-as.data.frame(BombergerPAO1$counts)
BombergerPAO1$Loci_ID<-row.names(BombergerPAO1)
BombergerPAO1<-BombergerPAO1 %>% relocate("Loci_ID")
#remove all but genes that appear in 95% of human samples.
row.names(ClinicalPAO1)<-ClinicalPAO1$Loci_ID
ClinicalPAO1<-subset(ClinicalPAO1, select=-c(1,2))
ClinicalPAO1$count <- apply(ClinicalPAO1, 1, function(x) length(which(x=="0")))
#95% of 25 is 23.75, so there can only be a "0" in <=1 sample.
Clinical<-ClinicalPAO1 %>% filter(ClinicalPAO1$count <= 1)
Clinical$Loci_ID<-row.names(Clinical)
Clinical<-Clinical %>% relocate("Loci_ID")

rawallcounts<-inner_join(BombergerPAO1, Clinical, by="Loci_ID")
rawallcountsfinal<-inner_join(rawallcounts, HunterPAO1, by="Loci_ID")
write.csv(rawallcountsfinal, file ="CHBfinalrawcountsbombergerexon.csv")




#####Actual AS2 score prep work.
#load TIGR annotations
TIGR_ann <- read.csv("AnnotationTIGRFAM (1).txt", sep = "\t")
SunburstInput <- sunburst_data("SunburstCategories (1).txt", sep = "\t", type = "node_parent", node_attributes = "color")

#make your metadata file
metadata_file<-as.data.frame(colnames(CHBcombinedcorenorm))
metadata_file$ind_condition<-(colnames(CHBcombinedcorenorm))
write.csv(metadata_file, "~/LewinFiles/SunburstMetadata.csv")
##open the metadata file in a data editor (like excel or notepad)
##add in a short group descriptor in column 2, see my file for reference
##you will ultimately need to add your own for your data
#set metadata: change list names and filter condition to match datasets you want to compare
metadata_file <- read.table(file = "SunburstMetadata.csv", sep=",",header = T)



#setsession
invitro_dataset <- "AEC" #set based on model set you want to analyze
outputversion <- "run1"

human_list <- metadata_file %>% filter(type1 == "Human") %>% .$filename %>% str_replace_all("-", "_")
invitro_list <- metadata_file %>% filter(type1 == invitro_dataset) %>% .$filename %>% str_replace_all("-", "_") #adjust filter (currently "type1") to match correct column of metadata_file

#########################################################################################
#1. calculate AS2s and plot sunbursts for groups of samples, based on parameters set above (could also work for individual samples depending on parameters)
#run functions below first

AS_2 <- draw_sunburst_STANDARD_PA(CHBcombinedcorenorm, TIGR_ann, SunburstInput, human_list, invitro_list, 2, TRUE)

#output AS2
penalty_output <- AS_2[[2]]
AS_2_score <- penalty_output %>% mutate(score_2 = ifelse((penalty > 2 | penalty < -2), 0, 1)) %>% mutate(AS2 = mean(score_2))
write.table(AS_2_score, paste(outputversion, invitro_dataset, "AS2penalty.txt", sep = "."), sep = "\t", row.names = FALSE) 
head(AS_2_score)

#plot sunburst
AS_2[[1]] 
#pdf(file = paste(invitro_dataset, "AS2sunburst.pdf", sep = "."), width = 8.5, height = 11)
#View(AS_2[[1]])

dev.off()
#we finish editing this sunburst image manually in Adobe Illustrator



#######################################################################################
#2. to calculate AS2 from subsampling X replicates at once (set below), doing all possible combinations.
#run functions below first

#set session
reps <- 2 #number of replicates to run at once if subsampling
outputversion <- paste("run1.sub", reps, sep = "")

#JUST OUTPUT AS2: outputs average AS2 value across iterations
c(1:choose(length(invitro_list), reps)) %>% map(function(.x) {dummy1 <- score_target_vs_model_PA(CHBcombinedcorenorm, human_list, combn(invitro_list, reps)[,.x]) %>% .$penalty
dummy2 <- abs(dummy1) < 2
dummy3 <- mean(dummy2)
return(dummy3)}) %>% purrr::reduce(c) %>% mean

#outputs table of each iteration with avg zscores across iterations. Can edit then use as input into sunburst to directly graph this version of the analysis 
zscores <- c(1:choose(length(invitro_list), reps)) %>% map(function(.x) {dummy1 <- score_target_vs_model_PA(CHBcombinedcorenorm, human_list, combn(invitro_list, reps)[,.x]) %>% .$penalty
return(dummy1)})
zscores2 <- as.data.frame(do.call(cbind, zscores))
zscores3 <- cbind(CHBcombinedcorenorm$locus_tag, zscores2)
colnames(zscores3)[1] <- "locus_tag"
zscores4 <- zscores3 %>% mutate(xmean = rowMeans(across(starts_with("V"))))
write.table(zscores4, paste(outputversion, invitro_dataset, "zscores.txt", sep = "."), sep = "\t", row.names = FALSE) 



#Individual counts of interest
#You can edit these to output your similar genes as you go. 
#These are meant to be run as you go, so as you do the first zscore calculation
#the output is zscores4, which will be the zscores for the first group (in this case ana4h)
##Once the first csv is written, return to line __ and select a new group and rerun lines ___ to ____
ana4hAS2bygene<-zscores4
ana4hAS2bygene<-subset(ana4hAS2bygene, select=-c(2,3,4))
ana4hAS2bygene %>% count(xmean >= 2 | xmean<=-2)
count(ana4hAS2bygene$xmean>=2)
write.csv(ana4hAS2bygene,"ana4hAS2bygene.csv")

ana16hAS2bygene<-zscores4
ana16hAS2bygene<-subset(ana16hAS2bygene, select=-c(2,3,4))
write.csv(ana16hAS2bygene,"ana16hAS2bygene.csv")

hne16hAS2bygene<-zscores4
hne16hAS2bygene<-subset(hne16hAS2bygene, select=-c(2,3,4))
write.csv(hne16hAS2bygene,"hne16hAS2bygene.csv")

hne4hAS2bygene<-zscores4
hne4hAS2bygene<-subset(hne4hAS2bygene, select=-c(2,3,4))
write.csv(hne4hAS2bygene,"hne4hAS2bygene.csv")

unt4hAS2bygene<-zscores4
unt4hAS2bygene<-subset(unt4hAS2bygene, select=-c(2,3,4))
write.csv(unt4hAS2bygene,"unt4hAS2bygene.csv")

unt16hAS2bygene<-zscores4
unt16hAS2bygene<-subset(unt16hAS2bygene, select=-c(2,3,4))
write.csv(unt16hAS2bygene,"unt16hAS2bygene.csv")

AECAS2bygene<-zscores4
AECAS2bygene<-subset(AECAS2bygene, select=-c(2,3,4))
write.csv(AECAS2bygene,"AECAS2bygene.csv")

EPIScfm2AS2bygene<-zscores4
EPIScfm2AS2bygene<-subset(EPIScfm2AS2bygene, select=-c(2,3,4))
write.csv(EPIScfm2AS2bygene,"EPIScfm2AS2bygene.csv")

#########################################################################################
#3. to calculate zscores for individual replicates
#run functions below first
##read in a previously created file OR
datasets <- (read.csv("PA_individual_input_all.txt", header = FALSE)) #text file with individual sample names listed in single column

#create one here.
datasets<-colnames(CHBcombinedcorenorm)
datasets<-as.data.frame(datasets)
#error trouble shooting
##these next lines also reference a previously created metadata file that was formed earlier
#if you get errors, make sure that file is read in.
#Additionally, check the names of column 1 in the metadata file, and ensure that they are the same as the list of file names in "datasets"
#if not, fix whichever one is incorrect (and then rename the column name in the CHBcombinedcorenormalized)

AS2_bygene <- data.frame(locus_tag=c(0))

for(i in 1:length(datasets)) {
  #invitro_list <- metadata_file %>% filter(ind_condition == datasets[i])# %>% .$filename %>% str_replace_all("-", "_")
  invitro_list <- datasets[i]
  dummy1 <- score_target_vs_model_PA(CHBcombinedcorenorm, human_list, invitro_list) #%>% .$penalty
  colnames(dummy1) <- c("locus_tag", invitro_list)
  AS2_bygene <- merge(AS2_bygene, dummy1, all=T)
}


AS2_bygene <- tail(AS2_bygene, n=-1)
write.table(AS2_bygene, paste("CHBcombined", "AS2_bygene_individual.csv", sep = "."), sep = ",", row.names = FALSE) #input into excel and calculate AS2 from penalties




#########################################################################################

library(ggplot2)
library(vegan)
library(dplyr)

#improting data
parent_info <-read.table("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos",sep=";")
FlexDepth <- read.table("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Data_lists/FlexSeq_Depth.txt")
ParentDepth<- read.table("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Data_lists/parent3000_depth.txt")
Flex_site_sum <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/FlexSeq_site_summary.txt",header = TRUE, sep = "\t")
Parent_site_sum <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parents3000_all_sorted_site_summary.txt",header = TRUE, sep = "\t")
Flex_hmp_names <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/progeny_hmp_names.txt", sep = "\t")
Parent_hmp_names <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parent_hmp_names.txt", sep = "\t")

#This is to get the MAF from tassel site summary data
FlexMAF <- subset(Flex_site_sum, select = c(Minor.Allele.Frequency))
ParentMAF <- subset(Parent_site_sum, select = c(Minor.Allele.Frequency))

#Renaming depth columns
FlexDepth <- FlexDepth %>% 
  rename(
    Depth = V1,
  )
ParentDepth <- ParentDepth %>% 
  rename(
    Depth = V1,
  )
# This function graphs the depth of our data sets
graphhistogramDepth <- function(dat){
  plot1 <- ggplot(dat, aes(x=log10(Depth))) + geom_density(colour="#1441ba") + theme_bw()
  return(plot1)
}
# This function graphs the MAF of our data sets
graphhistogramMAF <- function(dat){
  plot1 <- ggplot(dat, aes(x=Minor.Allele.Frequency)) + geom_density(colour="#f21818") + theme_bw()
  return(plot1)
}
graphhistogramDepth(FlexDepth)
graphhistogramDepth(ParentDepth)
graphhistogramMAF(Flex_site_sum)
graphhistogramMAF(Parent_site_sum)







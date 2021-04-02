#library(ggplot2)
#library(vegan)
#library(dplyr)
library(reshape2)
library(tidyverse)
library(vegan)


Flex_hmp_names <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/progeny_hmp_names.txt", sep = "\t")
Parent_hmp_names <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parent_hmp_names_clean.txt", sep = "\t")
Scores <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parent_progeny_score_table_normalized.txt", sep = " ")
Scores_not_normal <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parent_progeny_score_table.txt", sep = " ")
max_scores <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/max_scores.txt", sep = " ")

#Fixing up the name data to just have the names and not other stuff
Flex_hmp_names2 <- select(Flex_hmp_names, -1:-11)
Parent_hmp_names2 <- select(Parent_hmp_names, -1:-11)

Parent_hmp_names2 <- gather(Parent_hmp_names2, key)
Flex_hmp_names2 <- gather(Flex_hmp_names2, key)

typeof(Parent_hmp_names2)

# The Colnames functino renames columns with 1 for every column. Did this with 2 sets of data. 
new_col_names <- Parent_hmp_names2[, "value"]
colnames(Scores) <- new_col_names

new_row_names <- Flex_hmp_names2[, "value"]
rownames(Scores) <- new_row_names

Scores2 <- tibble::rownames_to_column(Scores, "Progeny")

Scores3 <- melt(Scores2, id = c("Progeny"))
head(Scores3)
###########################################  Plots 
# Scatter plot of all progeny. this works but the plots will not apper if made in the loop
MakeScatter <- function(parent_name){
  plot1 <- ggplot(Scores2, aes_(x=as.name('Progeny'), y=as.name(parent_name))) + geom_point(shape=1)
  return(plot1)
}

for (i in 2:18){
  Parent = colnames(Scores2[i])
  MakeScatter(Parent)
} 

#Histograms (Kinda useless tbh because it shows counts of progeny with certain scores)
MakeHis <- function(parent_name){
  plot1 <- ggplot(Scores2, aes_(x=as.name(parent_name))) + geom_histogram(binwidth=1) + theme_bw()
  return(plot1)
}

i=0
for (i in 2:18){
  Parent = colnames(Scores2[i])
  print(MakeHis(Parent))
} 

# Density plot showing all parents with their respective scores, interesting to see peak scores of all parents
ggplot(Scores3, aes_(x=as.name('value'), colour=as.name('variable'))) + geom_density() + theme_bw()

# Density plot of a single progeny compared to all parents
randomNum <- sample(1:3000, 1, replace=TRUE)
for (i in 1:10){
Progeny_name = Flex_hmp_names2[randomNum,2]
Scores3_single <- subset(Scores3, Progeny == Progeny_name)
print(ggplot(Scores3_single, aes_(x=as.name('value'), colour=as.name('Progeny'))) + geom_density() + theme_bw())
randomNum <- sample(1:3000, 1, replace=TRUE)

}


# scatter plot of individual progeny against all parents 
randomNum <- sample(1:3000, 1, replace=TRUE)
for (i in 1:20){
  Progeny_name = Flex_hmp_names2[randomNum,2]
  Scores3_single <- subset(Scores3, Progeny == Progeny_name)
  print(ggplot(Scores3_single, aes_(x=as.name('variable'), y=as.name('value'), colour=as.name('variable'))) + geom_point(size = 5) +
          theme(axis.text.x=element_blank()) + ylab(paste0("Score for ",Progeny_name)))
  randomNum <- sample(1:3000, 1, replace=TRUE)
}
  



#Other plot idea stuff
ggplot(dat, aes(x=rating, fill=cond)) + geom_histogram(binwidth=.5) 

ggplot(data=tips, aes(x=day)) +
  geom_bar(stat="count")

ggplot(dat, aes(x=xvar, y=yvar)) +
  geom_point(shape=1)

ggplot(df, aes(x=weight)) + 
  geom_histogram(binwidth=1)

ggplot(data=dat, aes(x=time, y=total_bill, fill=time)) +
  geom_bar(colour="black", stat="identity")
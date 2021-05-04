#library(ggplot2)
#library(vegan)
#library(dplyr)
library(reshape2)
library(tidyverse)
library(vegan)
library(Rfast)

Flex_hmp_names <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/progeny_hmp_names.txt", sep = "\t")
Parent_hmp_names <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parent_hmp_names_clean.txt", sep = "\t")
Scores <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parent_progeny_score_table_normalized.txt", sep = " ")
Scores_not_normal <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parent_progeny_score_table.txt", sep = " ")
max_scores <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/max_scores.txt", sep = " ")
Parents_Percent_Heterozygous <- read.table ("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/Parents_Percent_Heterozygous.txt", sep = "\t", header = TRUE)
parent_progeny_key <- read.table ("/home/drt83172/Documents/Tall_fescue/parent_projeny.txt", sep = "\t", header = FALSE)


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
Scores3 <- Scores3 %>% separate(variable, c("variable", "delete_me"), "_")
Scores3 <- subset(Scores3, select = -c(delete_me))

# Make a table with the standard deviations away from the mean things are
standardDevs <- apply(Scores,2,sd)
Means <- colMeans(Scores)

standardDevs[1]/100

ScoreDistances <- Scores
i = 1
j = 1
for (i in 1:ncol(ScoreDistances)){
  parentMean <- Means[i]
  parentSTD <- standardDevs[i]
  for (j in 1:nrow(ScoreDistances)){
    ScoreDistances[j,i] <- (Scores[j,i] - parentMean)/parentSTD
  }
}
ScoreDistances2 <- tibble::rownames_to_column(ScoreDistances, "Progeny")
ScoreDistances3 <- melt(ScoreDistances2, id = c("Progeny"))

# Finding parents with highest two psoitive stds away from the mean. Checks if others are close and if so adds two more
predicted_Parents <- as.data.frame(matrix(ncol = ncol(Scores), nrow = nrow(Scores)))
rownames(predicted_Parents) <- Scores2[,1]
predicted_Parents <- tibble::rownames_to_column(predicted_Parents, "Progeny")

for (i in 1:nrow(predicted_Parents)){ 
  Progeny_name <- predicted_Parents[i,1]
  single_Progeny <- subset(ScoreDistances3, Progeny == Progeny_name)
  single_Progeny_ordered <- single_Progeny[order(single_Progeny$value),c(1,2,3)] #Orders data set by chosen column
  predicted_Parents[i,2] <- single_Progeny_ordered[ncol(ScoreDistances),2] # For some reason this translates the names into numbers representing the names in alphabetical order
  predicted_Parents[i,3] <- single_Progeny_ordered[ncol(ScoreDistances)-1,2]
  if ((single_Progeny_ordered[ncol(ScoreDistances)-2,3] - single_Progeny_ordered[ncol(ScoreDistances)-3,3]) < .5){
    predicted_Parents[i,4] <- single_Progeny_ordered[ncol(ScoreDistances)-2,2]
    predicted_Parents[i,5] <- single_Progeny_ordered[ncol(ScoreDistances)-3,2]
    }
}

# Making a refrence key to know what progeny go with what parents
new_row_names <- parent_progeny_key[, "V1"]
rownames(parent_progeny_key) <- new_row_names
parent_progeny_key <- subset(parent_progeny_key, select = -c(V1))

# Adding key  to data that may need it
Scores2 <- merge(Scores2,parent_progeny_key, by.x = "Progeny", by.y = 0)


###########################################  Plots 
# Scatter plot of all progeny. this works but the plots will not apper if made in the loop
MakeScatter <- function(parent_name){
  plot1 <- ggplot(Scores2, aes_(x=as.name('Progeny'), y=as.name(parent_name))) + geom_point(shape=1)
  return(plot1)
}

for (i in 2:18){
  Parent = colnames(Scores2[i])
  print(MakeScatter(Parent))
} 



#Histograms of 1 parent vs all progeny
MakeHis <- function(parent_name){
  plot1 <- ggplot(Scores2, aes_(x=as.name(parent_name))) + geom_histogram(binwidth=.01) + theme_bw() +
    xlim(.6, 1) + ylim(0, 1800)
  return(plot1)
}

i=0
x = 301
for (i in 2:18){
  Parent = colnames(Scores2[i])
  KnownParent = subset(Scores2, V2 == x)
  cols <- c(i, 19)
  KnownParent <- KnownParent[,cols] 
  print(ggplot(Scores2, aes_(x=as.name(Parent))) + geom_histogram(binwidth=.01)  + geom_histogram(data = KnownParent, fill = "darkblue", binwidth=.01)+ theme_bw() +
    xlim(.6, 1) + ylim(-1, 1800)) 
  if (x == 308 | x == 310 | x == 316){
    x = x+2}
}


# Density plot showing all parents with their respective scores, interesting to see peak scores of all parents
ggplot(Scores3, aes_(x=as.name('value'), colour=as.name('variable'))) + geom_density() + theme_bw()
# Using standard deviations away from the mean
ggplot(ScoreDistances3, aes_(x=as.name('value'), colour=as.name('variable'))) + geom_density() + theme_bw()

# Density plot of a single progeny compared to all parents
randomNum <- sample(1:3000, 1, replace=TRUE)
for (i in 1:10){
Progeny_name = Flex_hmp_names2[randomNum,2]
Scores3_single <- subset(Scores3, Progeny == Progeny_name)
print(ggplot(Scores3_single, aes_(x=as.name('value'), colour=as.name('Progeny'))) + geom_density() + theme_bw())
randomNum <- sample(1:3000, 1, replace=TRUE)
}

# Density plot of a single parent to all progeny
Scores3_single_parent <- subset(Scores3, variable == 307)
ggplot(Scores3_single_parent, aes_(x=as.name('value'), colour=as.name('variable'))) + geom_density() + theme_bw() +
  scale_color_manual(values=c("#0018e6"))


# scatter plot of individual progeny against all parents 
randomNum <- sample(1:3000, 1, replace=TRUE)
for (i in 1:10){
  Progeny_name = Flex_hmp_names2[randomNum,2]
  Name <- which(rownames(parent_progeny_key)== Progeny_name)
  Parent_Name <- parent_progeny_key[Name,1]
  Scores3_single <- subset(Scores3, Progeny == Progeny_name)
  print(ggplot(Scores3_single, aes_(x=as.name('variable'), y=as.name('value'), colour=as.name('variable'))) + geom_point(size = 5) +
    theme_bw() + ylab(paste0("Score for ",Progeny_name)) + 
    annotate("text", x = 1, y = .96, label = Parent_Name, size = 5, colour = "red"))
  randomNum <- sample(1:3000, 1, replace=TRUE)
}

# annotate("text", label = Parent_Name, x = Parent_Name, y = .96, size = 5, colour = "red")
# Using standard deviations away from the mean
randomNum <- sample(1:3000, 1, replace=TRUE)
for (i in 1:20){
  Progeny_name = Flex_hmp_names2[randomNum,2]
  Scores3_single <- subset(ScoreDistances3, Progeny == Progeny_name)
  print(ggplot(Scores3_single, aes_(x=as.name('variable'), y=as.name('value'), colour=as.name('variable'))) + geom_point(size = 5) +
          theme(axis.text.x=element_blank()) + ylab(paste0("Score for ",Progeny_name)))
  randomNum <- sample(1:3000, 1, replace=TRUE)
}  

# Graphing the mean score of each parent against their heterozygosity
Parents_Percent_Heterozygous
Parent_Means <- as.data.frame(colMeans(Scores))
colnames(Parent_Means) <- "Mean_Score"
Parent_Data <- merge(Parent_Means,Parents_Percent_Heterozygous, by.x = 0, by.y = "Taxa.Name")
colnames(Parent_Data)[1] <- "Parent"

lm_eqn <- function(df){
  m <- lm(Proportion.Heterozygous ~ Mean_Score, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

ggplot(Parent_Data, aes_(x=as.name('Mean_Score'), y=as.name('Proportion.Heterozygous'))) + geom_point(size = 5) +
  theme_bw() + geom_text(aes(label=Parent),hjust=0, vjust=0) + 
  geom_text(x = .97, y = .48, label = lm_eqn(Parent_Data), parse = TRUE)




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
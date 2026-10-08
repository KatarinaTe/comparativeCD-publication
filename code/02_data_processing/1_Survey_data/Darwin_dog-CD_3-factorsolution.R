### R-script used for 14 items from the Darwin's Ark survey:
### 2024-02-13 ###


# LOAD LIBRARIES
require(tidyverse)
library(dplyr)
library(phateR)

# plotting and visualization
require(ggpubr)
require(colorspace)
require(viridis)
require("wesanderson")
require(scales)
require(ggrepel)

# data management
require(data.table)
require(zoo)

# dates and times
require(anytime)
require(lubridate)

# free text and strings
require(stringr)
require(stringdist)
require(english)
require(tm)
require(topicmodels)
require(lattice)
require(tidytext)
data(stop_words)

# locations
require(zipcode)

# correlation and regression
require(pspearman)
require(jtools)
require(caret)

# factor analysis
require(psych)
require(nFactors)
require(FactoMineR)
require(factoextra)
library(corrplot)
library(car)
library(BBmisc)
library(mice)

library(mirt)
library(lavaan)
#library(rstan)

library(tidyverse)

#install.packages(ggmirt)
#devtools::install_github("masurp/ggmirt")
library(ggmirt)


#### load data ######

setwd("/path/to/")

questions <- read.csv("DarwinsArk_20220315_questions.csv", sep= ";")
dog <- read.csv("DarwinsArk_20220715_dogs.csv", sep= ",")
answers <- read.csv("DarwinsArk_20220715_answers.csv")

table(questions$string, questions$id)

df <- data.frame(dog=answers$dog, question=answers$question, answer=answers$answer)
head(df)
unique(df$question)
####select all dogCD questions #####
selected_questions<- c(7,93,95,145,146,147,148,149,150,151,152,153,154,155) 

### Filter the data frame to include only the selected questions
filtered_df <- df %>%
  filter(question %in% selected_questions)

# Pivot the filtered table using a list-column
pivoted_df <- filtered_df %>%
  pivot_wider(names_from = question, values_from = answer, values_fn = list)

# Convert the list-column to a data frame
response_df <- as.data.frame(pivoted_df)
nrow(response_df)


# Handle missing column names
colnames(response_df) <- make.unique(colnames(response_df))

# Convert numeric columns to numeric, handling NULL values
numeric_columns <- setdiff(colnames(response_df), "dog")
response_df[numeric_columns] <- lapply(response_df[numeric_columns], function(x) {
  x <- lapply(x, function(y) ifelse(is.null(y), NA, y))
  as.numeric(x)
})

summary(response_df)

# Switch the values in the response_df table (OBS define which questions need to be switched!)
response_df1 <- response_df %>%
  mutate(across(c("7","93", "153","154"), ~ ifelse(is.na(.), NA, ifelse(. == 3, 1, ifelse(. == 4, 0, ifelse(. == 1, 3, ifelse(. == 0, 4, .)))))))

#create the right variables (1-5 instead of 0-4) for all questions, 
#this means that questions with response 0-5 becomes 1-5 with group 5 including the top two scores.

response_df2 <- response_df1 %>%
  mutate(across(c("7","93","95","145","146","147","148","149","150","151","152", "153","154","155"), ~ ifelse(. %in% 0:5, . + 1, .)))

###print response_df2:
#write.table(response_df2, file= "response_df_dogCDitems.txt", row.names =F, sep = "\t", col.names= T, quote = F ) 
####load factor 1 scores####
response_df2 <- read.table("response_df_dogCDitems.txt", header= T, sep= "\t", check.names = FALSE)

# Remove the "dog" column from response_df
response_df3 <- response_df2[-1]

#### FACTOR ANALYSIS ####
#### run: ocd data ####

# perform nfactors test
fa.ocd <- nfactors(response_df3, n=10)

# remove NAs to continue:
response_df4 <- na.omit(response_df3)
fa.ocd <- nfactors(response_df4, n=10)

nrow(response_df3) #25302
nrow(response_df4) #12989 no NA
colnames(response_df3)

# perform parallel analysis #nFactors::parallel
fa.ocd.ap <- parallel(subject=nrow(response_df4),
                      var=ncol(response_df4),
                      rep=100,
                      cent=.05)
plotParallel(fa.ocd.ap)

fa.ocd.nf <- nfactors(response_df4, n=25)
fa.ocd.ev <- eigen(cor(as.matrix(response_df4))) 
fa.ocd.ns <- nScree(x=fa.ocd.ev$values, aparallel=fa.ocd.ap$eigen$qevpea)

plotnScree(fa.ocd.ns) #4 factors? tried but SSloadings below 1 for the fourth factor. Try 3.

#plotuScree(fa.ocd.ev$values)

mat_ocd.cor <- cor(as.matrix(response_df4))

#try the rest of the script by using these:
mat_cor <- mat_ocd.cor
fa.ap <- fa.ocd.ap 
fa.nf <- fa.ocd.nf 
fa.ev <- fa.ocd.ev 
fa.ns <- fa.ocd.ns 
mat <- response_df4
qN = ncol(mat)
dN = nrow(mat)
freezeDate = format(Sys.Date(), "%a %b %d")
mat_dogs = rownames(mat)
mat_ques = colnames(mat)

####KMO####
#testing whether a correlation matrix is factorable (i.e., the correlations differ from 0) is the Kaiser-Meyer-Olkin KMO test of factorial adequacy.
#Henry Kaiser introduced a Measure of Sampling Adequacy (MSA) of factor analytic data matrices in 1970.[2] Kaiser and Rice then modified it in 1974.[3]
# The Kaiser-Meyer-Olkin (KMO) used to measure sampling adequacy is a better measure of factorability. According to Kaiser’s (1974) guidelines, 
#a suggested cutoff for determining the factorability of the sample data is KMO ≥ 60. The total KMO is 0.94, indicating that, 
#based on this test, we can probably conduct a factor analysis.
fa.kmo = KMO(r=mat_cor) # Overall MSA = 0.94 (DARWIN). Overall MSA = 0.77 (dog-cd questions)
fa.kmo
fa.kmo$MSAi



#####Bartlett’s Test#####
# Bartlett’s Test of Sphericity: Small values (8.84e-290 < 0.05) of the significance level indicate that a factor analysis may be 
#useful with our data.
fa.bart = cortest.bartlett(mat_cor, n = nrow(mat)) # p = 0 for both Darwin and ASD and ocd
fa.bart

#Bartlett (1951) proposed that -ln(det(R)*(N-1 - (2p+5)/6) was distributed as chi square if R were an identity matrix. 
#A useful test that residuals correlations are all zero. Contrast to the Kaiser-Meyer-Olkin test.


####determinant####
#'det' calculates the determinant of a matrix. 
#determinant is a generic function that returns separately the modulus of the determinant, 
#optionally on the logarithm scale, and the sign of the determinant.

# Positive determinates will run
fa.det = det(cor(mat_cor)) # positive, but very small
fa.det #dog-cd=6.917628e-19

####cronbach's alpha####
library(ltm)
cronbach.alpha(response_df4, CI=TRUE) #alpha: 0.754 (ocd)

####run factor analysis #####
nFac = 3 

#### Exploratory Factor Analysis with Varimax Rotation and Principal Axes ####
#write.table(mat_cor, file= "mat_cor.txt", row.names =T, sep = " ", col.names= T, quote = F ) #used to check pairwise correlations
fa.varimax <- fa(r=mat_cor,cor="poly",
                 nfactors = nFac,
                 fm="pa",
                 max.iter=200,
                 rotate="varimax", SMC=FALSE) #

fa.varimax$loadings
fa.varimax$uniquenesses

#If, however, a solution fails to be achieved, it is useful to try again using ones (SMC =FALSE). 
#for ocd -data no difference was seen between SMC = True or false.

#fm="pa" will do the principal factor solution, - this was used for Darwin I think
#fm="ml" will do a maximum likelihood factor analysis. 

#impute	 -"median" or "mean" values are used to replace missing values - try THIS??

# Item Uniqueness: uniqueness, sometimes referred to as noise, corresponds to the proportion of variability, 
#which can not be explained by a linear combination of the factors. 
#A high uniqueness for a variable indicates that the factors do not account well for its variance. 
#(what exists in the unique variance instead of common variance)
fa.varimax.uniqueness = data.table(question = names(fa.varimax$uniquenesses),
                                   uniqueness = fa.varimax$uniquenesses) %>% merge((questions %>% arrange(id) %>% 
                                                                                      dplyr::select(id,string,source=tags,options,survey=title) %>% 
                                                                                      mutate(id = as.character(id))), 
                                                                                   by.x = "question", by.y = "id") %>% arrange(as.numeric(question))
plot(fa.varimax.uniqueness$question, fa.varimax.uniqueness$uniqueness)
low_uniq <- fa.varimax.uniqueness[fa.varimax.uniqueness$uniqueness <0.5,]
high_uniq <- fa.varimax.uniqueness[fa.varimax.uniqueness$uniqueness >=0.5,]


#### run IRT for the three-factor solution ####

#### F1: repetition severity ####
selected_questions<- c(7,153,154,155) 

### Filter the data frame to include only the selected questions
filtered_df <- df %>%
  filter(question %in% selected_questions)

# Pivot the filtered table using a list-column
pivoted_df <- filtered_df %>%
  pivot_wider(names_from = question, values_from = answer, values_fn = list)

# Convert the list-column to a data frame
response_df <- as.data.frame(pivoted_df)

# Handle missing column names
colnames(response_df) <- make.unique(colnames(response_df))

# Convert numeric columns to numeric, handling NULL values
numeric_columns <- setdiff(colnames(response_df), "dog")
response_df[numeric_columns] <- lapply(response_df[numeric_columns], function(x) {
  x <- lapply(x, function(y) ifelse(is.null(y), NA, y))
  as.numeric(x)
})


# Switch the values in the response_df table (OBS define which questions need to be switched!)
response_df1 <- response_df %>%
  mutate(across(c("7","153","154"), ~ ifelse(is.na(.), NA, ifelse(. == 3, 1, ifelse(. == 4, 0, ifelse(. == 1, 3, ifelse(. == 0, 4, .)))))))

#create the right variables (1-5 instead of 0-4) for all questions
response_df2 <- response_df1 %>%
  mutate(across(c("7","153","154","155"), ~ ifelse(. %in% 0:5, . + 1, .)))

# Remove the "dog" column from response_df
response_df3 <- response_df2[-1]
head(response_df2)

fitGraded <- mirt(response_df3, 1, itemtype = "graded")
fitGraded


##get IRT estimates & plots
setnames(response_df3, old = c("7","153","154","155"), 
         new = c("q7","q153","q154","q155"))

model <- "f1 =~ q7+q153+q154+q155"
fitCTT <- cfa(model, data = response_df3) 

# IRT solution
summary(fitGraded)

# CTT solution
standardizedsolution(fitCTT) %>%
  filter(op == "=~") %>% dplyr::select(rhs, F1 = est.std)

#IRT parameters
params <- coef(fitGraded, IRTpars = TRUE, simplify = TRUE)
round(params$items, 2) 
M2(fitGraded, type = "M2*", calcNULL = FALSE, na.rm=TRUE) #Error: M2() statistic cannot be calculated due to too few degrees of freedom
itemfit(fitGraded,  na.rm=TRUE)

tracePlot(fitGraded, theta_range = c(-3, 8)) +
  labs(color = "Answer Options")

itemInfoPlot(fitGraded, c(1:4), facet = T)
testInfoPlot(fitGraded,theta_range = c(-3, 3), adj_factor = .5)

itemInfoPlot(fitGraded, c(1:4), theta_range = c(-3, 3), facet = F, legend = T) + 
  scale_color_brewer(palette = "Set3") 

#create the theta values for each factor:
summary(fitGraded)
item.weights <- c(0.81,4.86,25.77,1.44)
tabscores1 <- fscores(fitGraded,method='EAP',full.scores.SE=TRUE, item_weights = item.weights) #
tabscores1 <- as.data.frame(tabscores1)
head(tabscores1)
tabscores3 <- fscores(fitGraded,method='MAP',full.scores.SE=TRUE, item_weights = item.weights) #There were 50 or more warnings (use warnings() to see the first 50)
tabscores3 <- as.data.frame(tabscores3)

## Plot everything on top
plot(tabscores1$F1, tabscores1$SE_F1, xlab="Factor score", ylab="Standard error", title("F1 (CCD_3F-solution): Rep severity"), xlim=c(-3, 3), ylim=c(0,1))
points(tabscores3$F1, tabscores3$SE_F1, col=c("blue"), pch=4, cex=.8)
legend("topright",  legend=c("EAP", "MAP"), pch=19, bty="n", col=c("green", "blue"))

#combine with dog ID - using MAP
tabscores3$number <- rownames(tabscores3)
response_df2$number <- rownames(response_df2)
nrow(tabscores3)
nrow(response_df2)
F1 <- left_join(response_df2, tabscores3, by ="number")
head(F1)
nrow(F1)
F1_CCD3F <-data.frame(dog = F1$dog, F1 = F1$F1, F1_SE = F1$SE_F1,
                      item7 = F1$"7",item153 = F1$"153",item154 = F1$"154",item155 = F1$"155" )
head(F1_CCD3F)
head(F1)
nrow(F1_CCD3F)

x <- F1
x$na_count <- apply(x[,2:5], 1, function(x) {sum(is.na(x))})
length(x$na_count)
length(x$F1)
plot(x$SE_F1, x$na_count, xlab="SE",pch=19, ylab="NA count", title("F1_CCD3F"))
plot(x$F1, x$na_count, xlab="Factor score",pch=19, ylab="NA count", title("F1_CCD3F"))
plot(x$F1, x$SE_F1, xlab="Factor score",pch=19, ylab="SE", title("F1_CCD3F"))

table(F1_CCD3F$F1_SE, exclude=NULL)
plot(F1_CCD3F$F1_SE, F1_CCD3F$item7, ylab="item7",pch=19, xlab="SE", title("F1_CCD3F"))

#save factor scores
#write.table(F1_CCD3F, file= "F1_CCD3F_240213.txt", row.names =F, sep = "\t", col.names= T, quote = F ) 
####load factor 1 scores####
F1_CCD3F <- read.table("F1_CCD3F_240213.txt",header= T,   sep= "\t")
head(F1_CCD3F)
nrow(F1_CCD3F)

### F2: compulsive staring/trancing/pacing  ###

selected_questions<- c(93,95,150) 

### Filter the data frame to include only the selected questions
filtered_df <- df %>%
  filter(question %in% selected_questions)

# Pivot the filtered table using a list-column
pivoted_df <- filtered_df %>%
  pivot_wider(names_from = question, values_from = answer, values_fn = list)

# Convert the list-column to a data frame
response_df <- as.data.frame(pivoted_df)

# Handle missing column names
colnames(response_df) <- make.unique(colnames(response_df))

# Convert numeric columns to numeric, handling NULL values
numeric_columns <- setdiff(colnames(response_df), "dog")
response_df[numeric_columns] <- lapply(response_df[numeric_columns], function(x) {
  x <- lapply(x, function(y) ifelse(is.null(y), NA, y))
  as.numeric(x)
})


# Switch the values in the response_df table (OBS define which questions need to be switched!)
response_df1 <- response_df %>%
  mutate(across(c("93"), ~ ifelse(is.na(.), NA, ifelse(. == 3, 1, ifelse(. == 4, 0, ifelse(. == 1, 3, ifelse(. == 0, 4, .)))))))

#create the right variables (1-5 instead of 0-4) for all questions
response_df2 <- response_df1 %>%
  mutate(across(c("93","95","150"), ~ ifelse(. %in% 0:5, . + 1, .)))

# Remove the "dog" column from response_df
response_df3 <- response_df2[-1]

#response_df2 <- response_df1
head(response_df2)
fitGraded <- mirt(response_df3, 1, itemtype = "graded")
fitGraded

####F2 get IRT estimates & plots ####
setnames(response_df3, old = c("93","95","150"), 
         new = c("q93","q95","q150"))

model <- "f1 =~ q93+q95+q150"
fitCTT <- cfa(model, data = response_df3) 

# IRT solution
summary(fitGraded)

# CTT solution
standardizedsolution(fitCTT) %>%
  filter(op == "=~") %>% dplyr::select(rhs, F1 = est.std)

#IRT parameters
params <- coef(fitGraded, IRTpars = TRUE, simplify = TRUE)
round(params$items, 2) 
M2(fitGraded, type = "C2", calcNULL = FALSE, na.rm=TRUE) #Error: M2() statistic cannot be calculated due to too few degrees of freedom
itemfit(fitGraded,  na.rm=TRUE)

tracePlot(fitGraded, theta_range = c(-3, 8)) +
  labs(color = "Answer Options")

itemInfoPlot(fitGraded, c(1:3), facet = T)
testInfoPlot(fitGraded,theta_range = c(-4, 4), adj_factor = .5)

itemInfoPlot(fitGraded, c(1:3), theta_range = c(-4, 6), facet = F, legend = T) + 
  scale_color_brewer(palette = "Set2") 

#create the theta values for each factor:
summary(fitGraded)
item.weights <- c(1.01, 2.69,1.58)
tabscores1 <- fscores(fitGraded,method='EAP',full.scores.SE=TRUE, item_weights = item.weights) #
tabscores1 <- as.data.frame(tabscores1)
tabscores3 <- fscores(fitGraded,method='MAP',full.scores.SE=TRUE, item_weights = item.weights) #
tabscores3 <- as.data.frame(tabscores3)

#warnings() #A/Inf replaced by maximum positive value

## Plot everything on top
plot(tabscores1$F1, tabscores1$SE_F1, xlab="Factor score", ylab="Standard error", title("F2_CCD3F"), xlim=c(-2.5, 5), ylim=c(0,1))
points(tabscores3$F1, tabscores3$SE_F1, col=c("blue"), pch=4, cex=.8)
legend("topright",  legend=c("EAP", "MAP"), pch=19, bty="n", col=c("green", "blue"))

#combine with dog ID - using EAP
tabscores1$number <- rownames(tabscores1)
response_df2$number <- rownames(response_df2)
nrow(tabscores1)
nrow(response_df2)
F2 <- left_join(response_df2, tabscores1, by ="number")
head(F2)
F2_CCD3F <-data.frame(dog = F2$dog, F2 = F2$F1, F2_SE = F2$SE_F1, item93 = F2$"93", 
                      item95 = F2$"95", item150 = F2$"150")
tail(F2_CCD3F)

x <- F2
x$na_count <- apply(x[,2:4], 1, function(x) {sum(is.na(x))})
length(x$na_count)
length(x$F1)
plot(x$SE_F1, x$na_count, xlab="SE",pch=19, ylab="NA count", title("F2_CCD3F"))
plot(x$F1, x$na_count, xlab="Factor score",pch=19, ylab="NA count", title("F2_CCD3F"))
plot(x$F1, x$SE_F1, xlab="Factor score",pch=19, ylab="SE", title("F2_CCD3F"))

#save factor scores
#write.table(F2_CCD3F, file= "F2_CCD3F_240213.txt", row.names =F, sep = "\t", col.names= T, quote = F ) 
####load factor 2 scores####
F2_CCD3F <- read.table("F2_CCD3F_240213.txt",header= T,   sep= "\t")
head(F2_CCD3F)
F2_CCD3F[!is.na(F2_CCD3F$F2),] #no NA
nrow(F2_CCD3F)

#### F3 - repetition frequency CCD3F  ####

selected_questions<- c(145,146,147,148,149,151,152)

### Filter the data frame to include only the selected questions
filtered_df <- df %>%
  filter(question %in% selected_questions)

# Pivot the filtered table using a list-column
pivoted_df <- filtered_df %>%
  pivot_wider(names_from = question, values_from = answer, values_fn = list)

# Convert the list-column to a data frame
response_df <- as.data.frame(pivoted_df)

# Handle missing column names
colnames(response_df) <- make.unique(colnames(response_df))

# Convert numeric columns to numeric, handling NULL values
numeric_columns <- setdiff(colnames(response_df), "dog")
response_df[numeric_columns] <- lapply(response_df[numeric_columns], function(x) {
  x <- lapply(x, function(y) ifelse(is.null(y), NA, y))
  as.numeric(x)
})

#create the right variables (1-5 instead of 0-4) for all questions
response_df1 <- response_df

response_df2 <- response_df1 %>%
  mutate(across(c("145","146","147","148","149","151","152"), ~ ifelse(. %in% 0:5, . + 1, .)))

# Remove the "dog" column from response_df
response_df3 <- response_df2[-1]

#response_df2 <- response_df1
head(response_df2)
fitGraded <- mirt(response_df3, 1, itemtype = "graded")
fitGraded

##get IRT estimates & plots
setnames(response_df3, old = c("145","146","147","148","149","151","152"), 
         new = c("q145","q146","q147","q148","q149","q151","q152"))

model <- "f1 =~ q145+q146+q147+q148+q149+q151+q152"
fitCTT <- cfa(model, data = response_df3) 

# IRT solution
summary(fitGraded)

# CTT solution
standardizedsolution(fitCTT) %>%
  filter(op == "=~") %>% dplyr::select(rhs, F1 = est.std)

#IRT parameters
params <- coef(fitGraded, IRTpars = TRUE, simplify = TRUE)
round(params$items, 2) 
M2(fitGraded, type = "M2*", calcNULL = FALSE, na.rm=TRUE) #Error: M2() statistic cannot be calculated due to too few degrees of freedom

itemfit(fitGraded,  na.rm=TRUE)

tracePlot(fitGraded, theta_range = c(-3, 12)) +
  labs(color = "Answer Options")

itemInfoPlot(fitGraded, c(1:7), facet = T)
testInfoPlot(fitGraded,theta_range = c(-5, 13), adj_factor = .5)

itemInfoPlot(fitGraded, c(1:7), theta_range = c(-5, 13), facet = F, legend = T) + 
  scale_color_brewer(palette = "Set3") 

#create the theta values for each factor:
summary(fitGraded)
item.weights <- c(0.83,1.12,1.21,1.07,1.12,1.10,1.12)
tabscores1 <- fscores(fitGraded,method='EAP',full.scores.SE=TRUE, item_weights = item.weights) #
tabscores1 <- as.data.frame(tabscores1)
tabscores3 <- fscores(fitGraded,method='MAP',full.scores.SE=TRUE, item_weights = item.weights) #
tabscores3 <- as.data.frame(tabscores3)

## Plot everything on top
plot(tabscores1$F1, tabscores1$SE_F1, xlab="Factor score", ylab="Standard error", title("F3_CCD3F"), xlim=c(-2.5, 5), ylim=c(0,1))
points(tabscores3$F1, tabscores3$SE_F1, col=c("black"), pch=4, cex=.8)
legend("topright",  legend=c("EAP", "MAP"), pch=19, bty="n", col=c("green", "black"))

#combine with dog ID - using MAP
tabscores3$number <- rownames(tabscores3)
response_df2$number <- rownames(response_df2)
nrow(tabscores3)
nrow(response_df2)
F3 <- left_join(response_df2, tabscores3, by ="number")
head(F3)
F3_CCD3F <-data.frame(dog = F3$dog, F3 = F3$F1, F3_SE = F3$SE_F1, item145 = F3$"145",item146 = F3$"146",
                      item147 = F3$"147",item148 = F3$"148",item149 = F3$"149",item151 = F3$"151",item152 = F3$"152")
head(F3_CCD3F)


x <- F3
x$na_count <- apply(x[,2:8], 1, function(x) {sum(is.na(x))})
length(x$na_count)
length(x$F1)
plot(x$SE_F1, x$na_count, xlab="SE",pch=19, ylab="NA count", title("F3_CCD3F"))
plot(x$F1, x$na_count, xlab="Factor score",pch=19, ylab="NA count", title("F3_CCD3F"))
plot(x$F1, x$SE_F1, xlab="Factor score",pch=19, ylab="SE", title("F3_CCD3F"))
nrow(F3_CCD3F)
#save factor scores
#write.table(F3_CCD3F, file= "F3_CCD3F_240213.txt", row.names =F, sep = "\t", col.names= T, quote = F ) 
####load factor 3 scores####
F3_CCD3F <- read.table("F3_CCD3F_240213.txt",header= T,   sep= "\t")
head(F3_CCD3F)


####plot distribution ####

(summary(F1_CCD3F$F1)[6])/(sd(F1_CCD3F$F1))

SD2VÄRDE <- sd(F1_CCD3F$F1)*2+summary(F1_CCD3F$F1)[4]
Q3VÄRDE <- summary(F1_CCD3F$F1)[5]
MEAN <- summary(F1_CCD3F$F1)[4]
MEDIAN <- summary(F1_CCD3F$F1)[3]
ggplot(F1_CCD3F, aes(x=F1)) + geom_histogram(color = "black", fill = "green", alpha = 0.5, binwidth = 0.1) +
  geom_vline(aes(xintercept = SD2VÄRDE), color = "red", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEAN), color = "purple", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEDIAN), color = "orange", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = Q3VÄRDE), color = "blue", linetype = "dashed", size = 0.5) + theme_minimal()

(summary(F2_CCD3F$F2)[6])/(sd(F2_CCD3F$F2))

SD2VÄRDE <- sd(F2_CCD3F$F2)*2+summary(F2_CCD3F$F2)[4]
Q3VÄRDE <- summary(F2_CCD3F$F2)[5]
MEAN <- summary(F2_CCD3F$F2)[4]
MEDIAN <- summary(F2_CCD3F$F2)[3]
ggplot(F2_CCD3F, aes(x=F2)) + geom_histogram(color = "black", fill = "green", alpha = 0.5, binwidth = 0.1) +
  geom_vline(aes(xintercept = SD2VÄRDE), color = "red", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEAN), color = "purple", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEDIAN), color = "orange", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = Q3VÄRDE), color = "blue", linetype = "dashed", size = 0.5) + theme_minimal()

(summary(F3_CCD3F$F3)[6])/(sd(F3_CCD3F$F3))

SD2VÄRDE <- sd(F3_CCD3F$F3)*2+summary(F3_CCD3F$F3)[4]
Q3VÄRDE <- summary(F3_CCD3F$F3)[5]
MEAN <- summary(F3_CCD3F$F3)[4]
MEDIAN <- summary(F3_CCD3F$F3)[3]
ggplot(F3_CCD3F, aes(x=F3)) + geom_histogram(color = "black", fill = "green", alpha = 0.5, binwidth = 0.1) +
  geom_vline(aes(xintercept = SD2VÄRDE), color = "red", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEAN), color = "purple", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEDIAN), color = "orange", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = Q3VÄRDE), color = "blue", linetype = "dashed", size = 0.5) + theme_minimal()

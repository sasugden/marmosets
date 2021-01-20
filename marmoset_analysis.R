#### PREPARE WORKSPACE ####
library(phyloseq)
library(reshape2)
library(lme4)
library(ggplot2)
library(magrittr)
library(MuMIn)
library(gridExtra)
library(scales)
library(vegan)

plot_theme <- theme(axis.text=element_text(size=10, color="black"),
                    axis.title=element_text(size=10, color="black"),
                    legend.text=element_text(size=10, color="black"),
                    legend.title=element_text(size=10, color="black"),
                    panel.grid=element_blank())

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

################################################# PART ONE: DATA PROCESSING ####################################################
#### DATA PROCESSING #####
# Read data files ####
sample_data <- read.csv("raw_data/metadata_monkeys.csv", fileEncoding="UTF-8-BOM", row.names="SampleID",
                        colClasses=rep("factor", 8))

levels(sample_data$Diet)[levels(sample_data$Diet)=="WBS"] <- 0
levels(sample_data$Diet)[levels(sample_data$Diet)=="YBS"] <- 1
sample_data$Diet <- as.integer(as.character(sample_data$Diet))

levels(sample_data$Pathology)[levels(sample_data$Pathology)=="No"] <- 0
levels(sample_data$Pathology)[levels(sample_data$Pathology)=="Yes"] <- 1
sample_data$Pathology <- as.integer(as.character(sample_data$Pathology))

levels(sample_data$Time)[levels(sample_data$Time)=="t1"] <- 1
levels(sample_data$Time)[levels(sample_data$Time)=="t2"] <- 2
levels(sample_data$Time)[levels(sample_data$Time)=="t3"] <- 3
levels(sample_data$Time)[levels(sample_data$Time)=="t4"] <- 4
sample_data$Time <- as.integer(as.character(sample_data$Time))

otu_tables <- list(
  phylum = read.csv("raw_data/abundance_table_phylum.csv", fileEncoding="UTF-8-BOM", row.names="SampleID"),
  family = read.csv("raw_data/abundance_table_family.csv", fileEncoding="UTF-8-BOM", row.names="SampleID"),
  genus = read.csv("raw_data/abundance_table_genus.csv", fileEncoding="UTF-8-BOM", row.names="SampleID"),
  otu = read.csv("raw_data/abundance_table_otu.csv", fileEncoding="UTF-8-BOM", row.names="SampleID")
)

tax_tables <- list(
  phylum = read.csv("raw_data/tax_table_phylum.csv", fileEncoding="UTF-8-BOM", row.names="TaxID"),
  family = read.csv("raw_data/tax_table_family.csv", fileEncoding="UTF-8-BOM", row.names="TaxID"),
  genus = read.csv("raw_data/tax_table_genus.csv", fileEncoding="UTF-8-BOM", row.names="TaxID"),
  otu = read.csv("raw_data/tax_table_otu.csv", fileEncoding="UTF-8-BOM", row.names="TaxID")
)

for(i in 1:length(tax_tables)){
  tax_tables[[i]] <- subset(tax_tables[[i]], rownames(tax_tables[[i]]) %in% colnames(otu_tables[[i]]))
  otu_tables[[i]] <- as.data.frame(t(otu_tables[[i]]))
}

# Generate clr-transformed, log-transformed, and raw relative abundance OTU tables. Calculate alpha diversity. ####
phyloseq_objects <- list()
otu_tables_clr <- list()
otu_tables_log <- list()
otu_tables_rel <- list()

for(i in 1:length(otu_tables)){
  temp <- otu_table(as.matrix(otu_tables[[i]]), taxa_are_rows=TRUE)
  temp <- phyloseq(otu_table(temp), sample_data(sample_data))
  
  # Calculate alpha diversity.
  phyloseq_objects[[i]] <- temp
  if(i==4){alpha_diversity <- microbiome::alpha(temp)}
  
  # CLR-transformed.
  clrdf <- microbiome::transform(temp, "clr")
  clrdf <- as.data.frame(t(otu_table(clrdf)))
  clrdf <- transform(merge(clrdf, sample_data, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
  otu_tables_clr[[i]] <- clrdf
  
  # Log-transformed.
  logdf <- microbiome::transform(temp, "log10p")
  logdf <- as.data.frame(t(otu_table(logdf)))
  logdf <- transform(merge(logdf, sample_data, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
  otu_tables_log[[i]] <- logdf
  
  # Straight relative abundance.
  reldf <- transform_sample_counts(temp, function(x) 100*x/sum(x))
  reldf <- as.data.frame(t(otu_table(reldf)))
  reldf <- transform(merge(reldf, sample_data, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
  otu_tables_rel[[i]] <- reldf
  
  rm(temp, clrdf, logdf, reldf)
}

names(otu_tables_clr) <- names(otu_tables)
names(otu_tables_log) <- names(otu_tables_log)
names(otu_tables_rel) <- names(otu_tables_rel)

# Add alpha diversity data to sample metadata object.
alpha_diversity$SampleID <- rownames(alpha_diversity)
sample_data$SampleID <- rownames(sample_data)
sample_data <- merge(sample_data, alpha_diversity, by="SampleID", all=TRUE)
rownames(sample_data) <- sample_data$SampleID
rm(alpha_diversity)

# Create a new grouping variable for plotting relative abundances (diet_time_pathology) ####
# First do this for marmosets based on their individual group.
sample_data$Group <- paste0(sample_data$Diet, "_",
                               sample_data$Time, "_",
                               sample_data$Pathology)
sample_data$Group <- factor(sample_data$Group,
                               levels=c("0_1_0", "0_2_0", "0_3_0", "0_4_0", "0_4_1",
                                        "1_1_0", "1_2_0", "1_3_0", "1_4_0", "1_4_1"))

# Also separate classifications at t3 based on whether they get sick later.
temp <- subset(sample_data, Pathology==1)
temp$Subject <- factor(temp$Subject)
temp <- temp$Subject

t1 <- subset(sample_data, Subject %in% temp & Time==3)
t2 <- subset(sample_data, !(SampleID) %in% t1$SampleID)

t1$PathologyAnticipate <- rep(1)
t2$PathologyAnticipate <- rep(0)
sample_data_anticipate <- rbind(t1, t2)

sample_data_anticipate$Group.Anticipate <- paste0(sample_data_anticipate$Diet, "_",
                                          sample_data_anticipate$Time, "_",
                                          sample_data_anticipate$PathologyAnticipate)
sample_data_anticipate$Group.Anticipate <- factor(sample_data_anticipate$Group.Anticipate,
                                          levels=c("0_1_0", "0_2_0", "0_3_0", "0_3_1", "0_4_0", "0_4_1",
                                                   "1_1_0", "1_2_0", "1_3_0", "1_3_1", "1_4_0", "1_4_1"))

sample_data <- merge(sample_data, sample_data_anticipate[,c("SampleID", "PathologyAnticipate", "Group.Anticipate")],
                     by="SampleID", all=TRUE)

rm(t1, t2, temp, sample_data_anticipate)

# Calculate Bray-Curtis inter-sample difference at each time for each diet ####
# Obtain data frame of relative OTU abundances.
colnames(otu_tables_rel[[4]])
data <- otu_tables_rel[[4]]
data[,c(74:81)] <- NULL
dist_matrix_bray <- vegdist(data, method="bray")

# Convert distance matrix to a data frame.
dist_temp <- as.data.frame(as.matrix(dist_matrix_bray))

# Null distances between identical samples.
dist_temp[dist_temp==0] <- NA

# Create strings of sample names in eahc comparison.
tests <- list(
  wbst1 = as.character(subset(sample_data, Diet==0 & Time==1)$SampleID),
  wbst2 = as.character(subset(sample_data, Diet==0 & Time==2)$SampleID),
  wbst3 = as.character(subset(sample_data, Diet==0 & Time==3)$SampleID),
  wbst4 = as.character(subset(sample_data, Diet==0 & Time==4)$SampleID),
  ybst1 = as.character(subset(sample_data, Diet==1 & Time==1)$SampleID),
  ybst2 = as.character(subset(sample_data, Diet==1 & Time==2)$SampleID),
  ybst3 = as.character(subset(sample_data, Diet==1 & Time==3)$SampleID),
  ybst4 = as.character(subset(sample_data, Diet==1 & Time==4)$SampleID)
)

# Create a place to store the data.
bray.data <- data.frame(matrix(ncol=2))
colnames(bray.data) <- c("SampleID", "BraySmall")

# Calculate the mean Bray-Curtis distance among samples in each diet/time group.
for(i in c(1:8)){
  temp <- subset(dist_temp, rownames(dist_temp) %in% tests[[i]])
  temp <- as.data.frame(t(temp))
  temp <- subset(temp, rownames(temp) %in% tests[[i]])
  
  temp <- data.frame(
    SampleID = rownames(temp),
    BraySmall = rowMeans(temp, na.rm=TRUE)
  )
  
  bray.data <- rbind(bray.data, temp)
}

# Add data to final sample data frame.
bray.data <- subset(bray.data, is.na(BraySmall)=="FALSE")
sample_data <- merge(sample_data, bray.data, by="SampleID", all=TRUE)

rm(temp, dist_temp, bray.data, data)

#### ABUNDANCE SUMMARY INFORMATION ####
# Calculate Z-score standardized abundances in a "ggplot-able" format ####
taxa_numbers <- c(nrow(otu_tables[[1]]),
                  nrow(otu_tables[[2]]),
                  nrow(otu_tables[[3]]),
                  nrow(otu_tables[[4]]))

heatmaps_z_raw_abundance <- otu_tables_rel

for(i in 1:length(heatmaps_z_raw_abundance)){
  # Subset data frame to only taxa, no metadata.
  heatmaps_z_raw_abundance[[i]] <- heatmaps_z_raw_abundance[[i]][,c(1:taxa_numbers[i])]
  
  # Calculate Z scores.
  for(j in 1:ncol(heatmaps_z_raw_abundance[[i]])){
    heatmaps_z_raw_abundance[[i]][,j] <- scale(heatmaps_z_raw_abundance[[i]][,j], center=TRUE, scale=TRUE)
  }
  heatmaps_z_raw_abundance[[i]] <- as.data.frame(t(heatmaps_z_raw_abundance[[i]]))
  
  # Add sample metadata.
  heatmaps_z_raw_abundance[[i]]$TaxID <- rownames(heatmaps_z_raw_abundance[[i]])
  heatmaps_z_raw_abundance[[i]] <- melt(heatmaps_z_raw_abundance[[i]])
  colnames(heatmaps_z_raw_abundance[[i]]) <- c("TaxID", "SampleID", "Abundance")
  
  heatmaps_z_raw_abundance[[i]] <- merge(heatmaps_z_raw_abundance[[i]],
                                       sample_data[,c("Subject","Diet","Time","Cage","Pathology",
                                                      "Sex","TwinPair","SampleID","Group","Group.Anticipate")],
                                       by="SampleID", all=TRUE)
  
  # Calculate mean relative abundance in each group.
  heatmaps_z_raw_abundance[[i]] <- heatmaps_z_raw_abundance[[i]] %>% 
    dplyr::group_by(TaxID, Group) %>% 
    dplyr::summarise(Abund=mean(Abundance))
  
  # Add taxon metadata.
  tax_tables[[i]]$TaxID <- rownames(tax_tables[[i]])
  heatmaps_z_raw_abundance[[i]] <- merge(heatmaps_z_raw_abundance[[i]], tax_tables[[i]], by="TaxID", all=TRUE)
  tax_tables[[i]]$TaxID <- NULL
  
  heatmaps_z_raw_abundance[[i]] <- subset(heatmaps_z_raw_abundance[[i]], is.na(Abund)=="FALSE")
  
}

heatmaps_z_raw_abundance[[1]]$Phylum <- factor(heatmaps_z_raw_abundance[[1]]$Phylum)
heatmaps_z_raw_abundance[[2]]$Family <- factor(heatmaps_z_raw_abundance[[2]]$Family)
heatmaps_z_raw_abundance[[3]]$Genus <- factor(heatmaps_z_raw_abundance[[3]]$Genus)
heatmaps_z_raw_abundance[[4]]$TaxID <- factor(heatmaps_z_raw_abundance[[4]]$TaxID)

# Save a table showing the mean abundances per group ####
summaries_raw_abundances <- otu_tables_rel

for(i in 1:length(summaries_raw_abundances)){
  summaries_raw_abundances[[i]] <- summaries_raw_abundances[[i]][,c(1:taxa_numbers[i])]
  summaries_raw_abundances[[i]]$SampleID <- rownames(summaries_raw_abundances[[i]])
  
  # Add sample data.
  summaries_raw_abundances[[i]] <- merge(summaries_raw_abundances[[i]],
                                         sample_data[,c("SampleID", "Group")],
                                         by="SampleID", all=TRUE)
  
  # Calculate the mean abundance per group.
  summaries_raw_abundances[[i]] <- summaries_raw_abundances[[i]] %>% 
    dplyr::group_by(Group) %>% dplyr::summarise_if(is.numeric, mean) %>%
    as.data.frame()
  
  rownames(summaries_raw_abundances[[i]]) <- summaries_raw_abundances[[i]]$Group
  summaries_raw_abundances[[i]]$Group <- NULL
  summaries_raw_abundances[[i]] <- as.data.frame(t(summaries_raw_abundances[[i]]))
  
  # Add taxon information.
  summaries_raw_abundances[[i]] <- merge(tax_tables[[i]],
                                         summaries_raw_abundances[[i]], 
                                         by=0, all=TRUE)
  
  summaries_raw_abundances[[i]] <- subset(summaries_raw_abundances[[i]], is.na(`1_3_0`)=="FALSE")
}

#### CREATE A DATA FRAME OF DELTA ABUNDANCES AND 'LOG' DELTA ABUNDANCES ####
otu_tables_delta <- otu_tables_rel
otu_tables_delta_log <- list()

# For each taxonomic level
for(q in 1:length(otu_tables_delta)){
  # Clean the data frame.
  otu_tables_delta[[q]]$SampleID <- rownames(otu_tables_delta[[q]])
  otu_tables_delta[[q]][,c("Day","Cage","Pathology","Sex","TwinPair")] <- NULL
  
  # Calculate the mean relative abundance per subject, diet, and timepoint.
  # (This is just a way to reorganize the data frame).
  temp <- otu_tables_delta[[q]] %>%
    dplyr::group_by(Subject, Diet, Time) %>% 
    dplyr::summarise_if(is.numeric, mean) %>%
    as.data.frame()
  
  # Calculate differences in relative abundances per timepoint.
  temp2 <- data.frame(matrix())
  for(i in c(4:ncol(temp))){
    # Create a vector of numbers indicating where data should be stored.
    vector <- seq(from=1, to=46, by=3)
      
    for(j in 1:nlevels(temp$Subject)){
      # Subtract time 1 from time 2.
      temp2[vector[j],i-3] <- subset(temp, Subject==levels(temp$Subject)[j] & Time %in% c(1,2))[2,i]-
        subset(temp, Subject==levels(temp$Subject)[j] & Time %in% c(1,2))[1,i]
      
      # Subtract time 2 from time 3.
      temp2[vector[j]+1,i-3] <- subset(temp, Subject==levels(temp$Subject)[j] & Time %in% c(2,3))[2,i]-
        subset(temp, Subject==levels(temp$Subject)[j] & Time %in% c(2,3))[1,i]
      
      # Subtract time 3 from time 4.
      temp2[vector[j]+2,i-3] <- subset(temp, Subject==levels(temp$Subject)[j] & Time %in% c(3,4))[2,i]-
        subset(temp, Subject==levels(temp$Subject)[j] & Time %in% c(3,4))[1,i]}
  }
  
  # Rename columns and clean data.
  colnames(temp2) <- colnames(temp)[c(4:ncol(temp))]
  temp2$Subject <- rep(levels(temp$Subject), each=3)
  temp2$Difference <- rep(c("1_2","2_3","3_4"), nlevels(temp$Subject))
  
  # Label the differences between times 1-2 and 2-3 as "not pathology." 
  temp2a <- subset(temp2, Difference %in% c("1_2", "2_3"))
  temp2a <- merge(temp2a,
                  subset(sample_data, Time==4)[,c("Subject","Diet","TwinPair","Cage")], by="Subject", all=FALSE)
  temp2a$Pathology <- rep(0)
  
  # Separate the differences between times 3-4 and add the appropriate pathology label.
  temp2b <- subset(temp2, Difference %in% c("3_4"))
  temp2b <- merge(temp2b,
                  subset(sample_data, Time==4)[,c("Subject","Diet","TwinPair","Cage","Pathology")], by="Subject", all=FALSE)
  
  # Put the data back together.
  temp2 <- rbind(temp2a, temp2b)
 
  # Calculate the log-transform of the absolute value of each "delta abundance."
  temp3 <- temp2
  temp3[,c(2:(taxa_numbers[q]+1))] <- abs(temp3[,c(2:(taxa_numbers[q]+1))])
  temp3[,c(2:(ncol(temp)-2))] <- log(temp3[,c(2:(ncol(temp)-2))]+0.01, base=10)
  
  otu_tables_delta[[q]] <- temp2
  otu_tables_delta_log[[q]] <- temp3
  
  rm(temp, temp2, temp2a, temp2b, temp3, vector)
}

# Table with delta abundances per group #
summaries_delta_abundances <- otu_tables_delta

for(i in 1:length(summaries_delta_abundances)){
  # Create grouping variable.
  summaries_delta_abundances[[i]]$Variable <- paste0(summaries_delta_abundances[[i]]$Diet, ".",
                                                     summaries_delta_abundances[[i]]$Difference, ".",
                                                     summaries_delta_abundances[[i]]$Pathology)
  summaries_delta_abundances[[i]][,c("Subject","Diet","Difference","Pathology","TwinPair","Cage")] <- NULL
  
  # Calculate mean abundance.
  summaries_delta_abundances[[i]] <- summaries_delta_abundances[[i]] %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarise_if(is.numeric, mean) %>%
    as.data.frame()
  
  # Rename grouping variables for clarity.
  summaries_delta_abundances[[i]]$Variable[summaries_delta_abundances[[i]]$Variable=="0.1_2.0"] <- "0_1v2_0"
  summaries_delta_abundances[[i]]$Variable[summaries_delta_abundances[[i]]$Variable=="0.2_3.0"] <- "0_2v3_0"
  summaries_delta_abundances[[i]]$Variable[summaries_delta_abundances[[i]]$Variable=="0.3_4.0"] <- "0_3v4_0"
  summaries_delta_abundances[[i]]$Variable[summaries_delta_abundances[[i]]$Variable=="0.3_4.1"] <- "0_3v4_1"
  summaries_delta_abundances[[i]]$Variable[summaries_delta_abundances[[i]]$Variable=="1.1_2.0"] <- "1_1v2_0"
  summaries_delta_abundances[[i]]$Variable[summaries_delta_abundances[[i]]$Variable=="1.2_3.0"] <- "1_2v3_0"
  summaries_delta_abundances[[i]]$Variable[summaries_delta_abundances[[i]]$Variable=="1.3_4.0"] <- "1_3v4_0"
  summaries_delta_abundances[[i]]$Variable[summaries_delta_abundances[[i]]$Variable=="1.3_4.1"] <- "1_3v4_1"
  
  # Clean data
  rownames(summaries_delta_abundances[[i]]) <- summaries_delta_abundances[[i]]$Variable
  summaries_delta_abundances[[i]]$Variable <- NULL
  
  summaries_delta_abundances[[i]] <- as.data.frame(t(summaries_delta_abundances[[i]]))
  
  # Add taxonomic information
  summaries_delta_abundances[[i]] <- merge(tax_tables[[i]],
                                           summaries_delta_abundances[[i]],
                                           by=0, all=TRUE)
  summaries_delta_abundances[[i]] <- subset(summaries_delta_abundances[[i]], is.na(`0_1v2_0`)=="FALSE")
}



################################################# PART TWO: MODELS AND HEAT MAPS ################################################
#### LOG-TRANSFORMED ABUNDANCE MODELS ####
# Run models that predict log-transformed abundances as a function of diet, pathology, and time.
otu_names <- forcats::fct_inorder(tax_tables[[4]]$PlotName)

models_actual_abundance_log <- list()

for(i in c(1:4)){ 
  model_data <- otu_tables_log[[i]]
  
  # Rescale data so coefficients are comparable.
  model_data[,c(1:(ncol(model_data)-8))] <- scale(model_data[,c(1:(ncol(model_data)-8))])
  
  # Create a data frame to store coefficients.
  models_actual_abundance_log[[i]] <- data.frame(matrix(ncol=5))
  colnames(models_actual_abundance_log[[i]]) <- c("TaxID", "Variable", "Lower", "Coefficient", "Upper")
  
  for(p in c(1:(ncol(model_data)-8))){ # For each taxon in the OTU table
    
    model <- lmer(model_data[,p] ~ Diet + Time + Pathology + (1|TwinPair/Subject), # Construct a mixed-effect model
                  data=model_data, na.action="na.fail", REML=FALSE)
    
    model_dredge <- dredge(model) # Dredge all subsets.
    
    # If there is only one top model, no coefficient-averaging is required.
    # Just identify the top model, extract the coefficients, and calculate confidence intervals.
    if(nrow(subset(as.data.frame(model_dredge), delta < 2))==1){
      x <- as.data.frame(model_dredge)
      x <- subset(x, delta < 2)
      x <- x[,c(2:4)] 
      x <- as.data.frame(t(x))
      x <- subset(x, is.na(x[,1])=="FALSE")
      vars <- rownames(x) # Identify the variables present in the top model.
      
      if(length(vars)>0){ # Recreate the model.
        if(length(vars)==1){model <- lmer(model_data[,p] ~ model_data[,vars] + (1|TwinPair/Subject),
                                          data=model_data, na.action="na.fail", REML=FALSE)}
        if(length(vars)==2){model <- lmer(model_data[,p] ~ model_data[,vars[1]] + model_data[,vars[2]] + (1|TwinPair/Subject),
                                          data=model_data, na.action="na.fail", REML=FALSE)}
        if(length(vars)==3){model <- lmer(model_data[,p] ~ model_data[,vars[1]] + model_data[,vars[2]] + model_data[,vars[3]] + (1|TwinPair/Subject),
                                          data=model_data, na.action="na.fail", REML=FALSE)}
        temp <- confint.merMod(model, level=0.95, method="Wald")
        
        fill <- data.frame(
          Lower=temp[,1][c(5:(length(vars)+4))],
          Coefficient=model@beta[c(2:(length(vars)+1))],
          Upper=temp[,2][c(5:(length(vars)+4))])
        fill$Variable <- vars
        
        if(length(vars)==1){fill$TaxID <- rep(colnames(model_data)[p], 1)}
        if(length(vars)==2){fill$TaxID <- rep(colnames(model_data)[p], 2)}
        if(length(vars)==3){fill$TaxID <- rep(colnames(model_data)[p], 3)}
      }
      
      if(length(vars)==0){
        fill <- data.frame(
          Lower=rep(NA, 3),
          Coefficient=rep(NA, 3),
          Upper=rep(NA, 3),
          Variable=c("Diet", "Pathology", "Time"),
          TaxID=rep(colnames(model_data)[p], 3)
        )
      }
      
      model_coeff <- fill[,c("TaxID", "Variable", "Lower", "Coefficient", "Upper")]
    }
    
    # If there is more than one top model - average coefficients
    if(nrow(subset(as.data.frame(model_dredge), delta < 2)) > 1){
      model_coeff <- model.avg(model_dredge, subset = delta < 2)
      
      conf <- confint(model_coeff)
      conf <- subset(conf, rownames(conf) !="(Intercept)")
      
      model_coeff <- as.data.frame(model_coeff$coefficients)
      model_coeff <- subset(model_coeff, rownames(model_coeff)=="full")
      model_coeff <- as.data.frame(t(model_coeff))
      model_coeff <- subset(model_coeff, rownames(model_coeff) != "(Intercept)")
      
      model_coeff <- data.frame(
        TaxID=colnames(model_data)[p],
        Variable=rownames(model_coeff),
        Lower=conf[,1],
        Coefficient=model_coeff[,1],
        Upper=conf[,2]
      )
    }
    
    models_actual_abundance_log[[i]] <- rbind(models_actual_abundance_log[[i]],
                                              model_coeff)
    
  }
  
  # Use this hack to label coefficients as "significant" (confidence intervals do not overlap zero) or not.
  models_actual_abundance_log[[i]]$Significance <- models_actual_abundance_log[[i]]$Lower*models_actual_abundance_log[[i]]$Upper
  models_actual_abundance_log[[i]]$Significance[models_actual_abundance_log[[i]]$Significance > 0] <- 0.05
  models_actual_abundance_log[[i]]$Significance[models_actual_abundance_log[[i]]$Significance < 0] <- 1
  models_actual_abundance_log[[i]]$Confidence <- cut(models_actual_abundance_log[[i]]$Significance,
                                                     breaks=c(0,0.07,Inf),
                                                     labels=c("**", " "))
  
  models_actual_abundance_log[[i]] <- subset(models_actual_abundance_log[[i]],
                                             is.na(TaxID)=="FALSE")
  
  rm(model_data, model_dredge, model, fill, conf, x, temp)
}

for(i in 1:length(models_actual_abundance_log)){
  # Add taxon information to each model result.
  tax_tables[[i]]$TaxID <- rownames(tax_tables[[i]])
  models_actual_abundance_log[[i]] <- merge(models_actual_abundance_log[[i]], tax_tables[[i]], by="TaxID", all=TRUE)
  tax_tables[[i]]$TaxID <- NULL
  
  models_actual_abundance_log[[i]] <- subset(models_actual_abundance_log[[i]],
                                             is.na(Coefficient)=="FALSE")
  
  # Complete the data frame for missing coefficients (i.e, coefficients dropped from AIC model selection)
  if(i==1){
    models_actual_abundance_log[[1]]$Phylum <- factor(models_actual_abundance_log[[1]]$Phylum)
    models_actual_abundance_log[[i]] = models_actual_abundance_log[[i]] %>%
      dplyr::group_by(Phylum) %>%
      tidyr::complete(Variable = c("Diet", "Time", "Pathology"))
    models_actual_abundance_log[[i]]$Level <- rep("Phylum")
  }
  if(i==2){
    models_actual_abundance_log[[2]]$Family <- factor(models_actual_abundance_log[[2]]$Family)
    models_actual_abundance_log[[i]] = models_actual_abundance_log[[i]] %>%
      dplyr::group_by(Family) %>%
      tidyr::complete(Variable = c("Diet", "Time", "Pathology"))
    models_actual_abundance_log[[i]]$Level <- rep("Family")
  }
  
  if(i==3){
    models_actual_abundance_log[[3]]$Genus <- factor(models_actual_abundance_log[[3]]$Genus)
    models_actual_abundance_log[[i]] = models_actual_abundance_log[[i]] %>%
      dplyr::group_by(Genus) %>%
      tidyr::complete(Variable = c("Diet", "Time", "Pathology"))
    models_actual_abundance_log[[i]]$Level <- rep("Genus")
  }
  if(i==4){
    models_actual_abundance_log[[4]]$PlotName <- factor(models_actual_abundance_log[[4]]$PlotName)
    models_actual_abundance_log[[i]] = models_actual_abundance_log[[i]] %>%
      dplyr::group_by(PlotName) %>%
      tidyr::complete(Variable = c("Diet", "Time", "Pathology"))
    models_actual_abundance_log[[i]]$Level <- rep("OTU")
  }
  
}

# Take a look at the final plot to ensure it worked.
ggplot(models_actual_abundance_log[[3]], aes(x=Variable, y=Genus)) + 
  geom_tile(aes(fill=Coefficient)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="white") +
  geom_text(aes(label=Confidence)) +
  scale_y_discrete(limits=rev(levels(models_actual_abundance_log[[3]]$Genus))) +
  labs(x="Variable") +
  plot_theme 

#### DELTA ABUNDANCE MODELS ####
# This code does the same thing as above, using the "delta abundance" data frames.
models_delta_abundance <- list()

for(i in 1:length(otu_tables_delta_log)){ # For each OTU table...
  model_data <- otu_tables_delta_log[[i]]
  
  model_data$Difference[model_data$Difference=="1_2"] <- 0
  model_data$Difference[model_data$Difference=="2_3"] <- 1
  model_data$Difference[model_data$Difference=="3_4"] <- 2
  model_data$Difference <- as.integer(model_data$Difference)
  
  model_data[,c(2:(ncol(model_data)-5))] <- scale(model_data[,c(2:(ncol(model_data)-5))])
  
  models_delta_abundance[[i]] <- data.frame(matrix(ncol=5))
  colnames(models_delta_abundance[[i]]) <- c("TaxID", "Variable", "Lower", "Coefficient", "Upper")
  
  for(p in c(2:(ncol(model_data)-5))){ # For each taxon in that OTU table...
    
    model <- lmer(model_data[,p] ~ Diet + Difference + Pathology + (1|TwinPair/Subject),
                  data=model_data, na.action="na.fail", REML=FALSE)
    
    model_dredge <- dredge(model)
    
    if(nrow(subset(as.data.frame(model_dredge), delta < 2))==1){
      x <- as.data.frame(model_dredge)
      x <- subset(x, delta < 2)
      x <- x[,c(2:4)] 
      x <- as.data.frame(t(x))
      x <- subset(x, is.na(x[,1])=="FALSE")
      vars <- rownames(x) # Identify variables present in the top model.
      
      if(length(vars)>0){ # Recreate the model.
        if(length(vars)==1){model <- lmer(model_data[,p] ~ model_data[,vars] + (1|TwinPair/Subject),
                                          data=model_data, na.action="na.fail", REML=FALSE)}
        if(length(vars)==2){model <- lmer(model_data[,p] ~ model_data[,vars[1]] + model_data[,vars[2]] + (1|TwinPair/Subject),
                                          data=model_data, na.action="na.fail", REML=FALSE)}
        if(length(vars)==3){model <- lmer(model_data[,p] ~ model_data[,vars[1]] + model_data[,vars[2]] + model_data[,vars[3]] + (1|TwinPair/Subject),
                                          data=model_data, na.action="na.fail", REML=FALSE)}
        temp <- confint.merMod(model, level=0.95, method="Wald")
        
        fill <- data.frame(
          Lower=temp[,1][c(5:(length(vars)+4))],
          Coefficient=model@beta[c(2:(length(vars)+1))],
          Upper=temp[,2][c(5:(length(vars)+4))])
        fill$Variable <- vars
        
        if(length(vars)==1){fill$TaxID <- rep(colnames(model_data)[p], 1)}
        if(length(vars)==2){fill$TaxID <- rep(colnames(model_data)[p], 2)}
        if(length(vars)==3){fill$TaxID <- rep(colnames(model_data)[p], 3)}
      }
      
      if(length(vars)==0){
        fill <- data.frame(
          Lower=rep(NA, 3),
          Coefficient=rep(NA, 3),
          Upper=rep(NA, 3),
          Variable=c("Diet", "Difference", "Pathology"),
          TaxID=rep(colnames(model_data)[p], 3)
        )
      }
      
      model_coeff <- fill[,c("TaxID", "Variable", "Lower", "Coefficient", "Upper")]
    }
    
    # If there is more than one top model - average coefficients
    if(nrow(subset(as.data.frame(model_dredge), delta < 2)) > 1){
      model_coeff <- model.avg(model_dredge, subset = delta < 2)
      
      conf <- confint(model_coeff)
      conf <- subset(conf, rownames(conf) !="(Intercept)")
      
      model_coeff <- as.data.frame(model_coeff$coefficients)
      model_coeff <- subset(model_coeff, rownames(model_coeff)=="full")
      model_coeff <- as.data.frame(t(model_coeff))
      model_coeff <- subset(model_coeff, rownames(model_coeff) != "(Intercept)")
      
      model_coeff <- data.frame(
        TaxID=colnames(model_data)[p],
        Variable=rownames(model_coeff),
        Lower=conf[,1],
        Coefficient=model_coeff[,1],
        Upper=conf[,2]
      )
    }
    
    models_delta_abundance[[i]] <- rbind(models_delta_abundance[[i]],
                                         model_coeff)
    
  }
  
  
  models_delta_abundance[[i]]$Significance <- models_delta_abundance[[i]]$Lower*models_delta_abundance[[i]]$Upper
  models_delta_abundance[[i]]$Significance[models_delta_abundance[[i]]$Significance > 0] <- 0.05
  models_delta_abundance[[i]]$Significance[models_delta_abundance[[i]]$Significance < 0] <- 1
  models_delta_abundance[[i]]$Confidence <- cut(models_delta_abundance[[i]]$Significance,
                                                breaks=c(0,0.07,Inf),
                                                labels=c("**", " "))
  
  models_delta_abundance[[i]] <- subset(models_delta_abundance[[i]],
                                        is.na(TaxID)=="FALSE")
  
  rm(model_data, model_dredge, model, fill, conf, x, temp)
}

for(i in 1:length(models_delta_abundance)){
  tax_tables[[i]]$TaxID <- rownames(tax_tables[[i]])
  models_delta_abundance[[i]] <- merge(models_delta_abundance[[i]], tax_tables[[i]], by="TaxID", all=TRUE)
  tax_tables[[i]]$TaxID <- NULL
  
  models_delta_abundance[[i]] <- subset(models_delta_abundance[[i]],
                                        is.na(Coefficient)=="FALSE")
  
  if(i==1){
    models_delta_abundance[[1]]$Phylum <- factor(models_delta_abundance[[1]]$Phylum)
    models_delta_abundance[[i]] = models_delta_abundance[[i]] %>%
      dplyr::group_by(Phylum) %>%
      tidyr::complete(Variable = c("Diet", "Difference", "Pathology"))
    models_delta_abundance[[i]]$Level <- rep("Phylum")
  }
  if(i==2){
    models_delta_abundance[[2]]$Family <- factor(models_delta_abundance[[2]]$Family)
    models_delta_abundance[[i]] = models_delta_abundance[[i]] %>%
      dplyr::group_by(Family) %>%
      tidyr::complete(Variable = c("Diet", "Difference", "Pathology"))
    models_delta_abundance[[i]]$Level <- rep("Family")
  }
  
  if(i==3){
    models_delta_abundance[[3]]$Genus <- factor(models_delta_abundance[[3]]$Genus)
    models_delta_abundance[[i]] = models_delta_abundance[[i]] %>%
      dplyr::group_by(Genus) %>%
      tidyr::complete(Variable = c("Diet", "Difference", "Pathology"))
    models_delta_abundance[[i]]$Level <- rep("Genus")
  }
  if(i==4){
    models_delta_abundance[[4]]$PlotName <- factor(models_delta_abundance[[4]]$PlotName)
    models_delta_abundance[[i]] = models_delta_abundance[[i]] %>%
      dplyr::group_by(PlotName) %>%
      tidyr::complete(Variable = c("Diet", "Difference", "Pathology"))
    models_delta_abundance[[i]]$Level <- rep("OTU")
  }
}

# Plot to make sure it worked.
ggplot(models_delta_abundance[[1]], aes(x=Variable, y=Phylum)) + 
  geom_tile(aes(fill=Coefficient)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="white") +
  geom_text(aes(label=Confidence)) +
  scale_y_discrete(limits=rev(levels(models_delta_abundance[[2]]$Phylum))) +
  labs(x="Variable") +
  plot_theme 

#### CONVERT LOG-TRANSFORMED MODEL RESULTS INTO A SINGLE DATA FRAME FOR PLOTTING ####
megaplot_actual_abundance <- list()

for(i in 1:length(models_actual_abundance_log)){
  # First, remove "Uncl" because they show up at higher taxonomic levels.
  df <- models_actual_abundance_log[[i]]
  if(i==2){df <- df[-grep("Uncl", df$Family),]}
  if(i==3){df <- df[-grep("Uncl", df$Genus),]}
  
  # Then, remove taxa where no variable are significant.
  if(i==1){
    temp <- subset(df, Phylum %in% c("Actinobacteria","Bacteroidetes","Firmicutes","Fusobacteria","Proteobacteria","Tenericutes"))
    temp$Phylum <- factor(temp$Phylum)
    df <- subset(df, Phylum %in% temp$Phylum)
    df$PlotName <- df$Phylum
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), " "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]}
  if(i==2){
    temp <- subset(df, Significance==0.05)
    temp$Family <- factor(temp$Family)
    df <- subset(df, Family %in% temp$Family)
    df$PlotName <- df$Family
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "  "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
  }
  if(i==3){
    temp <- subset(df, Significance==0.05)
    temp$Genus <- factor(temp$Genus)
    df <- subset(df, Genus %in% temp$Genus)
    df$PlotName <- df$Genus
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "   "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
  }
  
  if(i==4){
    temp <- subset(df, Significance==0.05)
    temp$PlotName <- factor(temp$PlotName)
    
    x <- as.data.frame(as.character(otu_names))
    x <- subset(x, x[,1] %in% temp$PlotName)
    x[,1] <- factor(x[,1])
    x <- forcats::fct_inorder(x[,1])
    
    df <- subset(df, PlotName %in% temp$PlotName)
    df$PlotName <- factor(df$PlotName, levels=x)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "    "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
    
  }
  
  # Save it!
  megaplot_actual_abundance[[i]] <- df
}

megaplot_actual_abundance <- rbind(as.data.frame(megaplot_actual_abundance[[1]]),
                                   as.data.frame(megaplot_actual_abundance[[2]]),
                                   as.data.frame(megaplot_actual_abundance[[3]]),
                                   as.data.frame(megaplot_actual_abundance[[4]]))

#### CONVERT DELTA ABUNDANCE RESULTS INTO A SINGLE DATA FRAME FOR PLOTTING #####
megaplot_delta_abundance <- list()

for(i in 1:length(models_delta_abundance)){
  # First, remove "Uncl" because they show up at higher taxonomic levels.
  df <- models_delta_abundance[[i]]
  if(i==2){df <- df[-grep("Uncl", df$Family),]}
  if(i==3){df <- df[-grep("Uncl", df$Genus),]}
  
  # Then, remove taxa where no variable are significant.
  if(i==1){
    temp <- subset(df, Phylum %in% c("Actinobacteria","Bacteroidetes","Firmicutes","Fusobacteria","Proteobacteria","Tenericutes"))
    temp$Phylum <- factor(temp$Phylum)
    df <- subset(df, Phylum %in% temp$Phylum)
    df$PlotName <- df$Phylum
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), " "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]}
  if(i==2){
    temp <- subset(df, Significance==0.05)
    temp$Family <- factor(temp$Family)
    df <- subset(df, Family %in% temp$Family)
    df$PlotName <- df$Family
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "  "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
  }
  if(i==3){
    temp <- subset(df, Significance==0.05)
    temp$Genus <- factor(temp$Genus)
    df <- subset(df, Genus %in% temp$Genus)
    df$PlotName <- df$Genus
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "   "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
  }
  
  if(i==4){
    temp <- subset(df, PlotName %in% megaplot_actual_abundance$PlotName)
    temp$PlotName <- factor(temp$PlotName)
    
    x <- as.data.frame(as.character(otu_names))
    x <- subset(x, x[,1] %in% temp$PlotName)
    x[,1] <- factor(x[,1])
    x <- forcats::fct_inorder(x[,1])
    
    df <- subset(df, PlotName %in% temp$PlotName)
    df$PlotName <- factor(df$PlotName, levels=x)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "    "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
    
  }
  
  # Save it!
  megaplot_delta_abundance[[i]] <- df
}

megaplot_delta_abundance <- rbind(as.data.frame(megaplot_delta_abundance[[1]]),
                                  as.data.frame(megaplot_delta_abundance[[2]]),
                                  as.data.frame(megaplot_delta_abundance[[3]]),
                                  as.data.frame(megaplot_delta_abundance[[4]]))

#### PRODUCE PLOTS ####
plot1 <- ggplot(megaplot_actual_abundance, aes(x=Variable, y=PlotName)) + 
  geom_tile(aes(fill=Coefficient)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="white") +
  geom_text(aes(label=Confidence), vjust=0.75, hjust=0.5) +
  scale_y_discrete(limits=rev(levels(megaplot_actual_abundance$PlotName)), expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  labs(x="Variable", y="Taxon") +
  theme_bw() + plot_theme + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                                  axis.text.y=element_text(size=8),
                                  panel.grid=element_blank())

plot2 <- ggplot(megaplot_delta_abundance, aes(x=Variable, y=PlotName)) + 
  geom_tile(aes(fill=Coefficient)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="white") +
  geom_text(aes(label=Confidence), vjust=0.75, hjust=0.5) +
  scale_y_discrete(limits=rev(levels(megaplot_delta_abundance$PlotName)), expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  labs(x="Variable", y="Taxon") +
  theme_bw() + plot_theme + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                                  axis.text.y=element_text(size=8),
                                  panel.grid=element_blank())

grid.arrange(plot1, plot2, ncol=2)

# Create a data frame that identifies which taxa with significant "delta" coefficients
# do not appear in the "log" coefficients data frame.
df <- data.frame(
  delta_model_ID = levels(megaplot_delta_abundance$PlotName),
  in_log = levels(megaplot_delta_abundance$PlotName) %in% levels(megaplot_actual_abundance$PlotName)
)

# Use this data frame to identify taxa in the 'delta' models that don't appear in the 'log' models.
# Then, go back and reconstruct the two 'megaplot' data frames, based on these taxa.
# In this case, the taxa are: Bacillales_Incertae Sedis XI, Coriobacteriaceae, Porphyromonadaceae, and Gemella
# Manually ensure these are retained.
rm(df)

#### RECONSTRUCT LOG ABUNDANCE MODEL MEGA DATA FRAME ####
megaplot_actual_abundance <- list()

for(i in 1:length(models_actual_abundance_log)){
  # First, remove "Uncl" because they show up at higher taxonomic levels.
  df <- models_actual_abundance_log[[i]]
  if(i==2){df <- df[-grep("Uncl", df$Family),]}
  if(i==3){df <- df[-grep("Uncl", df$Genus),]}
  
  # Then, remove taxa where no variable are significant.
  if(i==1){
    temp <- subset(df, Phylum %in% c("Actinobacteria","Bacteroidetes","Firmicutes","Fusobacteria","Proteobacteria","Tenericutes"))
    temp$Phylum <- factor(temp$Phylum)
    df <- subset(df, Phylum %in% temp$Phylum)
    df$PlotName <- df$Phylum
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), " "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]}
  if(i==2){
    temp <- subset(df, Significance==0.05 | Family %in% c("Bacillales_Incertae Sedis XI",
                                                          "Coriobacteriaceae",
                                                          "Porphyromonadaceae"))
    temp$Family <- factor(temp$Family)
    df <- subset(df, Family %in% temp$Family)
    df$PlotName <- df$Family
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "  "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
  }
  if(i==3){
    temp <- subset(df, Significance==0.05 | Genus %in% c("Gemella"))
    temp$Genus <- factor(temp$Genus)
    df <- subset(df, Genus %in% temp$Genus)
    df$PlotName <- df$Genus
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "   "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
  }
  
  if(i==4){
    temp <- subset(df, Significance==0.05)
    temp$PlotName <- factor(temp$PlotName)
    
    x <- as.data.frame(as.character(otu_names))
    x <- subset(x, x[,1] %in% temp$PlotName)
    x[,1] <- factor(x[,1])
    x <- forcats::fct_inorder(x[,1])
    
    df <- subset(df, PlotName %in% temp$PlotName)
    df$PlotName <- factor(df$PlotName, levels=x)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "    "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
    
  }
  
  # Save it!
  megaplot_actual_abundance[[i]] <- df
}

megaplot_actual_abundance <- rbind(as.data.frame(megaplot_actual_abundance[[1]]),
                                   as.data.frame(megaplot_actual_abundance[[2]]),
                                   as.data.frame(megaplot_actual_abundance[[3]]),
                                   as.data.frame(megaplot_actual_abundance[[4]]))

#### RECONSTRUCT 'DELTA' MEGA DATA FRAME ####
megaplot_delta_abundance <- list()

for(i in 1:length(models_delta_abundance)){
  # First, remove "Uncl" because they show up at higher taxonomic levels.
  df <- models_delta_abundance[[i]]
  if(i==2){df <- df[-grep("Uncl", df$Family),]}
  if(i==3){df <- df[-grep("Uncl", df$Genus),]}
  
  # Then, remove taxa where no variable are significant.
  if(i==1){
    temp <- subset(df, Phylum %in% c("Actinobacteria","Bacteroidetes","Firmicutes","Fusobacteria","Proteobacteria","Tenericutes"))
    temp$Phylum <- factor(temp$Phylum)
    df <- subset(df, Phylum %in% temp$Phylum)
    df$PlotName <- df$Phylum
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), " "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]}
  if(i==2){
    temp <- subset(df, Family %in% megaplot_actual_abundance$PlotName)
    temp$Family <- factor(temp$Family)
    df <- subset(df, Family %in% temp$Family)
    df$PlotName <- df$Family
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "  "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
  }
  if(i==3){
    temp <- subset(df, Genus %in% megaplot_actual_abundance$PlotName)
    temp$Genus <- factor(temp$Genus)
    df <- subset(df, Genus %in% temp$Genus)
    df$PlotName <- df$Genus
    df$PlotName <- factor(df$PlotName)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "   "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
  }
  
  if(i==4){
    temp <- subset(df, PlotName %in% megaplot_actual_abundance$PlotName)
    temp$PlotName <- factor(temp$PlotName)
    
    x <- as.data.frame(as.character(otu_names))
    x <- subset(x, x[,1] %in% temp$PlotName)
    x[,1] <- factor(x[,1])
    x <- forcats::fct_inorder(x[,1])
    
    df <- subset(df, PlotName %in% temp$PlotName)
    df$PlotName <- factor(df$PlotName, levels=x)
    df$PlotName <- factor(df$PlotName, levels=c(levels(df$PlotName), "    "))
    df <- df[,c("TaxID","Level","PlotName","Variable","Lower","Coefficient","Upper","Significance","Confidence")]
    df <- df[order(df$PlotName),]
    
  }
  
  # Save it!
  megaplot_delta_abundance[[i]] <- df
}

megaplot_delta_abundance <- rbind(as.data.frame(megaplot_delta_abundance[[1]]),
                                  as.data.frame(megaplot_delta_abundance[[2]]),
                                  as.data.frame(megaplot_delta_abundance[[3]]),
                                  as.data.frame(megaplot_delta_abundance[[4]]))

#### COMPARE PLOTS AGAIN ####
# The same taxa should now appear in both plots.
plot1 <- ggplot(megaplot_actual_abundance, aes(x=Variable, y=PlotName)) + 
  geom_tile(aes(fill=Coefficient)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="white") +
  geom_text(aes(label=Confidence), vjust=0.75, hjust=0.5) +
  scale_y_discrete(limits=rev(levels(megaplot_actual_abundance$PlotName)), expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  labs(x="Variable", y="Taxon") +
  theme_bw() + plot_theme + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                                  axis.text.y=element_text(size=8),
                                  panel.grid=element_blank())

plot2 <- ggplot(megaplot_delta_abundance, aes(x=Variable, y=PlotName)) + 
  geom_tile(aes(fill=Coefficient)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="white") +
  geom_text(aes(label=Confidence), vjust=0.75, hjust=0.5) +
  scale_y_discrete(limits=rev(levels(megaplot_delta_abundance$PlotName)), expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  labs(x="Variable", y="Taxon") +
  theme_bw() + plot_theme + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                                  axis.text.y=element_text(size=8),
                                  panel.grid=element_blank())

grid.arrange(plot1, plot2, ncol=2)

rm(plot1, plot2)

#### SUBSET TAXA WITH THE SMALLEST COEFFICIENTS, ESPECIALLY TAXA AT THE OTU LEVEL ####
# Keep all the phyla.
log1 <- subset(megaplot_actual_abundance, Level=="Phylum")
log1$PlotName <- factor(log1$PlotName)
log1$PlotName <- forcats::fct_inorder(log1$PlotName)
log1$PlotName <- factor(log1$PlotName, levels=c(levels(log1$PlotName), " "))

delta1 <- subset(megaplot_delta_abundance, Level=="Phylum")
delta1$PlotName <- factor(delta1$PlotName)
delta1$PlotName <- forcats::fct_inorder(delta1$PlotName)
delta1$PlotName <- factor(delta1$PlotName, levels=c(levels(delta1$PlotName), " "))

# Check on families that have low coefficients in either column.
log2 <- subset(megaplot_actual_abundance, Level=="Family")
delta2 <- subset(megaplot_delta_abundance, Level=="Family")

log2x <- subset(log2, abs(Coefficient) > 0.5 | PlotName=="Bifidobacteriaceae")
log2x$TaxID <- factor(log2x$TaxID)

delta2x <- subset(delta2, abs(Coefficient) > 0.5 | PlotName=="Bifidobacteriaceae")
delta2x$TaxID <- factor(delta2x$TaxID)

log2b <- subset(log2, TaxID %in% log2x$TaxID | TaxID %in% delta2x$TaxID)
log2b$PlotName <- factor(log2b$PlotName)
log2b$PlotName <- forcats::fct_inorder(log2b$PlotName)
log2b$PlotName <- factor(log2b$PlotName, levels=c(levels(log2b$PlotName), "  "))

delta2b <- subset(delta2, PlotName %in% log2b$PlotName)
delta2b$PlotName <- factor(delta2b$PlotName)
delta2b$PlotName <- forcats::fct_inorder(delta2b$PlotName)
delta2b$PlotName <- factor(delta2b$PlotName, levels=c(levels(delta2b$PlotName), "  "))

rm(log2, delta2, log2x, delta2x)

# Check on genera that have low coefficients in either column.
log3 <- subset(megaplot_actual_abundance, Level=="Genus")
delta3 <- subset(megaplot_delta_abundance, Level=="Genus")

log3x <- subset(log3, abs(Coefficient) > 0.5 | PlotName=="Bifidobacterium")
log3x$TaxID <- factor(log3x$TaxID)

delta3x <- subset(delta3, abs(Coefficient) > 0.5 | PlotName=="Bifidobacterium")
delta3x$TaxID <- factor(delta3x$TaxID)

log3b <- subset(log3, TaxID %in% log3x$TaxID | TaxID %in% delta3x$TaxID)
log3b$PlotName <- factor(log3b$PlotName)
log3b$PlotName <- forcats::fct_inorder(log3b$PlotName)
log3b$PlotName <- factor(log3b$PlotName, levels=c(levels(log3b$PlotName), "   "))

delta3b <- subset(delta3, PlotName %in% log3b$PlotName)
delta3b$PlotName <- factor(delta3b$PlotName)
delta3b$PlotName <- forcats::fct_inorder(delta3b$PlotName)
delta3b$PlotName <- factor(delta3b$PlotName, levels=c(levels(delta3b$PlotName), "   "))

rm(log3, delta3, log3x, delta3x)

# Check on OTUs that have low coefficients in either model.
log4 <- subset(megaplot_actual_abundance, Level=="OTU")
delta4 <- subset(megaplot_delta_abundance, Level=="OTU")

log4x <- subset(log4, abs(Coefficient) > 0.5)
log4x$TaxID <- factor(log4x$TaxID)

delta4x <- subset(delta4, abs(Coefficient) > 0.5)
delta4x$TaxID <- factor(delta4x$TaxID)

log4b <- subset(log4, TaxID %in% log4x$TaxID | TaxID %in% delta4x$TaxID)
log4b$PlotName <- factor(log4b$PlotName)
log4b$PlotName <- forcats::fct_inorder(log4b$PlotName)
log4b$PlotName <- factor(log4b$PlotName, levels=c(levels(log4b$PlotName), "   "))

delta4b <- subset(delta4, PlotName %in% log4b$PlotName)
delta4b$PlotName <- factor(delta4b$PlotName)
delta4b$PlotName <- forcats::fct_inorder(delta4b$PlotName)
delta4b$PlotName <- factor(delta4b$PlotName, levels=c(levels(delta4b$PlotName), "   "))

rm(log4, delta4, log4x, delta4x)

megaplot_actual_abundance <- rbind(log1, log2b, log3b, log4b)
rm(log1, log2b, log3b, log4b)

megaplot_delta_abundance <- rbind(delta1, delta2b, delta3b, delta4b)
rm(delta1, delta2b, delta3b, delta4b)

#### PLOT HEAT MAPS AGAIN ####
levels(megaplot_actual_abundance$PlotName)[levels(megaplot_actual_abundance$PlotName)=="Peptococcaceae 1"] <- "Peptococcaceae"
levels(megaplot_delta_abundance$PlotName)[levels(megaplot_delta_abundance$PlotName)=="Peptococcaceae 1"] <- "Peptococcaceae"
megaplot_delta_abundance$Variable <- factor(megaplot_delta_abundance$Variable,
                                            levels=c("Diet", "Pathology", "Difference"))

plot1 <- ggplot(megaplot_actual_abundance, aes(x=Variable, y=PlotName)) + 
  geom_tile(aes(fill=Coefficient)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="white", guide=guide_colorbar(title.position="top")) +
  geom_text(aes(label=Confidence), vjust=0.75, hjust=0.5) +
  scale_y_discrete(limits=rev(levels(megaplot_actual_abundance$PlotName)), expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  labs(x="Variable", y="Taxon") +
  theme_bw() + 
  plot_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.text.y=element_text(size=7),
        legend.position="bottom",
        legend.title=element_text(hjust=0.5),
        legend.text=element_text(angle=90, hjust=0.3, vjust=0.5),
        panel.grid=element_blank())

plot2 <- ggplot(megaplot_delta_abundance, aes(x=Variable, y=PlotName)) + 
  geom_tile(aes(fill=Coefficient)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="white", guide=guide_colorbar(title.position="top")) +
  geom_text(aes(label=Confidence), vjust=0.75, hjust=0.5) +
  scale_y_discrete(limits=rev(levels(megaplot_delta_abundance$PlotName)), expand=c(0,0)) +
  scale_x_discrete(labels=c("Diet", "Pathology", "\u0394 Time"), 
                   expand=c(0,0)) +
  labs(x="Variable", y="Taxon") +
  theme_bw() + plot_theme + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                                  axis.text.y=element_text(size=7),
                                  legend.position="bottom",
                                  legend.title=element_text(hjust=0.5),
                                  legend.text=element_text(angle=90, hjust=0.3, vjust=0.5),
                                  panel.grid=element_blank())

grid.arrange(plot1, plot2, ncol=2)
rm(plot1, plot2)

#### SUBSET TAXA TO REMOVE TAXA WITH LOW RELATIVE ABUNDANCES ####
# Create a list of taxa with mean abundances > 0.04% #
rel_sum <- list()
for(i in 1:length(otu_tables_rel)){
  df <- otu_tables_rel[[i]]
  df <- as.data.frame(t(df[,c(1:(ncol(df)-8))]))
  df <- data.frame(
    TaxID=rownames(df),
    Mean=rowMeans(df))
  
  tax_tables[[i]]$TaxID <- rownames(tax_tables[[i]])
  df <- merge(df, tax_tables[[i]], by="TaxID", all=FALSE)
  tax_tables[[i]]$TaxID <- NULL
  
  df <- subset(df, Mean > 0.041)
  
  if(i==1){df <- df[-grep("Uncl", df$Phylum),]}
  if(i==2){df <- df[-grep("Uncl", df$Family),]}
  if(i==3){df <- df[-grep("Uncl", df$Genus),]}
  
  if(i==1){
    df$Level <- rep("Phylum")
    df$Name <- df$Phylum
    df$Name <- as.character(df$Name)
    df <- df[order(df$Name),]
    df$Name <- forcats::fct_inorder(df$Name)
    df$Name <- factor(df$Name, levels=c(levels(df$Name), " "))
  }  
  if(i==2){
    df$Level <- rep("Family")
    df$Name <- df$Family
    df$Name <- as.character(df$Name)
    df <- df[order(df$Name),]
    df$Name <- forcats::fct_inorder(df$Name)
    df$Name <- factor(df$Name, levels=c(levels(df$Name), "  "))
  }  
  if(i==3){
    df$Level <- rep("Genus")
    df$Name <- df$Genus
    df$Name <- as.character(df$Name)
    df <- df[order(df$Name),]
    df$Name <- forcats::fct_inorder(df$Name)
    df$Name <- factor(df$Name, levels=c(levels(df$Name), "   "))
  }  
  if(i==4){
    df$Level <- rep("OTU")
    df$Name <- df$PlotName
    df$Name <- factor(df$Name, levels=otu_names)
    df <- df[order(df$Name),]
    df$Name <- as.character(df$Name)
    df$Name <- forcats::fct_inorder(df$Name)
  }  
  
  df <- df[,c("TaxID", "Level","Name", "Mean")]
  rel_sum[[i]] <- df
  rm(df)
}

rel_sum <- rbind(rel_sum[[1]],
                 rel_sum[[2]],
                 rel_sum[[3]],
                 rel_sum[[4]])

# Then trim to relative abundance. #
megaplot_actual_abundance <- subset(megaplot_actual_abundance,
                                    PlotName %in% rel_sum$Name & is.na(Coefficient)=="FALSE")

t1 <- subset(megaplot_actual_abundance, Level=="Phylum")
t1$PlotName <- factor(t1$PlotName)
t1$PlotName <- forcats::fct_inorder(t1$PlotName)
t1$PlotName <- factor(t1$PlotName, levels=c(levels(t1$PlotName), " "))

t2 <- subset(megaplot_actual_abundance, Level=="Family")
t2$PlotName <- factor(t2$PlotName)
t2$PlotName <- forcats::fct_inorder(t2$PlotName)
t2$PlotName <- factor(t2$PlotName, levels=c(levels(t2$PlotName), "  "))

t3 <- subset(megaplot_actual_abundance, Level=="Genus")
t3$PlotName <- factor(t3$PlotName)
t3$PlotName <- forcats::fct_inorder(t3$PlotName)
t3$PlotName <- factor(t3$PlotName, levels=c(levels(t3$PlotName), "   "))

t4 <- subset(megaplot_actual_abundance, Level=="OTU")
t4$PlotName <- factor(t4$PlotName)
t4$PlotName <- forcats::fct_inorder(t4$PlotName)

megaplot_actual_abundance <- rbind(t1, t2, t3, t4)
rm(t1, t2, t3, t4)

megaplot_delta_abundance <- subset(megaplot_delta_abundance,
                                   PlotName %in% rel_sum$Name & is.na(Coefficient)=="FALSE")

t1 <- subset(megaplot_delta_abundance, Level=="Phylum")
t1$PlotName <- factor(t1$PlotName)
t1$PlotName <- forcats::fct_inorder(t1$PlotName)
t1$PlotName <- factor(t1$PlotName, levels=c(levels(t1$PlotName), " "))

t2 <- subset(megaplot_delta_abundance, Level=="Family")
t2$PlotName <- factor(t2$PlotName)
t2$PlotName <- forcats::fct_inorder(t2$PlotName)
t2$PlotName <- factor(t2$PlotName, levels=c(levels(t2$PlotName), "  "))

t3 <- subset(megaplot_delta_abundance, Level=="Genus")
t3$PlotName <- factor(t3$PlotName)
t3$PlotName <- forcats::fct_inorder(t3$PlotName)
t3$PlotName <- factor(t3$PlotName, levels=c(levels(t3$PlotName), "   "))

t4 <- subset(megaplot_delta_abundance, Level=="OTU")
t4$PlotName <- factor(t4$PlotName)
t4$PlotName <- forcats::fct_inorder(t4$PlotName)

megaplot_delta_abundance <- rbind(t1, t2, t3, t4)
megaplot_delta_abundance$PlotName <- factor(megaplot_delta_abundance$PlotName,
                                            levels=levels(megaplot_actual_abundance$PlotName))
rm(t1, t2, t3, t4)

#### FINAL MODEL COEFFICIENT PLOTS #####
plot1 <- ggplot(megaplot_actual_abundance, aes(x=Variable, y=PlotName)) + 
  geom_tile(aes(fill=Coefficient)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="white", guide=guide_colorbar(title.position="top")) +
  geom_text(aes(label=Confidence), vjust=0.75, hjust=0.5) +
  scale_y_discrete(limits=rev(levels(megaplot_actual_abundance$PlotName)), expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  labs(x="Variable", y="Taxon", fill="Coefficient\n(Effect of Predictor)") +
  theme_bw() + 
  plot_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.text.y=element_text(size=7),
        legend.position="bottom",
        legend.title=element_text(hjust=0.5),
        legend.text=element_text(angle=90, hjust=0.3, vjust=0.5),
        panel.grid=element_blank())

plot2 <- ggplot(megaplot_delta_abundance, aes(x=Variable, y=PlotName)) + 
  geom_tile(aes(fill=Coefficient)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="white", guide=guide_colorbar(title.position="top")) +
  geom_text(aes(label=Confidence), vjust=0.75, hjust=0.5) +
  scale_y_discrete(limits=rev(levels(megaplot_delta_abundance$PlotName)), expand=c(0,0)) +
  scale_x_discrete(labels=c("Diet", "Pathology", "\u0394 Time"), 
                   expand=c(0,0)) +
  labs(x="Variable", y="Taxon", fill="Coefficient\n(Effect of Predictor)") +
  theme_bw() + plot_theme + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                                  axis.text.y=element_text(size=7),
                                  legend.position="bottom",
                                  legend.title=element_text(hjust=0.5),
                                  legend.text=element_text(angle=90, hjust=0.3, vjust=0.5),
                                  panel.grid=element_blank())

grid.arrange(plot1, plot2, ncol=2)

ggsave("log_abundance_mixed_model_coefficients.jpg", plot1, width=4, height=7, units="in", dpi=300)
ggsave("delta_abundance_mixed_model_coefficients.jpg", plot2, width=4, height=7, units="in", dpi=300)

write.csv(megaplot_actual_abundance, "log_abundance_mixed_model_coefficients.csv")
write.csv(megaplot_delta_abundance, "delta_abundance_mixed_model_coefficients.csv")

rm(plot1, plot2)

#### ADD RAW ABUNDANCE HEAT MAP DATA AND SAVE PLOTS ########
# Create a single data frame for plotting Z-score transformed abundances.
temp1 <- heatmaps_z_raw_abundance[[1]]
temp1$PlotName <- temp1$Phylum
temp1$Level <- rep("Phylum")
temp1[,c("Domain", "Phylum")] <- NULL

temp2 <- heatmaps_z_raw_abundance[[2]]
temp2$PlotName <- temp2$Family
temp2$Level <- rep("Family")
temp2[,c("Domain", "Phylum", "Class", "Order", "Family")] <- NULL

temp3 <- heatmaps_z_raw_abundance[[3]]
temp3$PlotName <- temp3$Genus
temp3$Level <- rep("Genus")
temp3[,c("Domain", "Phylum", "Class", "Order", "Family", "Genus")] <- NULL

temp4 <- heatmaps_z_raw_abundance[[4]]
temp4$Level <- rep("OTU")
temp4[,c("Species", "Identity")] <- NULL

heatmap_master <- rbind(temp1, temp2, temp3, temp4)
rm(temp1, temp2, temp3, temp4)
heatmap_master <- subset(heatmap_master, PlotName %in% megaplot_actual_abundance$PlotName)
heatmap_master$PlotName <- factor(heatmap_master$PlotName,
                                  levels=levels(megaplot_actual_abundance$PlotName))

# Rename levels for plotting.
levels(heatmap_master$Group)[levels(heatmap_master$Group)=="0_1_0"] <- "WBS Day -56"
levels(heatmap_master$Group)[levels(heatmap_master$Group)=="0_2_0"] <- "WBS Day -7"
levels(heatmap_master$Group)[levels(heatmap_master$Group)=="0_3_0"] <- "WBS Day 21"
levels(heatmap_master$Group)[levels(heatmap_master$Group)=="0_4_0"] <- "WBS Day 49\n(healthy)"
levels(heatmap_master$Group)[levels(heatmap_master$Group)=="0_4_1"] <- "WBS Day 49\n(pathology)"
levels(heatmap_master$Group)[levels(heatmap_master$Group)=="1_1_0"] <- "YBS Day -56"
levels(heatmap_master$Group)[levels(heatmap_master$Group)=="1_2_0"] <- "YBS Day -7"
levels(heatmap_master$Group)[levels(heatmap_master$Group)=="1_3_0"] <- "YBS Day 21"
levels(heatmap_master$Group)[levels(heatmap_master$Group)=="1_4_0"] <- "YBS Day 49\n(healthy)"
levels(heatmap_master$Group)[levels(heatmap_master$Group)=="1_4_1"] <- "YBS Day 49\n(pathology)"

# Export file and add annotations from significance tests (see Part IV)
write.csv(heatmap_master, "~/heatmap_master.csv")
heatmap_master_tagged <- read.csv("raw_data/tagged_abundance_heatmap.csv")
heatmap_master_tagged$Group <- factor(heatmap_master_tagged$Group, levels=c("WBS Day -56",
                                                                            "WBS Day -7",
                                                                            "WBS Day 21",
                                                                            "WBS Day 49\n(healthy)",
                                                                            "WBS Day 49\n(pathology)",
                                                                            "YBS Day -56",
                                                                            "YBS Day -7",
                                                                            "YBS Day 21",
                                                                            "YBS Day 49\n(healthy)",
                                                                            "YBS Day 49\n(pathology)"))

abundance_heatmap <- ggplot(heatmap_master_tagged, aes(x=Group, y=PlotName, fill=Abund)) + geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", guide=guide_colorbar(title.position="top")) +
  geom_text(aes(label=Symbol), vjust=0.4, hjust=0.5, fontface="bold") +
  scale_y_discrete(limits=rev(levels(heatmap_master$PlotName)), expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  labs(x="Diet / Time / Pathology", y="Taxon", fill="Z-score standardized\nabundance") +
  geom_vline(xintercept=5.5, color="black") +
  theme_bw() + plot_theme + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=7),
                                  axis.text.y=element_text(size=7),
                                  panel.grid=element_blank(),
                                  panel.border=element_rect(color="black"),
                                  legend.position="bottom",
                                  legend.title=element_text(hjust=0.5, size=10))

ggsave("~/Documents/marmosets/figures_2020.11.18/abundance_heatmap.jpg", abundance_heatmap, width=6.5, height=7, units="in", dpi=300)

plot1 <- ggplotGrob(plot1 + theme(axis.text.x=element_text(size=7)))
abundance_heatmap <- ggplotGrob(abundance_heatmap + theme(axis.text.y=element_blank(), axis.title.y=element_blank()))

plot1$heights <- abundance_heatmap$heights

grid.arrange(plot1, abundance_heatmap, ncol=2)
figure_1_arranged <- arrangeGrob(plot1, abundance_heatmap, ncol=2)
ggsave("~/Documents/marmosets/figures_2020.11.26/figure_1_arranged.jpg", figure_1_arranged,
       width=7, height=8, units="in", dpi=300)


#### ADD DELTA ABUNDANCE HEAT MAP DATA BY SCALING FROM -10 to 10 #####
temp <- summaries_delta_abundances[[1]]
temp1 <- data.frame(
  TaxID = temp$Row.names,
  Level = rep("Phylum"),
  PlotName = temp$Phylum
)
temp <- as.data.frame(t(temp[,c(4:ncol(temp))]))
for(i in 1:ncol(temp)){
  temp[,i] <- rescale_mid(temp[,i], to=c(-10,10), mid=0)
}
temp <- as.data.frame(t(temp))
temp1 <- cbind(temp1, temp)
rm(temp)
temp1 <- subset(temp1, PlotName %in% megaplot_delta_abundance$PlotName)
temp1$PlotName <- factor(temp1$PlotName)
temp1$PlotName <- forcats::fct_inorder(temp1$PlotName)
temp1$PlotName <- factor(temp1$PlotName, levels=c(levels(temp1$PlotName), " "))

# Family level #
temp <- summaries_delta_abundances[[2]]
temp2 <- data.frame(
  TaxID = temp$Row.names,
  Level = rep("Family"),
  PlotName = temp$Family
)
temp <- as.data.frame(t(temp[,c(7:ncol(temp))]))
for(i in 1:ncol(temp)){
  temp[,i] <- rescale_mid(temp[,i], to=c(-10,10), mid=0)
}
temp <- as.data.frame(t(temp))
temp2 <- cbind(temp2, temp)
rm(temp)
temp2 <- subset(temp2, PlotName %in% megaplot_delta_abundance$PlotName)
temp2$PlotName <- factor(temp2$PlotName)
temp2$PlotName <- forcats::fct_inorder(temp2$PlotName)
temp2$PlotName <- factor(temp2$PlotName, levels=c(levels(temp2$PlotName), "  "))

# Genus level
temp <- summaries_delta_abundances[[3]]
temp3 <- data.frame(
  TaxID = temp$Row.names,
  Level = rep("Genus"),
  PlotName = temp$Genus
)
temp <- as.data.frame(t(temp[,c(8:ncol(temp))]))
for(i in 1:ncol(temp)){
  temp[,i] <- rescale_mid(temp[,i], to=c(-10,10), mid=0)
}
temp <- as.data.frame(t(temp))
temp3 <- cbind(temp3, temp)
rm(temp)
temp3 <- subset(temp3, PlotName %in% megaplot_delta_abundance$PlotName)
temp3$PlotName <- factor(temp3$PlotName)
temp3$PlotName <- forcats::fct_inorder(temp3$PlotName)
temp3$PlotName <- factor(temp3$PlotName, levels=c(levels(temp3$PlotName), "   "))

# OTU level
temp <- summaries_delta_abundances[[4]]
temp4 <- data.frame(
  TaxID = temp$Row.names,
  Level = rep("OTU"),
  PlotName = temp$PlotName
)
temp <- as.data.frame(t(temp[,c(5:ncol(temp))]))
for(i in 1:ncol(temp)){
  temp[,i] <- rescale_mid(temp[,i], to=c(-10,10), mid=0)
}
temp <- as.data.frame(t(temp))
temp4 <- cbind(temp4, temp)
rm(temp)
temp4 <- subset(temp4, PlotName %in% megaplot_delta_abundance$PlotName)
temp4$PlotName <- factor(temp4$PlotName)
temp4$PlotName <- forcats::fct_inorder(temp4$PlotName)

megaplot_delta_map <- rbind(temp1, temp2, temp3, temp4)
megaplot_delta_map <- melt(megaplot_delta_map)

levels(megaplot_delta_map$variable)[levels(megaplot_delta_map$variable)=="0_1v2_0"] <- "WBS\nDays -56 to -7"
levels(megaplot_delta_map$variable)[levels(megaplot_delta_map$variable)=="0_2v3_0"] <- "WBS\nDays -7 to 21"
levels(megaplot_delta_map$variable)[levels(megaplot_delta_map$variable)=="0_3v4_0"] <- "WBS\nDays 21 to 49\n(healthy)"
levels(megaplot_delta_map$variable)[levels(megaplot_delta_map$variable)=="0_3v4_1"] <- "WBS\nDays 21 to 49\n(pathology)"
levels(megaplot_delta_map$variable)[levels(megaplot_delta_map$variable)=="1_1v2_0"] <- "YBS\nDays -56 to -7"
levels(megaplot_delta_map$variable)[levels(megaplot_delta_map$variable)=="1_2v3_0"] <- "YBS\nDays -7 to 21"
levels(megaplot_delta_map$variable)[levels(megaplot_delta_map$variable)=="1_3v4_0"] <- "YBS\nDays 21 to 49\n(healthy)"
levels(megaplot_delta_map$variable)[levels(megaplot_delta_map$variable)=="1_3v4_1"] <- "YBS\nDays 21 to 49\n(pathology)"

megaplot_delta_map$PlotName <- factor(megaplot_delta_map$PlotName, levels=levels(megaplot_delta_abundance$PlotName))

megaplot_delta_map <- merge(megaplot_delta_map, 
                            megaplot_delta_map_tagged[,c("TaxID", "variable", "Symbol")],
                            by=c("TaxID", "variable"),
                            all=TRUE)

# Export to add annotations from significance tests.
write.csv(megaplot_delta_map, "megaplot_delta_map.csv")
megaplot_delta_map_tagged <- read.csv("raw_data/tagged_delta_heatmap.csv")

megaplot_delta_map_tagged$variable <- factor(megaplot_delta_map_tagged$variable,
                                             levels=c( "WBS\nDays -56 to -7",
                                                       "WBS\nDays -7 to 21",
                                                       "WBS\nDays 21 to 49\n(healthy)",
                                                       "WBS\nDays 21 to 49\n(pathology)",
                                                       "YBS\nDays -56 to -7",
                                                       "YBS\nDays -7 to 21",
                                                       "YBS\nDays 21 to 49\n(healthy)",
                                                       "YBS\nDays 21 to 49\n(pathology)"))

delta_heatmap <- ggplot(megaplot_delta_map, aes(x=variable, y=PlotName, fill=value)) + geom_tile() +
  scale_fill_gradient2(low="#6640f7", mid="white", high="#fb0200", midpoint=0, 
                       guide=guide_colorbar(title.position="top")) +
  scale_y_discrete(limits=rev(levels(megaplot_delta_map$PlotName)), expand=c(0,0)) +
  geom_text(aes(label=Symbol), vjust=0.4, hjust=0.5, fontface="bold") +
  scale_x_discrete(expand=c(0,0)) +
  labs(x=paste("Diet / ", expression("\u0394 Time")," / Pathology"),
       y="Taxon",
       fill=paste("Standardized\n", expression("\u0394 abundance"))) +
  geom_vline(xintercept=4.5, color="black") +
  theme_bw() + plot_theme + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=7),
                                  axis.text.y=element_text(size=7),
                                  panel.grid=element_blank(),
                                  panel.border=element_rect(color="black"),
                                  legend.position="bottom",
                                  legend.title=element_text(hjust=0.5, size=10))

ggsave("delta_heatmap.jpg", delta_heatmap, width=6.5, height=7, units="in", dpi=300)

plot2 <- ggplotGrob(plot2 + theme(axis.text.x=element_text(size=7)))
delta_heatmap <- ggplotGrob(delta_heatmap + theme(axis.text.y=element_blank(), axis.title.y=element_blank()))

plot2$heights <- delta_heatmap$heights

grid.arrange(plot2, delta_heatmap, ncol=2)
figure_2_arranged <- arrangeGrob(plot2, delta_heatmap, ncol=2)
ggsave("figure_2_arranged.jpg", figure_2_arranged,
       width=7, height=8, units="in", dpi=300)

#### MIXED-EFFECT MODELS FOR ALPHA DIVERSITY ####
model_alpha <- data.frame(matrix(nrow=0, ncol=5))
colnames(model_alpha) <- c("Measure", "Variable", "Lower", "Coefficient", "Upper")

# Convert inverse simpson (1/D) to modified Simpson (1-D)
sample_data$diversity_inverse_simpson <- 1/sample_data$diversity_inverse_simpson
sample_data$diversity_inverse_simpson <- 1-sample_data$diversity_inverse_simpson

# Calcualte mixed-effect models, as before.
model_data <- sample_data

for(p in c("observed","chao1","diversity_shannon", "diversity_inverse_simpson", "evenness_pielou",
           "BraySmall")){ # For each taxon in that OTU table...
  
  model_data[,p] <- scale(model_data[,p])
  
  model <- lmer(model_data[,p] ~ Diet + Time + Pathology + (1|TwinPair/Subject),
                data=model_data, na.action="na.fail", REML=FALSE)
  
  model_dredge <- dredge(model)
  
  # If there is only one top model, no coefficient-averaging is required.
  # Just identify the top model, extract the coefficients, and calculate confidence intervals.
  if(nrow(subset(as.data.frame(model_dredge), delta < 2))==1){
    x <- as.data.frame(model_dredge)
    x <- subset(x, delta < 2)
    x <- x[,c(2:4)] 
    x <- as.data.frame(t(x))
    x <- subset(x, is.na(x[,1])=="FALSE")
    vars <- rownames(x) # Identify variables present in the top model.
    
    if(length(vars)>0){ # Recreate the model.
      if(length(vars)==1){model <- lmer(model_data[,p] ~ model_data[,vars] + (1|TwinPair/Subject),
                                        data=model_data, na.action="na.fail", REML=FALSE)}
      if(length(vars)==2){model <- lmer(model_data[,p] ~ model_data[,vars[1]] + model_data[,vars[2]] + (1|TwinPair/Subject),
                                        data=model_data, na.action="na.fail", REML=FALSE)}
      if(length(vars)==3){model <- lmer(model_data[,p] ~ model_data[,vars[1]] + model_data[,vars[2]] + model_data[,vars[3]] + (1|TwinPair/Subject),
                                        data=model_data, na.action="na.fail", REML=FALSE)}
      temp <- confint.merMod(model, level=0.95, method="Wald")
      
      fill <- data.frame(
        Lower=temp[,1][c(5:(length(vars)+4))],
        Coefficient=model@beta[c(2:(length(vars)+1))],
        Upper=temp[,2][c(5:(length(vars)+4))])
      fill$Variable <- vars
      
      if(length(vars)==1){fill$Measure <- rep(p, 1)}
      if(length(vars)==2){fill$Measure <- rep(p, 2)}
      if(length(vars)==3){fill$Measure <- rep(p, 3)}
    }
    
    if(length(vars)==0){
      fill <- data.frame(
        Lower=rep(NA, 3),
        Coefficient=rep(NA, 3),
        Upper=rep(NA, 3),
        Variable=c("Diet", "Pathology", "Time"),
        Measure=rep(p, 3)
      )
    }
    
    model_coeff <- fill[,c("Measure", "Variable", "Lower", "Coefficient", "Upper")]
  }
  
  # If there is more than one top model - average coefficients
  if(nrow(subset(as.data.frame(model_dredge), delta < 2)) > 1){
    model_coeff <- model.avg(model_dredge, subset = delta < 2)
    
    conf <- confint(model_coeff)
    conf <- subset(conf, rownames(conf) !="(Intercept)")
    
    model_coeff <- as.data.frame(model_coeff$coefficients)
    model_coeff <- subset(model_coeff, rownames(model_coeff)=="full")
    model_coeff <- as.data.frame(t(model_coeff))
    model_coeff <- subset(model_coeff, rownames(model_coeff) != "(Intercept)")
    
    model_coeff <- data.frame(
      Measure=p,
      Variable=rownames(model_coeff),
      Lower=conf[,1],
      Coefficient=model_coeff[,1],
      Upper=conf[,2]
    )
  }
  
  model_alpha <- rbind(model_alpha, model_coeff)
  
}

rm(model_data, fill, model, model_dredge)

# Define significance based on whether the coefficient overlaps zero.
model_alpha$Significance <- model_alpha$Lower*model_alpha$Upper
model_alpha$Significance[model_alpha$Significance > 0] <- 0.05
model_alpha$Significance[model_alpha$Significance < 0] <- 1
model_alpha$Confidence <- cut(model_alpha$Significance,
                              breaks=c(0,0.07,Inf),
                              labels=c("**", " "))

# Clean data for plotting.
model_alpha <- subset(model_alpha, is.na(Measure)=="FALSE")

model_alpha <- model_alpha %>%
  dplyr::group_by(Measure) %>%
  tidyr::complete(Variable = c("Diet", "Time", "Pathology"))

model_alpha$Measure <- factor(model_alpha$Measure)

levels(model_alpha$Measure)[levels(model_alpha$Measure)=="observed"] <- "Species\nrichness"
levels(model_alpha$Measure)[levels(model_alpha$Measure)=="diversity_inverse_simpson"] <- "Simpson\ndiversity"
levels(model_alpha$Measure)[levels(model_alpha$Measure)=="evenness_pielou"] <- "Pielou\nevenness"
levels(model_alpha$Measure)[levels(model_alpha$Measure)=="BraySmall"] <- "Bray-Curtis\ndissimilarity"

model_alpha$Measure <- factor(model_alpha$Measure, levels=c("Bray-Curtis\ndissimilarity", "Pielou\nevenness", "Species\nrichness", "Simpson\ndiversity"))

# Generate plot.
alpha_div_mixed_models <- ggplot(subset(model_alpha,
                                        Measure %in% c("Species\nrichness", "Simpson\ndiversity", "Pielou\nevenness", "Bray-Curtis\ndissimilarity")),
                                 aes(x=Variable, y=Measure)) +
  geom_tile(aes(fill=Coefficient)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="white") +
  geom_text(aes(label=Confidence), vjust=0.75, hjust=0.5) +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  labs(x="Variable", y="Measure") +
  theme_bw() + plot_theme + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                                  axis.text.y=element_text(size=9),
                                  panel.grid=element_blank())

ggsave("diversity_mixed_models.jpg", alpha_div_mixed_models, width=5, height=3, units="in", dpi=300)

rm(alpha_div_mixed_models)

#### CALCULATE MARGINAL AND CONDITIONAL R2 FOR EACH MODEL ####
model_data <- sample_data
r2c_alpha <- data.frame()

# Alpha diversity models.
for(p in c("observed","chao1","diversity_shannon", "diversity_inverse_simpson", "evenness_pielou",
           "BraySmall")){ # For each measure...
  
  model_data[,p] <- scale(model_data[,p])
  
  model <- lmer(model_data[,p] ~ Diet + Time + Pathology + (1|TwinPair/Subject),
                data=model_data, na.action="na.fail", REML=FALSE)
  
  x <- 100*r.squaredGLMM(model)
  colnames(x) <- c("fixed", "total")
  x <- as.data.frame(x)
  x$random <- x$total - x$fixed
  x$var <- p

  x <- x[,c("var", "fixed", "random", "total")]
  
  r2c_alpha <- rbind(r2c_alpha, x)
  rm(x, model)
}

# Log-transformed models.
r2c_log <- data.frame()

for(i in c(1:4)){ 
  model_data <- otu_tables_log[[i]]
  model_data[,c(1:(ncol(model_data)-8))] <- scale(model_data[,c(1:(ncol(model_data)-8))])
  
  for(p in c(1:(ncol(model_data)-8))){ # For each taxon in the OTU table
    
    model <- lmer(model_data[,p] ~ Diet + Time + Pathology + (1|TwinPair/Subject), # Construct a mixed-effect model
                  data=model_data, na.action="na.fail", REML=FALSE)
    
    x <- 100*r.squaredGLMM(model)
    colnames(x) <- c("fixed", "total")
    x <- as.data.frame(x)
    x$random <- x$total - x$fixed
    x$taxa <- colnames(model_data)[p]
    
    x <- x[,c("taxa", "fixed", "random", "total")]
    
    r2c_log <- rbind(r2c_log, x)
  }
  
  rm(x, model)
}

# Determine whether fixed or random effects explain more variance.
r2c_log$winner <- rep("random")
for(i in 1:nrow(r2c_log)){
  if(r2c_log[i,"fixed"] > r2c_log[i, "random"]){
    r2c_log[i,"winner"] <- "fixed"
  }
}

r2c_log %>% dplyr::group_by(winner) %>% dplyr::count()

# Delta abundance models.
r2c_delta <- data.frame()

for(i in c(1:4)){ 
  model_data <- otu_tables_delta_log[[i]]
  
  model_data$Difference[model_data$Difference=="1_2"] <- 0
  model_data$Difference[model_data$Difference=="2_3"] <- 1
  model_data$Difference[model_data$Difference=="3_4"] <- 2
  model_data$Difference <- as.integer(model_data$Difference)
  
  model_data[,c(2:(ncol(model_data)-5))] <- scale(model_data[,c(2:(ncol(model_data)-5))])
  
  for(p in c(2:(ncol(model_data)-5))){ # For each taxon in the OTU table
    
    model <- lmer(model_data[,p] ~ Diet + Difference + Pathology + (1|TwinPair/Subject), # Construct a mixed-effect model
                  data=model_data, na.action="na.fail", REML=FALSE)
    
    x <- 100*r.squaredGLMM(model)
    colnames(x) <- c("fixed", "total")
    x <- as.data.frame(x)
    x$random <- x$total - x$fixed
    x$taxa <- colnames(model_data)[p]
    
    x <- x[,c("taxa", "fixed", "random", "total")]
    
    r2c_delta <- rbind(r2c_delta, x)
  }
  
  rm(x, model)
}

# Determine whether fixed or random effects explain more variance.
r2c_delta$winner <- rep("random")
for(i in 1:nrow(r2c_delta)){
  if(r2c_delta[i,"fixed"] > r2c_delta[i, "random"]){
    r2c_delta[i,"winner"] <- "fixed"
  }
}

# Count the number of models where fixed vs. random effects explain more variance.
r2c_log %>% dplyr::group_by(winner) %>% dplyr::count()
r2c_delta %>% dplyr::group_by(winner) %>% dplyr::count()

# Calculate the mean variance explained.
colMeans(r2c_log[,c("fixed", "random")])
colMeans(r2c_delta[,c("fixed", "random")])
colMeans(r2c_alpha[,c("fixed", "random")])

# Add taxon information to r2c_log and save.
r2clog_otu <- r2c_log[c(60:nrow(r2c_log)), ]
tax_tables[[4]]$taxa <- rownames(tax_tables[[4]])
r2clog_otu <- merge(r2clog_otu, tax_tables[[4]], by="taxa", all=FALSE)
write.csv(r2clog_otu, "r2c_log_otu.csv")

# Add taxon information to r2c_delta.
r2c_delta1 <- r2c_delta[c(1:7), ]
r2c_delta2 <- r2c_delta[c(8:28), ]
r2c_delta3 <- r2c_delta[c(29:59), ]
r2c_delta4 <- r2c_delta[c(60:nrow(r2c_delta)), ]

tax_tables[[1]]$taxa <- rownames(tax_tables[[1]])
tax_tables[[2]]$taxa <- rownames(tax_tables[[2]])
tax_tables[[3]]$taxa <- rownames(tax_tables[[3]])

r2c_delta1 <- merge(r2c_delta1, tax_tables[[1]][,c("taxa", "Phylum")], by="taxa", all=FALSE)
r2c_delta2 <- merge(r2c_delta2, tax_tables[[2]][,c("taxa", "Family")], by="taxa", all=FALSE)
r2c_delta3 <- merge(r2c_delta3, tax_tables[[3]][,c("taxa", "Genus")], by="taxa", all=FALSE)
r2c_delta4 <- merge(r2c_delta4, tax_tables[[4]][,c("taxa", "PlotName")], by="taxa", all=FALSE)

colnames(r2c_delta1)[6] <- "Name"
colnames(r2c_delta2)[6] <- "Name"
colnames(r2c_delta3)[6] <- "Name"
colnames(r2c_delta4)[6] <- "Name"

r2c_delta <- rbind(r2c_delta1, r2c_delta2, r2c_delta3, r2c_delta4)
write.csv(r2c_delta, "r2c_delta.csv")

################################################# PART THREE: ORDINATIONS ########################################################
#### DISTANCE MATRIX AND PERMANOVA RESULTS ####
colnames(otu_tables_rel[[4]])
data <- otu_tables_rel[[4]]
data[,c(74:81)] <- NULL

dist_matrix_bray <- vegdist(data, method="bray")
adonis(dist_matrix_bray ~ Diet*Time + Diet/Cage + Pathology + TwinPair/Subject, sample_data)

adonis(dist_matrix_bray ~ TwinPair, sample_data)
adonis(dist_matrix_bray ~ Diet/Cage, sample_data)
adonis(dist_matrix_bray ~ Pathology, sample_data)
adonis(dist_matrix_bray ~ Time, sample_data)
adonis(dist_matrix_bray ~ Diet, sample_data)

#### NMDS SERIES: COMPARING TWO DIETS AT EACH TIME ####
# The following code sections all follow this same format.
# Only this section is annotated.

nmds_time_series <- list()
nmds_time_series[["data"]] <- list()
nmds_time_series[["ellipses"]] <- list()
nmds_time_series[["means"]] <- list()
nmds_time_series[["models"]] <- list()
nmds_time_series[["plots"]] <- list()
permanova_time_series <- data.frame(matrix(nrow=4, ncol=7))
colnames(permanova_time_series) <- c("R2", "F", "df", "p", "disp_F", "disp_df", "disp_p")

for(i in c(1:4)){
  data <- otu_tables_rel[[4]]
  data[,c(74:81)] <- NULL
  
  # Subset the data to the appropriate comparison.
  if(i==1){data <- subset(data, rownames(data) %in% subset(sample_data, Time==1)$SampleID)}
  if(i==2){data <- subset(data, rownames(data) %in% subset(sample_data, Time==2)$SampleID)}
  if(i==3){data <- subset(data, rownames(data) %in% subset(sample_data, Time==3)$SampleID)}
  if(i==4){data <- subset(data, rownames(data) %in% subset(sample_data, Time==4)$SampleID)}
  
  # Calculate NMDS
  nmds_model <- metaMDS(data, k=2, trymax=100)
  
  # Extract axis loadings and add sample metadata.
  nmds_model_data <- as.data.frame(scores(nmds_model, display="sites"))
  nmds_model_data$SampleID <- rownames(nmds_model_data)
  nmds_model_data <- merge(nmds_model_data, sample_data[,c("SampleID", "TwinPair", "Subject", "Diet", "Time", "Pathology", "Cage")],
                           by="SampleID", all=FALSE)
  nmds_model_data$Diet <- as.factor(nmds_model_data$Diet)
  nmds_model_data$Pathology <- as.factor(nmds_model_data$Pathology)
  nmds_model_data$Time <- as.factor(nmds_model_data$Time)
  
  # Calculate PERMANOVA for each comparison.
  dist <- vegdist(data, method="bray")
  temp <- adonis(dist ~ Diet, nmds_model_data)
  
  permanova_time_series[i,1] <- temp$aov.tab[1,5] # R2
  permanova_time_series[i,2] <- temp$aov.tab[1,4] # F
  permanova_time_series[i,3] <- temp$aov.tab[1,1] # df
  permanova_time_series[i,4] <- temp$aov.tab[1,6] # p
  
  # Calculate dispersion for each comparison.
  temp <- betadisper(dist, nmds_model_data$Diet)
  temp <- permutest(temp)
  
  permanova_time_series[i,5] <- temp$tab[1,4] # F
  permanova_time_series[i,6] <- temp$tab[1,1] # df
  permanova_time_series[i,7] <- temp$tab[1,6] # p
  
  rm(temp)
  
  # Calculate ellipses.
  df_ell <- data.frame()
  for(g in levels(nmds_model_data$Diet)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(nmds_model_data[nmds_model_data$Diet==g,],
                                                     veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                  ,group=g))
  }
  
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==0] <- "WBS"
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  
  df_ell$group <- factor(df_ell$group)
  levels(df_ell$group)[levels(df_ell$group)==0] <- "WBS"
  levels(df_ell$group)[levels(df_ell$group)==1] <- "YBS"
  
  # Calculate the centroid of each group.
  NMDS.mean=aggregate(nmds_model_data[,c("NMDS1","NMDS2")], list(group=nmds_model_data$Diet), mean)  
  
  # Plot each timepoint, separately colored by diet.
  nmds_time_series$data[[i]] <- nmds_model_data
  nmds_time_series$ellipses[[i]] <- df_ell
  nmds_time_series$means[[i]] <- NMDS.mean
  nmds_time_series$models[[i]] <- nmds_model
  
  plot <- ggplot(data = nmds_model_data, aes(NMDS1, NMDS2)) + geom_point(aes(color = Diet)) + 
    scale_color_manual(values=c("deepskyblue3", "orangered1")) +
    geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean$NMDS1, y=NMDS.mean$NMDS2,label=NMDS.mean$group, size=3.5) +
    theme_bw() +
    theme(
      legend.justification = c(1, 0),
      legend.position = "none",
      legend.direction = "horizontal",
      plot.title=element_text(size=10),
      panel.grid=element_blank(),
      legend.text = element_text(colour = "black",size = 12),
      axis.text.x = element_text(size = 8,colour = "black"),
      axis.text.y = element_text(size = 8,colour = "black"),
      axis.title.x =element_text(size = 10),
      axis.title.y = element_text(size = 10)) +
    guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                 title.position = "top", title.hjust = 0.5))
  
  nmds_time_series$plots[[i]] <- plot
}

rownames(permanova_time_series) <- c("WBSvYBS_t1", "WBSvYBS_t2", "WBSvYBS_t3", "WBSvYBS_t4")
grid.arrange(nmds_time_series$plots[[1]] + labs(title="Day -56"),
             nmds_time_series$plots[[2]] + labs(title="Day -7"),
             nmds_time_series$plots[[3]] + labs(title="Day 21"),
             nmds_time_series$plots[[4]] + labs(title="Day 49"))

#### NMDS SERIES: COMPARING CONSECUTIVE DAYS, YBS DIET ####
nmds_ybs_series <- list()
nmds_ybs_series[["data"]] <- list()
nmds_ybs_series[["ellipses"]] <- list()
nmds_ybs_series[["means"]] <- list()
nmds_ybs_series[["models"]] <- list()
nmds_ybs_series[["plots"]] <- list()
permanova_ybs_series <- data.frame(matrix(nrow=3, ncol=7))
colnames(permanova_ybs_series) <- c("R2", "F", "df", "p", "disp_F", "disp_df", "disp_p")

for(i in c(1:3)){
  data <- otu_tables_rel[[4]]
  data[,c(74:81)] <- NULL
  
  if(i==1){data <- subset(data, rownames(data) %in% subset(sample_data, Time==1 | Time==2 & Diet==1)$SampleID)}
  if(i==2){data <- subset(data, rownames(data) %in% subset(sample_data, Time==2 | Time==3 & Diet==1)$SampleID)}
  if(i==3){data <- subset(data, rownames(data) %in% subset(sample_data, Time==3 | Time==4 & Diet==1)$SampleID)}
  
  nmds_model <- metaMDS(data, k=2, trymax=100)
  nmds_model_data <- as.data.frame(scores(nmds_model, display="sites"))
  nmds_model_data$SampleID <- rownames(nmds_model_data)
  nmds_model_data <- merge(nmds_model_data, sample_data[,c("SampleID", "TwinPair", "Subject", "Diet", "Time", "Pathology", "Cage")],
                           by="SampleID", all=FALSE)
  nmds_model_data$Diet <- as.factor(nmds_model_data$Diet)
  nmds_model_data$Pathology <- as.factor(nmds_model_data$Pathology)
  nmds_model_data$Time <- as.factor(nmds_model_data$Time)
  
  dist <- vegdist(data, method="bray")
  temp <- adonis(dist ~ Time, nmds_model_data)
  
  permanova_ybs_series[i,1] <- temp$aov.tab[1,5] # R2
  permanova_ybs_series[i,2] <- temp$aov.tab[1,4] # F
  permanova_ybs_series[i,3] <- temp$aov.tab[1,1] # df
  permanova_ybs_series[i,4] <- temp$aov.tab[1,6] # p
  
  temp <- betadisper(dist, nmds_model_data$Time)
  temp <- permutest(temp)
  
  permanova_ybs_series[i,5] <- temp$tab[1,4] # F
  permanova_ybs_series[i,6] <- temp$tab[1,1] # df
  permanova_ybs_series[i,7] <- temp$tab[1,6] # p
  
  rm(temp)
  
  df_ell <- data.frame()
  for(g in levels(nmds_model_data$Time)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(nmds_model_data[nmds_model_data$Time==g,],
                                                     veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                  ,group=g))
  }
  
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==0] <- "WBS"
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  
  df_ell$group <- factor(df_ell$group)
  levels(df_ell$group)[levels(df_ell$group)==1] <- "Day -56"
  levels(df_ell$group)[levels(df_ell$group)==2] <- "Day -7"
  levels(df_ell$group)[levels(df_ell$group)==3] <- "Day 21"
  levels(df_ell$group)[levels(df_ell$group)==4] <- "Day 49"
  
  
  levels(nmds_model_data$Time)[levels(nmds_model_data$Time)==1] <- "Day -56"
  levels(nmds_model_data$Time)[levels(nmds_model_data$Time)==2] <- "Day -7"
  levels(nmds_model_data$Time)[levels(nmds_model_data$Time)==3] <- "Day 21"
  levels(nmds_model_data$Time)[levels(nmds_model_data$Time)==4] <- "Day 49"
  
  NMDS.mean=aggregate(nmds_model_data[,c("NMDS1","NMDS2")], list(group=nmds_model_data$Time), mean)  
  
  # PLOT EACH TWO-DAY COMPARISON
  nmds_ybs_series$data[[i]] <- nmds_model_data
  nmds_ybs_series$ellipses[[i]] <- df_ell
  nmds_ybs_series$means[[i]] <- NMDS.mean
  nmds_ybs_series$models[[i]] <- nmds_model
  
  plot <- ggplot(data = nmds_model_data, aes(NMDS1, NMDS2)) + geom_point(aes(color = Time)) + 
    #scale_color_manual(values=c("magenta", "darkorange")) +
    geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, colour=group), size=1, linetype=2)+
    #annotate("text",x=NMDS.mean$NMDS1, y=NMDS.mean$NMDS2,label=NMDS.mean$group, size=3.5) +
    theme_bw() +
    theme(
      legend.justification = c(1, 0),
      legend.position = "none",
      legend.direction = "horizontal",
      plot.title=element_text(size=10),
      panel.grid=element_blank(),
      legend.text = element_text(colour = "black",size = 12),
      axis.text.x = element_text(size = 8,colour = "black"),
      axis.text.y = element_text(size = 8,colour = "black"),
      axis.title.x =element_text(size = 10),
      axis.title.y = element_text(size = 10)) +
    guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                 title.position = "top", title.hjust = 0.5))
  
  nmds_ybs_series$plots[[i]] <- plot
}

rownames(permanova_ybs_series) <- c("t1vt2_YBS", "t2vt3_YBS", "t3vt4_YBS")
grid.arrange(nmds_ybs_series$plots[[1]] + labs(title="YBS: Day -56 vs. Day -7"),
             nmds_ybs_series$plots[[2]] + labs(title="YBS: Day -7 vs. Day 21"),
             nmds_ybs_series$plots[[3]] + labs(title="YBS: Day 21 vs. Day 49"))

#### NMDS SERIES: COMPARING CONSECUTIVE DAYS, WBS DIET ####
nmds_wbs_series <- list()
nmds_wbs_series[["data"]] <- list()
nmds_wbs_series[["ellipses"]] <- list()
nmds_wbs_series[["means"]] <- list()
nmds_wbs_series[["models"]] <- list()
nmds_wbs_series[["plots"]] <- list()
permanova_wbs_series <- data.frame(matrix(nrow=3, ncol=7))
colnames(permanova_wbs_series) <- c("R2", "F", "df", "p", "disp_F", "disp_df", "disp_p")

for(i in c(1:3)){
  data <- otu_tables_rel[[4]]
  data[,c(74:81)] <- NULL
  
  if(i==1){data <- subset(data, rownames(data) %in% subset(sample_data, Time==1 | Time==2 & Diet==0)$SampleID)}
  if(i==2){data <- subset(data, rownames(data) %in% subset(sample_data, Time==2 | Time==3 & Diet==0)$SampleID)}
  if(i==3){data <- subset(data, rownames(data) %in% subset(sample_data, Time==3 | Time==4 & Diet==0)$SampleID)}
  
  nmds_model <- metaMDS(data, k=2, trymax=100)
  nmds_model_data <- as.data.frame(scores(nmds_model, display="sites"))
  nmds_model_data$SampleID <- rownames(nmds_model_data)
  nmds_model_data <- merge(nmds_model_data, sample_data[,c("SampleID", "TwinPair", "Subject", "Diet", "Time", "Pathology", "Cage")],
                           by="SampleID", all=FALSE)
  nmds_model_data$Diet <- as.factor(nmds_model_data$Diet)
  nmds_model_data$Pathology <- as.factor(nmds_model_data$Pathology)
  nmds_model_data$Time <- as.factor(nmds_model_data$Time)
  
  dist <- vegdist(data, method="bray")
  temp <- adonis(dist ~ Time, nmds_model_data)
  
  permanova_wbs_series[i,1] <- temp$aov.tab[1,5] # R2
  permanova_wbs_series[i,2] <- temp$aov.tab[1,4] # F
  permanova_wbs_series[i,3] <- temp$aov.tab[1,1] # df
  permanova_wbs_series[i,4] <- temp$aov.tab[1,6] # p
  
  temp <- betadisper(dist, nmds_model_data$Time)
  temp <- permutest(temp)
  
  permanova_wbs_series[i,5] <- temp$tab[1,4] # F
  permanova_wbs_series[i,6] <- temp$tab[1,1] # df
  permanova_wbs_series[i,7] <- temp$tab[1,6] # p
  
  rm(temp)
  
  
  df_ell <- data.frame()
  for(g in levels(nmds_model_data$Time)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(nmds_model_data[nmds_model_data$Time==g,],
                                                     veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                  ,group=g))
  }
  
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==0] <- "WBS"
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  
  df_ell$group <- factor(df_ell$group)
  levels(df_ell$group)[levels(df_ell$group)==1] <- "Day -56"
  levels(df_ell$group)[levels(df_ell$group)==2] <- "Day -7"
  levels(df_ell$group)[levels(df_ell$group)==3] <- "Day 21"
  levels(df_ell$group)[levels(df_ell$group)==4] <- "Day 49"
  
  
  levels(nmds_model_data$Time)[levels(nmds_model_data$Time)==1] <- "Day -56"
  levels(nmds_model_data$Time)[levels(nmds_model_data$Time)==2] <- "Day -7"
  levels(nmds_model_data$Time)[levels(nmds_model_data$Time)==3] <- "Day 21"
  levels(nmds_model_data$Time)[levels(nmds_model_data$Time)==4] <- "Day 49"
  
  NMDS.mean=aggregate(nmds_model_data[,c("NMDS1","NMDS2")], list(group=nmds_model_data$Time), mean)  
  
  # PLOT by TREATMENT (DIET)
  nmds_wbs_series$data[[i]] <- nmds_model_data
  nmds_wbs_series$ellipses[[i]] <- df_ell
  nmds_wbs_series$means[[i]] <- NMDS.mean
  nmds_wbs_series$models[[i]] <- nmds_model
  
  plot <- ggplot(data = nmds_model_data, aes(NMDS1, NMDS2)) + geom_point(aes(color = Time)) + 
    #scale_color_manual(values=c("forestgreen", "purple")) +
    geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean$NMDS1, y=NMDS.mean$NMDS2,label=NMDS.mean$group, size=3.5) +
    theme_bw() +
    theme(
      legend.justification = c(1, 0),
      legend.position = "none",
      legend.direction = "horizontal",
      plot.title=element_text(size=10),
      panel.grid=element_blank(),
      legend.text = element_text(colour = "black",size = 12),
      axis.text.x = element_text(size = 8,colour = "black"),
      axis.text.y = element_text(size = 8,colour = "black"),
      axis.title.x =element_text(size = 10),
      axis.title.y = element_text(size = 10)) +
    guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                 title.position = "top", title.hjust = 0.5))
  
  nmds_wbs_series$plots[[i]] <- plot
}

rownames(permanova_wbs_series) <- c("t1vt2_WBS", "t2vt3_WBS", "t3vt4_WBS")
grid.arrange(nmds_wbs_series$plots[[1]] + labs(title="WBS: Day -56 vs. Day -7"),
             nmds_wbs_series$plots[[2]] + labs(title="WBS: Day -7 vs. Day 21"),
             nmds_wbs_series$plots[[3]] + labs(title="WBS: Day 21 vs. Day 49"))

#### NMDS SERIES: PATHOLOGY BY DIET ####
nmds_pathology_series <- list()
nmds_pathology_series[["data"]] <- list()
nmds_pathology_series[["ellipses"]] <- list()
nmds_pathology_series[["means"]] <- list()
nmds_pathology_series[["models"]] <- list()
nmds_pathology_series[["plots"]] <- list()
permanova_pathology_series <- data.frame(matrix(nrow=3, ncol=7))
colnames(permanova_pathology_series) <- c("R2", "F", "df", "p", "disp_F", "disp_df", "disp_p")

for(i in c(1:3)){
  data <- otu_tables_rel[[4]]
  data[,c(74:81)] <- NULL
  
  if(i==1){data <- subset(data, rownames(data) %in% subset(sample_data, Time %in% c(3,4) & Diet==0)$SampleID)}
  if(i==2){data <- subset(data, rownames(data) %in% subset(sample_data, Time %in% c(3,4) & Diet==1)$SampleID)}
  if(i==3){data <- subset(data, rownames(data) %in% subset(sample_data, Time %in% c(3,4))$SampleID)}
  
  nmds_model <- metaMDS(data, k=2, trymax=100)
  nmds_model_data <- as.data.frame(scores(nmds_model, display="sites"))
  nmds_model_data$SampleID <- rownames(nmds_model_data)
  nmds_model_data <- merge(nmds_model_data, sample_data[,c("SampleID", "TwinPair", "Subject", "Diet", "Time",
                                                           "PathologyAnticipate", "Cage")],
                           by="SampleID", all=FALSE)
  nmds_model_data$Diet <- as.factor(nmds_model_data$Diet)
  nmds_model_data$PathologyAnticipate <- as.factor(nmds_model_data$PathologyAnticipate)
  nmds_model_data$Time <- as.factor(nmds_model_data$Time)
  
  dist <- vegdist(data, method="bray")
  temp <- adonis(dist ~ PathologyAnticipate, nmds_model_data)
  
  permanova_pathology_series[i,1] <- temp$aov.tab[1,5] # R2
  permanova_pathology_series[i,2] <- temp$aov.tab[1,4] # F
  permanova_pathology_series[i,3] <- temp$aov.tab[1,1] # df
  permanova_pathology_series[i,4] <- temp$aov.tab[1,6] # p
  
  temp <- betadisper(dist, nmds_model_data$PathologyAnticipate)
  
  temp <- permutest(temp)
  
  permanova_pathology_series[i,5] <- temp$tab[1,4] # F
  permanova_pathology_series[i,6] <- temp$tab[1,1] # df
  permanova_pathology_series[i,7] <- temp$tab[1,6] # p
  
  rm(temp)
  
  df_ell <- data.frame()
  for(g in levels(nmds_model_data$PathologyAnticipate)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(nmds_model_data[nmds_model_data$PathologyAnticipate==g,],
                                                     veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                  ,group=g))
    }
  
  
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==0] <- "WBS"
  levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
  
  df_ell$group <- factor(df_ell$group)
  levels(df_ell$group)[levels(df_ell$group)==0] <- "healthy"
  levels(df_ell$group)[levels(df_ell$group)==1] <- "pathology"
  
  levels(nmds_model_data$PathologyAnticipate)[levels(nmds_model_data$PathologyAnticipate)==0] <- "healthy"
  levels(nmds_model_data$PathologyAnticipate)[levels(nmds_model_data$PathologyAnticipate)==1] <- "pathology"
  
  NMDS.mean=aggregate(nmds_model_data[,c("NMDS1","NMDS2")], list(group=nmds_model_data$PathologyAnticipate), mean)  
  
  # PLOT by TREATMENT (DIET)
  nmds_pathology_series$data[[i]] <- nmds_model_data
  nmds_pathology_series$ellipses[[i]] <- df_ell
  nmds_pathology_series$means[[i]] <- NMDS.mean
  nmds_pathology_series$models[[i]] <- nmds_model
  
  if(i >1){
  plot <- ggplot(data = nmds_model_data, aes(NMDS1, NMDS2)) + geom_point(aes(color = PathologyAnticipate)) + 
    scale_color_manual(values=c("black", "goldenrod")) +
    geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean$NMDS1, y=NMDS.mean$NMDS2,label=NMDS.mean$group, size=3.5) +
    theme_bw() +
    theme(
      legend.justification = c(1, 0),
      legend.position = "none",
      legend.direction = "horizontal",
      plot.title=element_text(size=10),
      panel.grid=element_blank(),
      legend.text = element_text(colour = "black",size = 12),
      axis.text.x = element_text(size = 8,colour = "black"),
      axis.text.y = element_text(size = 8,colour = "black"),
      axis.title.x =element_text(size = 10),
      axis.title.y = element_text(size = 10)) +
    guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                 title.position = "top", title.hjust = 0.5))}
  
  if(i==1){
    plot <- ggplot(data = nmds_model_data, aes(NMDS1, NMDS2)) + geom_point(aes(color = PathologyAnticipate)) + 
      scale_color_manual(values=c("black", "goldenrod")) +
      geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, colour=group), size=1, linetype=2)+
      annotate("text",
               x=subset(NMDS.mean, group=="healthy")$NMDS1,
               y=subset(NMDS.mean, group=="healthy")$NMDS2-0.05,
               label=subset(NMDS.mean, group=="healthy")$group, size=3.5) +
      annotate("text",
               x=subset(NMDS.mean, group=="pathology")$NMDS1,
               y=subset(NMDS.mean, group=="pathology")$NMDS2+0.05,
               label=subset(NMDS.mean, group=="pathology")$group, size=3.5) +
      theme_bw() +
      theme(
        legend.justification = c(1, 0),
        legend.position = "none",
        legend.direction = "horizontal",
        panel.grid=element_blank(),
        plot.title=element_text(size=10),
        legend.text = element_text(colour = "black",size = 12),
        axis.text.x = element_text(size = 8,colour = "black"),
        axis.text.y = element_text(size = 8,colour = "black"),
        axis.title.x =element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
      guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                   title.position = "top", title.hjust = 0.5))}
  
  nmds_pathology_series$plots[[i]] <- plot
}

rownames(permanova_pathology_series) <- c("p0vp1_WBS", "p0vp1_YBS", "p0vp1_all")
grid.arrange(nmds_pathology_series$plots[[1]] + labs(title="WBS: Pathology"),
             nmds_pathology_series$plots[[2]] + labs(title="YBS: Pathology"),
             nmds_pathology_series$plots[[3]] + labs(title="All animals: Pathology"))

# -Grouped by diet/pathlogy (four groups) ####
data <- otu_tables_rel[[4]]
data[,c(74:81)] <- NULL
data <- subset(data, rownames(data) %in% subset(sample_data, Time %in% c(3,4))$SampleID)
nmds_model <- metaMDS(data, k=2, trymax=100)
nmds_model_data <- as.data.frame(scores(nmds_model, display="sites"))
nmds_model_data$SampleID <- rownames(nmds_model_data)
nmds_model_data <- merge(nmds_model_data, sample_data[,c("SampleID", "TwinPair", "Subject", "Diet", "Time", "PathologyAnticipate", "Cage")],
                         by="SampleID", all=FALSE)
nmds_model_data$Diet <- as.factor(nmds_model_data$Diet)
nmds_model_data$PathologyAnticipate <- as.factor(nmds_model_data$PathologyAnticipate)
nmds_model_data$Time <- as.factor(nmds_model_data$Time)
nmds_model_data$Diet_Path <- paste0(nmds_model_data$Diet, "_", nmds_model_data$PathologyAnticipate)
nmds_model_data$Diet_Path <- factor(nmds_model_data$Diet_Path)

dist <- vegdist(data, method="bray")
adonis(dist ~ Diet_Path, nmds_model_data)
permutest(betadisper(dist, nmds_model_data$Diet_Path), pairwise=TRUE)
pairwiseAdonis::pairwise.adonis2(dist ~ Diet_Path, nmds_model_data)

df_ell <- data.frame()
for(g in levels(nmds_model_data$Diet_Path)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(nmds_model_data[nmds_model_data$Diet_Path==g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                ,group=g))
}

levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==0] <- "WBS"
levels(nmds_model_data$Diet)[levels(nmds_model_data$Diet)==1] <- "YBS"
levels(nmds_model_data$Diet_Path)[levels(nmds_model_data$Diet_Path)=="0_0"] <- "WBS, healthy"
levels(nmds_model_data$Diet_Path)[levels(nmds_model_data$Diet_Path)=="0_1"] <- "WBS, pathology"
levels(nmds_model_data$Diet_Path)[levels(nmds_model_data$Diet_Path)=="1_0"] <- "YBS, healthy"
levels(nmds_model_data$Diet_Path)[levels(nmds_model_data$Diet_Path)=="1_1"] <- "YBS, pathology"

df_ell$group <- factor(df_ell$group)
levels(df_ell$group)[levels(df_ell$group)=="0_0"] <- "WBS, healthy"
levels(df_ell$group)[levels(df_ell$group)=="0_1"] <- "WBS, pathology"
levels(df_ell$group)[levels(df_ell$group)=="1_0"] <- "YBS, healthy"
levels(df_ell$group)[levels(df_ell$group)=="1_1"] <- "YBS, pathology"

NMDS.mean=aggregate(nmds_model_data[,c("NMDS1","NMDS2")], list(group=nmds_model_data$Diet_Path), mean)  
levels(NMDS.mean$group)[levels(NMDS.mean$group)=="WBS, healthy"] <- "WBS\nhealthy"
levels(NMDS.mean$group)[levels(NMDS.mean$group)=="WBS, pathology"] <- "WBS\npathology"
levels(NMDS.mean$group)[levels(NMDS.mean$group)=="YBS, healthy"] <- "YBS\nhealthy"
levels(NMDS.mean$group)[levels(NMDS.mean$group)=="YBS, pathology"] <- "YBS\npathology"

# PLOT by TREATMENT (DIET)
nmds_pathology_series$data[[4]] <- nmds_model_data
nmds_pathology_series$ellipses[[4]] <- df_ell
nmds_pathology_series$means[[4]] <- NMDS.mean
nmds_pathology_series$models[[4]] <- nmds_model

ggplot(data = nmds_model_data, aes(NMDS1, NMDS2)) + geom_point(aes(color = Diet_Path)) + 
    scale_color_manual(values=c("black", "gold2", "darkgrey", "darkgoldenrod4")) +
    geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1, linetype=2)+
  #annotate("text",x=NMDS.mean[1,"NMDS1"]+0.05, y=NMDS.mean[1,"NMDS2"],label=NMDS.mean[1,"group"], size=3) +
  #annotate("text",x=NMDS.mean[2,"NMDS1"]+0.05, y=NMDS.mean[1,"NMDS2"]+0.2,label=NMDS.mean[1,"group"], size=3) +
  #annotate("text",x=NMDS.mean[3,"NMDS1"], y=NMDS.mean[3,"NMDS2"],label=NMDS.mean[3,"group"], size=3) +
  #annotate("text",x=NMDS.mean[4,"NMDS1"]-0.05, y=NMDS.mean[4,"NMDS2"]+0.05,label=NMDS.mean[4,"group"], size=3) +
    theme_bw() +
    theme(
      
      legend.position = "bottom",
      plot.title=element_text(size=10),
      panel.grid=element_blank(),
      legend.text = element_text(colour = "black",size = 10),
      legend.title = element_blank(),
      axis.text.x = element_text(size = 8,colour = "black"),
      axis.text.y = element_text(size = 8,colour = "black"),
      axis.title.x =element_text(size = 10),
      axis.title.y = element_text(size = 10)) +
  labs(title="Pathology vs. healthy") +
    guides(color = guide_legend(ncol=1))

rm(nmds_model_data, plot, df_ell, nmds_model, NMDS.mean, data)

#### NMDS: SAVE NMDS SERIES ####
nmds_time_arrange <- arrangeGrob(nmds_time_series$plots[[1]] + labs(title="Day -56"),
                                  nmds_time_series$plots[[2]] + labs(title="Day -7"),
                                  nmds_time_series$plots[[3]] + labs(title="Day 21"),
                                  nmds_time_series$plots[[4]] + labs(title="Day 49"), ncol=4)
nmds_ybs_arrange <- arrangeGrob(nmds_ybs_series$plots[[1]] + labs(title="YBS: Day -56 vs. Day -7") + scale_color_manual(values=c("#66a61f", "#e6ab01")),
             nmds_ybs_series$plots[[2]] + labs(title="YBS: Day -7 vs. Day 21") + scale_color_manual(values=c("#e6ab01", "#7570b3")),
             nmds_ybs_series$plots[[3]] + labs(title="YBS: Day 21 vs. Day 49") + scale_color_manual(values=c("#7570b3", "#a65629")), ncol=3)
nmds_wbs_arrange <- arrangeGrob(nmds_wbs_series$plots[[1]] + labs(title="WBS: Day -56 vs. Day -7") + scale_color_manual(values=c("#66a61f", "#e6ab01")),
             nmds_wbs_series$plots[[2]] + labs(title="WBS: Day -7 vs. Day 21")+ scale_color_manual(values=c("#e6ab01", "#7570b3")),
             nmds_wbs_series$plots[[3]] + labs(title="WBS: Day 21 vs. Day 49")+ scale_color_manual(values=c("#7570b3", "#a65629")), ncol=3)
nmds_pathology_arrange <- arrangeGrob(nmds_pathology_series$plots[[1]] + labs(title="WBS: Pathology"),
             nmds_pathology_series$plots[[2]] + labs(title="YBS: Pathology"),
             nmds_pathology_series$plots[[3]] + labs(title="All animals: Pathology"), ncol=3)

ggsave("nmds_time_arrangement.jpg", nmds_time_arrange, width=8.67, height=2.5, units="in", dpi=300)
ggsave("nmds_ybs_arrangement.jpg", nmds_ybs_arrange, width=6.5, height=2.5, units="in", dpi=300)
ggsave("nmds_wbs_arrangement.jpg", nmds_wbs_arrange, width=6.5, height=2.5, units="in", dpi=300)
ggsave("nmds_pathology_arrangement.jpg", nmds_pathology_arrange, width=6.5, height=2.5, units="in", dpi=300)

permanova_series <- rbind(permanova_time_series,
                          permanova_wbs_series,
                          permanova_ybs_series,
                          permanova_pathology_series)
rm(permanova_time_series,
   permanova_wbs_series,
   permanova_ybs_series,
   permanova_pathology_series)

write.csv(permanova_series, "permanova_results.csv")

#### DISTANCE-BASED REDUNDANCY ANALYSIS ####
data <- as.data.frame(otu_tables_rel[[4]])
data[,c(74:81)] <- NULL
dist_matrix_bray <- vegdist(data, method="bray")

# Calculate the dbRDA model.
dbrda.model <- capscale(dist_matrix_bray ~ Diet/Cage + Diet*Time + Pathology + TwinPair/Subject, data = sample_data)
summary(dbrda.model)
anova(dbrda.model, by="term")
RsquareAdj(dbrda.model)$adj.r

# Check adjusted R2 values in a forward selection process
null.model <- capscale(dist_matrix_bray ~ 1, sample_data)
full.model <- capscale(dist_matrix_bray ~ Diet/Cage + Diet*Time + Pathology + TwinPair/Subject, data = sample_data)
selection <- ordiR2step (null.model, scope = formula(full.model), R2scope = adjR2, direction = 'forward', permutations = 999)
selection$anova
anova(null.model, full.model)
rm(null.model, full.model, selection)

# Extract model scores for plotting.
dbrda.model_fort <- scores(dbrda.model, display="sites")
dbrda.model_fort <- cbind(sample_data, dbrda.model_fort)
dbrda.model_fort$Time <- as.factor(as.integer(dbrda.model_fort$Time))
dbrda.model_fort$Diet <- as.factor(as.integer(dbrda.model_fort$Diet))

dbrda.model_fort$DietTime <- paste0(dbrda.model_fort$Diet, "_", dbrda.model_fort$Time)

levels(dbrda.model_fort$Diet)[levels(dbrda.model_fort$Diet)=="0"] <- "WBS"
levels(dbrda.model_fort$Diet)[levels(dbrda.model_fort$Diet)=="1"] <- "YBS"
levels(dbrda.model_fort$Time)[levels(dbrda.model_fort$Time)=="1"] <- "Day -56"
levels(dbrda.model_fort$Time)[levels(dbrda.model_fort$Time)=="2"] <- "Day -7"
levels(dbrda.model_fort$Time)[levels(dbrda.model_fort$Time)=="3"] <- "Day 21"
levels(dbrda.model_fort$Time)[levels(dbrda.model_fort$Time)=="4"] <- "Day 49"

# Plot: Sibling pair
plot1 <- ggplot(dbrda.model_fort, aes(x=CAP1, y=CAP2, color=TwinPair)) + 
  geom_point(size=2) +
  stat_ellipse(aes(group=TwinPair), level=0.5) +
  scale_color_discrete(guide=guide_legend(ncol=2, title.position="top", title="Sibling pair")) +
  labs(title="Sibling pair", x="CAP1 (12.3%)", y="CAP2 (5.8%)") +
  theme_bw() + 
  plot_theme +
  theme(plot.title=element_text(size=10),
        legend.position="bottom",
        legend.text=element_text(size=9),
        legend.title=element_text(hjust=0.5),
        legend.justification=c(0.5,1))

# Plot: Cage
plot2 <- ggplot(subset(dbrda.model_fort, Cage %in% c("C1", "C2", "C3", "C7", "C8", "C9")),
              aes(x=CAP1, y=CAP2, color=Cage)) + geom_point(size=2) +
  stat_ellipse(aes(group=Cage), level=0.5) +   
  scale_color_discrete(guide=guide_legend(ncol=2, title.position="top", title="Cage")) +
  labs(title="Cage", x="CAP1 (12.3%)", y="CAP2 (5.8%)") +
  theme_bw() + 
  plot_theme +
  theme(plot.title=element_text(size=10),
        legend.position="bottom",
        legend.text=element_text(size=9),
        legend.title=element_text(hjust=0.5),
        legend.justification=c(0.5,1))

# Plot: Time
plot3 <- ggplot(dbrda.model_fort, aes(x=CAP1, y=CAP2, color=Time)) + geom_point(size=2) +
  stat_ellipse(aes(group=Time), level=0.5) + 
  scale_color_manual(values=c("#66a61f", "#e6ab01", "#7570b3", "#a65629"),
                     guide=guide_legend(ncol=2, title.position="top", title="Time")) +
  labs(title="Time", x="CAP1 (12.3%)", y="CAP2 (5.8%)") +
  theme_bw() + 
  plot_theme +
  theme(plot.title=element_text(size=10),
        legend.position="bottom",
        legend.text=element_text(size=9),
        legend.title=element_text(hjust=0.5),
        legend.justification=c(0.5,1))

# Plot: Pathology (include "anticipatory pathology" at t3)
subjects <- as.character(subset(dbrda.model_fort, Pathology=="1")$Subject)
dbrda.model_fort$Pathology2 <- dbrda.model_fort$Pathology
dbrda.model_fort$Pathology2[dbrda.model_fort$Subject %in% subjects] <- 1
dbrda.model_fort$Pathology2 <- as.factor(as.character(dbrda.model_fort$Pathology2))

levels(dbrda.model_fort$Pathology2)[levels(dbrda.model_fort$Pathology2)=="0"] <- "healthy"
levels(dbrda.model_fort$Pathology2)[levels(dbrda.model_fort$Pathology2)=="1"] <- "pathology"
dbrda.model_fort$Pathology <- dbrda.model_fort$Pathology2

plot4 <- ggplot(subset(dbrda.model_fort, Time %in% c("Day 21","Day 49")),
       aes(x=CAP1, y=CAP2, color=Pathology)) + geom_point(size=2) +
  stat_ellipse(aes(group=Pathology), level=0.5) + 
  scale_color_manual(values=c("black", "goldenrod"),
                     guide=guide_legend(ncol=1, title.position="top", title="Pathology")) +
  labs(title="Pathology", x="CAP1 (12.3%)", y="CAP2 (5.8%)") +
  theme_bw() + 
  plot_theme +
  theme(plot.title=element_text(size=10),
        legend.position="bottom",
        legend.text=element_text(size=9),
        legend.title=element_text(hjust=0.5),
        legend.justification=c(0.5,1))

# Plot: Diet
plot5 <- ggplot(subset(dbrda.model_fort, Time %in% c("Day -7", "Day 21", "Day 49")),
       aes(x=CAP1, y=CAP2, color=Diet)) + geom_point(size=2) +
  stat_ellipse(aes(group=Diet), level=0.5) +
  scale_color_manual(values=c("deepskyblue3", "orangered1"),
                     guide=guide_legend(ncol=1, title.position="top", title="Diet")) +
  labs(title="Diet", x="CAP1 (12.3%)", y="CAP2 (5.8%)") +
  theme_bw() + 
  plot_theme +
  theme(plot.title=element_text(size=10),
        legend.position="bottom",
        legend.text=element_text(size=9),
        legend.title=element_text(hjust=0.5),
        legend.justification=c(0.5,1))

grid.arrange(plot1 + labs(title="Sibling pair"),
             plot3 + labs(title="Time"),
             plot4 + labs(title="Pathology"),
             plot2 + labs(title="Cage"),
             plot5 + labs(title="Diet"), layout_matrix=rbind(c(1,1,2,2,3,3),
                                                               c(NA,4,4,5,5,NA)))

# Standardize plot margins and sizes.
plot1 <- ggplotGrob(plot1)
plot2 <- ggplotGrob(plot2)
plot3 <- ggplotGrob(plot3)
plot4 <- ggplotGrob(plot4)
plot5 <- ggplotGrob(plot5)

plot3$heights <- plot1$heights
plot2$heights <- plot1$heights
plot4$heights <- plot1$heights
plot5$heights <- plot1$heights

grid.arrange(plot1, plot3, plot4, plot2, plot5,
             layout_matrix=rbind(c(1,1,2,2,3,3),
                                 c(NA,4,4,5,5,NA)))

dbrda_broad_arrange <- arrangeGrob(plot1,
                                   plot3,
                                   plot4,
                                   plot2,
                                   plot5,
                                   layout_matrix=rbind(c(1,1,2,2,3,3),
                                                       c(NA,4,4,5,5,NA)))


ggsave("dbrda_broad_arrangement.jpg", dbrda_broad_arrange,
       width=6.5, height=7, units="in", dpi=300)

rm(plot1, plot2, plot3, plot4, plot5)

################################################# PART FOUR: STATISTICAL TESTS ###################################################
#### ALPHA DIVERSITY ####
# Significant differences between each diet, at each timepoint.
alpha_div_significance <- data.frame(matrix(nrow=1, ncol=12))
colnames(alpha_div_significance) <- c("Time", "Diet", "Measure",
                                      "LeveneF", "Levenedf", "Levenep",
                                      "Student", "Studentdf", "Studentp",
                                      "Welcht", "Welchdf", "Welchp")

vector <- c("diversity_inverse_simpson", "observed", "evenness_pielou", "BraySmall")

# Comparisons between diets at each time
for(p in c(1:4)){ # For each time
  for(x in c(1:4)){ # For each measurement
    model_data <- subset(sample_data, Time==p)
    levene <- car::leveneTest(model_data[,vector[x]], as.factor(as.character(model_data$Diet)))
    student <- t.test(model_data[,vector[x]] ~ as.factor(as.character(model_data$Diet)), var.equal=TRUE)
    welch <- t.test(model_data[,vector[x]] ~ as.factor(as.character(model_data$Diet)), var.equal=FALSE)
    
    if(p==1){m <- 0}
    if(p==2){m <- 4}
    if(p==3){m <- 8}
    if(p==4){m <- 12}
    
    alpha_div_significance[x+m,1] <- x # Time
    alpha_div_significance[x+m,2] <- NA # Skip
    alpha_div_significance[x+m,3] <- vector[x] # Measure
    alpha_div_significance[x+m,4] <- levene[1,2] # Levene F
    alpha_div_significance[x+m,5] <- levene[1,1] # Levene df
    alpha_div_significance[x+m,6] <- levene[1,3] # Levene p
    alpha_div_significance[x+m,7] <- student[["statistic"]][["t"]] # Student t
    alpha_div_significance[x+m,8] <- student[["parameter"]][["df"]] # Student df
    alpha_div_significance[x+m,9] <- student$p.value # Student p
    alpha_div_significance[x+m,10] <- welch[["statistic"]][["t"]]
    alpha_div_significance[x+m,11] <- welch[["parameter"]][["df"]]
    alpha_div_significance[x+m,12] <- welch$p.value
    
    rm(model_data, student, welch, levene)
  }
}
  
# Comparisons to the next time for each diet
temp_storage1 <- data.frame(matrix(nrow=1, ncol=12))
temp_storage2 <- data.frame(matrix(nrow=1, ncol=12))

for(p in c(0:1)){ # For each diet
  
  for(m in c(1:3)){ # For each set of two timepoints
    for(x in c(1:4)){ # For each measurement
      model_data <- subset(sample_data, Diet==p & Time %in% c(m, m+1))
      
      levene <- car::leveneTest(model_data[,vector[x]], as.factor(as.character(model_data$Time)))
      student <- t.test(model_data[,vector[x]] ~ as.factor(as.character(model_data$Time)), var.equal=TRUE)
      welch <- t.test(model_data[,vector[x]] ~ as.factor(as.character(model_data$Time)), var.equal=FALSE)
      
      if(m==1){y <- 0}
      if(m==2){y <- 4}
      if(m==3){y <- 8}
      
      if(p==0){
        temp_storage1[x + m-1 + y, 1] <- paste0(m, ",", m+1)
        temp_storage1[x + m-1 + y, 2] <- p
        temp_storage1[x + m-1 + y, 3] <- vector[x] # Measure
        temp_storage1[x + m-1 + y, 4] <- levene[1,2] # Levene F
        temp_storage1[x + m-1 + y, 5] <- levene[1,1] # Levene df
        temp_storage1[x + m-1 + y, 6] <- levene[1,3] # Levene p
        temp_storage1[x + m-1 + y, 7] <- student[["statistic"]][["t"]] # Student t
        temp_storage1[x + m-1 + y, 8] <- student[["parameter"]][["df"]] # Student df
        temp_storage1[x + m-1 + y, 9] <- student$p.value # Student p
        temp_storage1[x + m-1 + y, 10] <- welch[["statistic"]][["t"]]
        temp_storage1[x + m-1 + y, 11] <- welch[["parameter"]][["df"]]
        temp_storage1[x + m-1 + y, 12] <- welch$p.value
      }
      
      if(p==1){
        temp_storage2[x + m-1 + y, 1] <- paste0(m, ",", m+1)
        temp_storage2[x + m-1 + y, 2] <- p
        temp_storage2[x + m-1 + y, 3] <- vector[x] # Measure
        temp_storage2[x + m-1 + y, 4] <- levene[1,2] # Levene F
        temp_storage2[x + m-1 + y, 5] <- levene[1,1] # Levene df
        temp_storage2[x + m-1 + y, 6] <- levene[1,3] # Levene p
        temp_storage2[x + m-1 + y, 7] <- student[["statistic"]][["t"]] # Student t
        temp_storage2[x + m-1 + y, 8] <- student[["parameter"]][["df"]] # Student df
        temp_storage2[x + m-1 + y, 9] <- student$p.value # Student p
        temp_storage2[x + m-1 + y, 10] <- welch[["statistic"]][["t"]]
        temp_storage2[x + m-1 + y, 11] <- welch[["parameter"]][["df"]]
        temp_storage2[x + m-1 + y, 12] <- welch$p.value
      }
      
      rm(model_data, student, welch, levene)
    }
  }
}

temp_storage <- rbind(temp_storage1, temp_storage2)
rm(temp_storage1, temp_storage2)
colnames(temp_storage) <- colnames(alpha_div_significance)

alpha_div_significance <- rbind(alpha_div_significance,
                                temp_storage)

openxlsx::write.xlsx(alpha_div_significance, "~/Documents/marmosets/signif_tests_adiv.xlsx", rownames=FALSE)

rm(temp_storage)

#### TAXON RELATIVE ABUNDANCE ####
taxa_numbers <- c(7, 21, 31, 73)

# -WBS v. YBS at each time point (separated by pathology) ####
taxon_diff_abund_by_diet_at_each_time <- list()

for(j in c(1:4)){ # For each taxonomic level...
  df_master <- data.frame(matrix(nrow=taxa_numbers[j], ncol=0))
  
  for(p in c(1:5)){ # For each of the five tests between diets
    model_data <- otu_tables_rel[[j]] # Extract the OTU table used for testing
    
    df <- data.frame(matrix(nrow=taxa_numbers, ncol=3))
    
    # Subset to the items being compared
    if(p==1){model_data <- subset(model_data, Time==1)}
    if(p==2){model_data <- subset(model_data, Time==2)}
    if(p==3){model_data <- subset(model_data, Time==3)}
    if(p==4){model_data <- subset(model_data, Time==4 & Pathology==0)}
    if(p==5){model_data <- subset(model_data, Time==4 & Pathology==1)}
    
    for(m in c(1:taxa_numbers[j])){ # For each taxon in the OTU table, calculate differential abundance
      levene <- car::leveneTest(model_data[,m] ~ as.factor(as.character(model_data$Diet)))
      student <- t.test(model_data[,m] ~ as.factor(as.character(model_data$Diet)), var.equal=TRUE)
      
      df[m,1] <- levene[1,3]
      df[m,2] <- student$p.value
    }
    
    df[,3] <- p.adjust(df[,2], method="BH")
    df_master <- cbind(df_master, df)
    rm(df)
  }
  rownames(df_master) <- colnames(model_data)[c(1:taxa_numbers[j])]
  colnames(df_master) <- c("lev_t1_p0", "stu_t1_p0", "bh_t1_p0",
                           "lev_t2_p0", "stu_t2_p0", "bh_t2_p0",
                           "lev_t3_p0", "stu_t3_p0", "bh_t3_p0",
                           "lev_t4_p0", "stu_t4_p0", "bh_t4_p0",
                           "lev_t4_p1", "stu_t4_p1", "bh_t4_p1")
                           
  
  df_master <- merge(tax_tables[[j]], df_master, by=0, all=FALSE)
  
  taxon_diff_abund_by_diet_at_each_time[[j]] <- df_master
  rm(df_master)
}

# -Healthy vs. sick for each diet ####
taxon_diff_abund_by_pathology_within_diets <- list()

for(j in c(1:4)){ # For each taxonomic level...
  df_master <- data.frame(matrix(nrow=taxa_numbers[j], ncol=0))
  
  for(p in c(1:2)){ # For each of the five tests between diets
    model_data <- otu_tables_rel[[j]] # Extract the OTU table used for testing
    
    df <- data.frame(matrix(nrow=taxa_numbers, ncol=3))
    
    # Subset to the items being compared
    if(p==1){model_data <- subset(model_data, Time==4 & Diet==0)}
    if(p==2){model_data <- subset(model_data, Time==4 & Diet==1)}
    
    for(m in c(1:taxa_numbers[j])){ # For each taxon in the OTU table, calculate differential abundance
      levene <- car::leveneTest(model_data[,m] ~ as.factor(as.character(model_data$Pathology)))
      student <- t.test(model_data[,m] ~ as.factor(as.character(model_data$Pathology)), var.equal=TRUE)
      
      df[m,1] <- levene[1,3]
      df[m,2] <- student$p.value
    }
    
    df[,3] <- p.adjust(df[,2], method="BH")
    df_master <- cbind(df_master, df)
    rm(df)
  }
  rownames(df_master) <- colnames(model_data)[c(1:taxa_numbers[j])]
  colnames(df_master) <- c("lev_d0_t4", "stu_d0_t4", "bh_d0_t4",
                           "lev_d1_t4", "stu_d1_t4", "bh_d1_t4")
  
  df_master <- merge(tax_tables[[j]], df_master, by=0, all=FALSE)
  
  taxon_diff_abund_by_pathology_within_diets[[j]] <- df_master
  rm(df_master)
}

# -Time vs. previous time for each diet ####
taxon_diff_abund_by_time_within_diets <- list()

for(j in c(1:4)){ # For each taxonomic level...
  df_master <- data.frame(matrix(nrow=taxa_numbers[j], ncol=0))
  
  for(p in c(1:8)){ # For each of the five tests between diets
    model_data <- otu_tables_rel[[j]] # Extract the OTU table used for testing
    
    df <- data.frame(matrix(nrow=taxa_numbers, ncol=3))
    
    # Subset to the items being compared
    if(p==1){model_data <- subset(model_data, Diet==0 & Time %in% c(1,2))}
    if(p==2){model_data <- subset(model_data, Diet==0 & Time %in% c(2,3))}
    if(p==3){model_data <- subset(model_data, Diet==0 & Time %in% c(3,4) & Pathology==0)}
    if(p==4){model_data <- subset(model_data, (Diet==0 & Time==3) | (Diet==0 & Time==4 & Pathology==1))}
    if(p==5){model_data <- subset(model_data, Diet==1 & Time %in% c(1,2))}
    if(p==6){model_data <- subset(model_data, Diet==1 & Time %in% c(2,3))}
    if(p==7){model_data <- subset(model_data, Diet==1 & Time %in% c(3,4) & Pathology==0)}
    if(p==8){model_data <- subset(model_data, (Diet==1 & Time==3) | (Diet==1 & Time==4 & Pathology==1))}
    
    for(m in c(1:taxa_numbers[j])){ # For each taxon in the OTU table, calculate differential abundance
      levene <- car::leveneTest(model_data[,m] ~ as.factor(as.character(model_data$Time)))
      student <- t.test(model_data[,m] ~ as.factor(as.character(model_data$Time)), var.equal=TRUE)
      
      df[m,1] <- levene[1,3]
      df[m,2] <- student$p.value
    }
    
    df[,3] <- p.adjust(df[,2], method="BH")
    df_master <- cbind(df_master, df)
    rm(df)
  }
  rownames(df_master) <- colnames(model_data)[c(1:taxa_numbers[j])]
  colnames(df_master) <- c("lev_d0_t1v2", "stu_d0_t1v2", "bh_d0_t1v2",
                           "lev_d0_t2v3", "stu_d0_t2v3", "bh_d0_t2v3",
                           "lev_d0_t3v4_p0", "stu_d0_t3v4_p0", "bh_d0_t3v4_p0",
                           "lev_d0_t3v4_p1", "stu_d0_t3v4_p1", "bh_d0_t3v4_p1",
                           "lev_d1_t1v2", "stu_d1_t1v2", "bh_d1_t1v2",
                           "lev_d1_t2v3", "stu_d1_t2v3", "bh_d1_t2v3",
                           "lev_d1_t3v4_p0", "stu_d1_t3v4_p0", "bh_d1_t3v4_p0",
                           "lev_d1_t3v4_p1", "stu_d1_t3v4_p1", "bh_d1_t3v4_p1")
  
  df_master <- merge(tax_tables[[j]], df_master, by=0, all=FALSE)
  
  taxon_diff_abund_by_time_within_diets[[j]] <- df_master
  rm(df_master)
}



names(taxon_diff_abund_by_diet_at_each_time) <- c("Phylum", "Family", "Genus", "OTU")
names(taxon_diff_abund_by_pathology_within_diets) <- c("Phylum", "Family", "Genus", "OTU")
names(taxon_diff_abund_by_time_within_diets) <- c("Phylum", "Family", "Genus", "OTU")

openxlsx::write.xlsx(taxon_diff_abund_by_diet_at_each_time,
                     "~/Documents/marmosets/signif_tests_by_diet_at_each_time.xlsx", rownames=FALSE)
openxlsx::write.xlsx(taxon_diff_abund_by_pathology_within_diets,
                     "~/Documents/marmosets/signif_tests_by_pathology_within_diets.xlsx", rownames=FALSE)
openxlsx::write.xlsx(taxon_diff_abund_by_time_within_diets,
                     "~/Documents/marmosets/signif_tests_by_time_within_diets.xlsx", rownames=FALSE)

#### SIGNIFICANTLY DIFFERENT DELTAS (?) ####
# -WBS v. YBS at each delta (separated by pathology) ####
delta_diff_abund_by_diet_at_each_time <- list()

for(j in c(1:4)){ # For each taxonomic level...
  df_master <- data.frame(matrix(nrow=taxa_numbers[j], ncol=0))
  
  for(p in c(1:4)){ # For each of the five tests between diets
    model_data <- otu_tables_delta[[j]] # Extract the OTU table used for testing
    model_data$Subject <- NULL
    
    df <- data.frame(matrix(nrow=taxa_numbers, ncol=3))
    
    # Subset to the items being compared
    if(p==1){model_data <- subset(model_data, Difference=="1_2")}
    if(p==2){model_data <- subset(model_data, Difference=="2_3")}
    if(p==3){model_data <- subset(model_data, Difference=="3_4" & Pathology==0)}
    if(p==4){model_data <- subset(model_data, Difference=="3_4" & Pathology==1)}
    
    for(m in c(1:taxa_numbers[j])){ # For each taxon in the OTU table, calculate differential abundance
      levene <- car::leveneTest(model_data[,m] ~ as.factor(as.character(model_data$Diet)))
      student <- t.test(model_data[,m] ~ as.factor(as.character(model_data$Diet)), var.equal=TRUE)
      
      df[m,1] <- levene[1,3]
      df[m,2] <- student$p.value
    }
    
    df[,3] <- p.adjust(df[,2], method="BH")
    df_master <- cbind(df_master, df)
    rm(df)
  }
  rownames(df_master) <- colnames(model_data)[c(1:taxa_numbers[j])]
  colnames(df_master) <- c("lev_t12_p0", "stu_t12_p0", "bh_t12_p0",
                           "lev_t23_p0", "stu_t23_p0", "bh_t23_p0",
                           "lev_t34_p0", "stu_t34_p0", "bh_t34_p0",
                           "lev_t34_p1", "stu_t34_p1", "bh_t34_p1")
  
  
  df_master <- merge(tax_tables[[j]], df_master, by=0, all=FALSE)
  
  delta_diff_abund_by_diet_at_each_time[[j]] <- df_master
  rm(df_master)
}

# -Healthy vs. sick for each diet ####
delta_diff_abund_by_pathology_within_diets <- list()

for(j in c(1:4)){ # For each taxonomic level...
  df_master <- data.frame(matrix(nrow=taxa_numbers[j], ncol=0))
  
  for(p in c(1:2)){ # For each of the five tests between diets
    model_data <- otu_tables_delta[[j]] # Extract the OTU table used for testing
    model_data$Subject <- NULL
    
    df <- data.frame(matrix(nrow=taxa_numbers, ncol=3))
    
    # Subset to the items being compared
    if(p==1){model_data <- subset(model_data, Difference=="3_4" & Diet==0)}
    if(p==2){model_data <- subset(model_data, Difference=="3_4" & Diet==1)}
    
    for(m in c(1:taxa_numbers[j])){ # For each taxon in the OTU table, calculate differential abundance
      levene <- car::leveneTest(model_data[,m] ~ as.factor(as.character(model_data$Pathology)))
      student <- t.test(model_data[,m] ~ as.factor(as.character(model_data$Pathology)), var.equal=TRUE)
      
      df[m,1] <- levene[1,3]
      df[m,2] <- student$p.value
    }
    
    df[,3] <- p.adjust(df[,2], method="BH")
    df_master <- cbind(df_master, df)
    rm(df)
  }
  rownames(df_master) <- colnames(model_data)[c(1:taxa_numbers[j])]
  colnames(df_master) <- c("lev_d0_t34", "stu_d0_t34", "bh_d0_t34",
                           "lev_d1_t34", "stu_d1_t34", "bh_d1_t34")
  
  df_master <- merge(tax_tables[[j]], df_master, by=0, all=FALSE)
  
  delta_diff_abund_by_pathology_within_diets[[j]] <- df_master
  rm(df_master)
}


# -Time vs. previous time for each diet ####
delta_diff_abund_by_time_within_diets <- list()

for(j in c(1:4)){ # For each taxonomic level...
  df_master <- data.frame(matrix(nrow=taxa_numbers[j], ncol=0))
  
  for(p in c(1:6)){ # For each of the five tests between diets
    model_data <- otu_tables_delta[[j]] # Extract the OTU table used for testing
    model_data$Subject <- NULL
    
    df <- data.frame(matrix(nrow=taxa_numbers, ncol=3))
    
    # Subset to the items being compared
    if(p==1){model_data <- subset(model_data, Diet==0 & Difference %in% c("1_2","2_3"))}
    if(p==2){model_data <- subset(model_data, (Diet==0 & Difference=="2_3") | (Diet==0 & Difference=="3_4" & Pathology==0))}
    if(p==3){model_data <- subset(model_data, (Diet==0 & Difference=="2_3") | (Diet==0 & Difference=="3_4" & Pathology==1))}
    if(p==4){model_data <- subset(model_data, Diet==1 & Difference %in% c("1_2","2_3"))}
    if(p==5){model_data <- subset(model_data, (Diet==1 & Difference=="2_3") | (Diet==1 & Difference=="3_4" & Pathology==0))}
    if(p==6){model_data <- subset(model_data, (Diet==1 & Difference=="2_3") | (Diet==1 & Difference=="3_4" & Pathology==1))}
    
    model_data$Difference <- factor(model_data$Difference)
    
    for(m in c(1:taxa_numbers[j])){ # For each taxon in the OTU table, calculate differential abundance
      levene <- car::leveneTest(model_data[,m] ~ as.factor(as.character(model_data$Difference)))
      student <- t.test(model_data[,m] ~ as.factor(as.character(model_data$Difference)), var.equal=TRUE)
      
      df[m,1] <- levene[1,3]
      df[m,2] <- student$p.value
    }
    
    df[,3] <- p.adjust(df[,2], method="BH")
    df_master <- cbind(df_master, df)
    rm(df)
  }
  rownames(df_master) <- colnames(model_data)[c(1:taxa_numbers[j])]
  colnames(df_master) <- c("lev_d0_12v23", "stu_d0_12v23", "bh_d0_12v23",
                           "lev_d0_23v34_p0", "stu_d0_23v34_p0", "bh_d0_23v34_p0",
                           "lev_d0_23v34_p1", "stu_d0_23v34_p1", "bh_d0_23v34_p1",
                           "lev_d1_12v23", "stu_d1_12v23", "bh_d1_12v23",
                           "lev_d1_23v34_p0", "stu_d1_23v34_p0", "bh_d1_23v34_p0",
                           "lev_d1_23v34_p1", "stu_d1_23v34_p1", "bh_d1_23v34_p1")
  
  df_master <- merge(tax_tables[[j]], df_master, by=0, all=FALSE)
  
  delta_diff_abund_by_time_within_diets[[j]] <- df_master
  rm(df_master)
}



names(delta_diff_abund_by_diet_at_each_time) <- c("Phylum", "Family", "Genus", "OTU")
names(delta_diff_abund_by_pathology_within_diets) <- c("Phylum", "Family", "Genus", "OTU")
names(delta_diff_abund_by_time_within_diets) <- c("Phylum", "Family", "Genus", "OTU")

openxlsx::write.xlsx(delta_diff_abund_by_diet_at_each_time,
                     "~/Documents/marmosets/delta_signif_tests_by_diet_at_each_time.xlsx", rownames=FALSE)
openxlsx::write.xlsx(delta_diff_abund_by_pathology_within_diets,
                     "~/Documents/marmosets/delta_signif_tests_by_pathology_within_diets.xlsx", rownames=FALSE)
openxlsx::write.xlsx(delta_diff_abund_by_time_within_diets,
                     "~/Documents/marmosets/delta_signif_tests_by_time_within_diets.xlsx", rownames=FALSE)

################################################# PART FIVE: REMAINING FIGURES ###################################################
#### FIGURE 1: STACKED BAR CHART ####
fig1_data <- data.frame(
  Phylum=c("Actinobacteria","Firmicutes","Bacteroidetes","Fusobacteria","Proteobacteria","Other"),
  Kap=c(66,20,11,0,1,2),
  Current=c(69,13.7,14,1.8,1.1,0.4)
)

levels(fig1_data$Phylum) <- c("Actinobacteria","Firmicutes","Bacteroidetes","Fusobacteria","Proteobacteria","Other")

fig1_data <- melt(fig1_data)
fig1_data$axis <- rep(1)

levels(fig1_data$variable)[levels(fig1_data$variable)=="Kap"] <- "Kap et al.\n(2018)"
levels(fig1_data$variable)[levels(fig1_data$variable)=="Current"] <- "This\nreview"

rel_abund_plot <- ggplot(fig1_data, aes(x=variable, y=value, fill=Phylum)) + 
  geom_bar(stat="identity", position=position_fill(reverse=TRUE)) +
  scale_fill_manual(values=c("#005ea5", "#ababab", "#f68b2c", "#a65629", "#9969a7", "#00b1e3")) +
  #facet_grid(~variable) + 
  scale_y_continuous(expand=c(0,0)) +
  theme_bw() + 
  plot_theme + 
  theme(axis.text.x=element_text(size=10),
        axis.title.x=element_blank(),
        strip.text=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y="Relative abundance (%)")

ggsave("~/Fig1a_labels_on_bottom.jpg", rel_abund_plot, width=4, height=3, units="in", dpi=300)

#### FIGURE 4: ALPHA DIVERSITY BOX PLOTS ####
fig2_data <- reshape2::melt(sample_data, id.vars=c("SampleID", "Subject", "Diet","Time", "Pathology"),
                  measure.vars=c("observed", "diversity_inverse_simpson","evenness_pielou","BraySmall"))

fig2_data$Time <- as.factor(fig2_data$Time)
fig2_data$Diet <- as.factor(fig2_data$Diet)
fig2_data$Pathology <- factor(fig2_data$Pathology)
levels(fig2_data$Pathology)[levels(fig2_data$Pathology)=="0"] <- "healthy"
levels(fig2_data$Pathology)[levels(fig2_data$Pathology)=="1"] <- "pathology"
levels(fig2_data$Diet)[levels(fig2_data$Diet)=="0"] <- "WBS"
levels(fig2_data$Diet)[levels(fig2_data$Diet)=="1"] <- "YBS"

fig2_data <- subset(fig2_data, variable %in% c("observed","diversity_inverse_simpson","evenness_pielou", "BraySmall"))
levels(fig2_data$variable)[levels(fig2_data$variable)=="observed"] <- "Species\nrichness"
levels(fig2_data$variable)[levels(fig2_data$variable)=="diversity_inverse_simpson"] <- "Simpson\ndiversity"
levels(fig2_data$variable)[levels(fig2_data$variable)=="evenness_pielou"] <- "Pielou\nevenness"
levels(fig2_data$variable)[levels(fig2_data$variable)=="BraySmall"] <- "Bray-Curtis\ndissimilarity"
fig2_data$variable <- factor(fig2_data$variable, levels=c("Simpson\ndiversity",
                                                          "Species\nrichness",
                                                          "Pielou\nevenness",
                                                          "Bray-Curtis\ndissimilarity"))

levels(fig2_data$Time)[levels(fig2_data$Time)==1] <- "Day -56"
levels(fig2_data$Time)[levels(fig2_data$Time)==2] <- "Day -7"
levels(fig2_data$Time)[levels(fig2_data$Time)==3] <- "Day 21"
levels(fig2_data$Time)[levels(fig2_data$Time)==4] <- "Day 49"

fig2_data$value[fig2_data$variable=="Bray-Curtis\ndissimilarity" & fig2_data$value > 0.8] <- NA

fig2_data$PlotVar <- paste0(fig2_data$Diet, fig2_data$Pathology)
levels

fig2 <- ggplot(fig2_data, aes(x=Time, y=value)) +
  #geom_boxplot(aes(fill=Diet)) +
  geom_boxplot_pattern(aes(fill=Diet, pattern=Pathology)) +
  geom_point(aes(x=Time, y=1.05*value), color="transparent") +
  #scale_fill_manual(values=c("deepskyblue3", "orangered1")) +
  facet_wrap(~variable, scales="free") +
  scale_pattern_type_discrete(choices=c("a", "b")) +
  theme_bw() + plot_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        panel.grid=element_blank(),
        legend.title=element_text(hjust=0.5))

ggsave("~/Documents/marmosets/figures_2020.11.20/diversity_barplot.jpg", fig2, width=6.5, height=6, units="in")



#### SAVE DATA ####
mixed_model_coefficients <- list(LogAbundance = megaplot_actual_abundance,
                                 DeltaAbundance = megaplot_delta_abundance)
openxlsx::write.xlsx(mixed_model_coefficients, "~/Documents/marmosets/figures_2020.11.26/mixed_model_coefficients.xlsx")
rm(mixed_model_coefficients)

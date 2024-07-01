#Load library
library(tidyverse)
library(data.table)
options(scipen=999)


#Iterate through the summary files, compute SD, and output 
setwd("~/Documents/mPCR/data/simulations/results")
fnames = list.files(pattern="*ldPruned_classification.out", full.names=TRUE, recursive=FALSE) #input vcfs

df_list <- list() #store data frames

for (query in seq_along(fnames)) {
  
  #read file in
  data = read_delim(file = fnames[query], delim = "\t")
  
  #compute accuracy and SD
  data$accuracy <- NA
  data$SD <- NA
  data[data$Type == 'T', 'accuracy'] <- data[data$Type == 'T', 'Replicate 1']
  data[data$Type == 'T', 'SD'] <- 0
  replicate_cols <- grep("Replicate", names(data), value = TRUE)
  data[data$Type == 'R', 'accuracy'] <- rowMeans(data[data$Type == 'R', replicate_cols], na.rm = TRUE)
  data[data$Type == 'R', 'SD'] <- apply(data[data$Type == 'R', replicate_cols], 1, sd, na.rm = TRUE)
  first_three <- data[, 1:3]
  last_two <- data[, c('accuracy', 'SD')]
  data <- cbind(first_three, last_two) #make new df with accuracy and sd
  data$misclassification = 1-data$accuracy
  
  #add columns for fst and data type and jitter for plotting
  data$FST = gsub("_samp100each_ldPruned_classification.out", "", unlist(str_split(fnames[query], "FST"))[2], )
  
  #append to list of data frames
  df_list <- append(df_list, list(data))
}

#Combine for plotting
merge = rbindlist(df_list)
set.seed(120)  # reproducibility
merge$jittered_x <- jitter(as.numeric(merge$`Num.of.Markers`), amount = 0.3)

#write.table(merge, 
#            file = "merged_classification_plotResults.txt", 
#            sep = "\t", quote = FALSE, row.names = FALSE)


#Plot
ggplot(merge %>% filter(Type == "R"), aes(x = `Num.of.Markers`, y = `Num.of.Individuals`, fill = misclassification)) +
  geom_tile(color = "white") +
  geom_text(aes(label=round(misclassification, 2))) +
  facet_wrap(~ FST, nrow = 2) +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +  # Adjust breaks as needed
  scale_y_continuous(breaks = seq(10, 50, by = 10)) +  # Adjust breaks as needed
  scale_fill_gradientn(colours = colorspace::heat_hcl(10),
                       breaks = seq(0,0.4, by = 0.08), 
                       labels = scales::label_number(accuracy = 0.04), limits = c(0, 0.4)) +
  labs(x = "Number of Markers", y = "Number of Individuals", fill = "Misclassification\nRate") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) 


ggplot(merge %>% filter(Type == "R"), aes(x = `Num.of.Markers`, y = `Num.of.Individuals`, fill = misclassification)) +
  geom_tile(color = "white") +
  facet_wrap(~ FST, nrow = 2) +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +  # Adjust breaks as needed
  scale_y_continuous(breaks = seq(10, 50, by = 10)) +  # Adjust breaks as needed
  scale_fill_gradientn(colours = colorspace::heat_hcl(10),
                       breaks = seq(0,0.4, by = 0.08), 
                       labels = scales::label_number(accuracy = 0.08), limits = c(0, 0.4)) +
  labs(x = "Number of Markers", y = "Number of Individuals", fill = "Misclassification\nRate") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) 

#check out summary data
merge %>%
  filter(Type == "R") %>%
  group_by(FST, Num.of.Markers) %>%
  summarise_at(vars(misclassification), list(Mean = mean, Min = min, Max = max)) %>%
  mutate_if(is.numeric, round, digits=3)  %>% 
  arrange(-desc(FST)) %>% View()
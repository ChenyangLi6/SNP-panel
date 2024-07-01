#load libraries
options(scipen=999)
suppressPackageStartupMessages(library(optparse)) 
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(caret))


#split inputs from cluster
###Use input parameters from command line, if not found set to default
option_list <- list( 
  make_option("--numberSNPs", type="integer",default=4000, 
              help="number of available SNPs"),
  make_option("--numIndivs", type="integer",default=100, 
              help="number of samples total"),
  make_option("--ArrayID", type="integer", 
              help="ID for job array"),
  make_option("--inFilePath", type="character", default="/project/jazlynmo_738/Jazlyn/mpcr/sim_vcfs/",
              help="specify input file path", metavar="character")
)

#Parse the parameters and assign to variables
opt = parse_args(OptionParser(option_list=option_list)) #this must come first before we assign 
nsites = opt$numberSNPs
nsamps = opt$numIndivs
query = opt$ArrayID #set the array to end at the total number of files we need to iterate through, ArrayID = query 
infilePath = opt$inFilePath
fnames = list.files(path=infilePath, pattern="*_ldPruned\\.raw", full.names=TRUE, recursive=FALSE) #input vcfs only need to do for pop1 because we will use gsub to get pop2 dfs

#iterate through the file
#for (query in seq_along(fnames)){
  dfpopA = read_delim(file = fnames[query], delim = "\t") %>% 
    select(-c(1:6)) %>% #select the first hundred indivs and remove columns that aren't markers
    slice(1:100) #first 100 are pop1 
  
  dfpopB = read_delim(file = fnames[query], delim = "\t") %>% 
    select(-c(1:6)) %>%  #select the first hundred indivs and remove columns that aren't markers
    slice(101:200)  #second 100 are pop2 
  
  #Compute allele frequencies
  markers <- colnames(dfpopA)
  frequencies <- as.numeric(colSums(dfpopA) / (2 * nsamps))
  freq_popA <- data.frame(markers, frequencies)
  
  markers2 <- markers
  frequencies2 <- as.numeric(colSums(dfpopB) / (2 * nsamps))
  freq_popB <- data.frame(markers2, frequencies2)
  
  #Compute FST
  Fst_list <- vector("numeric", length = nsites)
  for(i in 1:nsites){
    p1 = freq_popA[i,2]
    q1 = 1-p1
    p2 = freq_popB[i,2]
    q2 = 1-p2
    pbar = (p1+p2)/2
    qbar = (q1+q2)/2
    varp = ((pbar-p1)^2 + (pbar-p2)^2)/2
    Fst <- varp/(pbar * qbar)
    Fst_list[i] <- Fst
  }
  
  #clean up FST
  nan_positions <- which(is.nan(Fst_list))
  Fst_list[nan_positions] <- 0
  top_indices<- order(unlist(Fst_list), decreasing = TRUE)[1:100]#top 100 markers
  
  #top marker data frame
  dfFST = as.data.frame(Fst_list) %>%
    rownames_to_column("numMarker") %>% 
    mutate(topMarker = ifelse(numMarker %in% top_indices, 1, 0)) 

  #write FST out
  outfileFST = gsub("ldPruned.raw", "ldPruned_classification.fst", fnames[query]) #outfile name
  write.table(dfFST, outfileFST, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  #Compute delta for forca
  deltapopA = as.data.frame(map2_dfc(dfpopA, rep(1, nsites), ~ ifelse(.x == 1, 0, 1)))
  deltapopB = as.data.frame(map2_dfc(dfpopB, rep(1, nsites), ~ ifelse(.x == 1, 0, 1)))
  
  #Run the classification with random and top markers
  dcontrol <- 20 #number of replicates
  k = 1 # control number of individuals
  num_individuals <- rep(c(10, 20, 30, 40, 50), each = 10) #random sample of individuals
  num_markers <- rep(c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), times = 5) #markers to select
  type <- rep(c("T", "R"), each = 50) #top or random markers
  accuracy <- rep(NA, 50)  #vector for accuracy
  result <- data.frame(`Num of Individuals` = num_individuals,`Num of Markers` = num_markers,`Type` = type, `Accuracy` = accuracy) #data frame to store results
  result$Accuracy <- NULL #set this to null
  
  replicate_columns <- paste0("Replicate ", 1:20)
  for(col in replicate_columns) {
    result[[col]] <- NA
  }
  
  for(w in 1:5){
    k = w*10 
    print(paste("num of individuals:", k))
    for(x in 1:10){
      m <- x*10
      print("Top")
      print(paste("num of Markers:", m))
      ppl1value = numeric(k)
      for(i in 1:k){
        ppl1value[i] = 1
      }
      
      for(i in 1:k){
        for(j in 1:m){
          if(dfpopA[i,top_indices[j]] == 2){
            ppl1value[i] = ppl1value[i]*(2-deltapopA[i,top_indices[j]])*freq_popA[top_indices[j], 2]*freq_popA[top_indices[j], 2]
          }
          else if(dfpopA[i,top_indices[j]]== 0){
            ppl1value[i] = ppl1value[i]*(2-deltapopA[i,top_indices[j]])*(1-freq_popA[top_indices[j], 2])*(1-freq_popA[top_indices[j], 2])
          }
          else {
            ppl1value[i] = ppl1value[i]*(2-deltapopA[i,top_indices[j]])*freq_popA[top_indices[j], 2]*(1-freq_popA[top_indices[j], 2])
          }
        }
      }
      ppl2value = numeric(k)
      for(i in 1:k){
        ppl2value[i] = 1
      }
      for(i in 1:k){
        for(j in 1:m){
          if(dfpopA[i,top_indices[j]] == 2){
            ppl2value[i] = ppl2value[i]*(2-deltapopA[i,top_indices[j]])*freq_popB[top_indices[j], 2]*freq_popB[top_indices[j], 2]
          }
          else if(dfpopA[i,top_indices[j]] == 0){
            ppl2value[i] = ppl2value[i]*(2-deltapopA[i,top_indices[j]])*(1-freq_popB[top_indices[j], 2])*(1-freq_popB[top_indices[j], 2])
          }
          else {
            ppl2value[i] = ppl2value[i]*(2-deltapopA[i,top_indices[j]])*freq_popB[top_indices[j], 2]*(1-freq_popB[top_indices[j], 2])
          }
        }
      }
      matrix1 <- matrix(1, nrow = 1, ncol = k)
      for(i in 1:k){
        if(ppl1value[i] >= ppl2value[i]){
          matrix1[i] = 1
        }else
          matrix1[i] = 2
      }
      ppl1value = numeric(k)
      for(i in 1:k){
        ppl1value[i] = 1
      }
      for(i in 1:k){
        for(j in 1:m){
          if(dfpopB[i,top_indices[j]] == 2){
            ppl1value[i] = ppl1value[i]*(2-deltapopB[i,top_indices[j]])*freq_popA[top_indices[j], 2]*freq_popA[top_indices[j], 2]
          }
          else if(dfpopB[i,top_indices[j]]== 0){
            ppl1value[i] = ppl1value[i]*(2-deltapopB[i,top_indices[j]])*(1-freq_popA[top_indices[j], 2])*(1-freq_popA[top_indices[j], 2])
          }
          else {
            ppl1value[i] = ppl1value[i]*(2-deltapopB[i,top_indices[j]])*freq_popA[top_indices[j], 2]*(1-freq_popA[top_indices[j], 2])
          }
        }
      }
      ppl2value = numeric(k)
      for(i in 1:k){
        ppl2value[i] = 1
      }
      for(i in 1:k){
        for(j in 1:m){
          if(dfpopB[i,top_indices[j]] == 2){
            ppl2value[i] = ppl2value[i]*(2-deltapopB[i,top_indices[j]])*freq_popB[top_indices[j], 2]*freq_popB[top_indices[j], 2]
          }
          else if(dfpopB[i,top_indices[j]] == 0){
            ppl2value[i] = ppl2value[i]*(2-deltapopB[i,top_indices[j]])*(1-freq_popB[top_indices[j], 2])*(1-freq_popB[top_indices[j], 2])
          }
          else {
            ppl2value[i] = ppl2value[i]*(2-deltapopB[i,top_indices[j]])*freq_popB[top_indices[j], 2]*(1-freq_popB[top_indices[j], 2])
          }
        }
      }
      matrix2 <- matrix(0, nrow = 1, ncol = k)
      for(i in 1:k){
        if(ppl1value[i]>=ppl2value[i]){
          matrix2[i] = 1
        }else
          matrix2[i] = 2
      }
      
      if (all(matrix1 == 2)) {
        accuracy <- 0.5
      }
      else if (all(matrix2 == 1)) {
        accuracy <- 0.5
      }
      else{
        actual_label <- matrix(c(rep(1, k), rep(2, k)), nrow = 2, ncol = k, byrow = TRUE)
        combined_matrix <- rbind(matrix1, matrix2)
        predicted_classes <- combined_matrix
        confusion_matrix<- confusionMatrix(table(predicted_classes, actual_label))
        accuracy <- confusion_matrix$overall["Accuracy"]
      }
      l <- 0
      l <- (w-1)*10 + x
      result[l,4] <- accuracy
    }
    
    ####Now Random Markers
    for(y in 1:10){
      m <- y*10
      print("Random")
      print(paste("num of Markers:", m))
      sum = 0
      for(d in 1:dcontrol){
        set.seed(d)
        random_indices <- sample(1:nsites, m, replace = FALSE)
        ppl1value = numeric(k)
        for(i in 1:k){
          ppl1value[i] = 1
        }
        
        for(i in 1:k){
          for(j in 1:m){
            if(dfpopA[i,random_indices[j]] == 2){
              ppl1value[i] = ppl1value[i]*(2-deltapopA[i,random_indices[j]])*freq_popA[random_indices[j], 2]*freq_popA[random_indices[j], 2]
            }
            else if(dfpopA[i,random_indices[j]] == 0){
              ppl1value[i] = ppl1value[i]*(2-deltapopA[i,random_indices[j]])*(1-freq_popA[random_indices[j], 2])*(1-freq_popA[random_indices[j], 2])
            }
            else {
              ppl1value[i] = ppl1value[i]*(2-deltapopA[i,random_indices[j]])*freq_popA[random_indices[j], 2]*(1-freq_popA[random_indices[j], 2])
            }
          }
        }
        ppl2value = numeric(k)
        for(i in 1:k){
          ppl2value[i] = 1
        }
        for(i in 1:k){
          for(j in 1:m){
            if(dfpopA[i,random_indices[j]] == 2){
              ppl2value[i] = ppl2value[i]*(2-deltapopA[i,random_indices[j]])*freq_popB[random_indices[j], 2]*freq_popB[random_indices[j], 2]
            }
            else if(dfpopA[i,random_indices[j]] == 0){
              ppl2value[i] = ppl2value[i]*(2-deltapopA[i,random_indices[j]])*(1-freq_popB[random_indices[j], 2])*(1-freq_popB[random_indices[j], 2])
            }
            else {
              ppl2value[i] = ppl2value[i]*(2-deltapopA[i,random_indices[j]])*freq_popB[random_indices[j], 2]*(1-freq_popB[random_indices[j], 2])
            }
          }
        }
        
        matrix1 <- matrix(1, nrow = 1, ncol = k)
        for(i in 1:k){
          if(ppl1value[i] > ppl2value[i]){
            matrix1[i] = 1
          }else
            matrix1[i] = 2
        }
        ppl1value = numeric(k)
        for(i in 1:k){
          ppl1value[i] = 1
        }
        for(i in 1:k){
          for(j in 1:m){
            if(dfpopB[i,random_indices[j]] == 2){
              ppl1value[i] = ppl1value[i]*(2-deltapopB[i,random_indices[j]])*freq_popA[random_indices[j], 2]*freq_popA[random_indices[j], 2]
            }
            else if(dfpopB[i,random_indices[j]]== 0){
              ppl1value[i] = ppl1value[i]*(2-deltapopB[i,random_indices[j]])*(1-freq_popA[random_indices[j], 2])*(1-freq_popA[random_indices[j], 2])
            }
            else {
              ppl1value[i] = ppl1value[i]*(2-deltapopB[i,random_indices[j]])*freq_popA[random_indices[j], 2]*(1-freq_popA[random_indices[j], 2])
            }
          }
        }
        ppl2value = numeric(k)
        for(i in 1:k){
          ppl2value[i] = 1
        }
        for(i in 1:k){
          for(j in 1:m){
            if(dfpopB[i,random_indices[j]] == 2){
              ppl2value[i] = ppl2value[i]*(2-deltapopB[i,random_indices[j]])*freq_popB[random_indices[j], 2]*freq_popB[random_indices[j], 2]
            }
            else if(dfpopB[i,random_indices[j]]== 0){
              ppl2value[i] = ppl2value[i]*(2-deltapopB[i,random_indices[j]])*(1-freq_popB[random_indices[j], 2])*(1-freq_popB[random_indices[j], 2])
            }
            else {
              ppl2value[i] = ppl2value[i]*(2-deltapopB[i,random_indices[j]])*freq_popB[random_indices[j], 2]*(1-freq_popB[random_indices[j], 2])
            }
          }
        }
        matrix2 <- matrix(0, nrow = 1, ncol = k)
        for(i in 1:k){
          if(ppl1value[i]>ppl2value[i]){
            matrix2[i] = 1
          }else
            matrix2[i] = 2
        }
        
        #Make confusion matrix with accuracy
        combined_matrix <- rbind(matrix1, matrix2)
        predicted_classes <- combined_matrix
        confusion_matrix<- confusionMatrix(table(predicted_classes, actual_label))
        accuracy <- confusion_matrix$overall["Accuracy"]
        l <- (w-1)*10 + 50 + y
        result[l,3+d] <- accuracy
      }
    }
  }
  
  #write the results out 
  outfileClass = gsub("ldPruned.raw", "ldPruned_classification.out", fnames[query]) #outfile name
  write.table(result, file = outfileClass, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  print(paste0("done ", fnames[query]))
#}

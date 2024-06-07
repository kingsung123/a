run_hla_scan_command <- function(scan_dir) {
#############
  bam_file <- list.files(paste(scan_dir, "1.HLAscan_work_processing/2.bam_input/", sep = ""))
  bam_file <- strsplit(bam_file, "\\.")
  bam_file <- sapply(bam_file, function(x) ifelse(length(x) > 1, substring(x[[1]], 1), ""))
  bam_file <- unique(bam_file)
  bam_file_l<-length(bam_file)


  setwd(paste(scan_dir,sep = "")) 
   library("data.table")
   
   HLA_data <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1")
   
   split_data <- strsplit(HLA_data, "-")
   
   HLA_data_str <- sapply(split_data, function(x) ifelse(length(x) > 1, substring(x[[2]], 1), ""))
   HLA_data_str <- paste(HLA_data_str, "*", sep = "")
   HLA_data_str_rep <- rep(HLA_data_str, each = 2)
   
   HLA_data_rep <- rep(HLA_data, each = 2)
   
   HLA_data_ex <- c("_1", "_2")
   
   HLA_names <- paste(HLA_data_rep, HLA_data_ex, sep = "")

  
  ###############
  
  command<-c()
  
  for(rr in 1:bam_file_l)
  {
    command[[rr]] <- paste('for gene in'," ", paste(HLA_data, collapse = " "),'; do ',scan_dir,'2.HLAscan_program/2.hla_scan_program/hla_scan_r_v2.1.4 -b',paste(scan_dir, "1.HLAscan_work_processing/2.bam_input/", sep = ""), bam_file[[rr]], '.bam -v 38 -d', scan_dir,'5.HLAdb/2.hla_db/HLA-ALL.IMGT -g $gene -t 4; done >', paste(scan_dir, "1.HLAscan_work_processing/3.scan_file/",bam_file[[rr]],'.txt', sep = ""), sep = "")
    
    
    system(command[[rr]])
  }
  
  
  
}
run_hla_scan_combind <- function(scan_dir,scan_file123)
{
  #############
  bam_file <- list.files(paste(scan_dir, "1.HLAscan_work_processing/2.bam_input/", sep = ""))
  bam_file <- strsplit(bam_file, "\\.")
  bam_file <- sapply(bam_file, function(x) ifelse(length(x) > 1, substring(x[[1]], 1), ""))
  bam_file <- unique(bam_file)
  bam_file_l<-length(bam_file)


  setwd(paste(scan_dir,sep = "")) 
   library("data.table")
   
   HLA_data <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1")
   
   split_data <- strsplit(HLA_data, "-")
   
   HLA_data_str <- sapply(split_data, function(x) ifelse(length(x) > 1, substring(x[[2]], 1), ""))
   HLA_data_str <- paste(HLA_data_str, "*", sep = "")
   HLA_data_str_rep <- rep(HLA_data_str, each = 2)
   
   HLA_data_rep <- rep(HLA_data, each = 2)
   
   HLA_data_ex <- c("_1", "_2")
   
   HLA_names <- paste(HLA_data_rep, HLA_data_ex, sep = "")

  
  ###############
  
  
  output_lines <- readLines(paste(scan_dir, "1.HLAscan_work_processing/3.scan_file/",scan_file123,sep = ""))
  output_lines_l <- length(output_lines)
  output <- character()
  
  for (i in 1:output_lines_l) {
    line <- gsub("\t", "", output_lines[i])  
    output <- c(output, line)
  }
  
  hla_types_index <- which(output == "----------- HLA-Types -----------")
  hla_types_index_l <- length(hla_types_index)
  extracted_string1 <- c()
  extracted_string2 <- c()
  
  HLA_data_ex_l <- length(HLA_data_ex)
  extracted_blank <- c()
  extracted_blank1 <- c()
  hla_type1 <- c()
  hla_type2 <- c()
  original_list <- c()
  extracted_vector <- c()
  hla_types_index_l_num <- hla_types_index_l - 1
  
  for (j in 1:hla_types_index_l) {
    for (k in 1:HLA_data_ex_l) {
      hla_type1[[j]] <- output[hla_types_index[[j]] + 1]
      extracted_string1[[j]] <- substr(hla_type1[[j]], start = 9, stop = 13)
      
      hla_type2[[j]] <- output[hla_types_index[[j]] + 2]
      extracted_string2[[j]] <- substr(hla_type2[[j]], start = 9, stop = 13)
      
      extracted_blank[[j]] <- c(extracted_string1[[j]], extracted_string2[[j]])
      
 
        extracted_vector <-unlist(extracted_blank)
      
    }
  }
  
  empty_dt <- data.table(matrix(vector(), ncol = length(HLA_names)))
  setnames(empty_dt, HLA_names)
  
  data1 <- paste(HLA_data_str_rep, extracted_vector, sep = "")
  col_name <- character()
  param <- character()
  dd_balnak <- data.table()
  
  for (kk in seq_along(data1)) {
    col_name <- HLA_names[kk]
    param <- data1[[kk]]
    dd_balnak[, (col_name) := param]
  }
  dd_balnak1<-dd_balnak
  return(dd_balnak1)
  
}


run_hla_scan_all_bam <- function(scan_dir)
{  
bam_file <- list.files(paste(scan_dir, "1.HLAscan_work_processing/2.bam_input/", sep = ""))
bam_file <- strsplit(bam_file, "\\.")
bam_file <- sapply(bam_file, function(x) ifelse(length(x) > 1, substring(x[[1]], 1), ""))
bam_file <- unique(bam_file)
bam_file_l<-length(bam_file)


setwd(paste(scan_dir,sep = ""))
library("data.table")

HLA_data <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1")

split_data <- strsplit(HLA_data, "-")

HLA_data_str <- sapply(split_data, function(x) ifelse(length(x) > 1, substring(x[[2]], 1), ""))
HLA_data_str <- paste(HLA_data_str, "*", sep = "")
HLA_data_str_rep <- rep(HLA_data_str, each = 2)

HLA_data_rep <- rep(HLA_data, each = 2)

HLA_data_ex <- c("_1", "_2")

HLA_names <- paste(HLA_data_rep, HLA_data_ex, sep = "")

bam_file_l<-length(bam_file)
command<-c()

run_hla_scan_command(scan_dir)

scan_file <- list.files(paste(scan_dir, "1.HLAscan_work_processing/3.scan_file/", sep = ""))

scan_file_muti<-c()
dd_balnak2<-c()
for(uu in  1:bam_file_l)
{  
  scan_file_muti[[uu]]<-run_hla_scan_combind(scan_dir,scan_file[[uu]]) 
  
  dd_balnak2<-rbind(dd_balnak2, scan_file_muti[[uu]])
  
  
}
dd_balnak2<-cbind("sample_bam"=bam_file,dd_balnak2)



output_file <- paste(scan_dir, "3.HLAscan_result/combind_scan_", 'bam_file', ".csv", sep = "")
fwrite(dd_balnak2, output_file)
}



  

scan_dir <- "/mnt/d/google雲端資料備份/生物資訊分析程式碼庫/基因型分析/HLAscan/"

run_hla_scan_all_bam(scan_dir)

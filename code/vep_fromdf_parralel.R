vep_fromdf_parralel=function(input_df,chr,enst){
  # input_df=snps_sum[c(1:100),]
  # chr=9
  # enst="ENST00000275493"
  # enst="ENST00000318560"
  ##########This function parallelizes ENSEBML annotations, making them 10x faster
  snps_sum=input_df #please make all instances of snps_sum input df
  library(foreach)
  library(doParallel)
  library(tictoc)
  cores=detectCores()
  cl= makeCluster(cores[1]-1)
  registerDoParallel(cl)

  #I tried to iterate this function so that each core would get 10 rows of variants and not, but that didn't actually save much time.
  # snps_sum_ann=foreach(i=seq(1,nrow(a),10),combine=rbind) %dopar% {
  snps_sum_ann=foreach(i=seq(1,nrow(snps_sum[,1])),combine=rbind,.packages = c("dplyr",
                                                                               "httr",
                                                                               "xml2")) %dopar% {
source("code/variantcaller/vep_fromdf.R") #I don't know how to make it so that it doesn't have to initialize the function for each row
source("code/variantcaller/vep_fromquery.R")
    # seq(1,10,2)
     #for some reason, each of these packages needs to be re-initialized upon
    # library(jsonlite)
    # library(httr)
    # library(xml2)
    # snps_sum_ann=vep_fromdf(a[c(i:(i+9)),])
    snps_sum_ann=vep_fromdf(snps_sum[i,],chr=chr,enst=enst)
    snps_sum_ann
  }
  stopCluster(cl)

  snps_sum_ann=as.data.frame(do.call(rbind,snps_sum_ann))
  snps_sum_ann
}

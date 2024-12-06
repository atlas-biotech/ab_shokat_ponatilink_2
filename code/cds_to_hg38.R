cds_to_hg38=function(input_df){
  ####To convert coordinates from ref CDS e.g. (796-797A/G) to 1300000-130001A/G#####
  # Inputs: input_df dataframe containing alt_start_pos, alt_end_pos
  # Outputs: input_df with alt_start_hg38 and alt_end_hg38
  #Dependencies: ref_genomic_coordinates dataframe
#   #The function below converts cDNA start and end positions to HG38 positions for that transcript. Why do you need a function? Because at genomic splice sites for multi-nucleotide variants, you can't simply take the alt position and add the number of nucelotides in the MNV. Because that way, it might seem like an MNV starts at an exon and ends in a splice-site, wheras in reality, the MNV started in one exon and ended at the other. Ensembl's rest API does not
# like it when this happens and throws an error. This is a very fringe-case exception, I know, so this cds_to_hg38 function can probably be deprecated at some point.
  if(sum(as.numeric(grep("alt_end_pos",colnames(input_df))))==0){
  input_df$alt_end_pos=input_df$alt_start_pos+1
  #^basically saying that if you don't find a column called alt_end_hg38, make it. And initialize end positions as start positions + 1
}
  input_df$alt_start_hg38=NA
  input_df$alt_end_hg38=NA
  for(i in 1:nrow(input_df)){
    ref_i=ref_genomic_coordinates
    ref_i$diff_start=input_df$alt_start_pos[i]-ref_i$start
    ref_i$diff_end=input_df$alt_end_pos[i]-ref_i$start
    ref_i=ref_i%>%filter(diff_start>=0)%>%arrange(diff_start)
    ref_i=ref_i[1,]
#     #The above lines are basically doing, for a given start position (aka unique row in the input dataframe, find the codon that is closest to that position. That is the codon that your reference is pointing to (this assumes that
# the multinucleotide variants are in-frame MNVs, which is the case for the MNVs that we work with)
    input_df$alt_start_hg38[i]=ref_i$chr_start+ref_i$diff_start
    input_df$alt_end_hg38[i]=ref_i$chr_start+ref_i$diff_end
  }
  return(input_df)
}

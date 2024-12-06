vep_fromdf=function(input_df,chr,enst){
  # input_df=mnv_sum
  # chr=7
  # enst="ENST00000275493"
  #Version: August 2022
  #Inputs: Dataframe with these fields: alt_start_hg38, alt
  #Optional field in dataframe: alt_end_hg38 (for mnvs because their end position is different than their starting position)
  #This function adds ensembl variant effect predictor annotations to dataframes.
  #This needs the vep_fromquery function
  #Requires the dplyr package
  input_df=input_df%>%filter(!ref==alt)

  input_df=input_df%>%mutate(protein_start=NaN,protein_end=NaN,amino_acids="NaN",codons="NaN",impact="NaN",polyphen_prediction="NaN",consequence_terms="NaN")
  # input_df=snps_sum
  # input_df=mnv_sum
  # i=1
  ###Can probably get rid of this clause because it appears in the cds_to_hg38 funciton.
  if(sum(as.numeric(grep("alt_end_hg38",colnames(input_df))))==0){
    input_df$alt_end_hg38=input_df$alt_start_hg38+1
    #^basically saying that if you don't find a column called alt_end_hg38, make it
  }
  # i=4
  # for(i in 276:nrow(input_df)){
  for(i in 1:nrow(input_df)){
    # i=1
    transcript_consequences=vep_fromquery(paste(input_df$alt_start_hg38[i],":",input_df$alt_end_hg38[i]-1,"/",input_df$alt[i],sep = ""),chr=chr,enst=enst)

    # b=vep_fromquery("130862936:130862938/TAA")
    # vep_fromquery("130874932:130874932/A")
    # transcript_consequences=vep_fromquery(paste(130884390,":",130884390,"/",input_df$alt[i],sep = ""))
    # transcript_consequences=vep_fromquery(paste(130862864,":",130862864,"/","C",sep = ""))
    # transcript_consequences$protein_start

    if("consequence_terms"%in%names(transcript_consequences)&&transcript_consequences$consequence_terms=="bad_request"){
      input_df$consequence_terms[i]="bad_request"
    } else {
      if(grepl("intron|stop",as.character(transcript_consequences$consequence_terms))%in%TRUE){
        input_df$consequence_terms[i]=transcript_consequences$consequence_terms[[1]][1]
      } else {
        if("protein_start"%in%names(transcript_consequences)){
          input_df$protein_start[i]=transcript_consequences$protein_start[[1]][1]
        }
        if("protein_end"%in%names(transcript_consequences)){
          input_df$protein_end[i]=transcript_consequences$protein_end[[1]][1]
        }
        if("amino_acids"%in%names(transcript_consequences)){
          input_df$amino_acids[i]=transcript_consequences$amino_acids[[1]][1]
        }
        if("codons"%in%names(transcript_consequences)){
          input_df$codons[i]=transcript_consequences$codons[[1]][1]
        }
        if("impact"%in%names(transcript_consequences)){
          input_df$impact[i]=transcript_consequences$impact[[1]][1]
        }
        if("polyphen_score"%in%names(transcript_consequences)){
          input_df$polyphen_prediction[i]=transcript_consequences$polyphen_prediction[[1]][1]
        }

        # if(grepl("synonymous|frameshift",as.character(transcript_consequences$consequence_terms))%in%TRUE){
        #   input_df$polyphen_prediction[i]=NA
        # }
        # if(sum(as.numeric(grepl("polyphen_score",names(transcript_consequences))))==0){
        #   input_df$polyphen_prediction[i]=NA
        # }
        # else{
        #   input_df$polyphen_prediction[i]=transcript_consequences$polyphen_prediction[[1]][1]
        # }

        input_df$consequence_terms[i]=transcript_consequences$consequence_terms[[1]][1]
        # transcript_consequences$consequence_terms[[1]][1]
      }
    }
  }
  input_df
  # a=snps_ann%>%filter(alt_start_hg38%in%"130874932")
}

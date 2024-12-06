vep_fromquery=function(region,chr,enst){
  ####Functions to use Ensembl's VEP to find protein consequence of mutations.
  #Version: August 2022
  # The ensembl server where you can download multiple things from
  # https://rest.ensembl.org/
  # https://rest.ensembl.org/vep/human/region/9:130884090:130884091/A
  library(httr)
  library(jsonlite)
  library(xml2)
  # region="130862976:130862976/A"
  # region="130854096:130854096/A" #splice variant
  # region="130854240:130854240/G" #throws bad http error
  # region="130854096:130854096/A" #first missense in dataset
  # region="130872215:130872215/T"
  # region="130863033:130863036/TAC"
  # region="130872213..130872215/CAC"
  # region="130872213..130872861/CAC"
  # a=vep_fromquery(region)

  #This function fetches ensembl vep for queries that look like this
  #Needs the packages httr, jsonlite, xml2
  server <- "https://rest.ensembl.org"
  ext=paste("/vep/human/region/",chr,":",gsub("\"","",region),"?",sep = "")
  # ext <- "/vep/human/region/9:130862976:130862976/A?"

  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  #A 503 error occurs when the server is busy. In that case, just redo the request till you don't get the 503 error
  while(r$status_code%in%503){
    r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  }
  #If we give ensembl a bad request, returns an error. These exceptions will be noted down for the respective SNPs/MNVs. Ideally, our code won't be sending ensembl bad requests, but this should deal with it for now. Note that most of these "bad requests" are when the alternate allele is the same as the variant allele. This happens often times at splice junctions, especially for MNVs, where the variant starts at one exon and ends at another.
  if(r$status_code==400){
    transcript_consequences=data.frame(consequence_terms="bad_request")
  }
  else{
    stop_for_status(r)
    rest_results=(fromJSON(toJSON(content(r))))
    transcript_consequences=(rest_results$transcript_consequences)[[1]]
    transcript_consequences=transcript_consequences%>%filter(transcript_id%in%enst)
    transcript_consequences
  }
}

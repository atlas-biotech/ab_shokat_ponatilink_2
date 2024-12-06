transcript_verifier=function(pos,enst){
  #Transcript verifying function:
  #The following is code to verify that the user's CDS to Hg38 mapping coordinates are the same as the coordinates predicted by Ensemble's rest API:
  #Source: https://rest.ensembl.org/documentation/info/assembly_cds
  #pos is the position on the coding sequence
  #enst is the ensembl transcript ID such as "ENST00000318560"
  #output is the hg38 genomic coordinate of the coding sequence, according to ensembl
  server <- "https://rest.ensembl.org"
  # ext <- "/map/cds/ENST00000318560/907..909?"
  # ext <- "/map/cds/ENST00000275493/1..20?"
  ext=paste("/map/cds/",enst,"/",pos,"..",pos,"?",sep = "")

  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  #sometimes the server is busy and it throws a busy server error code 503. I'm using the while loop to keep trying the request until the 503 request changes
  count=1
  while(r$status_code%in%503){
    r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
    count=count+1
  }
  stop_for_status(r)

  rest_results=(fromJSON(toJSON(content(r))))
  tcrpt=rest_results$mappings
  return(list(tcrpt$start[[1]][1],count))
}

setOldClass("igraph")

STRINGdb <- setRefClass("STRINGdb",
    fields = list(
      annotations="data.frame",
      annotations_description="data.frame",
      graph="igraph",
      homology_graph="igraph",
      proteins="data.frame",
      aliases_tf="data.frame",
      aliases_type="character",
      #spe    ciesList="data.frame",
      species="numeric",
      version="character",    
      input_directory="character",
      backgroundV = "vector",
      score_threshold = "numeric",
      pathways_benchmark_blackList="data.frame",
      stable_url = "character"
    ),

    methods = list(
      
      #########################################
      ## add_diff_exp_color
      #########################################
      
      add_diff_exp_color = function(screen, logFcColStr="logFC" ){

'
Description:
  Take in input a dataframe containing a logFC column that reports the logarithm of the difference in expression level.
  Add a "color" column to the data frame such that strongly downregulated genes are colored in green and strong upregulated genes are in red.
  When the down or up-regulation is instead weak the intensity of the color gets weaker as well, accordingly.

Author(s):
   Andrea Franceschini
'

          screen_pval05_pos = subset(screen, as.matrix(screen[logFcColStr])[, 1] > 0)
          k=exp(screen_pval05_pos[,logFcColStr])/exp(max(screen_pval05_pos[,logFcColStr]))
          screen_pval05_pos_col = data.frame(screen_pval05_pos, color = color.scale(k,c(1,1),c(1,0),c(1,0)) )
          screen_pval05_neg = subset(screen, as.matrix(screen[logFcColStr])[, 1] <= 0)
          k=exp(-screen_pval05_neg[,logFcColStr])/exp(-min(screen_pval05_neg[,logFcColStr]))
          screen_pval05_neg_col = data.frame(screen_pval05_neg, color = color.scale(k,c(1,0),c(1,1),c(1,0)) )
          return(rbind(screen_pval05_pos_col, screen_pval05_neg_col))

      },
      
      
      #########################################
      ## add_proteins_description
      #########################################
      
      
      add_proteins_description = function(screen){

'
Description:
  Add description coluns to the proteins that are present in the data frame given in input. 
  The data frame must contain a column named "STRING_id".

Input parameters:
      "screen"  a data frame having a "STRING_id" column 

Author(s):
   Andrea Franceschini
'
        
          proteinsDf2 = get_proteins()
          proteinsDf3 = merge(screen, proteinsDf2, by.x="STRING_id", by.y="protein_external_id", all.x=TRUE, sort=FALSE)

          return(proteinsDf3)

      },
      
      #########################################
      ## initialize
      #########################################
      
      initialize = function(...) {

'
Description:
    Initiliazes STRING with a given version and species.
    Checks the current version against the API and sets the "take_first" parameter. 
    
Author(s):
   Andrea Franceschini
'
 
          ## TODO: allow different versions (v11 will be the lowest one available) 

          callSuper(...)

          if(length(species)==0) {

            cat("WARNING: You didn't specify a species. Hence we will set 9606 (Homo Sapiens) as your species.\n")
            species <<- 9606                        

          }

          if(length(score_threshold)==0) {
              cat("WARNING: Score threshold is not specified. We will be using medium stringency cut-off of 400.\n")
              score_threshold <<- 400
          }
          
          

          server_version = read.table(url("https://string-db.org/api/tsv-no-header/version"), colClasses = "character")$V1
          # stable_url = read.table(url("https://string-db.org/api/tsv-no-header/version"), colClasses = "character")$V2
          
          valid_versions = server_version ## for now we can only use the current version
          
          if(length(version)==0) {
            
            cat("WARNING: You didn't specify a version of the STRING database to use. Hence we will use STRING", server_version,"\n")
            version <<- server_version
            
          } else {

            version_dotted = paste(version, ".0", sep="")

            if(! (version %in% valid_versions || version_dotted %in% valid_versions)) {
              
              cat("ERROR: Currently STRINGdb only supports the most recent version of STRING: ", server_version)
              stop()
              
            } 
          }
          if (version_dotted %in% valid_versions) version <<- version_dotted


          version_hyphenated = gsub("\\.", "-", version)
          stable_url <<- paste("https://version-", version_hyphenated, ".string-db.org", sep="")          

          if(input_directory=="" || is.null(input_directory) || length(input_directory)==0) input_directory<<-tempdir()
          if(input_directory=="" || is.null(input_directory) || length(score_threshold)==0 || score_threshold<1) score_threshold <<- 1

          aliases_type <<- "take_first"

      },


      #########################################
      ## load_all
      #########################################
 

      load_all = function(){

'
Description:
  Force download and loading of all the files (so that you can later store the object on the hard disk if you like)
  It makes use of the variables:
  "backgroundV"       vector containing STRING identifiers to be used as background 
                      (i.e. the STRING network loaded will contain only the proteins that are present also in this vector)
  "score_threshold"   STRING combined score threshold (the network loaded contains only interactions having a combined score greater than this threshold)

Author(s):
   Andrea Franceschini
'
       
          ## TODO: DS: We have removed get_annotations.
          ## We leave the function in the code 
          ## in case at some point we dump these in the downloads
 
          #x1 = get_annotations()
          #x2 = get_annotations_desc()  
          x3 = get_proteins()
          #x4 = get_species()


          ## TODO: DS: Why we are running it like that twice?

          get_aliases(takeFirst=FALSE)
           get_aliases(takeFirst=TRUE)
           load()
       },
       
 
 
       #########################################
       ## get_proteins
       #########################################
  
 
       get_proteins = function(){

'
Description:
  Returns the STRING proteins data frame.
  (it downloads and caches the information the first time that is called).

Author(s):
   Andrea Franceschini
'
        
        if(nrow(proteins)==0){
          temp = downloadAbsentFile(paste("https://stringdb-static.org/download/protein.info.v", version, "/", species, ".protein.info.v", version, ".txt.gz", sep=""), oD=input_directory)
          proteinsDf <- read.table(temp, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE, quote="")
          proteinsDf2 = subset(proteinsDf, select=c("protein_external_id",  "preferred_name", "protein_size", "annotation"))
          proteins <<- proteinsDf2
        }
        return(proteins)
      },
 

      #########################################
      ## get_aliases
      #########################################
      
      get_aliases = function(takeFirst=TRUE){

'
Description:
  Loads and returns STRING alias table. 

Author(s):
   Andrea Franceschini
'        

          if(takeFirst &  aliases_type == "take_first" & nrow(aliases_tf)>0) return(aliases_tf)
          if(!takeFirst &  aliases_type == "all" & nrow(aliases_tf)>0) return(aliases_tf)

          
          ## TODO: DS: Test take first
          ## or better: implement it nicer
          
          temp = downloadAbsentFile(paste("https://stringdb-static.org/download/protein.aliases.v", version, "/", species, ".protein.aliases.v", version, ".txt.gz", sep=""), oD=input_directory)        
          if(!takeFirst){ 
            aliases_type <<- "all"
          } else {
            aliases_type <<- "take_first"
          }
          
          
          if(takeFirst &  aliases_type == "take_first" & nrow(aliases_tf)>0) return(aliases_tf)
          if(!takeFirst &  aliases_type == "all" & nrow(aliases_tf)>0) return(aliases_tf)
          
          proteins <<- get_proteins()
          
          aliasDf <- read.table(temp, skip=1, sep = "\t", header=FALSE, quote="", stringsAsFactors=FALSE, fill = TRUE)
          colnames(aliasDf) <- c("STRING_id", "alias", "sources")
          aliasDf = subset(aliasDf, select=c("STRING_id", "alias"))
          
          pr1=data.frame(STRING_id=proteins$protein_external_id, alias=proteins$preferred_name, stringsAsFactors=FALSE)        
          pr2=data.frame(STRING_id=proteins$protein_external_id, alias=proteins$protein_external_id, stringsAsFactors=FALSE)
          pr3=data.frame(STRING_id=proteins$protein_external_id, alias=unlist(strsplit(proteins$protein_external_id, "\\."))[seq(from=2, to=2*nrow(proteins), by=2)], stringsAsFactors=FALSE)
          
          #if(takeFirst){aliasDf = subset(aliasDf, !(alias %in% proteins$preferred_name) & !(alias %in% proteins$protein_external_id) )  }
          if(takeFirst){aliasDf = subset(aliasDf, !(toupper(iconv(alias, "WINDOWS-1252", "UTF-8")) %in% toupper(proteins$preferred_name)) & 
                                           !(toupper(iconv(alias, "WINDOWS-1252", "UTF-8")) %in% toupper(proteins$protein_external_id))  &
                                            !(toupper(iconv(alias, "WINDOWS-1252", "UTF-8")) %in% toupper(unlist(strsplit(proteins$protein_external_id, "\\."))[seq(from=2, to=2*nrow(proteins), by=2)])) )
          }
          
          aliasDf2=rbind(pr1,pr2,pr3, aliasDf)
          aliases_tf <<- aliasDf2

          return(aliasDf2)
      },
 

      #########################################
      ## get_bioc_graph
      #########################################
      
      get_bioc_graph = function(){

'
Description:
  Returns the interaction graph as an object of the graph package in Bioconductor

Author(s):
   Andrea Franceschini
'

        print(".... please wait about 20 minutes to create the graph ....")
        print(".... it may take up to 12GB for full network (initiate the package higher score_threshold to filter the interactions) ....")

        if("graph" %in% installed.packages()){

          library(graph)
          if(is.null(graph)) load()
         
          url <- paste("https://stringdb-static.org/download/protein.links.v", version, "/", species, ".protein.links.v", version, ".txt.gz", sep="")
          temp = downloadAbsentFile(url, oD=input_directory)
         
          PPI <- read.table(temp, sep = " ", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
          
          PPIselected = PPI

          if(length(score_threshold)!=0) PPIselected <- PPI[PPI$combined_score >= score_threshold,]

          if(!is.null(backgroundV)){
            PPIselected <- PPIselected[PPIselected$protein1 %in% backgroundV,]
            PPIselected <- PPIselected[PPIselected$protein2 %in% backgroundV,]
          }
           
          nel_graph <- ftM2graphNEL(cbind(PPIselected$protein1, PPIselected$protein2), W=PPIselected$combined_score)
          nel_graph2 <- ugraph(nel_graph)

          return(nel_graph2)

          #return(igraph.to.graphNEL(graph))

        } else {

          cat("ERROR: In order to run this function you must install the \"graph\" package from CRAN. \n Please install the package and run this function again.\n")

        }
      },
 
      #########################################
      ## get_clusters
      #########################################
 
      
      get_clusters = function(string_ids, algorithm="fastgreedy"){

'
Description:
  Returns a list of clusters of interacting proteins.

References:
  See the iGraph (http://igraph.sourceforge.net/) documentation for additional information on the algorithms.
  Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. 
  http://igraph.sf.net

Input parameters:
    "string_ids"    a vector of STRING identifiers.
    "algorithm"     algorithm to use for the clustering (fastgreedy, walktrap,  edge.betweenness)

Author(s):
   Andrea Franceschini
'
        
        if(is.null(graph)) load()

        string_ids = unique(string_ids)
        string_ids = string_ids[!is.na(string_ids)]
        
        if(algorithm=="fastgreedy") fgreedy<-fastgreedy.community(simplify(get_subnetwork(string_ids)),merges=TRUE, modularity=TRUE)
        if(algorithm=="walktrap") fgreedy<-walktrap.community(get_subnetwork(string_ids),merges=TRUE, modularity=TRUE)
        if(algorithm=="edge.betweenness") fgreedy<-edge.betweenness.community(get_subnetwork(string_ids),merges=TRUE, modularity=TRUE)

        memb = membership(fgreedy)
        clusters = NULL
        for(i in 1:max(memb)){
          clusters[[i]] = names(membership(fgreedy)[membership(fgreedy)==i])
        }
        return(clusters)
      },
      
      #########################################
      ## get_enrichment
      #########################################
     

      get_enrichment = function(string_ids, category='All', methodMT=NULL, iea=NULL, minScore=NULL){

'
Description:
  Returns the enrichment for various annotation categories (among others Gene Ontology, KEGG, domains, Reference Publications) for a vector of STRING proteins.

Input parameters:
  "string_ids"    a vector of STRING identifiers.
  "category" enrichment term category (All, Process, Component, Function, Keyword, KEGG, RCTM, Pfam, SMART, InterPro), default:ALL 

Author(s):
   Andrea Franceschini
'

         if (!(is.null(methodMT))) {
             warning("methodMT parameter is depecated. Only FDR correction is available.")
         }

         if (!(is.null(iea))) {
             warning("iea parameter is deprecated.")
         }

         if (!(is.null(minScore))) {
             warning("minScore parameter is deprecated.")
         }


         urlStr <- paste(stable_url, "/api/tsv/enrichment?", sep="") 


         string_ids = unique(string_ids)
         string_ids = string_ids[!is.na(string_ids)]

         identifiers = ""

         for(id in string_ids) {
           identifiers = paste(identifiers, id, "%0d", sep="")
         }
 
         if(length(backgroundV)==0) {

             params = list(species=species, identifiers=identifiers)

         } else {

             background = ""
             for(id in backgroundV) {
               background = paste(background, id, "%0d", sep="")
             }
 
             params = list(species=species, identifiers=identifiers, background_string_identifiers=background)
         }
 
         tempDfv =postFormSmart(urlStr, .params=params)

         ann = read.table(text=tempDfv, sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE, header=TRUE)

         # ann = renameColDf(ann, "V1", "category")
         # ann = renameColDf(ann, "V2", "term_id")
         # ann = renameColDf(ann, "V3", "number_of_genes")
         # ann = renameColDf(ann, "V4", "number_of_genes_in_background")
         # ann = renameColDf(ann, "V5", "species")
         # ann = renameColDf(ann, "V6", "STRING_ids")
         # ann = renameColDf(ann, "V7", "preferred_names")
         # ann = renameColDf(ann, "V8", "pvalue")
         # ann = renameColDf(ann, "V9", "fdr")
         # ann = renameColDf(ann, "V10", "description")

         requestedCategory = category

         if ( (requestedCategory !=  'All')  && (!is.null(requestedCategory))) {

             ann = subset(ann, category==requestedCategory)

         }

         return(ann)

     },

      #########################################
      ## plot_network
      #########################################
 
      ## DS: We have removed get_link functionality

      plot_network = function(string_ids, payload_id=NULL, required_score=NULL, add_link=TRUE, network_flavor='evidence', add_summary=TRUE) {

'
Description:
  Plots an image of the STRING network with the given proteins.

Input parameters:
  "string_ids"        a vector of STRING identifiers
  "payload_id"        an identifier of payload data on the STRING server (see method post_payload for additional informations)
  "required_score"   a threshold on the score that overrides the default score_threshold, that we use only for the picture
                      As default this option is active but we suggest to deactivate it in case one is generating many images (e.g. in a loop). 
                      Deactivating this option avoids to generate and store a lot of short-urls on our server.
  "add_summary"       parameter to specify whether you want to add a summary text to the picture. This summary includes a p-value and the number of proteins/interactions.

Author(s):
   Andrea Franceschini
'
        
        if(is.null(required_score) ) required_score = score_threshold
        img = get_png(string_ids, payload_id=payload_id, required_score=required_score, network_flavor=network_flavor)
        if(!is.null(img)){
          plot(1:(dim(img)[2]), type='n', xaxt='n', yaxt='n', xlab="", ylab="", ylim=c(1,dim(img)[1]), xlim=c(1,(dim(img)[2])), asp = 1 )
          if(add_summary) mtext(get_summary(string_ids, required_score), cex = 0.7)
          #if(add_link) mtext(get_link(string_ids, payload_id=payload_id, required_score=required_score), cex = 0.7, side=1)
          rasterImage(img, 1, 1, (dim(img)[2]), dim(img)[1])
        } 
      },


      #########################################
      ## get_summary
      #########################################


      get_summary = function(string_ids, required_score){

'
Description:
  Returns a summary of the STRING sub-network containing the identifiers provided in input.

Input parameters:
  "string_ids"      a vector of STRING identifiers.

Author(s):
   Andrea Franceschini
'
        string_ids = unique(string_ids)
        string_ids = string_ids[!is.na(string_ids)]

        enrichment = ppi_enrichment(string_ids, required_score)
        summaryText = paste("proteins: ", length(string_ids), '\n', 
                            "interactions: ", enrichment$edges, '\n',
                            "expected interactions: ", enrichment$lambda,' (',
                            "p-value: ",enrichment$enrichment, ')\n', sep=""
        )
        return(summaryText)
      },


      #########################################
      ## ppi_enrichment
      #########################################

      ppi_enrichment = function(string_ids, required_score=NULL){

'
Description:
  Queries STRING server for protein-protein interaction enrichment

Input parameters:
  "string_ids"      a vector of STRING identifiers.

Author(s):
   Damian Szklarczyk
'

          if (is.null(required_score)) required_score = score_threshold

          urlStr <- paste(stable_url, "/api/tsv-no-header/ppi_enrichment", sep="")

          string_ids = unique(string_ids)
          string_ids = string_ids[!is.na(string_ids)]

          identifiers = ""

          for(id in string_ids) {
            identifiers = paste(identifiers, id, sep="%0d")
          }


         if(length(backgroundV)==0) {

             params = list(species=species, identifiers=identifiers, required_score=required_score)

         } else {

             background = ""
             for(id in backgroundV) {
               background = paste(background, id, "%0d", sep="")
             }
             params = list(species=species, identifiers=identifiers, background_string_identifiers=background, required_score=required_score)
         }
 

          tempDfv=postFormSmart(urlStr, .params=params)

          answer = read.table(text=tempDfv, sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE, header=FALSE)

          result = list(enrichment = answer$V6, edges=answer$V2, lambda=answer$V5)
  
          return(result)
      },
 

 
      #########################################
      ## get_png
      #########################################
 

      get_png = function(string_ids, required_score=NULL, network_flavor="evidence", file=NULL, payload_id=NULL){

'
Description:
  Returns a png image of a STRING protein network with the given identifiers.

Input parameters:
  "string_ids"        a vector of STRING identifiers.
  "required_score"    minimum STRING combined score of the interactions 
                        (if left NULL we get the combined score of the object, which is 400 by default)
  "network_flavor"    specify the flavor of the network ("evidence", "confidence" or "actions".  default "evidence")
  "file"              file where to save the image (must have .png extension)

Author(s):
   Andrea Franceschini
'
        if(length(string_ids) > 2000) {
          cat("ERROR: We do not support lists with more than 2000 genes.\nPlease reduce the size of your input and rerun the analysis. \t")
          stop()
        }

        if(is.null(required_score) ) required_score = score_threshold

        string_ids = unique(string_ids)
        string_ids = string_ids[!is.na(string_ids)]

        #urlStr = paste(stable_url, "/api/image/network", sep="")
        urlStr = paste("http://version", version, ".string-db.org/api/image/network", sep="" )

        identifiers=""
        for(id in string_ids ){ identifiers = paste(identifiers, id, sep="%0d")}
        params = list(required_score=required_score, required_score=required_score, network_flavor=network_flavor, identifiers=identifiers, species=species, caller_identity='STRINGdb-package')

        if(!is.null(payload_id)) params["internal_payload_id"]= payload_id
        img <- readPNG(postFormSmart(  urlStr, .params=params) )
        if(!is.null(file))  writePNG(img,  file)
        
        return(img)
      },
 

      #########################################
      ## get_graph
      #########################################
      
      get_graph = function(){
'
Description:
  Return an igraph object with the entire STRING network. 
  We invite the user to use the functions of the iGraph package to conveniently 
  search/analyze the network.
  
References:
  Csardi G, Nepusz T: The igraph software package for complex network research, 
  InterJournal, Complex Systems 1695. 2006. 
  http://igraph.sf.net

See Also:
  In order to simplify the most common tasks, we do also provide convenient functions 
  that wrap some iGraph functions.
  get_interactions(string_ids)   # returns the interactions in between the input proteins
  get_neighbors(string_ids)      # Get the neighborhoods of a protein (or of a vector of proteins).
  get_subnetwork(string_ids)     # returns a subgraph from the given input proteins
  
Author(s):
   Andrea Franceschini
'
        if(is.null(graph)) load()
        return(graph)
      },
      
      #########################################
      ## post_payload
      #########################################
            
      post_payload = function(stringIds, colors=NULL, comments=NULL, links=NULL, iframe_urls=NULL, logo_imgF=NULL, legend_imgF=NULL ){

'
Description:
  Posts the input to STRING and returns an identifier that you can use to access the payload when you enter in our website.

Input parameters:
  "string_ids"        vector of STRING identifiers.
  "colors"            vector containing the colors to use for a every STRING identifier ( the order of the elements must match those in the string_ids vector)
  "comments"          vector containing the comments to use for every STRING identifier ( the order of the elements must match those in the string_ids vector)
  "links"             vector containing the links to use for every STRING identifier ( the order of the elements must match those in the string_ids vector)
  "iframe_urls"       vector containing the urls of the iframes to use for every STRING identifier 
                              ( the order of the elements must match those in the string_ids vector)
  "logo_imgF"         path to a file containing the logo image to be display in the STRING website
  "legend_imgF"       path to a file containing a legend image to be display in the STRING website

Author(s):
   Andrea Franceschini

'

          postFormParams = list(identifiers=paste(stringIds, collapse=" ") )
        
          if(!is.null(colors)) postFormParams = c(postFormParams, list(colors=paste(colors, collapse=" ")))
          if(!is.null(comments)) postFormParams = c(postFormParams, list(comments=paste(comments, collapse="____")))
          if(!is.null(links)) postFormParams = c(postFormParams, list(links=paste(links, collapse=" ")))
          if(!is.null(iframe_urls)) postFormParams = c(postFormParams, list(iframe_urls=paste(iframe_urls, collapse=" ")))
          if(!is.null(logo_imgF)) postFormParams = c(postFormParams, list(logo_img=fileUpload(logo_imgF)))
          if(!is.null(legend_imgF)) postFormParams = c(postFormParams, list(legend_img=fileUpload(legend_imgF)))
        
          postRs = postFormSmart(paste(stable_url, "/cgi/webservices/post_payload.pl", sep=""),  .params = postFormParams)
          
          return(postRs)
    },


     #########################################
     ## get_homology_graph
     #########################################
 

     get_homology_graph = function(min_homology_bitscore=60){

         temp = downloadAbsentFile(paste("https://stringdb-static.org/download/protein.homology.v", version, "/", species, ".protein.homology.v", version, ".txt.gz", sep=""), oD=input_directory)

         PPI <- read.table(temp, sep = "\t", header=FALSE, skip=1, stringsAsFactors=FALSE, fill = TRUE)
          
         PPIselected = PPI

         if(!is.null(min_homology_bitscore)) PPIselected <- PPI[PPI$V9 >= min_homology_bitscore,]

         if(!is.null(backgroundV)){
           PPIselected <- PPIselected[PPIselected$V1 %in% backgroundV,]
           PPIselected <- PPIselected[PPIselected$V2 %in% backgroundV,]
         }
         
         myg = graph.data.frame(PPIselected,FALSE)
         homology_graph <<- myg
         return(myg)

      },
 

     #########################################
     ## get_paralogs
     #########################################
      
      
     get_paralogs = function(string_ids) {

'
Description:
  Returns the list of paralogs of the given input in their species. 

Input parameters:
    "string_ids"          a vector of STRING identifiers.

Author(s):
   Andrea Franceschini

'

        string_ids = unique(string_ids)
        string_ids = string_ids[!is.na(string_ids)]

        
        if(length(string_ids) > 2000) {
          cat("ERROR: We support a maximum of 2000 STRING identifiers per call. Please reduce the size of the input and try again. \t")
          stop()
        }

        urlStr = paste(stable_url, "/api/tsv-no-header/homology", sep="")
        identifiers=""
        for(id in string_ids ){
          identifiers = paste(identifiers, id, "%0D", sep="")
        }

        params = list(identifiers=identifiers, species=species)

        tempDfv=postFormSmart(urlStr, .params=params)
        hhits <- read.table(text=tempDfv, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
        return(hhits)
      },
 


     #########################################
     ## get_homologs_besthits
     #########################################
      
      
     get_homologs_besthits = function(string_ids, target_species_id=NULL) {

'
Description:
  Returns the list of closest homologs (as measured by bitscore) of the given input identifiers in all STRING species or single target species.

Input parameters:
    "string_ids"          a vector of STRING identifiers.
    "target_species_id"   NCBI taxonomy identifier of the species to query for homologs (the species must be present in the STRING database); All spcies by default. 

Author(s):
   Andrea Franceschini

'

        string_ids = unique(string_ids)
        string_ids = string_ids[!is.na(string_ids)]

        
        if(length(string_ids) > 2000) {
          cat("ERROR: We support a maximum of 2000 STRING identifiers per call. Please reduce the size of the input and try again. \t")
          stop()
        }

        urlStr = paste(stable_url, "/api/tsv-no-header/homology_best", sep="")
        identifiers=""
        for(id in string_ids ){
          identifiers = paste(identifiers, id, "%0D", sep="")
        }

        if (is.null(target_species_id)) {
            params = list(identifiers=identifiers, species=species)
        } else { 
            params = list(species_b=target_species_id, identifiers=identifiers, species=species)
        }

        tempDfv=postFormSmart(urlStr, .params=params)
        hhits <- read.table(text=tempDfv, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
        return(hhits)
      },
      

      #########################################
      ## get_interactions
      #########################################

      get_interactions = function(string_ids){
'
Description:
  Shows the interactions in between the proteins that are given in input.

Input parameters:
    "string_ids"    a vector of STRING identifiers.

Author(s):
   Andrea Franceschini

'
        if(is.null(graph)) load()
        return(get.data.frame(get_subnetwork(string_ids),  c("edges")))
      },
      
      
      #########################################
      ## get_neighbors
      #########################################
    
      
      get_neighbors = function(string_ids){
'
Description:
Get the neighborhoods of a protein (or of a vector of proteins) that is given in input.

Input parameters:
  "string_ids" =  a vector of STRING identifiers.

Author(s):
   Andrea Franceschini
'
        if(is.null(graph)) load()

        vp <- intersect(string_ids, V(graph)$name) 

        ns <- c()

        for (vertex in vp) {
            ns <- append(ns, unique(V(graph)[neighbors(graph, vertex)]$name))
            
        }
        return(ns)
      },
      
      
      #########################################
      ## get_ppi_enrichment
      #########################################
     
      
      get_ppi_enrichment = function(string_ids){

'
Description:
  Returns a pvalue representing the enrichment in interactions of the list of proteins 
          (i.e. the probability to obtain such a number of interactions by chance).

Input parameters:
  "string_ids"    a vector of STRING identifiers, sorted.


Author(s):
   Andrea Franceschini
'


        if(is.null(graph)) load()
        return(ppi_enrichment(string_ids))
      },
      
      #########################################
      ## get_subnetwork
      #########################################
     
      
      get_subnetwork = function(string_ids){

'
Description:
  Returns the subgraph generated by the given input proteins.

Input parameters:
  "string_ids"      a vector of STRING identifiers.

Author(s):
   Andrea Franceschini
'

        if(is.null(graph)) load()
        return(induced.subgraph(graph, which(V(graph)$name %in% string_ids)))
      },
      

      #########################################
      ## get_annotations
      #########################################
      
      
      get_annotations = function(string_ids=NULL) {

'
Description:
  Returns the term annotation (GO, UniProt Keywords, PFAM etc.) of the given list of proteins.

Input parameters:
  "string_ids"      vector of STRING identifiers. If the variable is set, the method returns only the proteins that are present in this vector.

Author(s):
   Damian Szklarczyk
'


        string_ids = unique(string_ids)
        string_ids = string_ids[!is.na(string_ids)]

        if(length(string_ids) > 2000) {
          cat("ERROR: We do not support lists with more than 2000 genes.\nPlease reduce the size of your input and rerun the analysis. \t")
          stop()
        }

        urlStr = paste(stable_url, "/api/tsv-no-header/functional_annotation", sep="")
        identifiers=""
        for(id in string_ids ){ identifiers = paste(identifiers, id, sep="%0d")}
        params = list(identifiers=identifiers, species=species, caller_identity='STRINGdb-package')

        response <- postFormSmart(  urlStr, .params=params)

        ann = read.table(text=response, sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE, header=FALSE)

        ann = renameColDf(ann, "V1", "category")
        ann = renameColDf(ann, "V2", "term_id")
        ann = renameColDf(ann, "V3", "number_of_genes")
        ann = renameColDf(ann, "V4", "ratio_in_set")
        ann = renameColDf(ann, "V5", "species")
        ann = renameColDf(ann, "V6", "string_ids")
        ann = renameColDf(ann, "V7", "preferred_names")
        ann = renameColDf(ann, "V8", "description")

        return(ann)

     },
      
      #########################################
      ## load
      #########################################
 
      
      load = function() {

'
Description:
  Downloads and returns the STRING network (the network is set also in the graph variable of the STRING_db object).
  
It makes use of the variables:
    "backgroundV"         vector containing STRING identifiers to be used as background 
                            (i.e. the STRING network loaded will contain only the proteins that are present also in this vector)
    "score_threshold"     STRING combined score threshold (the network loaded contains only interactions having a combined score greater than this threshold)

Author(s):
   Andrea Franceschini
'
        

        url <- paste("https://stringdb-static.org/download/protein.links.v", version, "/", species, ".protein.links.v", version, ".txt.gz", sep="")
        temp = downloadAbsentFile(url, oD=input_directory)
        PPI <- read.table(temp, sep = " ", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
        
        PPIselected = PPI
        if(length(score_threshold)!=0) PPIselected <- PPI[PPI$combined_score >= score_threshold,]
        if(!is.null(backgroundV)){
          PPIselected <- PPIselected[PPIselected$protein1 %in% backgroundV,]
          PPIselected <- PPIselected[PPIselected$protein2 %in% backgroundV,]
        }
        
        myg = graph.data.frame(PPIselected,FALSE)
        graph <<- myg
        return(myg)
      },
      
      
      #########################################
      ## mp
      #########################################
 
      
      mp = function(protein_aliases){

'
Description:
    Maps the gene identifiers of the input vector to STRING identifiers (using a take first approach).
    It returns a vector with the STRING identifiers of the mapped proteins

Input parameters:
    "protein_aliases"       vector with the aliases of the proteins that we want to map to STRING

Author(s):
    Andrea Franceschini
'
        
        temp_df = data.frame(proteins=protein_aliases)
        temp_df_mapped = map(temp_df, "proteins", removeUnmappedRows=TRUE, quiet=TRUE)
        return(temp_df_mapped$STRING_id)
        
      },
      
      


      #########################################
      ## map
      #########################################


      map = function(my_data_frame,
                              my_data_frame_id_col_names,
                              takeFirst=TRUE, removeUnmappedRows=FALSE, quiet=FALSE
      ){

'
Description:
  Maps the gene identifiers of the input dataframe to STRING identifiers.
  It returns the input dataframe with the "STRING_id" additional column.

Input parameters:
  "my_data_frame"                 data frame provided as input. 
  "my_data_frame_id_col_names"    vector contatining the names of the columns of "my_data_frame" that have to be used for the mapping.
  "takeFirst"                     boolean indicating what to do in case of multiple STRING proteins that map to the same name. 
                                      If TRUE, only the first of those is taken. Otherwise all of them are used. (default TRUE)
  "removeUnmappedRows"            remove the rows that cannot be mapped to STRING 
                                      (by default those lines are left and their STRING_id is set to NA)
  "quiet"                         Setting this variable to TRUE we can avoid printing the warning relative to the unmapped values.

Author(s):
   Andrea Franceschini
'

        aliasDf2=get_aliases(takeFirst)

        tempDf = multi_map_df(my_data_frame, aliasDf2, my_data_frame_id_col_names, "alias", "STRING_id")
        naDf = subset(tempDf, is.na(STRING_id))
        if(nrow(naDf) > 0 & !quiet) cat(paste("Warning:  we couldn't map to STRING ", as.integer((nrow(naDf)/nrow(tempDf))*100), "% of your identifiers" , sep=""))
        if(removeUnmappedRows) tempDf = subset(tempDf, !is.na(STRING_id))

        return(tempDf)
        
      },
     



      #########################################
      ## remove_homologous_interactions
      #########################################




    remove_homologous_interactions = function(interactions_dataframe, bitscore_threshold = 60){

'
Description:
  With this method it is possible to remove the interactions that are composed by a pair of homologous/similar proteins, having a similarity bitscore between each other higher than a threshold.

Input parameters:
  "interactions_dataframe"  a data frame of protein-protein interactions, having the following column names: "proteinA", "proteinB", "score"
  "bitscore_threshold"       bitscore threshold used to detect the homology (two proteins are considered homologs, if their similarity bitscore is higher than this threshold)

Author(s):
   Andrea Franceschini
'

        if(!(names(interactions_dataframe)[1]=="proteinA" && names(interactions_dataframe)[2]=="proteinB" && names(interactions_dataframe)[3]=="score")){
          cat("ERROR: \nthe input dataframe should contain an header with the following names: proteinA, proteinB, score")
          stop()
        }
        interactions_dataframe=data.frame(interactions_dataframe, score_temp_sort_abc = seq(1:nrow(interactions_dataframe)))
        interactionsGraph = graph.data.frame(interactions_dataframe,FALSE)
        hg=get_homology_graph(bitscore_threshold)
        nonHomologousIntGraph = graph.difference(interactionsGraph, hg)
        nonHomologousIntDf=get.data.frame(nonHomologousIntGraph)
        names(nonHomologousIntDf)=names(interactions_dataframe)
        nonHomologousIntDf=arrange(nonHomologousIntDf, score_temp_sort_abc)
        nonHomologousIntDf=delColDf(nonHomologousIntDf, "score_temp_sort_abc")
        cat("Discarded", nrow(interactions_dataframe)-nrow(nonHomologousIntDf), "interactions (", 
            (nrow(interactions_dataframe)-nrow(nonHomologousIntDf))*100/nrow(interactions_dataframe),"% )", " between homologous proteins\n")

        return(nonHomologousIntDf)

    },


      #########################################
      ## set_background
      #########################################


      
      set_background = function(background_vector){

'
Description:
  With this method you can specify a vector of proteins to be used as background. 
  The network is reloaded and only the proteins that are present in the background vector are inserted in the graph.  
  Besides, the background is taken in consideration for all the enrichment statistics.

Input parameters:
  "background_vector"     vector of STRING protein identifiers

Author(s):
   Andrea Franceschini
'
        
        backgroundV <<- background_vector
        load()
      },
      
  
      show = function(){

        cat(paste(
          "***********  STRING - https://string-db.org   ***********", "\n",
          "(Search Tool for the Retrieval of Interacting Genes/Proteins)  ", "\n",
          "version: ", version, "\n",
          "species: ", species, "\n", 
          "............please wait............\n", sep=""
        ))
        
        cat(paste(
            "proteins: ", nrow(get_proteins()), "\n",
            "interactions: ", length(E(get_graph())), sep=""
                )
            )
        
      },



benchmark_ppi = function(interactions_dataframe, pathwayType = "KEGG", max_homology_bitscore = 60, precision_window=400, exclude_pathways="blacklist"){
    .Deprecated('Contact developers to request functionality')
 
##   'Description:  
##     benchmark a list of protein-protein interactions using pathways (e.g. KEGG). 
##     The function outputs a table where the interactions are mapped to KEGG and the number of TPs and FPs are counted.
## 
## Input parameters:
##       "interactions_dataframe"  a data frame of protein-protein interactions, having the following column names: "proteinA", "proteinB", "score"
##       "pathwayType"    the annotation category that we want to use for benchmarking
##       "max_homology_bitscore"  if this variable is set, we remove the homologous proteins from the benchmark that have a bitscore higher than this variable
##       "precision_window"   size of a window used to estimate the precision (i.e. the list of protein-protein interactions is scanned with a window
##                           and the precision is estimated only for the proteins in the window. The precision value that is reported in every row is the middle point of the window).
##                           At the beginning and at the end of the ppi list, the window is automatically expanded (and reduced).
##       "exclude_pathways"   Vector of terms (i.e. pathways) to exclude from benchmarking. By default, the vector is set to "blacklist" which means that a curated set 
##                                 of "dangerous" pathways is removed from the analysis.
## 
## Author(s):
##    Andrea Franceschini
## 
##   '
##   
##        if(!(names(interactions_dataframe)[1]=="proteinA" && names(interactions_dataframe)[2]=="proteinB" && names(interactions_dataframe)[3]=="score")){
##           cat("ERROR: \nthe input dataframe should contain an header with the following names: proteinA, proteinB, score")
##           stop()
##         }
##         if(!is.null(max_homology_bitscore)) {interactions_dataframe = remove_homologous_interactions(interactions_dataframe, max_homology_bitscore)}
##         
##         ann=get_annotations()
##         annCat=subset(ann, category==pathwayType)
##         if(!is.null(exclude_pathways)){
##           bbl=exclude_pathways
##           if("blacklist" %in% exclude_pathways) {
##             bbl=c(bbl, (get_pathways_benchmarking_blackList())$term_id)
##           }
##           annCat = subset(annCat, !(term_id %in% bbl))
##         }
##         protein_pathways = hash()
##         temp=apply(annCat, 1,
##                    function(x){
##                      paths=c()
##                      if(!is.null(protein_pathways[[x[1]]])){paths=protein_pathways[[x[1]]]}
##                      paths=c(paths, x[2])
##                      protein_pathways[[x[1]]]=paths
##                    }
##         )
##         
##         
##         mappedNonHomologousIntDf=subset(interactions_dataframe, proteinA %in% unique(annCat$STRING_id) & proteinB %in% unique(annCat$STRING_id))
##         cat("NOTE: ", nrow(mappedNonHomologousIntDf), "interactions (", round((nrow(mappedNonHomologousIntDf)/nrow(interactions_dataframe))*100, digits=1), "% )", "have been mapped to", pathwayType, "\n")
##         if(nrow(mappedNonHomologousIntDf)<precision_window){
##           precision_window=nrow(mappedNonHomologousIntDf)-1
##           cat("WARNING: The size of the precision window is larger than the number of mapped interactions. \n     We reduce the size of the precision window to ", precision_window)
##         }
##   
##         outcome=rep("FP", nrow(mappedNonHomologousIntDf))
##         numTP=rep(0, nrow(mappedNonHomologousIntDf))
##         numFP=rep(0, nrow(mappedNonHomologousIntDf))
##         pathways=rep(NA, nrow(mappedNonHomologousIntDf))
##         acc=rep(0, nrow(mappedNonHomologousIntDf))
##         counter=0
##         tps=0
##         fps=0
##         proteinsA=mappedNonHomologousIntDf$proteinA
##         proteinsB=mappedNonHomologousIntDf$proteinB
##         for(counter in seq(1:nrow(mappedNonHomologousIntDf))){
##           x=c(proteinsA[counter], proteinsB[counter])
##           paths=c()
##           if(!is.null(protein_pathways[[x[1]]]) && !is.null(protein_pathways[[x[1]]]) && length(intersect(protein_pathways[[x[1]]],protein_pathways[[x[2]]]))>0){
##             outcome[counter]="TP"
##             tps=tps+1
##             pathways[counter]=paste((intersect(protein_pathways[[x[1]]],protein_pathways[[x[2]]])), collapse=" ")
##           }else{
##             fps=fps+1
##           }
##           numTP[counter]=tps
##           numFP[counter]=fps
##           
##           if(counter<=precision_window){
##             acc[floor(counter/2)] = tps/(tps+fps)
##           }else{
##             acc[floor(counter - (precision_window/2))] = (tps-numTP[counter-precision_window])/((tps-numTP[counter-precision_window])+(fps-numFP[counter-precision_window]))
##           }
##           
##           if(counter>=(nrow(mappedNonHomologousIntDf)-precision_window/2)){
##             acc[counter] = (tps-numTP[ceiling(counter-((precision_window)/2))])/((tps-numTP[ceiling(counter-((precision_window)/2))])+(fps-numFP[ceiling(counter-((precision_window)/2))]))
##           }
##           
##         }
##         mappedNonHomologousIntDf2 = data.frame(nTP=numTP, nFP=numFP, mappedNonHomologousIntDf, outcome=outcome, precision=acc, pathways = pathways, stringsAsFactors=FALSE)
##         
##         return(mappedNonHomologousIntDf2)
},


benchmark_ppi_pathway_view = function(benchmark_ppi_data_frame, precision_threshold=0.2, pathwayType = "KEGG"){
    .Deprecated('Contact developers to request functionality')

## 'Description:  
##     Takes in input the results of the "benchmark_ppi" function, and constructs a new table that provides a view at the pathway level 
##     (i.e. it list all the pathways of the interactions)
## 
## Input parameters:
##       "benchmark_ppi_data_frame"  a data frame that comes out as result from the "benchmark_ppi" function
##       "precision_threshold"    precision threshold after which to stop reading the interactions.
##      
## Author(s):
##    Andrea Franceschini
## 
##   '
##   
##   ann=get_annotations()
##   annCat=subset(ann, category==pathwayType)
##   pathway_proteinsNum = hash()
##   temp=apply(annCat, 1,
##              function(x){
##                paths=0
##                if(!is.null(pathway_proteinsNum[[x[2]]])){paths=pathway_proteinsNum[[x[2]]]}
##                paths=paths+1
##                pathway_proteinsNum[[x[2]]]=paths
##              }
##   )
##   
##   benchmark_ppi_data_frame2=subset(benchmark_ppi_data_frame,!is.na(pathways))
##   totalInteractions = 0
##   pathway_interactions=hash()
##   inside=TRUE
##   precisions = benchmark_ppi_data_frame2$precision
##   pathways = as.character(benchmark_ppi_data_frame2$pathways)
##   for(i in seq(1:nrow(benchmark_ppi_data_frame2))){
##     if(precisions[i]<=precision_threshold){inside=FALSE}
##     if(inside){
##       totalInteractions=totalInteractions+1
##       splitted=strsplit(pathways[i], " ")[[1]]
##       sapply(splitted, function(y){
##         if(is.null(pathway_interactions[[y]])){pathway_interactions[[y]]=0}
##         pathway_interactions[[y]]=pathway_interactions[[y]]+1
##       }
##       )
##     }
##   }
##   
##   repdf=data.frame(coverage=numeric(), proteins=numeric(), interactions=numeric(), total_representation=numeric(), stringsAsFactors=FALSE)
##   i=0
##   for(ke in keys(pathway_interactions)){
##     i=i+1
##     interactions2 = pathway_interactions[[ke]]
##     proteins2 = pathway_proteinsNum[[ke]]
##     coverage = round(interactions2/((proteins2 * (proteins2-1))/2), digits=3)
##     total_representation = round(interactions2/totalInteractions, digits=6)
##     repdf[i,] = c(coverage, proteins2, interactions2, total_representation)
##   }
##   repdf2=data.frame(pathways=keys(pathway_interactions), repdf, stringsAsFactors=FALSE)
##   repdf3=merge(repdf2, get_annotations_desc(), by.x="pathways", by.y="term_id", all.x=TRUE)
##   return(arrange(repdf3, -total_representation))
},
 
get_homologs = function(string_ids, target_species_id, bitscore_threshold=NULL){
    .Deprecated('Contact developers to request functionality')
##'
##Description:
##  Returns the homologs of the given input identifiers that are present in the given target_species_id.
##
##Input parameters:
##    "string_ids"          a vector of STRING identifiers.
##    "target_species_id"   NCBI taxonomy identifier of the species to query for homologs (the species must be present in the STRING database)
##    "bitscore_threshold"  threshold on the bitscore of the blast alignment.
##
##Author(s):
##   Andrea Franceschini
##
##'
##        
##        if(length(string_ids) > 300) {
##          cat("ERROR: We support a maximum of 300 STRING identifiers per call. Please reduce the size of the input and try again. \t")
##          stop()
##        }
##        urlStr = paste("http://string-db.org/version_", version, "/newstring_cgi/webservices/homology_hits.pl", sep="")
##        identifiers=""
##        for(id in string_ids ){
##          identifiers = paste(identifiers, id, "%0D", sep="")
##        }
##        params = list(target_species_id=target_species_id, identifiers=identifiers)
##        if(!is.null(bitscore_threshold)) params["bitscore_threshold"] = bitscore_threshold
##        tempDfv=postFormSmart(urlStr, .params=params)
##        hhits <- read.table(text=tempDfv, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
##        return(hhits)
##      },
##      
##      
##      get_homology_graph = function(min_homology_bitscore=60){
##        temp = downloadAbsentFile(paste("http://string.uzh.ch/permanent/string/", version, "/homology_links/",
##                                        species, "_homology_links.tsv.gz", sep=""), oD=input_directory)
##        PPI <- read.table(temp, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
##        
##        PPIselected = PPI
##        if(!is.null(min_homology_bitscore)) PPIselected <- PPI[PPI$bitscore >= min_homology_bitscore,]
##        if(!is.null(backgroundV)){
##          PPIselected <- PPIselected[PPIselected$protein1 %in% backgroundV,]
##          PPIselected <- PPIselected[PPIselected$protein2 %in% backgroundV,]
##        }
##        
##        myg = graph.data.frame(PPIselected,FALSE)
##        homology_graph <<- myg
##        return(myg)
},
 

get_link = function(string_ids, required_score=NULL, network_flavor="evidence", payload_id = NULL){
    .Deprecated('Contact developers to request functionality')

##'
##Description:
##  Returns a short link to the network page of our STRING website that shows the protein interactions between the given identifiers.
##
##Input parameters:
##  "string_ids"        a vector of STRING identifiers.
##  "required_score"    minimum STRING combined score of the interactions 
##                        (if left NULL we get the combined score of the object, which is 400 by default)
##  "network_flavor"    specify the flavor of the network ("evidence", "confidence" or "actions".  default "evidence")
##
##Author(s):
##   Andrea Franceschini
##'
##        if(length(string_ids) > 400) {
##          cat("ERROR: We do not support lists with more than 400 genes.\nPlease reduce the size of your input and rerun the analysis. \t")
##          stop()
##        }
##        if(is.null(required_score) ) required_score = score_threshold
##        string_ids = unique(string_ids)
##        urlStr = paste("http://string-db.org/version_",version,"/newstring_cgi/link_to.pl", sep="" )
##        #urlStr = paste("http://130.60.240.90:8380/newstring_cgi/link_to.pl", sep="" )
##        identifiers=""
##        for(id in string_ids ){   identifiers = paste(identifiers, id, "%0D", sep="")}
##        params=list(required_score=required_score, limit=0, network_flavor=network_flavor, identifiers=identifiers)
##        if(!is.null(payload_id)) params["internal_payload_id"] = payload_id
##        tempDfv=postFormSmart(urlStr, .params=params)
##        df=read.table(text=tempDfv, stringsAsFactors=FALSE, fill = TRUE)
##        return(df$V1)
##

},

get_pathways_benchmarking_blackList = function(){
    .Deprecated('Contact developers to request functionality')
##        '
##Description:
##  Returns a list of pathways that we suggests to exclude for benchmarking (at the moment our list includes only KEGG maps)
##
##'
##        if(nrow(pathways_benchmark_blackList)==0){
##          temp = downloadAbsentFile(paste("http://string.uzh.ch/permanent/string/", version, "/benchmarking_blacklist.tsv", sep=""), oD=input_directory)
##          pathways_benchmark_blackList <<- read.table(temp, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE, quote="", colClasses=c("character", "character"))
##        }else{return(pathways_benchmark_blackList)}

},


get_ppi_enrichment_full = function(string_ids, sliceWindow = 20, edgeWindow  = 140, windowExtendedReferenceThreshold = 260, growingWindowLimit=NULL){
    .Deprecated('Contact developers to request functionality')
##'
##Description:
##  Returns a vector showing the enrichment in protein interactions in various positions of the list of genes in input.
##  In practice, a list of 3 vectors is returned:
##  1) enrichment  (i.e.  enrichment computed in the window from 1 to x)
##  2) enrichmentWindow (i.e. enrichment computed in a sliding window of size determined by the "edgeWindow" parameters 
##  and the sliding steps determined by the "sliceWindow" parameter)
##  3) enrichmentWindowExtended  (i.e. like the enrichmentWindow, 
##  but it also includes an initial window of size "windowExtendedReferenceThreshold" with respect to which to compute the enrichment )
##
##Input parameters:
##  "string_ids"                          vector of STRING identifiers, sorted.
##  "file"                                file where to save the graph as an image
##  "sliceWindow"                         defines the interval in proteins after which to compute the enrichment, scanning the list (i.e. the resolution)
##  "windowExtendedReferenceThreshold"    defines the size of a window at the beginning of the list. 
##                                              The enrichment will be computed always including the proteins in this window
##  "title"                               title of the graph.
##  "growingWindowLimit"                  threshold where to stop the computation of the enrichment
##
##Author(s):
##   Andrea Franceschini
##
##'  
##        if(is.null(graph)) load()
##        return(ppi_enrichment_full(string_ids, graph, sliceWindow = sliceWindow, edgeWindow  = edgeWindow, 
##                                   windowExtendedReferenceThreshold = windowExtendedReferenceThreshold, growingWindowLimit=growingWindowLimit))
},


get_pubmed = function(string_ids){
    .Deprecated('Contact developers to request functionality')

##'
##Description:
##  Returns vector with the PUBMED IDs of the publications that contain the names of the proteins in the input vector.
##
##Input parameters:
##  "string_ids"      a vector of STRING identifiers.
##
##Author(s):
##   Andrea Franceschini
##'        
##        if(length(string_ids) > 300) {
##          cat("ERROR: We support a maximum of 300 STRING identifiers per call. Please reduce the size of the input and try again. \t")
##          stop()
##        }
##        #urlStr = paste("http://string-db.org/version_", version, "/api/tsv/abstracts?limit=1000000&identifiers=", sep="")
##        urlStr = paste("http://string-db.org/version_", version, "/newstring_cgi/webservice_handler.pl", sep="")
##        identifiers=""
##        for(id in string_ids ){identifiers = paste(identifiers, id, "%0D", sep="")}
##        params = list(limit=1000000, identifiers=identifiers, output="tsv", request="abstracts")
##        tempDfv=postFormSmart(urlStr, .params=params)
##        pubmedIdsDf <- read.table(text=tempDfv, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
##        return(pubmedIdsDf$abstractId)

},
      
      
get_pubmed_interaction = function(STRING_id_a, STRING_id_b){
    .Deprecated('Contact developers to request functionality')

## '
## Description:
##   Returns vector with the PUBMED IDs of the publications that contain the names of both the input proteins.
## 
## Input parameters:
##   "STRING_id_a"      STRING identifier.
##   "STRING_id_b"      STRING identifier.
## 
## Author(s):
##    Andrea Franceschini
## '             
##         return(intersect(get_pubmed(STRING_id_a), get_pubmed(STRING_id_b) ) )

},
 
get_term_proteins = function(term_ids, string_ids=NULL, enableIEA = TRUE){
    .Deprecated('get_annotations')

##'
##Description:
##  Returns the proteins annotated to belong to a given term.
##
##Input parameters:
##  "term_ids"        vector of terms 
##  "string_ids"      vector of STRING identifiers. If the variable is set, the method returns only the proteins that are present in this vector.
##  "enableIEA"       whether to consider also Electronic Inferred Annotations
##
##Author(s):
##   Andrea Franceschini
##'
##        
##        annotations2 = get_annotations()
##        if(!enableIEA) { annotations2 = subset(annotations2, type!="IEA") }
##        annotations3 = subset(annotations2, term_id %in% term_ids, select=c("STRING_id", "term_id"))
##        if(!is.null(string_ids)){
##          annotations3 = subset(annotations3, STRING_id %in% string_ids)
##        }
##        annotations4 = add_proteins_description(annotations3)
##        return(annotations4)
},
 
plot_ppi_enrichment = function(string_ids, file=NULL, sliceWindow = 20, edgeWindow = 140, 
                                     windowExtendedReferenceThreshold = 260, minVal=0.0000000001, title="", quiet=FALSE){
    .Deprecated('Contact developers to request functionality')
## '
## Description:
##   Plots a graph showing the enrichment in protein interactions in various positions of the list of genes in input.
## 
## Input parameters:
##   "string_ids"                          vector of STRING identifiers, sorted.
##   "file"                                file where to save the graph as an image
##   "sliceWindow"                         defines the interval in proteins after which to compute the enrichment, scanning the list (i.e. the resolution)
##   "windowExtendedReferenceThreshold"    defines the size of a window at the beginning of the list. 
##                                           The enrichment will be computed always including the proteins in this window.
##   "minVal"                              minimum value of the pvalue (lower values with respect to this one will assume this minimum value.)
##   "title"                               title of the graph.
## 
## Author(s):
##    Andrea Franceschini
## 
## '
##           if(is.null(graph)) load()
##           plot_ppi_enrichment_graph(string_ids, graph, file=file, sliceWindow = sliceWindow, edgeWindow = edgeWindow, 
##                                     windowExtendedReferenceThreshold = windowExtendedReferenceThreshold, minVal=minVal, title=title, quiet=quiet)
},
      
enrichment_heatmap = function(genesVectors, vectorNames, output_file=NULL, title="", 
                                  enrichmentType="Process", limitMultiPicture=NULL, fdr_threshold=0.05, pvalue_threshold=NULL, 
                                  cexRow=NULL, cexCol=1, selectTermsVector=NULL, iea = TRUE, 
                                  sortingMethod="rowMeans", useInputAsBackground=FALSE, limit=NULL){

      .Deprecated('Contact developers to request functionality')
      
##       enrichList = list()
##       
##       if(useInputAsBackground) {
##         genes=genesVectors[[1]]
##         for(i in 2:length(genesVectors)){ genes=intersect(genes, genesVectors[[i]])  }
##         for(i in 1:length(genesVectors)){genesVectors[[i]] = intersect(genesVectors[[i]], genes)}
##         cat("INFO: we are benchmarking ",length(genes), "genes\n")
##       }else{
##         cat("WARNING: Your input has not been intersected; hence the datasets could not be perfectly comparable.\n")
##         for(i in 1:length(genesVectors)){genesVectors[[i]] = unique(genesVectors[[i]])}
##       }
##       
##       for(i in 1:length(genesVectors)){
##         genesInput = genesVectors[[i]]
##         if(!is.null(limit)) genesInput = genesVectors[[i]][0:limit]
##         stringGenes = mp(genesInput)
##         enrichList[[i]] = get_enrichment(stringGenes, category = enrichmentType, methodMT = "fdr", iea = iea )
##       }
##       enrichHash = hash()
##       for(i in 1:length(genesVectors)){
##         for(j in 1:nrow(enrichList[[i]])){
##           if((is.null(pvalue_threshold) || enrichList[[i]]$pvalue[j] <= pvalue_threshold) &&
##                (is.null(fdr_threshold) || enrichList[[i]]$pvalue_fdr[j] <= fdr_threshold) 
##           ){
##             enrichVect = rep(NA, length(genesVectors))
##             myterm = as.character(enrichList[[i]]$term_description[j])
##             
##             enter = FALSE
##             if(!is.null(selectTermsVector) ){ 
##               for(term in selectTermsVector){if(grepl(term, myterm)){enter=TRUE}}
##             }else{enter = TRUE}
##             
##             if(enter){
##               if(has.key(myterm, enrichHash)){  enrichVect = enrichHash[[myterm]] }
##               enrichVect[i] = -log(enrichList[[i]]$pvalue_fdr[j])
##               enrichHash[[myterm]] = enrichVect
##             }
##           }
##         }
##       }
##       
##       enrichMatr = matrix(NA, length(enrichHash), length(genesVectors))
##       rownames(enrichMatr) = rep(NA, length(keys(enrichHash)))
##       if(length(keys(enrichHash))>0){
##         for(i in 1:length(keys(enrichHash))){
##           k=keys(enrichHash)[i]
##           enrichMatr[i,] = enrichHash[[k]]
##           rownames(enrichMatr)[i] = k
##         }
##       }
##       colnames(enrichMatr) = vectorNames
##       
##       
##       if(is.null(limitMultiPicture)){limitMultiPicture= nrow(enrichMatr)+2}
##       if(!is.null(sortingMethod) && sortingMethod=="rowMeans" && nrow(enrichMatr)>1){
##         enrichMatr = enrichMatr[order(rowMeans(enrichMatr, na.rm=TRUE), decreasing=TRUE),]
##       }
##       
##       if(is.null(cexRow)){
##         if(nrow(enrichMatr) <= 60 ) {
##           cexRow=0.4 + 1/log(nrow(enrichMatr))
##         }else if(nrow(enrichMatr) <= 120 ) {
##           cexRow=0.2 + 1/log(nrow(enrichMatr))
##         }else if(nrow(enrichMatr) <= 200 ) {
##           cexRow=0.1 + 1/log2(nrow(enrichMatr))
##         }else{
##           cexRow=0.1
##         }
##         
##       }
##       
##       lmat = rbind(c(0,3),c(2,1),c(0,4))
##       lwid = c(0.7,5)
##       lhei = c(1,5,0.7)
##       
##       if(nrow(enrichMatr)%%limitMultiPicture == 1 && nrow(enrichMatr)!=1) limitMultiPicture=limitMultiPicture+1
##       if(nrow(enrichMatr)>1){
##         if(!is.null(output_file)) pdf(output_file, paper="a4", width=12, height=12, pagecentre=FALSE)
##         suppressWarnings(par(new=TRUE))
##         for(i in 1:ceiling(nrow(enrichMatr)/limitMultiPicture)){
##           enrichMatrTemp=enrichMatr[(limitMultiPicture*(i-1)+1):min(limitMultiPicture*(i), nrow(enrichMatr)),]
##           #heatmap(enrichMatrTemp, Rowv=NA, Colv = NA, col = brewer.pal(6,"Blues"), scale = "none",
##           #         margins = c(7,30),  cexRow=cexRow, cexCol=cexCol, main=paste("                                     ",title, sep=""))
##           tempV = c()
##           for(j in 1:ncol(enrichMatrTemp)){tempV = c(tempV, enrichMatrTemp[,j])}
##           tempV=unique(tempV[!is.na(tempV)])
##           
##           if(length(tempV)==1){
##             enrichMatrTemp[is.na(enrichMatrTemp)] <- 0 
##           }  
##           suppressWarnings(heatmap.2(enrichMatrTemp, density.info="none", trace="none", keysize=1,lmat = lmat, lhei=lhei, lwid=lwid,Rowv=NA, Colv = NA, col = brewer.pal(6,"Blues"), scale = "none",
##                                      margins = c(7,30),  cexRow=cexRow, cexCol=1, main=paste("                ",title, sep="")))
##         }  
##         if(!is.null(output_file)) dev.off()
##         suppressWarnings(par(new=FALSE))
##       }else if(nrow(enrichMatr)==1){
##         cat("Only one term has been found to be significantly enriched:\n")
##         print(enrichMatr)
##       }else{cat("No enriched terms below the p-value threshold.")}
##       return(enrichMatr)

}

)
)
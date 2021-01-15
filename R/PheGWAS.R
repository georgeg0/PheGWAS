# once added the roxygen syntax things then run following
# devtools::document()
BP <- CHR <- P <- PHENO <- SNP <- gene <- lab <- logp <- Entire_Val <- LDblock <- z <- bpdivision <-  NULL
flag <- FALSE
utils::globalVariables(c("phenos"))

## Function to add gene for repective rsid, if genes are not provided by user
addgene <- function(gwasmulti){

  ensembl <- useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")

  human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

  gwasmulti$gene <- "NA"
  for(i in rownames(gwasmulti)){
    xx <- getBM(attributes=c( "ensembl_gene_stable_id"),filters="snp_filter", values=gwasmulti[i, "SNP"],mart=ensembl, uniqueRows=TRUE)
    if (nrow(xx) == 0){
      gwasmulti[i,]$gene <- "NA"
    }else{
      gene <- getBM(attributes = c("hgnc_symbol"),
                    filters = "ensembl_gene_id", values = unlist(xx), mart = human)
      gwasmulti[i,]$gene <- if (nrow(gene) == 0) "NA" else  toString( unlist(gene))
    }}
  print("Finished mapping")
  return(gwasmulti)
}

#' Prepare the dataframe to pass to landscape function
#'
#' @import tidyverse
#' @import tidyr
#' @import dplyr
#' @importFrom stringr str_trim
#' @param x List of dataframes that need to do PheGWAS on. Arrange the dataframe in the order how the the phenotypes should align in y axis
#' @param phenos a vector of phenotypes that you are passing for in the summary statistics file. It should be in the order that you pass the List.
#' @param chromosome The chromosome in which the region lies
#' @param lab The region within the chromosme to export
#' @param LDblock If want to pass a custom LDblock file for division of BP groups (applicable only for chromosomal level)
#' @return A processed dataframe to pass to PheGWAS function
#' @details Make sure there are no duplicate rsid's in any of the dataframe, If there aremake sure to resolve it before passing it to this function.
#' @author George Gittu
#' @examples
#' \dontrun{
#' xprocess <- exportregion(x,phenos,chromosome = 19, lab = 11 )
#' }
#' @export
exportregion <- function(x,phenos, chromosome= NULL,lab = NULL,LDblock = NULL) {

  print(paste("Exporting the region ", lab, "in chromosome ", chromosome))

  for(i in 1:length(x))
  {
    x[[i]] <- x[[i]] %>% rename_at(vars(-(1:5)), ~ paste0(phenos[i],"_",.))
    x[[i]] <-     x[[i]][x[[i]]$CHR==chromosome,]

    if (is.null(LDblock)){
      x[[i]]$lab <- "NULL"
      x[[i]]$lab <- as.integer(x[[i]]$BP /bpdivision, 0)

    }else{
      bpfull = read.table(LDblock,header=TRUE,sep = "\t")
      bpp <- bpfull[str_trim(bpfull$chr)==paste0("chr",chromosome),]
      x[[i]]$lab <- "NULL"
      for (j in 1:nrow(bpp)) {
        print(bpp$start[j])
        x[[i]][x[[i]]$BP >= bpp$start[j] & x[[i]]$BP < bpp$stop[j],]$lab <- bpp$start[j]

      }
    }
    x[[i]] <- x[[i]][x[[i]]$lab==lab,]
  }
  x
}## end of the function export region


#' Prepare the dataframe to pass to landscape function
#'
#' @import tidyverse
#' @import tidyr
#' @import dplyr
#' @param x List of dataframes that need to do PheGWAS on. Arrange the dataframe in the order how the the phenotypes should align in y axis
#' @param phenos a vector of phenotypes that you are passing for in the summary statistics file. It should be in the order that you pass the List.
#' @return A processed dataframe to pass to PheGWAS function
#' @details Make sure there are no duplicate rsid's in any of the dataframe, If there aremake sure to resolve it before passing it to this function.
#' @author George Gittu
#' @examples
#' \dontrun{
#' x <- list(hdl,ldl,trig,tchol)
#' phenos <- c("HDL","LDL","TRIGS","TOTALCHOLESTROL")
#' ## y is ready to be passed to function landscape
#' y <- processphegwas(x, phenos)
#' }
#' @export
processphegwas <- function(x,phenos) {
  #rename the column for not to get duplicate column names
  for(i in 1:length(x))
  {
    x[[i]] <- x[[i]] %>% rename_at(vars(-(1:5)), ~ paste0(phenos[i],"xx_xx",.))
  }

  fulldf <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("CHR","BP","rsid","A1","A2"), all = TRUE),x)

  if(length(grep("gene", colnames(fulldf))) > 0){
    ## means it got gene added by the user
    print("Processing for gene provided by the user")
    fulldf$G <- apply(fulldf[, grep("gene", colnames(fulldf))], 1, function(x) paste(unique(unlist(strsplit(x[!is.na(x)], ";"))), collapse = ";"))
    fullinterm <- fulldf[, -grep("gene", colnames(fulldf))]
    ##Doing the merge of columns
    for(i in 1:length(x))
    {
      cccc <- colnames(fullinterm)[grepl( phenos[i], names( fullinterm ) ) ]
      fullinterm <- unite_(fullinterm,phenos[i], cccc, sep = "and", remove = FALSE)
    }
    fulldfff <- fullinterm[, -grep("xx_xx", colnames(fullinterm))]
    gwasmulti.meltF <- reshape2::melt(data = fulldfff, id.vars = c("CHR","BP","rsid","A1","A2","G"), variable.name = "color", value.name = "Entire_Val")
    gwasmulti.melt <- gwasmulti.meltF %>% separate(Entire_Val, c("BETA", "SE","p_value"), "and")
    d <- data.frame(CHR = gwasmulti.melt$CHR, BP = gwasmulti.melt$BP, A1 = gwasmulti.melt$A1, A2 = gwasmulti.melt$A2, SNP = gwasmulti.melt$rsid ,P = as.numeric(gwasmulti.melt$p_value),BETA = as.numeric(gwasmulti.melt$BETA),SE = as.numeric(gwasmulti.melt$SE),gene = gwasmulti.melt$G,PHENO = gwasmulti.melt$color)
  }else{
    ## Gene when not addded
    print("Gene will be added as part of BioMart module")
    fullinterm <- fulldf
    ##Doing the merge of columns
    for(i in 1:length(x))
    {
      cccc <- colnames(fullinterm)[grepl( phenos[i], names( fullinterm ) ) ]
      fullinterm <- unite_(fullinterm,phenos[i], cccc, sep = "and", remove = FALSE)
    }
    fulldfff <- fullinterm[, -grep("xx_xx", colnames(fullinterm))]
    gwasmulti.meltF <- reshape2::melt(data = fulldfff, id.vars = c("CHR","BP","rsid","A1","A2"), variable.name = "color", value.name = "Entire_Val")
    gwasmulti.melt <- gwasmulti.meltF %>% separate(Entire_Val, c("BETA", "SE","p_value"), "and")
    d <- data.frame(CHR = gwasmulti.melt$CHR, BP = gwasmulti.melt$BP, A1 = gwasmulti.melt$A1, A2 = gwasmulti.melt$A2, SNP = gwasmulti.melt$rsid ,P = as.numeric(gwasmulti.melt$p_value),BETA = as.numeric(gwasmulti.melt$BETA),SE = as.numeric(gwasmulti.melt$SE),PHENO = gwasmulti.melt$color)
  }
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  d$logp <- round(-log10(d$P),3)
  d$BETA <- round(d$BETA,3)
  d$logp[is.infinite(d$logp)] <- max(d$logp[!is.infinite(d$logp)],na.rm = TRUE) + 10
  d
}

#' Interactive 3-D association landscape for many phenotypes
#'
#' @import tidyverse
#' @import tidyr
#' @import plotly
#' @importFrom biomaRt useMart getBM
#' @import reshape2
#' @importFrom httr GET stop_for_status content_type content
#' @importFrom utils read.table
#' @import xml2
#' @import jsonlite
#' @importFrom stringr str_trim
#' @param d DataFrame output from processphegwas
#' @param sliceval Integer to indicate value of -log10(p) to do the sectionalcut. Usually value > -log10 6 is considered to be significant
#' @param chromosome Integer to indicate the chromosome number thats interested, If not given entire chromosome is given
#' @param bpdivision Integer value to indicate the base pair divisions. By default it is 100,000 bp. This vis only applciable on a single chromosome view. For entire
#' @param upperlimit Integer to indicate value of -log10(p) to do the upperlimit. If not given it will take the max -log10(p)
#' @param geneview This checks for the common genes across the section
#' @param betaplot This plot the beta instad of pvalue in the yaxis. The heat map layer is based on the pvalue.
#' @param LDblock If want to pass a custom LDblock file for division of BP groups (applicable only for chromosomal level)
#' @param calculateLD This shoudld be set to true if the calcualte LD logic needed to be added to the plot
#' @param pop The population to select for calculation the LD (default GB)
#' @param R2 The value to set to calculate LD
#' @param D THe value to set to calcualte LD
#' @param mutualLD Calcualte the mutual LD SNP between the phenotypes
#' @param levelsdown Used to find the independant signals
#' chromosoem view the max peak is selected
#' @author George Gittu
#' @examples
#' \dontrun{
#' x <- list(hdl,ldl,trig,tchol)
#' phenos <- c("HDL","LDL","TRIGS","TOTALCHOLESTROL")
#' y <- processphegwas(x, phenos)
#' ## pass the dataframe from the processphegwas
#'
#' # 3D landscape visualization of all the phenotypes across the base pair positions
#' landscape(y)
#'
#' # Showing only the crusts above a certain threshold
#' landscape(y,sliceval = 6)
#'
#' # 3D landscape visualization of chromosome number 16
#' landscape(y,sliceval = 6,chromosome = 16)
#'
#' # 3D landscape visualization of chromosome number 16, with bp division on 1,000,000
#' landscape(y,sliceval = 6,chromosome = 16,bpdivision = 1000000)
#' landscape(y,sliceval = 6,chromosome = 16,bpdivision = 1000000, betaplot = TRUE)
#' }
#' @export
landscape <- function(d, sliceval = 0, chromosome=FALSE, bpdivision= 1000000, upperlimit = NULL,geneview = FALSE,betaplot = FALSE,LDblock = NULL, calculateLD = FALSE, pop = "GBR", R2 = 0.75, D = 0.75, mutualLD = FALSE, levelsdown = 0) {

  if (chromosome == FALSE) {
    print("Processing for the entire chromosome")

    ### this is for the MAtrix to process the entire chromosome...
    gwasmultifull <- d %>%
      group_by(PHENO,CHR) %>%
      slice(which.max(logp))

    if(length(grep("gene", colnames(d))) == 0){
    print("Applying BioMArt module for matching gene to rsid")
    gwasmultifull <- addgene(gwasmultifull)
    }

    gwas_surface_full <- acast(gwasmultifull, gwasmultifull$PHENO ~ gwasmultifull$CHR, value.var = "logp")
    gwas_surface_full_use <- gwas_surface_full
    gwas_surface_prime_full_use <- gwas_surface_full_use
    upperlimit = max(gwas_surface_full_use)

    gwas_surface_prime_full_use[,]= upperlimit + 2
    gwas_surface_copy_full_use <- gwas_surface_full_use

    for (i in 1:nrow(gwas_surface_copy_full_use)) {
      for (j in 1:ncol(gwas_surface_copy_full_use))
        gwas_surface_copy_full_use[i, j] <-
          paste("Phenotype = ", rownames(gwas_surface_copy_full_use)[i], '\n',"Chromosome = ", colnames(gwas_surface_copy_full_use)[j],'\n',
                "BP =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$BP,"\n",
                "A1 =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$A1,' ',
                "A2 =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$A2,'\n',
                "Effect Size =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$BETA,' ',
                "SE =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$SE,'\n',
                "-log10 (P-value) = ", gwas_surface_copy_full_use[i, j],'\n',
                "SNPID =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$SNP,'\n',
                "Gene =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$gene)
    }

    p <- plot_ly(x=colnames(gwas_surface_full_use),y= rownames(gwas_surface_full_use),z = gwas_surface_full_use) %>%
      add_surface(cmin = sliceval,cmax = upperlimit,opacity=.95,surfacecolor = gwas_surface_full_use,text = gwas_surface_copy_full_use,hoverinfo = "text") %>%
      layout(title = 'PheGWAS',scene = list(camera = list(eye = list(x=2, y=2, z=2)),yaxis = list(title = "Phenotypes", tickmode = "array",nticks = 8),xaxis =list(title= "Chromosomes",autotick = F, dtick = 1),
                                              zaxis =list(title= "-log 10(P)",range=c(sliceval,upperlimit+2)),aspectmode = "manual",aspectratio = list( x = 2, y = .8, z = .8))) %>%
      add_trace(x=colnames(gwas_surface_prime_full_use),y= rownames(gwas_surface_prime_full_use),z = gwas_surface_prime_full_use, type = "surface",surfacecolor = gwas_surface_full_use,showscale = FALSE,
                text = gwas_surface_copy_full_use,hoverinfo = "text",cmin = sliceval,cmax = upperlimit) %>% colorbar(title = "-log10 (P-value)")
    p
  }
  else {
    paste0("Processing for the chromosome number ", chromosome)
    ## Pick the chromosome that is passed as param "chromosome"
    dchrom1initial <- d[d$CHR==chromosome,]
if (is.null(LDblock)){
  dchrom1initial$lab <- as.integer(dchrom1initial$BP /bpdivision, 0)
}else{
      bpfull = read.table(LDblock,header=TRUE,sep = "\t")
      bpp <- bpfull[str_trim(bpfull$chr)==paste0("chr",chromosome),]
      i = 1
      dchrom1initial$lab <- "NULL"
      for (i in 1:nrow(bpp)) {
        dchrom1initial[dchrom1initial$BP >= bpp$start[i] & dchrom1initial$BP < bpp$stop[i],]$lab <- bpp$start[i]
      }
    }
    whileivar = 0
    while (whileivar <= levelsdown) {
    gwasmulti <- dchrom1initial %>%
      group_by(lab) %>%
      filter(any(logp > sliceval)) %>%
      ungroup() %>%
      group_by(lab, PHENO) %>%
      slice(which.max(logp))
## addign for modifying gwasmulti logic
    if (calculateLD){
    BB <- dchrom1initial %>%
      group_by(lab) %>%
      filter(logp > sliceval)

    gwasmulti$snpsinld <- "NA"
    gwasmulti$snpsinldcount <- 0
    gwasmulti$snpsinldup <- "NA"
    z<-NULL

    for (i in 1:nrow(gwasmulti)) {

      stringcreate  <- sprintf("/1000GENOMES:phase_3:%s?",pop)
      server <- "http://grch37.rest.ensembl.org"
      ext <- paste0("/ld/human/",gwasmulti[i,]$SNP,stringcreate)
      r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
      tryCatch({
        stop_for_status(r)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      lddf <- fromJSON(toJSON(content(r)))
      ldlist <- unlist(lddf$variation2)
      if(!is.null(ldlist)){
      ldlist <- unlist(lddf[lddf$r2 >= R2 | lddf$d_prime >= D,]$variation2)
      hhh <- BB[BB$lab == gwasmulti[i,]$lab & BB$PHENO == gwasmulti[i,]$PHENO, ]
      snpsinld <- intersect(hhh$SNP,ldlist)
      snpsinldup <- length(snpsinld) / nrow(hhh)
      gwasmulti[i,]$snpsinld <- paste0(toString(snpsinld))
      gwasmulti[i,]$snpsinldcount <- length(snpsinld)
      gwasmulti[i,]$snpsinldup <- paste0(snpsinldup)
      }
    }
if(mutualLD){
  gwasmulti$common <- "NA"
    for(i in unique(gwasmulti$lab)){
      for (j in 1:length(phenos)) {
        for (k in phenos[-j]) {
        uuu <- strsplit(gwasmulti[gwasmulti$lab == i & gwasmulti$PHENO == phenos[j], ]$snpsinld,",")
        vfg <- strsplit(gwasmulti[gwasmulti$lab == i & gwasmulti$PHENO == k, ]$snpsinld,",")
        vfg <- lapply(vfg, str_trim)
        uuu <- lapply(uuu, str_trim)
        unld <- intersect(append(unlist(uuu),toString(gwasmulti[gwasmulti$lab == i & gwasmulti$PHENO == phenos[j], ]$SNP)), append(unlist(vfg), toString(gwasmulti[gwasmulti$lab == i & gwasmulti$PHENO == k, ]$SNP)))
        if(length(unld)== 0){
          unld <- "NA"
        }
        z <- c(z,paste(phenos[j],"with",k,":" ,paste0(toString(unld)), collapse = "," ))
        }
      gwasmulti[gwasmulti$lab == i & gwasmulti$PHENO == phenos[j],]$common <- paste(toString(z))
      z <- NULL
      }
    }
    }
}## the end of the calcualte LD loop, will apply only if the calculateLD is set to true
    ### the gwasmulti logic ends here so here we can put the loop above to go down the peaks
    if(levelsdown != 0 & calculateLD){
      dchrom1initialEXC <- dchrom1initial
    df_total = data.frame()
    for (ipheno in unique(gwasmulti$PHENO)) {
      bgt <- strsplit(gwasmulti[gwasmulti$PHENO == ipheno,]$snpsinld,",")
      xxxx <- lapply(bgt, str_trim)
      ppp <- union(unlist(xxxx),gwasmulti[gwasmulti$PHENO == ipheno,]$SNP)
      xxlist <- as.list(ppp)
      df <- dchrom1initial[dchrom1initial$PHENO == ipheno,][! dchrom1initial[dchrom1initial$PHENO == ipheno,]$SNP %in% xxlist , ]
      df_total <- rbind(df_total,df)
    }
    dchrom1initial <- df_total
    }
    whileivar = whileivar+1;
} ### end of the while is here
    if(nrow(gwasmulti) == 0){
      print("There are no SNP'S above this logp threshold, try decreasing the logp value")
    }
    if(length(grep("gene", colnames(d))) == 0){
      print("Applying BioMArt module for matching gene to rsid")
      gwasmulti <- addgene(gwasmulti)
    }
    if( length(unique(gwasmulti$lab)) == 1 & any(gwasmulti$logp>sliceval)){
      gwasmulticopybind1 <- gwasmulti
      gwasmulticopybind2 <- gwasmulti
      gwasmulticopybind1$lab  <- toString(as.numeric(unique(gwasmulti$lab)) + 1)
      gwasmulticopybind2$lab  <- toString(as.numeric(unique(gwasmulti$lab)) - 1)
      gwasmulticopybind1$logp  <- 0
      gwasmulticopybind2$logp  <- 0
      gwasmulticopybind1$gene <- as.character("NA")
      gwasmulticopybind2$gene <- as.character("NA")
      gwasmulti$gene <- as.character(gwasmulti$gene)
      gwasmulti <- rbind.data.frame(rbind.data.frame(gwasmulticopybind1,gwasmulticopybind2),gwasmulti)
      flag = TRUE
    }
    gwas_surface <- acast(gwasmulti, gwasmulti$PHENO ~ gwasmulti$lab, value.var = "logp")
    yy <- sort(as.numeric(colnames(gwas_surface)))
    gwas_surface <- gwas_surface[,as.character(yy)]
    gwas_surface_use <- gwas_surface
    if(is.null(upperlimit)){
      upperlimit = max(gwas_surface_use)
    }
    if( flag == FALSE ){
    ## take only the columns that have atleast one value which is greater than the slicevalue that is passed
    gwas_surface_use <- gwas_surface_use[,colSums(gwas_surface_use > sliceval) > 0]
    gwas_surface_use <- gwas_surface_use[,colSums(gwas_surface_use <  upperlimit +2) > nrow(gwas_surface_use)-1]
    }
    gwas_surface_prime_use <- gwas_surface_use
    gwas_surface_prime_use[,]=  upperlimit + 2
    gwas_surface_copy_use <- gwas_surface_use

    for (i in 1:nrow(gwas_surface_copy_use)) {
      for (j in 1:ncol(gwas_surface_copy_use))
        if(calculateLD){
          gwas_surface_copy_use[i, j] <-
            paste("Phenotype: ", rownames(gwas_surface_copy_use)[i], ' ',"kbp: ", colnames(gwas_surface_use)[j],' ',
                  "A1/A2: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$A1,'/',
                  gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$A2,
                  "Effect Size: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_use)[j],]$BETA,'\n',
                  "SE: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$SE,' ',
                  "-log10 (P-value): ", gwas_surface_copy_use[i, j],' ',
                  "SNPID: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$SNP,"\n",
                  gsub('(.{1,90})(\\s|$)', '\\1\n', paste("Gene: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$gene)),"\n",
                  gsub('(.{1,90})(\\s|$)', '\\1\n', paste("SNP's in LD: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$snpsinld)),"\n",
          "Linked SNP's ratio: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$snpsinldup,"\n",
gsub('(.{1,90})(\\s|$)', '\\1\n', paste("With other pheno: ",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$common)))

        }else{
        gwas_surface_copy_use[i, j] <-
          paste("Phenotype = ", rownames(gwas_surface_copy_use)[i], '\n',"kbp position range= ", colnames(gwas_surface_use)[j],'\n',
                "BP =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$BP,"\n",
                "A1 =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$A1,' ',
                "A2 =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$A2,'\n',
                "Effect Size =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_use)[j],]$BETA,' ',
                "SE =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$SE,'\n',
                "-log10 (P-value) = ", gwas_surface_copy_use[i, j],'\n',
                "SNPID =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$SNP,"\n",
                "Gene =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$gene)
        }
    }

    #### the new surfece for gene logic
    colorheatmap = NULL

if(geneview == TRUE){
  print("GENE View is active")
    gwas_surface_copy_gene_use <- gwas_surface_use
    for (i in 1:nrow(gwas_surface_copy_gene_use)) {
      for (j in 1:ncol(gwas_surface_copy_gene_use))
        gwas_surface_copy_gene_use[i, j] <-
          paste(gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_gene_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_gene_use)[j],]$gene)
    }

    xxx <- apply(gwas_surface_copy_gene_use, 2, function(x) unique(unlist(strsplit(x[!is.na(x)], ";"))[duplicated(unlist(strsplit(x[!is.na(x)], ";")))]))

    for (i in 1:nrow(gwas_surface_copy_gene_use)) {
      for (j in 1:ncol(gwas_surface_copy_gene_use)){
        x <- gwas_surface_copy_gene_use[i,j]
        xx <- unlist(strsplit(x[!is.na(x)], ";"))
        if(length( intersect(xx,unlist(xxx[j])) ) > 0 & gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_gene_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_gene_use)[j],]$logp > sliceval & !is.na(gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_gene_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_gene_use)[j],]$gene)){
          gwas_surface_copy_gene_use[i, j] <- max(gwas_surface_use)
        }else{
          gwas_surface_copy_gene_use[i, j] <- sliceval
        }
      }
    }
    surfacecolourheat = gwas_surface_copy_gene_use
    colorheatmap = "Viridis"
}else{
  surfacecolourheat = gwas_surface_use
}

    if(betaplot == TRUE){

      gwas_surface_beta <- gwas_surface_use
      for (i in 1:nrow(gwas_surface_copy_use)) {
        for (j in 1:ncol(gwas_surface_copy_use))
          gwas_surface_beta[i, j] <-
            gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_use)[j],]$BETA
      }

      gwas_surface_prime_use_beta <- gwas_surface_beta
      sliceval = min(gwas_surface_beta)
      upperlimit = max(gwas_surface_beta)
      gwas_surface_prime_use_beta[,]=  upperlimit
      gwas_surface_beta_use <- gwas_surface_beta

      for (i in 1:nrow(gwas_surface_beta_use)) {
        for (j in 1:ncol(gwas_surface_beta_use))
          gwas_surface_beta_use[i, j] <-
            paste("Phenotype = ", rownames(gwas_surface_beta_use)[i], '\n',"kbp position range= ", colnames(gwas_surface_beta)[j],'\n',
                  "A1 =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_beta_use)[i]& gwasmulti$lab==colnames(gwas_surface_beta_use)[j],]$A1,' ',
                  "A2 =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_beta_use)[i]& gwasmulti$lab==colnames(gwas_surface_beta_use)[j],]$A2,'\n',
                  "EffectSize =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_beta_use)[i]& gwasmulti$lab==colnames(gwas_surface_beta)[j],]$BETA,' ',
                  "SE =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_beta_use)[i]& gwasmulti$lab==colnames(gwas_surface_beta_use)[j],]$SE,'\n',
                  "-log10 (P-value) = ", gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_beta_use)[i]& gwasmulti$lab==colnames(gwas_surface_beta)[j],]$logp,'\n',
                  "SNPID =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_beta_use)[i]& gwasmulti$lab==colnames(gwas_surface_beta_use)[j],]$SNP,"\n",
                  "Gene =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_beta_use)[i]& gwasmulti$lab==colnames(gwas_surface_beta_use)[j],]$gene)
      }
      colorheatmap = NULL
      surfacecolourheat = gwas_surface_beta

      p <- plot_ly() %>%
        layout(title = 'PheGWAS',scene = list(camera = list(eye = list(x=2, y=2, z=2)),yaxis = list(title = "Phenotypes")
                                              ,xaxis =list(title= "100 kbp divisions" ,
                                                           type = 'category',tickmode = "array",tickvals = colnames(gwas_surface_beta), ticks = "outside"),
                                              zaxis =list(title= "Effect size",range=c(sliceval, upperlimit)),aspectmode = "manual",
                                              aspectratio = list( x = 2, y = .8, z = .8)))

      p <- p %>% add_surface(x=colnames(gwas_surface_beta),y= rownames(gwas_surface_beta),z = gwas_surface_beta,opacity=.98,
                             surfacecolor = gwas_surface_beta, cmin = sliceval ,cmax = upperlimit,
                             showscale = FALSE,text = gwas_surface_beta_use,hoverinfo = "text") %>%
        add_trace(x=colnames(gwas_surface_beta),y= rownames(gwas_surface_beta),z = gwas_surface_prime_use_beta,
                  type = "surface",surfacecolor = surfacecolourheat,showscale = FALSE,text = gwas_surface_beta_use,hoverinfo = "text",colorscale = colorheatmap)
      p

    } else{

    p <- plot_ly() %>%
      layout(title = 'PheGWAS',scene = list(camera = list(eye = list(x=2, y=2, z=2)),yaxis = list(title = "Phenotypes")
                                            ,xaxis =list(title= "100 kbp divisions" ,
                                                         type = 'category',tickmode = "array",tickvals = colnames(gwas_surface_use), ticks = "outside"),
                                            zaxis =list(title= "-log 10(P)",range=c(sliceval,upperlimit + 2)),aspectmode = "manual",
                                            aspectratio = list( x = 2, y = .8, z = .8)))


    p <- p %>% add_surface(x=colnames(gwas_surface_use),y= rownames(gwas_surface_use),z = gwas_surface_use,opacity=.98,
                           surfacecolor = gwas_surface_use, cmin = sliceval ,cmax = upperlimit,
                           showscale = FALSE,text = gwas_surface_copy_use,hoverinfo = "text") %>%
      add_trace(x=colnames(gwas_surface_use),y= rownames(gwas_surface_use),z = gwas_surface_prime_use,
                type = "surface",surfacecolor = surfacecolourheat,showscale = FALSE,text = gwas_surface_copy_use,hoverinfo = "text",colorscale = colorheatmap)
    p
  }}
  p
}


#' 3-D association table for many phenotypes
#'
#' @import tidyverse
#' @import tidyr
#' @import plotly
#' @importFrom biomaRt useMart getBM
#' @import reshape2
#' @param d DataFrame output from processphegwas
#' @param sliceval Integer to indicate value of -log10(p) to do the sectionalcut. Usually value > -log10 6 is considered to be significant
#' @param chromosome Integer to indicate the chromosome number thats interested, If not given entire chromosome is given
#' @param bpdivision Integer value to indicate the base pair divisions. By default it is 100,000 bp. This vis only applciable on a single chromosome view. For entire
#' chromosoem view the max peak is selected
#' @author George Gittu
#' @examples
#' \dontrun{
#' # table of chromosome number 16, with bp division on 1,000,000 (default is 1,000,000)
#' x <- list(hdl,ldl,trig,tchol)
#' phenos <- c("HDL","LDL","TRIGS","TOTALCHOLESTROL")
#' y <- processphegwas(x, phenos)
#' landscapetable(y,sliceval = 6,chromosome = 16,bpdivision = 1000000)
#' }
#' @export
landscapetable <- function(d, sliceval = 8, chromosome=1, bpdivision= 1000000) {
  paste0("Processing for the chromosome number ", chromosome)
  ## Pick the chromosome that is passed as param "chromosome"
  dchrom1initial <- d[d$CHR==chromosome,]
  ## to divide based on a single base pair
  dchrom1initial$lab <- as.integer(dchrom1initial$BP /bpdivision, 0)
  gwasmulti <- dchrom1initial %>%
    group_by(lab) %>%
    filter(any(logp > sliceval)) %>%
    ungroup() %>%
    group_by(lab, PHENO) %>%
    slice(which.max(logp))
if(nrow(gwasmulti) == 0){
  print("There are no SNP'S above this logp threshold, try decreasing the logp value")
}
  if(length(grep("gene", colnames(d))) == 0){
    print("Applying BioMArt module for matching gene to rsid")
    gwasmulti <- addgene(gwasmulti)
  }
  table <- gwasmulti

  CCCC <- table %>%
    group_by(lab) %>%
    filter(logp > sliceval) %>%
    summarise(CHR = paste(unique(CHR),collapse = " ") ,MarkerName=paste(SNP, collapse=" "),GX=paste(gene, collapse = ";"),AssociatedTraits=paste(PHENO, collapse=" "),P_value = paste(P,collapse = " "),logP_value= paste(round(logp,2),collapse=" "))

  vvv <- function(x) {
    if(length(unlist(strsplit(x[2][!is.na(x[2])],","))) >1){
      paste(unique(unlist(strsplit(x[1][!is.na(x[1])], ";"))[duplicated(unlist(strsplit(x[1][!is.na(x[1])], ";")))]),collapse=' ')
    }else{
      paste(unique(unlist(strsplit(x[1][!is.na(x[1])], ";"))), collapse = " ")
    }
  }


  CCCC$Genes <- apply(CCCC[,4:5], 1, vvv)
  drops <- c("GX")
  CCCC <- CCCC[ , !(names(CCCC) %in% drops)]

  values <- rbind(CCCC$CHR,CCCC$lab, CCCC$MarkerName,CCCC$AssociatedTraits,CCCC$Genes,CCCC$P_value,CCCC$logP_value)

  p <- plot_ly(
    type = 'table',
    columnorder = c(1,2,3,4,5,6,7),
    columnwidth = c(7,15,60,60,100,55,45),
    header = list(
      values = c('<b>Chr</b>','<b>Position group (Mb)</b>', '<b>Marker Name</b>','<b>Associated Traits</b>','<b>Genes</b>','<b>P-Values</b>','<b>-log10 (P-value)</b>'),
      line = list(color = '#506784'),
      fill = list(color = 'blue'),
      align = c('left','center'),
      font = list(color = 'white', size = 12),
      height = 40
    ),
    cells = list(
      values = values,
      line = list(color = '#506784'),
      fill = list(color = c( 'lightblue')),
      align = c('left', 'center'),
      font = list(color = c('black'), size = 12),
      height = 30
    ))
  p
}

#' 3-D association table for many phenotypes(entire genome)
#'
#' @import tidyverse
#' @import tidyr
#' @import plotly
#' @importFrom biomaRt useMart getBM
#' @import reshape2
#' @param d DataFrame output from processphegwas
#' @param sliceval Integer to indicate value of -log10(p) to do the sectionalcut. Usually value > -log10 6 is considered to be significant
#' @param bpdivision Integer value to indicate the base pair divisions. By default it is 100,000 bp. This vis only applciable on a single chromosome view. For entire
#' chromosoem view the max peak is selected
#' @author George Gittu
#' @examples
#' \dontrun{
#' # table of chromosome number 16, with bp division on 1,000,000 (default is 1,000,000)
#' x <- list(hdl,ldl,trig,tchol)
#' phenos <- c("HDL","LDL","TRIGS","TOTALCHOLESTROL")
#' y <- processphegwas(x, phenos)
#' landscapetable(y,sliceval = 6,chromosome = 16,bpdivision = 1000000)
#' }
#' @export
landscapefulltable <- function (d, sliceval = 6, bpdivision = 1e+06)
{
  values = list()
  for(chromosome in 1:22 )
  {
    paste0("Processing for the chromosome number ", chromosome)
    dchrom1initial <- d[d$CHR == chromosome, ]
    dchrom1initial$lab <- as.integer(dchrom1initial$BP/bpdivision,
                                     0)
    gwasmulti <- dchrom1initial %>% group_by(lab) %>% filter(any(logp >
                                                                   sliceval)) %>% ungroup() %>% group_by(lab, PHENO) %>%
      slice(which.max(logp))
    if (nrow(gwasmulti) == 0) {
      print("There are no SNP'S above this logp threshold, try decreasing the logp value")
    }
    if (length(grep("gene", colnames(d))) == 0) {
      print("Applying BioMArt module for matching gene to rsid")
      gwasmulti <- addgene(gwasmulti)
    }
    table <- gwasmulti
    CCCC <- table %>% group_by(lab) %>% filter(logp > sliceval) %>%
      summarise(CHR = paste(unique(CHR), collapse = " "), MarkerName = paste(SNP,
                                                                             collapse = " "), GX = paste(gene, collapse = ";"),
                AssociatedTraits = paste(PHENO, collapse = " "),
                P_value = paste(P, collapse = " "), logP_value = paste(round(logp,
                                                                             2), collapse = " "))
    vvv <- function(x) {
      if (length(unlist(strsplit(x[2][!is.na(x[2])], ","))) >
          1) {
        paste(unique(unlist(strsplit(x[1][!is.na(x[1])],
                                     ";"))[duplicated(unlist(strsplit(x[1][!is.na(x[1])],
                                                                      ";")))]), collapse = " ")
      }
      else {
        paste(unique(unlist(strsplit(x[1][!is.na(x[1])],
                                     ";"))), collapse = " ")
      }
    }
    if(nrow(CCCC) != 0) {
      CCCC$Genes <- apply(CCCC[, 4:5], 1, vvv)
      drops <- c("GX")
      CCCC <- CCCC[, !(names(CCCC) %in% drops)]
      values[[chromosome]] <- rbind(CCCC$CHR, CCCC$lab, CCCC$MarkerName, CCCC$AssociatedTraits,
                                    CCCC$Genes, CCCC$P_value, CCCC$logP_value)
    }
  }
  big_data = do.call(cbind, values)
  p <- plot_ly(type = "table", columnorder = c(1, 2, 3, 4,
                                               5, 6, 7), columnwidth = c(7, 15, 60, 60, 100, 55, 45),
               header = list(values = c("<b>Chr</b>", "<b>Position group (Mb)</b>",
                                        "<b>Marker Name</b>", "<b>Associated Traits</b>",
                                        "<b>Genes</b>", "<b>P-Values</b>", "<b>-log10 (P-value)</b>"),
                             line = list(color = "#506784"), fill = list(color = "blue"),
                             align = c("left", "center"), font = list(color = "white",
                                                                      size = 12), height = 40), cells = list(values = big_data,
                                                                                                             line = list(color = "#506784"), fill = list(color = c("lightblue")),
                                                                                                             align = c("left", "center"), font = list(color = c("black"),
                                                                                                                                                      size = 12), height = 30))
  p
}

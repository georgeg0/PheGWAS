
# once added the roxygen syntax things then run following
# devtools::document()
BP <- CHR <- P <- PHENO <- SNP <- gene <- lab <- logp <- NULL

## Function to add gene for repective rsid, if genes are not provided by user
addgene <- function(gwasmulti){

  ensembl <- useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")

  human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

  gwasmulti$gene <- NA
  for(i in rownames(gwasmulti)){
    xx <- getBM(attributes=c( "ensembl_gene_stable_id"),filters="snp_filter", values=gwasmulti[i, "SNP"],mart=ensembl, uniqueRows=TRUE)
    if (nrow(xx) == 0){
      gwasmulti[i,]$gene <- NA
    }else{
      gene <- getBM(attributes = c("hgnc_symbol"),
                    filters = "ensembl_gene_id", values = unlist(xx), mart = human)
      gwasmulti[i,]$gene <- if (nrow(gene) == 0) NA else  toString( unlist(gene))
    }}
  print("Finished mapping")
  gwasmulti
}

#' Prepare the dataframe to pass to landscape function
#'
#' @import tidyverse
#' @import tidyr
#' @param x List of dataframes that need to do PheGWAS on. Arrange the dataframe in the order how the the phenotypes should align in y axis
#' @return A processed dataframe to pass to PheGWAS function
#' @details Make sure there are no duplicate rsid's in any of the dataframe, If there aremake sure to resolve it before passing it to this function.
#' @author George Gittu
#' @examples
#' x <- list(bmimen,bmiwomen)
#'
#' ## y is ready to be passed to function landscape
#' y <- processphegwas(x)
#' @export
processphegwas <- function(x) {

  fulldf <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("CHR","BP","rsid"), all = TRUE),x)

  if(length(grep("gene", colnames(fulldf))) > 0){
    ## means it got gene added by the user
    print("Processing for gene provided by the user")
    fulldf$G <- apply(fulldf[, grep("gene", colnames(fulldf))], 1, function(x) paste(unique(unlist(strsplit(x[!is.na(x)], ";"))), collapse = ";"))
    fulldfff <- fulldf[, -grep("gene", colnames(fulldf))]
    gwasmulti.melt <- reshape2::melt(data = fulldfff, id.vars = c("CHR","BP","rsid","G"), variable.name = "color", value.name = "p_value")
    d <- data.frame(CHR = gwasmulti.melt$CHR, BP = gwasmulti.melt$BP, P = gwasmulti.melt$p_value,SNP = gwasmulti.melt$rsid,gene = gwasmulti.melt$G,PHENO = gwasmulti.melt$color)
  }else{
    ## Gene when not addded
    print("Gene will be added as part of BioMart module")
  gwasmulti.melt <- reshape2::melt(data = fulldf, id.vars = c("CHR","BP","rsid"), variable.name = "color", value.name = "p_value")
  d <- data.frame(CHR = gwasmulti.melt$CHR, BP = gwasmulti.melt$BP, P = gwasmulti.melt$p_value,SNP = gwasmulti.melt$rsid,PHENO = gwasmulti.melt$color)
  }
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  d$logp <- -log10(d$P)
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
#' @param d DataFrame output from processphegwas
#' @param sliceval Integer to indicate value of -log10(p) to do the sectionalcut. Usually value > -log10 6 is considered to be significant
#' @param chromosome Integer to indicate the chromosome number thats interested, If not given entire chromosome is given
#' @param bpdivision Integer value to indicate the base pair divisions. By default it is 100,000 bp. This vis only applciable on a single chromosome view. For entire
#' @param upperlimit Integer to indicate value of -log10(p) to do the upperlimit. If not given it will take the max -log10(p)
#' @param geneview This checks for the common genes across the section
#' chromosoem view the max peak is selected
#' @author George Gittuuuuu
#' @examples
#' x <- list(hdl,ldl,trig,tchol)
#' y <- processphegwas(x)
#' ## pass the dataframe from the processphegwas
#'
#' # 3D landscape visualization of all the phenotypes across the base pair positions
#' landscape(y)
#'
#' # 3D landscape visualization of all the phenotypes across the base pair positions, showing only the crusts above a certain threshold
#' landscape(y,sliceval = 6)
#'
#' # 3D landscape visualization of chromosome number 16
#' landscape(y,sliceval = 6,chromosome = 16)
#'
#' # 3D landscape visualization of chromosome number 16, with bp division on 1,000,000 (default is 1,000,000)
#' landscape(y,sliceval = 6,chromosome = 16,bpdivision = 1000000)
#' @export
landscape <- function(d, sliceval = 0, chromosome=FALSE, bpdivision= 1000000, upperlimit = NULL,geneview = FALSE) {

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
          paste("Phenotype = ", rownames(gwas_surface_copy_full_use)[i], '\n',"Chromosome = ", colnames(gwas_surface_copy_full_use)[j],'\n',"P value = ", gwas_surface_copy_full_use[i, j],'\n',
                "SNPID =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$SNP,'\n',
                "Gene =",gwasmultifull[gwasmultifull$PHENO==rownames(gwas_surface_copy_full_use)[i]& gwasmultifull$CHR==colnames(gwas_surface_copy_full_use)[j],]$gene)
    }

    p <- plot_ly(x=colnames(gwas_surface_full_use),y= rownames(gwas_surface_full_use),z = gwas_surface_full_use) %>%
      add_surface(cmin = sliceval,cmax = upperlimit,opacity=.95,surfacecolor = gwas_surface_full_use,text = gwas_surface_copy_full_use,hoverinfo = "text") %>%
      layout(title = 'FIGIWAS',scene = list(yaxis = list(title = "Phenotypes", tickmode = "array",nticks = 8),xaxis =list(title= "Chromosomes",autotick = F, dtick = 1),
                                              zaxis =list(title= "-log 10(P)",range=c(sliceval,upperlimit+2)),aspectmode = "manual",aspectratio = list( x = 1.8, y = 1, z = 1))) %>%
      add_trace(x=colnames(gwas_surface_prime_full_use),y= rownames(gwas_surface_prime_full_use),z = gwas_surface_prime_full_use, type = "surface",surfacecolor = gwas_surface_full_use,showscale = FALSE,
                text = gwas_surface_copy_full_use,hoverinfo = "text",cmin = sliceval,cmax = upperlimit) %>% colorbar(title = "P Value")
    p

  }
  else {
    paste0("Processing for the chromosome number ", chromosome)
    ## Pick the chromosome that is passed as param "chromosome"
    dchrom1initial <- d[d$CHR==chromosome,]
    ## to divide based on a single base pair
    # maxp = tail(sort(dchrom1initial$logp),3)[3]
    # maxt = tail(sort(dchrom1initial$logp),3)[2]
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
    gwas_surface <- acast(gwasmulti, gwasmulti$PHENO ~ gwasmulti$lab, value.var = "logp")
    gwas_surface_use <- gwas_surface
    if(is.null(upperlimit)){
      upperlimit = max(gwas_surface_use)
      print(upperlimit)
    }
    ## take only the columns that have atleast one value which is greater than the slicevalue that is passed
    gwas_surface_use <- gwas_surface_use[,colSums(gwas_surface_use > sliceval) > 0]
    gwas_surface_use <- gwas_surface_use[,colSums(gwas_surface_use <  upperlimit +2) > nrow(gwas_surface_use)-1]
    gwas_surface_prime_use <- gwas_surface_use

    gwas_surface_prime_use[,]=  upperlimit + 2
    gwas_surface_copy_use <- gwas_surface_use

    for (i in 1:nrow(gwas_surface_copy_use)) {
      for (j in 1:ncol(gwas_surface_copy_use))
        gwas_surface_copy_use[i, j] <-
          paste("Phenotype = ", rownames(gwas_surface_copy_use)[i], '\n',"kbp position range= ", colnames(gwas_surface_copy_use)[j],'\n',"P value = ", gwas_surface_copy_use[i, j],'\n',"SNPID =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$SNP,"\n","Gene =",gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_use)[j],]$gene)
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
        # print(paste("row =",i))
        # print(paste("column =",j))
        # unlistedstuff <- paste(unlist(xxx[j]))
        # print(unlistedstuff)

        if(length( intersect(xx,unlist(xxx[j])) ) > 0 & gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_gene_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_gene_use)[j],]$logp > sliceval & !is.na(gwasmulti[gwasmulti$PHENO==rownames(gwas_surface_copy_gene_use)[i]& gwasmulti$lab==colnames(gwas_surface_copy_gene_use)[j],]$gene)){
          # print(length( intersect(xx,unlist(xxx[j])) ))
          # print("Condition satisfied")
          gwas_surface_copy_gene_use[i, j] <- max(gwas_surface_use)
        }else{
          # print(length( intersect(xx,unlist(xxx[j])) ))
          # print("Condition NOT satisfied")
          gwas_surface_copy_gene_use[i, j] <- sliceval

        }

      }
    }
    surfacecolourheat = gwas_surface_copy_gene_use
    colorheatmap = "Viridis"
}else{
  surfacecolourheat = gwas_surface_use
}

    p <- plot_ly() %>%
      layout(title = 'FIGIWAS',scene = list(yaxis = list(title = "Phenotypes")
                                            ,xaxis =list(title= "100 kbp divisions" ,
                                                         type = 'category',tickmode = "array",tickvals = colnames(gwas_surface_use), ticks = "outside"),
                                            zaxis =list(title= "-log 10(P)",autotick = F, dtick = 5,range=c(sliceval,upperlimit + 2)),aspectmode = "manual",
                                            aspectratio = list( x = 2.5, y = 1, z = 1)))

    p <- p %>% add_surface(x=colnames(gwas_surface_use),y= rownames(gwas_surface_use),z = gwas_surface_use,opacity=.98,
                           surfacecolor = gwas_surface_use, cmin = sliceval ,cmax = upperlimit,
                           showscale = FALSE,text = gwas_surface_copy_use,hoverinfo = "text") %>%
      add_trace(x=colnames(gwas_surface_use),y= rownames(gwas_surface_use),z = gwas_surface_prime_use,
                type = "surface",surfacecolor = surfacecolourheat,showscale = FALSE,text = gwas_surface_copy_use,hoverinfo = "text",colorscale = colorheatmap)
    p
  }
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
#' # table of chromosome number 16, with bp division on 1,000,000 (default is 1,000,000)
#' x <- list(hdl,ldl,trig,tchol)
#' y <- processphegwas(x)
#' landscapetable(y,sliceval = 6,chromosome = 16,bpdivision = 1000000)
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

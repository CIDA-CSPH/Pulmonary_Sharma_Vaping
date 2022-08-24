library(sesame)

SigSetToRGChannel <- function(sset, manifest = NULL, controls = NULL) {
  
  if (is.null(manifest)) {
    dfAddress <- sesameDataGet(paste0(sset@platform,'.address'))
    manifest <- dfAddress$ordering
    controls <- dfAddress$controls
  }
  
  SSRed <- NULL
  SSGrn <- NULL
  
  IIdf <- manifest[
    manifest$COLOR_CHANNEL=='Both', c('Probe_ID','U')]
  SSRed <- c(SSRed, setNames(sset@II[match(
    IIdf$Probe_ID, rownames(sset@II)), 'U'],
    as.character(IIdf$U)))
  SSGrn <- c(SSGrn, setNames(sset@II[match(
    IIdf$Probe_ID, rownames(sset@II)), 'M'],
    as.character(IIdf$U)))
  
  IRdf <- manifest[
    manifest$COLOR_CHANNEL=='Red', c('Probe_ID','M','U')]
  SSRed <- c(SSRed, setNames(sset@IR[match(
    IRdf$Probe_ID, rownames(sset@IR)), 'M'],
    as.character(IRdf$M)))
  SSRed <- c(SSRed, setNames(sset@IR[match(
    IRdf$Probe_ID, rownames(sset@IR)), 'U'],
    as.character(IRdf$U)))
  ## OOB signals
  SSGrn <- c(SSGrn, setNames(sset@oobG[match(
    IRdf$Probe_ID, rownames(sset@oobG)), 'M'],
    as.character(IRdf$M)))
  SSGrn <- c(SSGrn, setNames(sset@oobG[match(
    IRdf$Probe_ID, rownames(sset@oobG)), 'U'],
    as.character(IRdf$U)))
  
  IGdf <- manifest[
    manifest$COLOR_CHANNEL=='Grn', c('Probe_ID','M','U')]
  SSGrn <- c(SSGrn, setNames(sset@IG[match(
    IGdf$Probe_ID, rownames(sset@IG)), 'M'], 
    as.character(IGdf$M)))
  SSGrn <- c(SSGrn, setNames(sset@IG[match(
    IGdf$Probe_ID, rownames(sset@IG)), 'U'], 
    as.character(IGdf$U)))
  ## OOB signals
  SSRed <- c(SSRed, setNames(sset@oobR[match(
    IGdf$Probe_ID, rownames(sset@oobR)), 'M'],
    as.character(IGdf$M)))
  SSRed <- c(SSRed, setNames(sset@oobR[match(
    IGdf$Probe_ID, rownames(sset@oobR)), 'U'],
    as.character(IGdf$U)))
  
  ## controls
  if (!is.null(controls)) {
    control.names <- make.names(controls$Name, unique = TRUE)
    SSGrn <- c(SSGrn, setNames(sset@ctl[match(
      control.names, rownames(sset@ctl)),'G'], 
      as.character(controls$Address)))
    SSRed <- c(SSRed, setNames(sset@ctl[match(
      control.names, rownames(sset@ctl)),'R'], 
      as.character(controls$Address)))
  } ## else TODO controls obtained from manifest
  
  list(grn=SSGrn, red=SSRed)
}

## annotation, if not given is guessed
guessMinfiAnnotation <- function(platform, annotation = NA) {
  if (is.na(annotation)) {
    if (platform %in% c("HM450", "HM27")) {
      'ilmn12.hg19'
    } else { # EPIC
      'ilm10b4.hg19'
    }
  } else {
    annotation
  }
}

#' Convert sesame::SigSet to minfi::RGChannelSet
#' 
#' @param ssets a list of sesame::SigSet
#' @param BPPARAM get parallel with MulticoreParam(n)
#' @param annotation the minfi annotation string, guessed if not given
#' @return a minfi::RGChannelSet
#' @import BiocParallel
#' @examples
#'
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' rgSet <- SigSetsToRGChannelSet(sset)
#'
#' @export 
SigSetsToRGChannelSet <- function(ssets, BPPARAM=SerialParam(), annotation=NA) {
  if (is(ssets, 'SigSet')) {
    ssets <- list(sample=ssets)
  }
  
  platform <- ssets[[1]]@platform
  annotation <- guessMinfiAnnotation(annotation)
  
  ss_all <- bplapply(ssets, SigSetToRGChannel, BPPARAM=BPPARAM)
  rgset <- minfi::RGChannelSet(
    Green=do.call(cbind, lapply(ss_all, function(ss) ss$grn)), 
    Red=do.call(cbind, lapply(ss_all, function(ss) ss$red)), 
    annotation=c(
      array=unname(platformSmToMinfi(platform)),
      annotation=annotation))
}

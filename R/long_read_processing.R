# functions related to long-read processing

#' Trim long reads to single nucleotide in their 3'ends
#'
#' @param GenomicRanges object from GenomicAlignments::readGAlignments() bam file
#'
#' @return trimmed genomic ranges object
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors
#' @examples
trim_single_nucleotide <- function(reads.gr) {
  pos <- reads.gr[ BiocGenerics::strand(reads.gr)== "+",]
  start( pos) <-  end(pos)
  neg <- reads.gr[BiocGenerics::strand(reads.gr)== "-",]
  end( neg) <-  start(neg)
  trim.gr <- c(pos,neg )
  return(trim.gr)
}

#' Create 3'end clusters
#'
#' @param nucleotide_regions
#' @param refAnnot
#' @param window
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors
#' @examples
create_clusters <- function(nucleotide_regions, refAnnot, window) {
  dataset_clusters <- nucleotide_regions %>%
    reduce(., min.gapwidth = window)
  # overlap to ensembl
  geneHitsExpanded <- IRanges::findOverlaps(dataset_clusters, refAnnot, maxgap = 2000)
  dataset_unassigned_subset <- dataset_clusters[ S4Vectors::queryHits(geneHitsExpanded)]
  # retrieve features from assigned gene
  dataset_unassigned_subset$gene_id <- refAnnot[ S4Vectors::subjectHits(geneHitsExpanded)]$gene_id
  dataset_unassigned_subset$gene_name <- refAnnot[ S4Vectors::subjectHits(geneHitsExpanded)]$gene_name
  dataset_unassigned_subset$gene_biotype <- refAnnot[ S4Vectors::subjectHits(geneHitsExpanded)]$gene_biotype
  # add gene positions
  dataset_unassigned_subset$geneEnd <- end(refAnnot[ S4Vectors::subjectHits(geneHitsExpanded)])
  dataset_unassigned_subset$geneStart <- start(refAnnot[ S4Vectors::subjectHits(geneHitsExpanded)])
  dataset_unassigned_subset$geneStrand <- as.character(BiocGenerics::strand(refAnnot[subjectHits(geneHitsExpanded)]))
  dataset_unassigned_subset$geneWidth <- width(refAnnot[ S4Vectors::subjectHits(geneHitsExpanded)])

  # calculate fraction of overlap
  dataset_unassigned_subset$fraction_overlap <- ifelse(dataset_unassigned_subset$geneStrand == "+",
                                                       c(c(end(dataset_unassigned_subset) - dataset_unassigned_subset$geneStart) / dataset_unassigned_subset$geneWidth),
                                                       c(c(start(dataset_unassigned_subset) - dataset_unassigned_subset$geneEnd) / dataset_unassigned_subset$geneWidth))
  dataset_unassigned_subset$fraction_overlap <- 1-abs(dataset_unassigned_subset$fraction_overlap)
  # sort by closest overlap to gene
  dataset_unassigned_subset <- dataset_unassigned_subset[order(dataset_unassigned_subset$fraction_overlap, decreasing = FALSE), ]
  # cluster_names
  dataset_unassigned_subset$cluster_name <- paste0(as.data.frame( dataset_unassigned_subset )$seqnames,
                                                   ":",
                                                   as.data.frame( dataset_unassigned_subset )$start,
                                                   ":",
                                                   as.data.frame( dataset_unassigned_subset )$end,
                                                   ":",
                                                   as.data.frame( dataset_unassigned_subset )$strand)
  # keep the assigment to most similar gene
  annotated_clusters  <- dataset_unassigned_subset[!duplicated(dataset_unassigned_subset$cluster_name),]
  #annotated_clusters <- dataset_unassigned_subset
  mcols(annotated_clusters) <- mcols(annotated_clusters)[c("gene_id", "gene_name", "gene_biotype") ]
  ## create cluster id with positional indexing
  # indexing for positive strand
  annotated_clusters.pos <- annotated_clusters[annotated_clusters@strand=="+",]
  # sort
  annotated_clusters.pos <- sort(annotated_clusters.pos, decreasing = FALSE)
  annotated_clusters.pos <- annotated_clusters.pos %>%
    as.data.frame(.) %>%
    group_by(gene_id) %>%
    mutate(cluster_name = paste0(gene_name,":","UTR",row_number(gene_name)),
           cluster_id = paste0(gene_id,":","UTR",row_number(gene_id))) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  # negative strand
  annotated_clusters.neg <- annotated_clusters[annotated_clusters@strand=="-",]
  annotated_clusters.neg <- sort(annotated_clusters.neg, decreasing = TRUE)
  annotated_clusters.neg <- annotated_clusters.neg %>%
    as.data.frame(.) %>%
    group_by(gene_id) %>%
    mutate(cluster_name = paste0(gene_name,":","UTR",row_number(gene_name)),
           cluster_id = paste0(gene_id,":","UTR",row_number(gene_id))) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  annotated_clusters <- c(annotated_clusters.pos,annotated_clusters.neg)
  return(annotated_clusters)
}


#' Extract features from clusters
#'
#' @param clusters3seq
#' @param genomeSeq
#' @param dwindow
#' @param upwindow
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors Biostrings
#' @examples
extract_features <- function(clusters3seq, genomeSeq, dwindow, upwindow){
  seqlevelsStyle(clusters3seq) <- "UCSC"
  clusters3seq <- keepStandardChromosomes(clusters3seq, pruning.mode = "coarse")
  # extend 3'coordinates for feature extraction
  message("Extending clusters in flanking windows")
  downstreamSite <-
    extend(clusters3seq, upstream = 0, downstream = dwindow)
  upstreamSite <-
    extend(clusters3seq, upstream = upwindow, downstream = 0)
  # make sequences
  seqsUpstream <- Biostrings::getSeq(genomeSeq, upstreamSite)
  names(seqsUpstream) <- upstreamSite$cluster_name
  seqsDownstream <- Biostrings::getSeq(genomeSeq, downstreamSite)
  names(seqsDownstream) <- downstreamSite$cluster_name
  message("retrieve nucleotide frequencies")
  seqsUpstream.dt <-
    as.data.frame(alphabetFrequency(seqsUpstream, baseOnly = TRUE, as.prob =
                                      TRUE)) %>% dplyr::select(-other)
  colnames(seqsUpstream.dt) <-
    paste0("nt.prob.upstream", colnames(seqsUpstream.dt))
  seqsUpstream.dt$cluster_id <- names(seqsUpstream)
  seqsDownstream.dt <-
    as.data.frame(alphabetFrequency(seqsDownstream, baseOnly = TRUE, as.prob =
                                      TRUE))  %>% dplyr::select(-other)
  colnames(seqsDownstream.dt) <-
    paste0("nt.prob.downstream", colnames(seqsDownstream.dt))
  seqsDownstream.dt$cluster_id <- names(seqsDownstream)
  clusterAContent <-
    left_join(seqsUpstream.dt, seqsDownstream.dt, by = "cluster_id")
  ## annotate poly(A) signals in upstream
  message("Retrieve poly(A) signals")
  message("Retrieve poly(A) signals Upstream")
  signalsUpStream <-
    getPolyAsignal(APAmotifs, hexamers , seqsUpstream)
  signalsUpStream <- lapply(signalsUpStream, as.data.frame)
  signalsUpStream <-
    bind_rows(signalsUpStream, .id = "column_label")
  message("Annotate poly(A) signals to clusters")
  signalsUpStream <-
    left_join(clusterAContent,
              signalsUpStream %>%
                dplyr::rename(cluster_id = names)) %>%
    dplyr::rename(polyA.Signal = column_label) %>%
    mutate(
      polyA.Signal = ifelse(is.na(polyA.Signal), "not detected", polyA.Signal),
      signalDistance = end - 50
    ) %>% dplyr::select(-c(width, start, end))
  ## remove duplicated clusters with multiple assignments keeping signal closer to cleavage site
  signalsUpStream <- signalsUpStream %>%
    arrange(abs(signalDistance), polyA.Signal) %>%
    distinct(cluster_id, .keep_all = TRUE)
  # retrieve counts
  signalsUpStream <- left_join(
    signalsUpStream ,
    mcols(clusters3seq) %>% as.data.frame(.) %>%
      dplyr::select(Counts, cluster_name, contains("gene")) %>%
      dplyr::rename(cluster_id = cluster_name),
    by = "cluster_id"
  )
  # add oligo kmers
  message("Search for kmers")
  oligoDT.kmers <-
    lapply(c(10, 15, 20), function (x) {
      strrep("A", x)
    })
  kmerClustersUpstream <-
    lapply(oligoDT.kmers, kmer_scan, seqsUpstream)
  kmerClustersDownstream <-
    lapply(oligoDT.kmers, kmer_scan, seqsDownstream)
  message("Prepare table")
  signalsUpStream <-
    signalsUpStream %>% mutate(
      Ant.8mer.upstream = ifelse(cluster_id %in% kmerClustersUpstream[1], TRUE, FALSE),
      Ant.10mer.upstream = ifelse(cluster_id %in% kmerClustersUpstream[2], TRUE, FALSE),
      Ant.10mer.upstream = ifelse(cluster_id %in% kmerClustersUpstream[3], TRUE, FALSE),
      Ant.8mer.downstream = ifelse(cluster_id %in% kmerClustersDownstream[1], TRUE, FALSE),
      Ant.10mer.downstream = ifelse(cluster_id %in% kmerClustersDownstream[2], TRUE, FALSE),
      Ant.10mer.downstream = ifelse(cluster_id %in% kmerClustersDownstream[3], TRUE, FALSE)
    ) %>%
    dplyr::select(c(
      cluster_id,
      Counts,
      gene_id,
      gene_name,
      gene_biotype,
      polyA.Signal,
      1 ,
      2,
      3,
      4,
      6,
      7,
      8,
      9,
      11,
      16,
      17,
      18,
      19
    )) %>%
    mutate(polyAsignalClass =
             ifelse(polyA.Signal %in% c("AATATA", "ATTAAA", "AATAAA"),
                    polyA.Signal,
                    ifelse(polyA.Signal == "not_detected",
                           "not_detected",
                           "non_canonical")))
  return(signalsUpStream)
}




APAmotifs <-
  c(
    "AATAAA",
    "ATTAAA",
    "AATATA",
    "AAGAAA",
    "AATACA",
    "AATAGA",
    "AATGAA",
    "ACTAAA",
    "CATAAA",
    "GATAAA",
    "TATAAA",
    "TTTAAA"
  )
hexamers <- list()

#' Title
#'
#' @param x
#' @param upstream
#' @param downstream
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors Biostrings
#' @examples
extend <- function(x, upstream=0, downstream=0)   {
  if (any(BiocGenerics::strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- BiocGenerics::strand(x) == "+" | BiocGenerics::strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

#' Title
#'
#' @param APAmotifs
#' @param hexamers
#' @param seqs
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors Biostrings
#' @examples
getPolyAsignal <- function(APAmotifs, hexamers, seqs){
  if(rlang::is_empty(APAmotifs)){ # no motifs remaining
    return(hexamers)
  } else {
    apaMotif = APAmotifs[1] #get first entry in list
    dnaseq <- DNAString(apaMotif) #canonical 1
    hexamers[[apaMotif]] <-unlist(vmatchPattern(dnaseq, seqs))
    seqs <- seqs[!names(seqs)  %in% names( hexamers[[apaMotif]] ),]
    APAmotifs <- APAmotifs[!APAmotifs %in%  apaMotif] # remove motif from list
    getPolyAsignal(APAmotifs, hexamers, seqs)
  }
}

#' Title
#'
#' @param clusters3seq
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors Biostrings
#' @examples
prepare_clusters <- function(clusters3seq){
  clusters3seq <-
    clusters3seq %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  end(clusters3seq) <- start(clusters3seq)
  seqlevelsStyle(clusters3seq) <- "UCSC"
  clusters3seq$cluster_name <-
    paste0(seqnames(clusters3seq),
           ":",
           start(clusters3seq),
           ":",
           BiocGenerics::strand(clusters3seq))
  names(clusters3seq) <- clusters3seq$cluster_name
  clusters3seq <- keepStandardChromosomes(clusters3seq, pruning.mode = "coarse")
  return(clusters3seq)
}

#' Title
#'
#' @param pattern
#' @param seqsUp
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors
#' @examples
kmer_scan <- function(pattern,seqsUp){
  upstreamA <- unlist(vmatchPattern("AAAAAAAA", seqsUp)) %>% as.data.frame(.) %>% dplyr::pull(names)
  return(upstreamA)
}

# plot nucleotides


#' Title
#'
#' @param cluster_regions
#' @param genome
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors Biostrings
#' @examples
getNucleotideProb <- function(cluster_regions, genome){
  seqlevelsStyle(cluster_regions) <- "UCSC"
  cluster_regions <- keepStandardChromosomes(cluster_regions, pruning.mode = "coarse")
  flank = 50
  pas.flank <- flank(cluster_regions , width = flank, both = TRUE)
  pas.flank_seq <- Biostrings::getSeq(genome, pas.flank)
  prob <- consensusMatrix(pas.flank_seq,as.prob = TRUE)
  df <- reshape2::melt(prob)
  prob1 <- prob
  prop3 <- list(prob1, prob)
  df <- reshape2::melt(prop3)
  levels(df$Var1)[4] <- 'U'
  df = subset(df, grepl("A|C|G|U", Var1))
  return(df)
}
  # plot
#' Title
#'
#' @param x
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors ggplot2
#' @examples
nucleotidePlot <- function(x) {
    nt.plot <- ggplot(x, aes(Var2 - 50, value)) +
      geom_line(aes(color = Var1), alpha = 1, size = 1) +
      scale_color_manual(
        values =
          c(
            "A" = "red",
            "C" = "blue",
            "G" = "#ECC54E",
            "U" = "#A76BCF"
          )
      ) +
      xlab("Distance to Poly(A) site") +
      ylab("Relative frequency") +
      theme_classic() +
      ggtitle("Nucleotide probabilities") +
      theme(
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)
      )
    return(nt.plot)
  }


#' Title
#'
#' @param three_prime_ends.gr
#' @param genome
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors
#' @examples
plotNucleotides <- function(three_prime_ends.gr, genome) {
  nucleotide.probs <- getNucleotideProb(three_prime_ends.gr, genome)
  nucleotidePlot(nucleotide.probs)

  }




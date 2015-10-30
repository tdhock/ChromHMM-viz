works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/namedCapture@a31a38be12d4dad4aec6257a47c0e2307679589c")

pattern <- paste0(
  "(?<chrom>.*?)",
  ":",
  "(?<regionStart>[0-9,]+)",
  "-",
  "(?<regionEnd>[0-9,]+)",
  " ",
  "(?<label>[^ ]+)")

label.lines <- readLines("one_iPS_sample_labels.txt")
toint <- function(x)as.integer(gsub(",", "", x))
label.df <- str_match_named(label.lines, pattern, list(
  regionStart=toint, regionEnd=toint))
label.df$chunk.id <- cumsum(is.na(label.df[,1]))+1
not.na <- subset(label.df, !is.na(chrom))

readBigWig <- function
### Read part of a bigWig file into R as a data.table (assumes
### bigWigToBedGraph is present on your PATH).
(bigwig.file,
### path or URL of bigwig file.
 chrom,
### chromosome to read.
 start,
### position before reading.
 end
### plain text file where coverage is saved before reading into R.
 ){
  stopifnot(length(bigwig.file) == 1)
  stopifnot(length(chrom) == 1)
  stopifnot(length(start) == 1)
  stopifnot(length(end) == 1)
  stopifnot(is.character(bigwig.file))
  stopifnot(is.character(chrom))
  start <- as.integer(start)
  end <- as.integer(end)
  stopifnot(0 <= start)
  stopifnot(start < end)
  stopifnot(end < Inf)
  cmd <-
    sprintf("bigWigToBedGraph -chrom=%s -start=%d -end=%d %s /dev/stdout",
            chrom, start, end,
            bigwig.file)
  bg <- fread(cmd, drop=1)
  if(is.null(bg)){
    data.table(chromStart=integer(),
               chromEnd=integer(),
               count=integer())
  }else{
    setnames(bg, c("chromStart", "chromEnd", "norm"))
    stopifnot(0 <= bg$norm)
    nonzero <- bg[0 < norm, ]
    min.nonzero.norm <- min(nonzero[, norm])
    nonzero[, count := as.integer(norm/min.nonzero.norm) ]
    nonzero[, .(
      chromStart,
      chromEnd,
      count
      )]
  }
### data.table with columns chrom chromStart chromEnd count.
}

bigwig.file.vec <- Sys.glob("inducedpluripotentstemcell/*.bigwig")
bigwig.pattern <- paste0(
  "(?<donor>.*?)",
  "[.]",
  "(?<cellType>.*?)",
  "[.]",
  "(?<experiment>.*?)",
  "[.]signal[.]bigwig")
bigwig.mat <- str_match_named(basename(bigwig.file.vec), bigwig.pattern)

labels.by.chunk <- split(not.na, not.na$chunk.id)
coverage.dt.list <- list()
for(chunk.str in names(labels.by.chunk)){
  chunk.labels <- labels.by.chunk[[chunk.str]]
  chrom <- chunk.labels$chrom[1]
  regionStart <- min(chunk.labels$regionStart)
  regionEnd <- max(chunk.labels$regionEnd)
  region.size <- regionEnd - regionStart
  zoom.bases <- as.integer(region.size/10)
  chunkStart <- regionStart - zoom.bases
  chunkEnd <- regionEnd + zoom.bases
  position <- unique(as.integer(seq(chunkStart, chunkEnd, l=1000)))
  for(bigwig.i in seq_along(bigwig.file.vec)){
    experiment <- bigwig.mat[bigwig.i, "experiment"]
    bigwig.file <- bigwig.file.vec[[bigwig.i]]
    counts <- readBigWig(bigwig.file, chrom, chunkStart, chunkEnd)
    coverage <- with(counts, {
      approx(chromStart+1, count, position, method="constant", rule=2)
    })$y
    coverage.dt.list[[paste(chrom, experiment)]] <- 
      data.table(chrom, experiment, position, coverage)
  }
}
coverage.dt <- do.call(rbind, coverage.dt.list)

one_iPS_sample_labels <- list(coverage=coverage.dt, labels=data.table(not.na))

save(one_iPS_sample_labels, file="one_iPS_sample_labels.RData")

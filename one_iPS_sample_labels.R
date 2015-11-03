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

label.lines <- readLines(file.path("labels", "one_iPS_sample_labels.txt"))
toint <- function(x)as.integer(gsub(",", "", x))
label.df <- str_match_named(label.lines, pattern, list(
  regionStart=toint, regionEnd=toint))
label.df$chunk.id <- cumsum(is.na(label.df[,1]))+1
not.na <- subset(label.df, !is.na(chrom))

### Run fread but do not stop for an error on an empty file.
fread.or.null <- function(...){
  tryCatch({
    fread(...)
  }, error=function(e){
    NULL
  })
}

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
  bg <- fread.or.null(cmd, drop=1)
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

data.dir <- file.path("labels", "inducedpluripotentstemcell")
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*.bigwig"))
bigBed.file.vec <- Sys.glob(file.path(data.dir, "*.bigBed"))
filename.pattern <- paste0(
  "(?<donor>.*?)",
  "[.]",
  "(?<cellType>.*?)",
  "[.]",
  "(?<experiment>.*?)",
  "[.]",
  "(?<dataType>.*?)",
  "[.]",
  "(?<fileType>.*)")
bigwig.mat <- str_match_named(basename(bigwig.file.vec), filename.pattern)
bigBed.mat <- str_match_named(basename(bigBed.file.vec), filename.pattern)
stopifnot(bigwig.mat[, "experiment"] == bigBed.mat[, "experiment"])

labels.by.chrom <- split(not.na, not.na$chrom)
coverage.dt.list <- list()
chunk.dt.list <- list()
peaks.dt.list <- list()
for(chrom in names(labels.by.chrom)){
  chunk.labels <- labels.by.chrom[[chrom]]
  regionStart <- min(chunk.labels$regionStart)
  regionEnd <- max(chunk.labels$regionEnd)
  region.size <- regionEnd - regionStart
  zoom.bases <- as.integer(region.size/10)
  chunkStart <- regionStart - zoom.bases
  chunkEnd <- regionEnd + zoom.bases
  position <- unique(as.integer(seq(chunkStart, chunkEnd, l=1000)))
  chunk.dt.list[[chrom]] <-
    data.table(chrom, regionStart, regionEnd, chunkStart, chunkEnd)
  for(bigwig.i in seq_along(bigwig.file.vec)){
    experiment <- bigwig.mat[bigwig.i, "experiment"]
    bigwig.file <- bigwig.file.vec[[bigwig.i]]
    counts <- readBigWig(bigwig.file, chrom, chunkStart, chunkEnd)
    coverage <- with(counts, {
      approx(chromStart+1, count, position, method="constant", rule=2)
    })$y
    coverage.dt.list[[paste(chrom, experiment)]] <- 
      data.table(chrom, experiment, position, coverage)
    bigBed.file <- bigBed.file.vec[[bigwig.i]]
    tmp.bed <- tempfile()
    cmd <- sprintf("bigBedToBed -chrom=%s -start=%d -end=%d %s %s",
                   chrom, chunkStart, chunkEnd,
                   bigBed.file, tmp.bed)
    system(cmd)
    peaks <- fread.or.null(tmp.bed, drop=4:5)
    if(is.data.table(peaks)){
      setnames(peaks, c("chrom", "peakStart", "peakEnd"))
      peaks.dt.list[[paste(chrom, experiment)]] <-
        data.table(chrom, experiment, peaks)
    }
  }
}
coverage.dt <- do.call(rbind, coverage.dt.list)
chunk.dt <- do.call(rbind, chunk.dt.list)
peaks.dt <- do.call(rbind, peaks.dt.list)
setkey(chunk.dt, chrom, chunkStart, chunkEnd)

expected.changes <- c(constant=0, oneChange=1)
out.dir.vec <- Sys.glob("ChromHMM/ipsStates2-32/out-*")
all.segs.list <- list()
all.errors.list <- list()
all.iterations.list <- list()
for(out.dir in out.dir.vec){
  iterations.txt <- file.path(out.dir, "iterations.txt")
  iterations <- read.table(iterations.txt, header=TRUE)
  last.iteration <- tail(iterations, 1)
  maxStates <- as.integer(sub(".*-", "", out.dir))
  all.iterations.list[[out.dir]] <-
    data.table(maxStates, last.iteration)
  segments.bed <- Sys.glob(paste0(out.dir, "/*_segments.bed"))
  stopifnot(length(segments.bed) == 1)
  seg.dt <- fread(segments.bed)
  setnames(seg.dt, c("chrom", "segStart", "segEnd", "state"))
  setkey(seg.dt, chrom, segStart, segEnd)
  over.dt <- foverlaps(seg.dt, chunk.dt, nomatch=0L)
  over.dt[, start := ifelse(segStart < chunkStart, chunkStart, segStart)]
  over.dt[, end := ifelse(chunkEnd < segEnd, chunkEnd, segEnd)]
  segs.by.chrom <- split(over.dt, over.dt$chrom)
  for(chrom in names(segs.by.chrom)){
    chrom.segs <- segs.by.chrom[[chrom]]
    all.segs.list[[paste(out.dir, chrom)]] <-
      data.table(maxStates, chrom.segs)
    chrom.labels <- data.table(labels.by.chrom[[chrom]])
    change.vec <- chrom.segs$end[-1]
    chrom.labels$changes <- NA
    for(label.i in 1:nrow(chrom.labels)){
      l <- chrom.labels[label.i, ]
      in.region <- l$regionStart < change.vec & change.vec < l$regionEnd
      chrom.labels$changes[label.i] <- sum(in.region)
    }
    chrom.labels[, expected := expected.changes[paste(label)] ]
    chrom.labels[, fp := expected < changes ]
    chrom.labels[, fn := changes < expected ]
    chrom.labels[, status := ifelse(fp, "false positive",
                            ifelse(fn, "false negative", "correct"))]
    all.errors.list[[paste(out.dir, chrom)]] <-
      data.table(maxStates, chrom.labels)
  }
}
all.errors <- do.call(rbind, all.errors.list)
all.segs <- do.call(rbind, all.segs.list)
all.iterations <- do.call(rbind, all.iterations.list)

bigWigInfo <- function
### Run bigWigInfo to find chrom sizes.
(bigwig.file
### path or URL of bigwig file.
 ){
  stopifnot(is.character(bigwig.file))
  stopifnot(length(bigwig.file) == 1)
  cmd <- paste("bigWigInfo", bigwig.file, "-chroms | grep '^\\s'")
  chroms <- fread(cmd, header=FALSE, sep=" ")
  setnames(chroms, c("chrom", "chromStart", "chromEnd"))
  chroms$chrom <- sub("\\s*", "", chroms$chrom)
  chroms
}

chroms <- bigWigInfo(bigwig.file)

one_iPS_sample_labels <- list(
  chroms=chroms,
  chunks=chunk.dt,
  peaks=peaks.dt,
  coverage=coverage.dt,
  labels=data.table(not.na),
  segments=all.segs,
  iterations=all.iterations,
  errors=all.errors)

save(one_iPS_sample_labels, file="one_iPS_sample_labels.RData")

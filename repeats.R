works_with_R("3.2.2",
             data.table="1.9.6")

emission.list <- list()
transition.list <- list()
metrics.list <- list()
data.dir.vec <- Sys.glob(file.path("ChromHMM", "ips-repeats-1-10", "repeat-*"))
for(data.dir.i in seq_along(data.dir.vec)){
  data.dir <- data.dir.vec[[data.dir.i]]
  cat(sprintf("%4d / %4d repeats %s\n", data.dir.i, length(data.dir.vec),
              data.dir))
  repeat.chr <- sub(".*-", "", data.dir)
  repeat.int <- as.integer(repeat.chr)
  repeat.fac <- factor(repeat.int, 1:10)
  emissions.txt.vec <- Sys.glob(file.path(data.dir, "emissions*.txt"))
  stopifnot(100 == length(emissions.txt.vec))

  for(emissions.txt in emissions.txt.vec){
    iteration.str <- sub("[.]txt$", "", sub(".*_i", "", emissions.txt))
    iteration <- as.integer(iteration.str)
    emission.mat <- fread(emissions.txt)
    state.fac <- factor(emission.mat[[1]])
    for(experiment in names(emission.mat)[-1]){
      frequency <- emission.mat[[experiment]]
      emission.list[[paste(repeat.int, emissions.txt, experiment)]] <-
        data.table(repeat.int, repeat.fac,
                   iteration, experiment,
                   state=state.fac, frequency)
    }
  }

  transition.txt.vec <- Sys.glob(file.path(data.dir, "transition*.txt"))
  stopifnot(100 == length(transition.txt.vec))
  for(transition.txt in transition.txt.vec){
    iteration.str <- sub("[.]txt$", "", sub(".*_i", "", transition.txt))
    iteration <- as.integer(iteration.str)
    transition.mat <- fread(transition.txt, header=TRUE)
    for(state in names(transition.mat)[-1]){
      probability <- transition.mat[[state]]
      state.to <- factor(state, levels(state.fac))
      state.from <- factor(transition.mat[[1]], levels(state.fac))
      transition.list[[paste(repeat.int, transition.txt, state)]] <-
        data.table(repeat.int, repeat.fac,
                   iteration, state.to, state.from, probability)
    }
  }  

  metric.mat <- read.table(
    file.path(data.dir, "iterations.txt"),
    header=TRUE, na.strings="-")
  names(metric.mat)[2] <- "log.lik"
  names(metric.mat)[4] <- "seconds"
  metric.mat$log.neg.log.lik <- log(-metric.mat$log.lik)
  for(metric.name in names(metric.mat)[-1]){
    sub.mat <- metric.mat#[-(1:2),]
    metrics.list[[paste(repeat.int, metric.name)]] <-
      data.table(repeat.int, repeat.fac, metric.name,
                 metric.value=sub.mat[[metric.name]],
                 iteration=sub.mat[[1]])
  }
}

repeats <- list(
  metrics=do.call(rbind, metrics.list),
  transition=do.call(rbind, transition.list),
  emission=do.call(rbind, emission.list))

ChromHMMinit <- repeats
save(ChromHMMinit, file="~/R/aninint-examples/data/ChromHMMinit.RData",
     compress="xz")

save(repeats, file="repeats.RData")

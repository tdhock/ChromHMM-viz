works_with_R("3.2.2",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             data.table="1.9.6")

emissions.txt.vec <- Sys.glob("toby/emissions.mat/*")
emission.list <- list()
for(emissions.txt in emissions.txt.vec){
  iteration.str <- sub(".*_i", "", emissions.txt)
  iteration <- as.integer(iteration.str)
  emission.mat <- fread(emissions.txt)
  state.fac <- factor(emission.mat[[1]])
  for(experiment in names(emission.mat)[-1]){
    frequency <- emission.mat[[experiment]]
    emission.list[[paste(emissions.txt, experiment)]] <-
      data.table(iteration, experiment, state=state.fac, frequency)
  }
}
emission <- do.call(rbind, emission.list)

ggplot()+
  coord_equal()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_wrap("iteration")+
  scale_fill_gradient(low="white", high="blue")+
  geom_tile(aes(state, experiment, fill=frequency),
            data=emission)

transition.txt.vec <- Sys.glob("toby/transitions.mat/*")
transition.list <- list()
for(transition.txt in transition.txt.vec){
  iteration.str <- sub(".*_i", "", transition.txt)
  iteration <- as.integer(iteration.str)
  transition.mat <- fread(transition.txt, header=TRUE)
  for(state in names(transition.mat)[-1]){
    probability <- transition.mat[[state]]
    state.to <- factor(state, levels(state.fac))
    state.from <- factor(transition.mat[[1]], levels(state.fac))
    transition.list[[paste(transition.txt, state)]] <-
      data.table(iteration, state.to, state.from, probability)
  }
}
transition <- do.call(rbind, transition.list)

ggplot()+
  coord_equal()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_wrap("iteration")+
  scale_fill_gradient(low="white", high="blue")+
  geom_tile(aes(state.to, state.from, fill=probability),
            data=transition)

metric.mat <- read.table(
  "toby/iterations/iterations.txt",
  header=TRUE, na.strings="-")
names(metric.mat)[2] <- "log.lik"
names(metric.mat)[4] <- "seconds"
metric.mat$log.neg.log.lik <- log(-metric.mat$log.lik)
metrics.list <- list()
for(metric.name in names(metric.mat)[-1]){
  sub.mat <- metric.mat#[-(1:2),]
  metrics.list[[metric.name]] <-
    data.table(metric.name,
               metric.value=sub.mat[[metric.name]],
               iteration=sub.mat[[1]])
}
metrics <- do.call(rbind, metrics.list)

ggplot()+
  geom_line(aes(iteration, metric.value),
            data=metrics)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(metric.name ~ ., scales="free_y")

iterations <- list(
  metrics=metrics,
  transition=transition,
  emission=emission)

save(iterations, file="iterations.RData")

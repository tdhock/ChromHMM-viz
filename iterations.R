works_with_R("3.2.2",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             data.table="1.9.6")

emission.mat <- fread("Toby/emissions_10_999.txt")
state.fac <- factor(emission.mat[[1]])

emission.list <- list()
for(experiment in names(emission.mat)[-1]){
  frequency <- emission.mat[[experiment]]
  emission.list[[experiment]] <-
    data.table(experiment, state=state.fac, frequency)
}
emission <- do.call(rbind, emission.list)
emission[, list(total=sum(frequency)), by=state]
emission[, list(total=sum(frequency)), by=experiment]

ggplot()+
  coord_equal()+
  scale_fill_gradient(low="white", high="blue")+
  geom_tile(aes(state, experiment, fill=frequency),
            data=emission)

transition.mat <- fread("Toby/transitions_10_999.txt", header=TRUE)
transition.list <- list()
for(state in names(transition.mat)[-1]){
  probability <- transition.mat[[state]]
  state.to <- factor(state, levels(state.fac))
  state.from <- factor(transition.mat[[1]], levels(state.fac))
  transition.list[[state]] <-
    data.table(state.to, state.from, probability)
}
transition <- do.call(rbind, transition.list)

ggplot()+
  coord_equal()+
  scale_fill_gradient(low="white", high="blue")+
  geom_tile(aes(state.to, state.from, fill=probability),
            data=transition)

metric.mat <- fread("Toby/iteration-loglike.txt", na.strings="-")
metrics.list <- list()
for(metric.name in names(metric.mat)[-1]){
  sub.mat <- metric.mat[-(1:2),]
  metrics.list[[metric.name]] <-
    data.table(metric.name,
               metric.value=sub.mat[[metric.name]],
               metric=sub.mat[[1]])
}
metrics <- do.call(rbind, metrics.list)

ggplot()+
  geom_line(aes(metric, metric.value),
            data=metrics)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(metric.name ~ ., scales="free_y")

iterations <- list(
  metrics=metrics,
  transition=transition,
  emission=emission)

save(iterations, file="iterations.RData")

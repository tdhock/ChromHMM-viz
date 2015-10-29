works_with_R("3.2.2",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@14835d5638d4cf234219a16ee313859c770164ff",
             data.table="1.9.6")

load("iterations.RData")

viz <- list(
  emission=ggplot()+
    coord_equal()+
    scale_fill_gradient(low="white", high="blue")+
    geom_tile(aes(state, experiment, fill=frequency,
                  key=paste(state, experiment),
                  showSelected=iteration),
              data=iterations$emission),
  transition=ggplot()+
    coord_equal()+
    scale_fill_gradient(low="white", high="red")+
    geom_tile(aes(state.to, state.from, fill=probability,
                  key=paste(state.from, state.to),
                  showSelected=iteration),
              data=iterations$transition),
  metrics=ggplot()+
    make_tallrect(metrics, "iteration")+
    geom_line(aes(iteration, metric.value),
              data=iterations$metrics)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(metric.name ~ ., scales="free_y"),
  time=list(variable="iteration", ms=500),
  duration=list(iteration=500),
  title="ChromHMM parameter fitting for one iPS sample")

animint2dir(viz, "figure-iterations")




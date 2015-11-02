works_with_R("3.2.2",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@4257e8cf76eb5021a98010b6629b779a4f383b24",
             data.table="1.9.6")

load("iterations.RData")

viz <- list(
  parameters=ggplot()+
    ggtitle("parameters at selected iteration")+
    scale_fill_gradient(low="white", high="blue")+
    geom_tile(aes(state, experiment, fill=frequency,
                  key=paste(state, experiment),
                  showSelected=iteration),
              data=data.table(iterations$emission, parameters="emission"))+
    scale_color_gradient(low="white", high="red")+
    theme_bw()+
    theme_animint(height=500, width=350)+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(parameters ~ .,
               space="free",
               scales="free_y")+
    scale_y_discrete(drop=FALSE)+
    geom_point(aes(state.to, state.from, color=probability,
                  key=paste(state.from, state.to),
                  showSelected=iteration),
               size=10,
               data=data.table(iterations$transition, parameters="transition")),
  metrics=ggplot()+
    ggtitle("convergence metrics, select iteration")+
    make_tallrect(iterations$metrics, "iteration")+
    geom_line(aes(iteration, metric.value),
              data=iterations$metrics)+
    theme_bw()+
    theme_animint(height=500)+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(metric.name ~ ., scales="free_y"),
  time=list(variable="iteration", ms=500),
  duration=list(iteration=500),
  title="ChromHMM parameter fitting for one iPS sample")

animint2dir(viz, "figure-iterations")




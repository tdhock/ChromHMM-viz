works_with_R("3.2.2",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@4257e8cf76eb5021a98010b6629b779a4f383b24",
             data.table="1.9.6")

load("repeats.RData")

active.marks <- repeats$emission[, {
  list(active.marks=sum(0.01 < frequency))
}, by=.(iteration, repeat.fac, state)]
active.states <- active.marks[, {
  list(active.states=sum(0 < active.marks))
}, by=.(iteration, repeat.fac)]
last.iteration <- repeats$metrics[iteration==100, ]

viz <- list(
  parameters=ggplot()+
    ggtitle("parameters at selected iteration")+
    scale_fill_gradient(low="white", high="blue")+
    geom_tile(aes(state, experiment, fill=frequency,
                  key=paste(state, experiment),
                  showSelected=repeat.fac,
                  showSelected2=iteration),
              chunk_vars=c("repeat.fac"),
              data=data.table(repeats$emission, parameters="emission"))+
    scale_color_gradient(low="white", high="red")+
    theme_bw()+
    theme_animint(height=500, width=400)+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(parameters ~ .,
               space="free",
               scales="free_y")+
    scale_y_discrete(drop=FALSE)+
    geom_point(aes(state.to, state.from, color=probability,
                   key=paste(state.from, state.to),
                   showSelected=repeat.fac,
                   showSelected2=iteration),
               size=8,
               chunk_vars=c("repeat.fac"),
               data=data.table(repeats$transition, parameters="transition")),
  metrics=ggplot()+
    ggtitle("convergence metrics, select iteration")+
    make_tallrect(repeats$metrics, "iteration")+
    geom_line(aes(iteration, active.states,
                   group=repeat.fac,
                   clickSelects=repeat.fac),
               size=3,
               alpha=0.6,
               data=data.table(metric.name="active.states",
                 active.states))+
    geom_line(aes(iteration, metric.value,
                  clickSelects=repeat.fac,
                  group=repeat.fac),
              size=3,
              alpha=0.6,
              data=repeats$metrics[metric.name != "Change", ])+
    theme_bw()+
    theme_animint(height=500)+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(metric.name ~ ., scales="free_y"),
  last=ggplot()+
    ggtitle("last iteration, select initialization")+
    theme_bw()+
    theme_animint(height=500, width=400)+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(metric.name ~ ., space="fixed", scales="free")+
    geom_point(aes(repeat.fac, active.states,
                   clickSelects=repeat.fac),
               size=5,
               data=data.table(metric.name="active.states",
                 active.states[iteration==100, ]))+
    geom_point(aes(repeat.fac, metric.value,
                   clickSelects=repeat.fac),
               size=5,
               data=last.iteration[metric.name != "Change", ])+
    scale_x_discrete("random initialization")+
    scale_y_continuous(""),
  duration=list(iteration=500),
  first=list(iteration=100),
  time=list(variable="iteration", ms=500),
  title="10 ChromHMM fits for 6 experiments on one iPS sample")

animint2dir(viz, "figure-repeats")




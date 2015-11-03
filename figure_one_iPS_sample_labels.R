works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@4257e8cf76eb5021a98010b6629b779a4f383b24")

load("one_iPS_sample_labels.RData")

fp.fn <- one_iPS_sample_labels$errors[, {
  list(fp=sum(fp), fn=sum(fn))
}, by=maxStates]
fp.fn[, errors := fp + fn]
fp.fn[order(maxStates), ]
error.curves.list <- list()
for(error.type in c("fp", "fn", "errors")){
  incorrect.regions <- fp.fn[[error.type]]
  maxStates <- fp.fn$maxStates
  error.curves.list[[error.type]] <-
    data.table(error.type, maxStates, incorrect.regions)
}
error.curves <- do.call(rbind, error.curves.list)

iterations <- one_iPS_sample_labels$iterations
iterations[, minutes := TotalTimeSec/60]
abbrev.vec <- c(
  Iteration="iterations",
  EstimatedLogLikelihood="log.lik",
  minutes="minutes")
metrics.list <- list()
for(full in names(abbrev.vec)){
  metric.name <- abbrev.vec[[full]]
  metric.value <- iterations[[full]]
  maxStates <- iterations$maxStates
  metrics.list[[metric.name]] <-
    data.table(metric.name, maxStates, metric.value)
}
metrics <- do.call(rbind, metrics.list)

max.log.lik <- metrics.list$log.lik[metric.value==max(metric.value), ]

fp.fn.colors <- c(FP="skyblue",
                  fp="skyblue",
                  fn="#E41A1C",
                  FN="#E41A1C",
                  tn="white",
                  tp="grey",
                  errors="black")
fp.fn.sizes <- c(errors=2, fp=1, fn=1)
ggplot()+
  geom_hline(aes(yintercept=metric.value),
             data=max.log.lik,
             color="grey50")+
  geom_point(aes(maxStates, metric.value),
             data=max.log.lik,
             color="grey50")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(metric.name ~ ., scales="free")+
  geom_line(aes(maxStates, metric.value),
            data=metrics)+
  geom_line(aes(maxStates, incorrect.regions,
                group=error.type,
                size=error.type,
                color=error.type),
            data=data.table(error.curves, metric.name="incorrect.regions"))+
  scale_color_manual(values=fp.fn.colors)+
  scale_size_manual(values=fp.fn.sizes)

chunks <- one_iPS_sample_labels$chunks
chroms <- one_iPS_sample_labels$chroms[chrom %in% chunks$chrom, ]

error.by.chrom <- one_iPS_sample_labels$errors[, {
  list(fp=sum(fp), fn=sum(fn))
}, by=.(maxStates, chrom)]
setkey(chunks, chrom)
setkey(error.by.chrom, chrom)
error.text <- error.by.chrom[chunks]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_wrap("maxStates")+  
  scale_x_continuous("position on chromosome (mega bases)")+
  geom_segment(aes(chromStart/1e6, chrom,
                   xend=chromEnd/1e6, yend=chrom),
               color="grey",
               data=chroms)+
  geom_point(aes(regionStart/1e6, chrom),
             size=5,
             data=chunks)+
  scale_color_manual(values=fp.fn.colors)+
  geom_text(aes(regionStart/1e6, chrom,
                label=paste(fp, "fp"),
                color=error.type),
            data=data.table(error.text[0 < fp, ], error.type="fp"),
            hjust=1)+
  geom_text(aes(regionStart/1e6, chrom,
                label=paste(fn, "fn"),
                color=error.type),
            data=data.table(error.text[0 < fn, ], error.type="fn"),
            hjust=0)

ann.colors <-
  c(noPeaks="#f6f4bf",
    constant="#f6f4bf",
    peakStart="#ffafaf",
    oneChange="#ff4c4c",
    peaks="#a445ee")

## Static plot of just the data.
(gg <- ggplot()+
  scale_fill_manual(values=ann.colors)+
  geom_tallrect(aes(xmin=regionStart/1e3, xmax=regionEnd/1e3,
                    fill=label),
                alpha=0.5,
                color="grey",
                data=one_iPS_sample_labels$labels)+
  geom_line(aes(position/1e3, coverage),
            color="grey50",
            data=one_iPS_sample_labels$coverage)+
  geom_segment(aes(peakStart/1e3, 0,
                   xend=peakEnd/1e3, yend=0),
               data=one_iPS_sample_labels$peaks,
               size=5,
               color="deepskyblue")+
  geom_point(aes(peakStart/1e3, 0),
             data=one_iPS_sample_labels$peaks,
             shape=1,
             color="deepskyblue")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(experiment ~ chrom, scales="free")+
  scale_x_continuous("position on chromosome (kilo bases = kb)")
 )

png("figure_one_iPS_sample_labels.png", w=10, h=7, res=100, units="in")
print(gg)
dev.off()

mergeChunk <- function(dt){
  setkey(dt, chrom)
  dt[chunks]
}
normalize <- function(pos.vec, chunkStart, chunkEnd){
  (pos.vec - chunkStart)/(chunkEnd-chunkStart)
}

viz <- list(
  title="Quantifying ChromHMM accuracy using labels",
  coverage=ggplot()+
    ggtitle("Coverage, labels, peaks, selected model")+
    scale_fill_manual(values=ann.colors)+
    scale_x_continuous("relative position on chromosome")+
    scale_y_continuous(breaks=function(lim.vec){
      floor(lim.vec[2])
    })+
    coord_cartesian(xlim=c(0, 1))+
    geom_tallrect(aes(xmin=normalize(regionStart, chunkStart, chunkEnd),
                      xmax=normalize(regionEnd, chunkStart, chunkEnd),
                      showSelected=chrom,
                      fill=label),
                  alpha=0.5,
                  color="grey",
                  data=mergeChunk(one_iPS_sample_labels$labels))+
    geom_line(aes(normalize(position, chunkStart, chunkEnd),
                  coverage,
                  showSelected=chrom),
              color="grey50",
              data=mergeChunk(one_iPS_sample_labels$coverage))+
    geom_segment(aes(normalize(peakStart, chunkStart, chunkEnd), 0,
                     showSelected=chrom,
                     xend=normalize(peakEnd, chunkStart, chunkEnd), yend=0),
                 data=mergeChunk(one_iPS_sample_labels$peaks),
                 size=5,
                 color="deepskyblue")+
    geom_point(aes(normalize(peakStart, chunkStart, chunkEnd), 0,
                   showSelected=chrom),
               data=mergeChunk(one_iPS_sample_labels$peaks),
               shape=1,
               color="deepskyblue")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    theme_animint(width=1200, height=500)+
    geom_tallrect(aes(xmin=normalize(start, chunkStart, chunkEnd),
                      xmax=normalize(end, chunkStart, chunkEnd),
                      showSelected=chrom,
                      showSelected2=maxStates),
                  color="green",
                  fill="transparent",
                  alpha=0.25,
                  data=one_iPS_sample_labels$segments)+
    geom_text(aes(normalize((start+end)/2, chunkStart, chunkEnd),
                  0,
                  label=sub("^E", "", state),
                  showSelected=chrom,
                  showSelected2=maxStates),
              data=data.table(one_iPS_sample_labels$segments,
                              experiment="H3K27ac"))+
    geom_tallrect(aes(xmin=normalize(regionStart, chunkStart, chunkEnd),
                      xmax=normalize(regionEnd, chunkStart, chunkEnd),
                      showSelected=maxStates,
                      showSelected2=chrom,
                      linetype=status),
                  color="black",
                  fill="transparent",
                  data=mergeChunk(one_iPS_sample_labels$errors))+
    scale_linetype_manual("error type",
                          limits=c("correct", 
                            "false negative",
                            "false positive"
                                   ),
                          values=c(correct=0,
                            "false negative"=3,
                            "false positive"=1))+
    facet_grid(experiment ~ ., scales="free"),
  chroms=ggplot()+
    ggtitle("Select labeled region")+
    scale_x_continuous("position on chromosome (mega bases)")+
    geom_segment(aes(chromStart/1e6, chrom,
                     xend=chromEnd/1e6, yend=chrom),
                 color="grey",
                 data=chroms)+
    geom_point(aes(regionStart/1e6, chrom,
                   clickSelects=chrom),
               size=5,
               data=chunks)+
    scale_color_manual(values=fp.fn.colors)+
    geom_text(aes(regionStart/1e6, chrom,
                  clickSelects=chrom,
                  showSelected=maxStates,
                  label=paste(fp, "fp"),
                  color=error.type),
              data=data.table(error.text[0 < fp, ], error.type="fp"),
              hjust=1)+
    geom_text(aes(regionStart/1e6, chrom,
                  clickSelects=chrom,
                  showSelected=maxStates,
                  label=paste(fn, "fn"),
                  color=error.type),
              data=data.table(error.text[0 < fn, ], error.type="fn"),
              hjust=0)+
   guides(color="none"),
  curves=ggplot()+
    ggtitle("Select number of states")+
    make_tallrect(error.curves, "maxStates")+
    geom_hline(aes(yintercept=metric.value),
               data=max.log.lik,
               alpha=0.5,
               color="black")+
    geom_point(aes(maxStates, metric.value),
               data=max.log.lik,
               size=4,
               alpha=0.5,
               color="black")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(metric.name ~ ., scales="free")+
    geom_line(aes(maxStates, metric.value),
              data=metrics)+
    geom_line(aes(maxStates, incorrect.regions,
                  group=error.type,
                  size=error.type,
                  color=error.type),
              data=data.table(error.curves, metric.name="incorrect.regions"))+
    scale_color_manual(values=fp.fn.colors)+
    scale_size_manual(values=fp.fn.sizes)
  )  

animint2dir(viz, "figure_one_iPS_sample_labels")

works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             "tdhock/animint@5f06aee037663c248297592fac6215388230d415")

load("one_iPS_sample_labels.RData")

ann.colors <-
  c(noPeaks="#f6f4bf",
    constant="#f6f4bf",
    peakStart="#ffafaf",
    oneChange="#ff4c4c",
    peaks="#a445ee")

gg <- ggplot()+
  scale_fill_manual(values=ann.colors)+
  geom_tallrect(aes(xmin=regionStart/1e3, xmax=regionEnd/1e3,
                    fill=label),
                alpha=0.5,
                data=one_iPS_sample_labels$labels)+
  geom_line(aes(position/1e3, coverage),
            color="grey50",
            data=one_iPS_sample_labels$coverage)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(experiment ~ chrom, scales="free")

png("figure_one_iPS_sample_labels.png", w=10, h=7, res=100, units="in")
print(gg)
dev.off()

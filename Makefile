figure_one_iPS_sample_labels/index.html: figure_one_iPS_sample_labels.R one_iPS_sample_labels.RData
	R --no-save < $<
one_iPS_sample_labels.RData: one_iPS_sample_labels.R
	R --no-save < $<
figure-iterations/index.html: figure-iterations.R iterations.RData
	R --no-save < $<
iterations.RData: iterations.R
	R --no-save < $<

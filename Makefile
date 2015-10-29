figure-iterations/index.html: figure-iterations.R iterations.RData
	R --no-save < $<
iterations.RData: iterations.R
	R --no-save < $<

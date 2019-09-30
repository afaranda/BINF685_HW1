
# Import Data Set and Plot Histogram
df <-read.csv("peptide_data.txt")
png("peptide_data_histogram.png")
hist(
	(0-df$Energy), freq = FALSE,
	xlab = "Energy (Sign inverted)",
	main = paste(
		"Distribution of Energy values \n over", 
		nrow(df), "random configurations"
	), 
	ylim = c(0, 0.7),
	col = "red"
)
dev.off()

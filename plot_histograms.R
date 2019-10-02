# Import Data Sets and Plot Histograms
for(f in list.files(pattern = "peptide_data")){
      temp<-gsub("peptide_data_", "", f)
      temp<-gsub(".txt","",temp)

      df <-read.csv(f)
      png(paste("peptide_data_",temp,".png"))
      hist(
	(0-df$Energy), freq = FALSE,
	xlab = "Energy (Sign inverted)",
	main = paste(
	     "Distribution of Energy values \n over", 
	     row(df), "random configurations at temperature", temp
	), 
	ylim = c(0, 0.7),
	col = "red"
	)
	dev.off()
}


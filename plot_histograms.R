# Import Data Sets and Plot Histograms
for(f in list.files(pattern = "peptide_data")){
      temp<-gsub("peptide_data_", "", f)
      temp<-gsub(".txt","",temp)

      df <-read.csv(f)
      png(paste("peptide_hist_",temp,".png", sep=""))
      hist(
		(0-df$Energy), freq = FALSE,
		xlab = "Energy (Sign inverted)",
		main = paste(
	    		"Distribution of Energy values \n over", 
	     	nrow(df), "random configurations at temperature", temp
		), 
		ylim = c(0, 0.8),
		col = "red"
	)
	dev.off()
}


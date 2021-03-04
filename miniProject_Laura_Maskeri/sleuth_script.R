#load sleuth package
library(sleuth)

#read in sleuth table, tab delimited
stab <- read.table("sleuth_table.txt",sep = "\t", header=TRUE,stringsAsFactors=FALSE)
#print(stab)

#initialize sleuth object
so <- sleuth_prep(stab)

##Perform differential expression analysis comparing 2dpi to 6dpi

#fit a model comparing the two conditions 
so <- sleuth_fit(so, ~condition, 'full')
#so

#fit the reduced model to compare in the likelihood ratio test 
so <- sleuth_fit(so, ~1, 'reduced')
#so

#likelihood ratio test for differential expression between conditions 
so <- sleuth_lrt(so, 'reduced', 'full')


##Extract the results and look at most significant

#load the dplyr package for data.frame filtering
library(dplyr)

#extract the test results from the sleuth object 
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 


#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval) 

#Write the following details for
#each significant transcript (FDR < 0.05) to your log file, include a header row, and tab-delimit each item:
#target_id test_stat pval qval

output <- dplyr::select(sleuth_significant,target_id, test_stat, pval, qval)
#output

#write a txt file to folder in order to write results to log file
write.table(output, file="sleuth_results.txt", quote = FALSE, row.names = FALSE)







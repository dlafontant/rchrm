# R-CHRM
Performs recursive clustering by heterogenerity of measures
The rchrmclus function returns a nested list. 
Suppose your object is called myclus then myclus[[1]] contains a number of data frames equivalent to the number of outcomes used. Each data frame contains the repeated measures in wide format and cluster variables assoicated with each time point. myclus[[2]] contains matrices equivalent to the number of outcomes used. If intermittent missing values were present in the data, these matrices would contain the imputed values unlike the raw data contained in the data frames in myclus[[1]]. 

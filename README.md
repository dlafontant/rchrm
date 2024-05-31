# R-CHRM
Performs recursive clustering by heterogenerity of measures
The rchrmclus function returns a nested list. 
Suppose your object is called myclus then myclus[[1]] contains a number of data frames equal to the number of outcomes used. Each data frame contains the repeated measures in wide format and cluster variables assoicated with each time point. myclus[[2]] contains matrices equal to the number of outcomes used. If intermittent missing values were present in the data, each matrix will display the imputed values not contained in the data frames in myclus[[1]]. Prior to running the algorithm, please make sure that your ID vector 'id' is sorted and that each outcome (in wide) is sorted accordingly. 

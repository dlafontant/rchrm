
suppressPackageStartupMessages({
  library(kml3d)
  library(kml)
  library(cluster)
  library(flexclust)
  library(whitening)
  library(Kmedians)
  library(Rcpp)
  library(Matrix)
  library(Rcpp)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(nlme)
  library(data.table)
  library(MASS)
  library(sas7bdat)
  library(boot)
})

#The following function returns a AR(1) matrix. It takes in a the length of the diagonal matrix and the correlation
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                    (1:n - 1))
  rho^exponent
}

#The following function returns a long matrix. It takes in an array or a list, concatenates each object and returns the matrix
# Function to reshape array of matrices
reshape_array <- function(input_array) {

  if(class(input_array)=="list"){
    k <- length(input_array) # Number of matrices in the array
    i <- nrow(input_array[[1]]) # Number of rows in each matrix
    j <- ncol(input_array[[1]]) # Number of columns in each matrix

    # Check if all matrices have the same number of rows and columns
    if (any(sapply(input_array, function(mat) nrow(mat) != i)) ||
        any(sapply(input_array, function(mat) ncol(mat) != j))) {
      stop("All matrices in the array must have the same dimensions.")
    }

    # Reshape the array of matrices into a long matrix
    long_matrix <- matrix(NA, nrow = k * i, ncol = j)
    for (m in 1:k) {
      long_matrix[((m - 1) * i + 1):(m * i), ] <- input_array[[m]]
    }
  }

  else if(class(input_array)=="array"){
    k <- dim(input_array)[3] # Number of matrices in the 3D array
    i <- dim(input_array)[1] # Number of rows in each matrix
    j <- dim(input_array)[2] # Number of columns in each matrix

    # Reshape the 3D array into a long matrix
    long_matrix <- matrix(NA, nrow = k * i, ncol = j)
    for (m in 1:k) {
      long_matrix[((m - 1) * i + 1):(m * i), ] <- input_array[, , m]
    }
  }

  return(long_matrix)
}

#The following function takes in a list or array of matrices of the same dimensions
#and turns it into a matrix.
#Type = wide returns a matrix of double the number of columns.
#type = long returns a matrix of double the number of rows
listtomat <- function(arr, type=c("wide","long")){
  if(class(arr)[1]=="array"){
    y1values_wide <- matrix(sapply(arr, cbind),nrow(arr[,,1]), ncol(arr[,,1])*dim(arr)[3])
    #y1values_long <- matrix(sapply(arr, rbind),nrow(arr[,,1])*dim(arr)[3], ncol(arr[,,1]))
    y1values_long <- reshape_array(arr)

  }
  #Stack data in long format as if we had one variable with two time points. This will allow us to calculate E the temporal correlation matrix.
  else if(class(arr)[1]=="list"){
    y1values_wide <- matrix(sapply(arr, cbind),nrow(arr[[1]]), ncol(arr[[1]])*length(arr))
    #y1values_long <- matrix(sapply(arr, rbind),nrow(arr[[1]])*length(arr), ncol(arr[[1]]))
    y1values_long <- reshape_array(arr)
  }
  if (type=="wide"){return(y1values_wide)}
  else if (type=="long"){return(y1values_long)}
}


#THe following function takes in a array and returns a covariance matrix that takes both
#variable and temporal correlations into account. You may provide your variable correlation
#matrix at baseline. Otherwise, it will compute the variable correlation using the first
#column of each array.
varcov <- function(arr, cor1=FALSE, var_a = TRUE, cor_o=TRUE, cor_t=TRUE, covmethod="pearson"){

  if(class(arr)[1]=="array"){
    y1values_wide <- matrix(sapply(arr, cbind),nrow(arr[,,1]), ncol(arr[,,1])*dim(arr)[3])
    y1values_long <- reshape_array(arr)
    blvalues <- matrix(0, nrow = dim(arr)[1], ncol=dim(arr[3]))
    for(e in 1:dim(arr)[3]){
      blvalues[,e] <- arr[,1,e]
    }
  }
  #Stack data in long format as if we had one variable with two time points. This will allow us to calculate E the temporal correlation matrix.
  else if(class(arr)[1]=="list"){
    y1values_wide <- matrix(sapply(arr, cbind),nrow(arr[[1]]), ncol(arr[[1]])*length(arr))
    y1values_long <- reshape_array(arr)
    blvalues <- matrix(0, nrow = nrow(arr[[1]]), ncol=length(arr))
    for(e in 1:length(arr)){
      blvalues[,e] <- arr[[e]][,1]
    }
  }

  if(class(cor1)[1]=="matrix"){}
  else if(class(cor1)[1]=="logical" & cor1==FALSE){cor1 <- cor(blvalues, method=covmethod)}

  cov2 <- cov(y1values_wide, method=covmethod)
  Std_dev_diag <- sqrt(diag(cov2))

  A <- diag(Std_dev_diag)
  P <- cor1
  E <- cor(y1values_long, method=covmethod)

  if(cor_o==FALSE){
    P <- diag(rep(1,nrow(P)))
  }
  if(cor_t==FALSE){
    E <- diag(rep(1,nrow(E)))
  }
  if(var_a==FALSE){
    A <- diag(rep(1,nrow(A)))
  }
  C <- P %x% E

  sigma <- A^(1/2) %*% C %*% A^(1/2)
  return(sigma)
}

#I created a whiten_cov function (borrowed from whitening package).
#whiten_cov function whitens a data matrix X using the empirical
#covariance matrix cov(X) or a provided covariance matrix S as basis
#for computing the whitening transformation. Whitening is a linear
#transformation z = Wx where the whitening matrix satisfies the
#constraint WTW = Σ−1 where Σ = Cov(x). This function implements various
#natural whitening transformations discussed in Kessy, Lewin, and
#Strimmer (2018).
whiten_cov <- function (X, S=FALSE, center = FALSE, method = c("ZCA", "ZCA-cor", "PCA",
                                                               "PCA-cor", "Cholesky"), covmethod = "pearson")
{

  if(class(S)[1] =="matrix"){}
  else if (class(S)=="logical"){
    if(S==FALSE){S=cov(X, method=covmethod)}
  }
  method = match.arg(method)
  W = whiteningMatrix(S, method = method)
  Z = tcrossprod(X, W)
  if (center)
    Z = sweep(Z, 2, colMeans(Z))
  colnames(Z) = paste0("L", 1:ncol(X))
  attr(Z, "method") = method
  return(Z)
}

#Use cholesky decomposition to change X to X*
trickkml <- function(x){
  traj2 <- whiten_cov(X=x, S=FALSE, center=FALSE, method=c("Cholesky"))
  return(traj2) # scale to trick kml
}

#The following function provides distances from many methods including mahalanobis
#c("euclidean","manhattan" , "mahalanobis","minkowski", "maximum", "canberra", "binary",
#"mahalanobis-ZCA", "mahalanobis-ZCA-cor","mahalanobis-PCA", "mahalanobis-PCA-cor",
#"mahalanobis-Cholesky"
distance <- function(y, center=TRUE, method="euclidean"){

  if(method %in% c("manhattan", "euclidean", "minkowski", "maximum", "canberra", "binary")){
    return(dist(y,method=method))
  }
  else if(method=="mahalanobis"){
    vc1 <- ginv(var(y))
    vc2 <- chol(vc1, pivot=TRUE)  # Cholesky decomposition: t(vc2) %*% vc2 = vc1
    traj2 <- y %*% t(vc2)  # scale
    #print(traj2)
    return(dist(traj2))
  }
  else if(method=="mahalanobis2"){
    #y <- (y - mean(y,na.rm=TRUE)) / sd(y,na.rm=TRUE)
    N <- nrow(y)
    vc1 <- ginv(var(y))
    P2 <- matrix(0, nrow = N, ncol = N)
    for(i1 in 1:N){
      t1 <- y[i1, , drop = FALSE]
      for(i2 in 1:N){
        if(i1 < i2){
          tmp1 <- t1 - y[i2, , drop = FALSE]
          tmp2 <- sqrt(tmp1 %*% vc1 %*% t(tmp1))  # Now there is vc1
          P2[i1, i2] <- tmp2
          P2[i2, i1] <- tmp2
        }
      }
    }
    return(as.dist(P2))
  }
  else if(method=="mahalanobis-ZCA"){
    traj2 <- whiten(X=y, center=center, method=c("ZCA"))
    return(dist(traj2))
  }
  else if(method=="mahalanobis-ZCA-cor"){
    traj2 <- whiten(X=y, center=center, method=c("ZCA-cor"))
    return(dist(traj2))
  }
  else if(method=="mahalanobis-PCA"){
    traj2 <- whiten(X=y, center=center, method=c("PCA"))
    return(dist(traj2))
  }
  else if(method=="mahalanobis-PCA-cor"){
    traj2 <- whiten(X=y, center=center, method=c("PCA-cor"))
    return(dist(traj2))
  }
  else if(method=="mahalanobis-Cholesky"){
    traj2 <- whiten(X=y, center=center, method=c("Cholesky"))
    return(dist(traj2))
  }
}


#The folloinwg function returns a vector of numeric values. It takes in any vector and
#repalces it's values with numbers 1,2,3.... regarless of order.
relabel_numbers <- function(vector) {
  unique_numbers <- unique(vector)
  relabeled_vector <- rep(NA, length(vector))
  for (i in 1:length(unique_numbers)) {
    relabeled_vector[vector == unique_numbers[i]] <- i
  }
  return(relabeled_vector)
}


#The following function performs kmeans 20 times starting with different centers
#each time and choosing the cluster vector resulting in the largest gap (i.e smallest Wk)
kmeans2 <- function(x,k,d.power=2, dis="euclidean"){

  W.kb <- function(X, clus) {
    n <- nrow(X)
    ii <- seq_len(n)
    0.5 * sum(vapply(split(ii, clus), function(I) {
      xs <- X[I, , drop = FALSE]
      #sum(dist(xs)^d.power/nrow(xs))
      sum(distance(xs, method=dis)^d.power/nrow(xs))
    }, 0))
  }
  #The following function computes the cluster dissimilarity mentioned in Tibshirani.
  #it is used to return the cluster vector yielding the largest GAP statistic.
  multikmeans <- function(num){
    kmeansalgo <- c("randomK", "maxDist", "kmeans+", "kmeans-","kmeans++","kmeans--")
    unique_rows <- unique(x)
    kmeanscenters <- initializePartition(nbClusters=k, lengthPart=nrow(unique_rows), method = kmeansalgo[round(runif(1, 1,length(kmeansalgo)))], data=unique_rows)
    # print(unique_rows[!is.na(kmeanscenters),])
    #kk <- kmeans(x, centers=unique_rows[!is.na(kmeanscenters),], nstart=round(runif(1, 5,20)), iter.max=100)
    ur <- unique_rows[!is.na(kmeanscenters),]
    ur2 <-ur[order(ur[,1]),]
    km <- kcca(x=x, k = ur2, family=kccaFamily("kmeans"), simple=TRUE)
    cluster <-  relabel_numbers(clusters(km))
    #print(cluster)
    return(list("clus"=cluster,"wk"= W.kb(X=x, clus=cluster)))
  }

  vecs <- lapply(1:20, multikmeans)
  allwk <- unlist(lapply(1:20, function(x){vecs[[x]][[2]]}))
  min_index <- which(allwk == min(allwk))[1]

  cluster <-vecs[[min_index]][[1]]
}

#The following function performs kmeans 20 times starting with different centers
#each time and choosing the cluster vector resulting in the largest gap (i.e smallest Wk)
kmedians2 <- function(x,k,d.power=2, dis="euclidean"){

  W.kb <- function(X, clus) {
    n <- nrow(X)
    ii <- seq_len(n)
    0.5 * sum(vapply(split(ii, clus), function(I) {
      xs <- X[I, , drop = FALSE]
      #sum(dist(xs)^d.power/nrow(xs))
      sum(distance(xs, method=dis)^d.power/nrow(xs))
    }, 0))
  }
  #The following function computes the cluster dissimilarity mentioned in Tibshirani.
  #it is used to return the cluster vector yielding the largest GAP statistic.
  multikmedians <- function(num){
    kmeansalgo <- c("randomK", "maxDist", "kmeans+", "kmeans-","kmeans++","kmeans--")
    unique_rows <- unique(x)
    kmeanscenters <- initializePartition(nbClusters=k, lengthPart=nrow(unique_rows), method = kmeansalgo[round(runif(1, 1,length(kmeansalgo)))], data=unique_rows)

    ur <- unique_rows[!is.na(kmeanscenters),]
    ur2 <-ur[order(ur[,1]),]
    #kk <- kmeans(x, centers=unique_rows[!is.na(kmeanscenters),], nstart=round(runif(1, 5,20)), iter.max=100)
    km <- kcca(x=x, k = ur2, family=kccaFamily("kmedians"), simple=TRUE)
    cluster <-  relabel_numbers(clusters(km))
    #print(cluster)
    return(list("clus"=cluster,"wk"= W.kb(X=x, clus=cluster)))
  }

  vecs <- lapply(1:20, multikmedians)
  allwk <- unlist(lapply(1:20, function(x){vecs[[x]][[2]]}))
  min_index <- which(allwk == min(allwk))[1]

  cluster <- vecs[[min_index]][[1]]
}

#The following function grabs an array of matrices,
#and performs kml3d on the matrices.
#It returns a vector of clusters.
KML2d <- function(mat, id, k){
  par_kml1 <- parALGO(
    saveFreq = Inf,  # turn off save
    imputationMethod = "copyMean"
  )
  kk <- cld(traj = mat, idAll = id)  # note the maxNA
  invisible(capture.output(kml(object = kk, parAlgo = par_kml1, nbClusters=k)))
  cluster <- relabel_numbers(as.numeric(factor(getClusters(kk, k))))
  return(cluster)
}


KMLmedian2d <- function(mat, id, k){
  #centerMethod = function(x)median(x,na.rm=TRUE)
  par_kml1 <- parALGO(
    saveFreq = Inf,  # turn off save
    imputationMethod = "copyMean",
    centerMethod =function(x)median(x,na.rm=TRUE)
  )
  kk <- cld(traj = mat, idAll = id)  # note the maxNA
  invisible(capture.output(kml(object = kk, parAlgo = par_kml1, nbClusters=k)))
  cluster <- relabel_numbers(as.numeric(factor(getClusters(kk, k))))

  return(cluster)
}




#The following function takes in a matrix (or a vector) and computes the Gap statistics
#It takes in as parameter the matrix or vector, the partitioning function (like kmeans),
#the vector of IDs ID for each individual
#the number of Clusters K
#The time vector to specify the time points associated with each visit
#The number of bootstrap stamples to draw for the Tibshirani Gap statisitic "B"
#The power of the distance d.power
#Whether to use scaledPCA or the original Tibshirani Gap statistic "spaceH0"
#Whether printing is needed every 50 samples drawn from the B replicates
#the type of distance to use gapdis: options include euclidean, Cholesky, etc...
gapfunction <- function (x, FUNcluster, id, K, timevector, B = 100, d.power = 2, spaceH0 = c("scaledPCA",
                                                                                             "original"), verbose = FALSE, gapdis="euclidean", ...)
{

  #ind <- 1:length(timevector)

  #x <- x[, ind]
  #Subset matrix to timevector
  if(length(timevector)==1){
    x <- t(t(x))
  }

  stopifnot(is.function(FUNcluster), length(dim(x)) == 2, K >=
              1, (n <- nrow(x)) >= 1, ncol(x) >= 1)
  if (B != (B. <- as.integer(B)) || (B <- B.) <= 0)
    stop("'B' has to be a positive integer")
  cl. <- match.call()
  if (is.data.frame(x))
    x <- as.matrix(x)
  ii <- seq_len(n)
  W.k <- function(X, kk) {
    clus <- if (kk > 1)
      FUNcluster(X, kk, ...)$cluster
    else rep.int(1L, nrow(X))
    list("w.k"=0.5 * sum(vapply(split(ii, clus), function(I) {
      xs <- X[I, , drop = FALSE]
      sum(distance(xs,gapdis)^d.power/nrow(xs))
    }, 0)),"clus"=clus)
  }
  logW <- E.logW <- SE.sim <- numeric(K)
  if (verbose)
    cat("Clustering k = ", K, "): .. ",
        sep = "")
  w.k=W.k(x, K)
  logW[K] <- log(w.k[[1]])
  if (verbose)
    cat("done\n")
  spaceH0 <- match.arg(spaceH0)
  xs <- scale(x, center = TRUE, scale = FALSE)
  m.x <- rep(attr(xs, "scaled:center"), each = n)
  switch(spaceH0, scaledPCA = {
    V.sx <- svd(xs, nu = 0)$v
    xs <- xs %*% V.sx
  }, original = {
  }, stop("invalid 'spaceH0':", spaceH0))
  rng.x1 <- apply(xs, 2L, range)
  logWks <- matrix(0, B, K)
  #if (verbose)
  #cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n",
  #sep = "")
  for (b in 1:B) {
    z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1],
                                                 max = M[2]), nn = n)
    z <- switch(spaceH0, scaledPCA = tcrossprod(z1, V.sx),
                original = z1) + m.x
    w.k2 = W.k(z, K)
    logWks[b, K] <- log(w.k2[[1]])

    if (verbose)
      cat(".", if (b%%50 == 0)
        paste(b, "\n"))
  }

  if (verbose && (B%%50 != 0))
    cat("", B, "\n")
  E.logW <- colMeans(logWks)
  SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))

  list("gap"=sum(E.logW) - sum(logW), "sdk"=sum(SE.sim), "k"=K, "logW"= sum(logW), "ElogW"=sum(E.logW), "clus"=w.k[[2]])
  #structure(class = "clusGap", list(Tab = cbind(logW, E.logW,
  #                                              gap = E.logW - logW, SE.sim), call = cl., spaceH0 = spaceH0,
  #                                  n = n, B = B, FUNcluster = FUNcluster))
}

#The following function takes in a wide matrix for two or more longitudinal variables merged
#together and returns the numeric vector of cluster partitions for each ID provided in the same order
#it takes in a wide matrix mat
#the number of clusters k
#the IDs associated with each cluster id
#the timevector associated with the visits for each ID
#the partitioning algorithm to use to obtain the cluster partitions: kmeans, kmeans2, kmedians, kmedians2, pam, kml, Kmedians, kmlmedian
kk <- NULL
obtainclusvectornew <- function(mat, k, id, timevector, algo="kml"){

  #Find min and max time intervals
  mintime <- min(timevector)
  maxtime <- max(timevector)

  #Create a vector with just subset of time vector
  #ind <- 1:length(timevector)
  #Subset matrix to timevector
  #mat1 <- mat[, ind]
  mat1 <- mat
  if(length(timevector)==1){
    mat1 <- t(t(mat1))
  }

  #If we only need 1 cluster, function returns a vector of 1s
  if(k==1){
    cluster <- rep(1, nrow(mat1))
  }

  #If number of requested clusters is greater than 1 then invoke clustering functions
  else if(k>1){
    if(algo=="kml"){
      #If only 1 time point, invoke kmeans and store the cluster vector
      if(length(timevector)==1){
        cluster <- kmeans2(mat1, k)
      }
      #If more than 1 time pointm invoke kml instead and store the cluster vector
      else if(length(timevector)>1){
        cluster <- relabel_numbers(KML2d(mat=mat1, id=id, k=k))
      }
    }
    else if(algo=="kmeans"){
      #If only 1 time point, invoke kmeans and store the cluster vector
      if(length(timevector)>=1){
        cluster <- relabel_numbers(kmeans(mat1, k)$cluster)
      }
    }
    else if(algo=="kmeans2"){
      #If only 1 time point, invoke kmeans and store the cluster vector
      if(length(timevector)>=1){
        cluster <- kmeans2(mat1, k)
      }
    }
    else if(algo=="pam"){
      if(length(timevector)>=1){
        cluster <- relabel_numbers(pam(x=mat1, k=k, metric="euclidean", cluster.only=TRUE))
      }
    }
    else if(algo=="kmlmedian"){
      if(length(timevector)==1){
        #cluster <- pam(x=mat1, k=k, metric="euclidean", cluster.only=TRUE)
        #cluster <- Kmedians(mat1, nclust=k,nit=40)$bestresult$cluster
        km <- kcca(x=mat1, k = k, family=kccaFamily("kmedians"),control=list(initcent="kmeanspp"), simple=TRUE)
        cluster <-  clusters(km)
      }
      else if(length(timevector)>1){
        cluster <-KMLmedian2d(mat=mat1, id=id, k=k)
      }
    }
    else if(algo=="Kmedians"){
      if(length(timevector)>=1){
        #cluster <- pam(x=mat1, k=k, metric="euclidean", cluster.only=TRUE)
        cluster <- relabel_numbers(Kmedians(mat1, nclust=k,nit=40)$bestresult$cluster)
      }
    }
    else if(algo=="kmedians"){
      if(length(timevector)>=1){
        #cluster <- pam(x=mat1, k=k, metric="euclidean", cluster.only=TRUE)
        km <- kcca(x=mat1, k = k, family=kccaFamily("kmedians"),control=list(initcent="kmeanspp"), simple=TRUE)
        cluster <- relabel_numbers(clusters(km))
      }
    }
    else if(algo=="kmedians2"){
      if(length(timevector)>=1){
        #cluster <- pam(x=mat1, k=k, metric="euclidean", cluster.only=TRUE)
        cluster <-  kmedians2(x=mat1,k=k)
      }
    }
  }
  #Return cluster vector
  return(list("cluster"=cluster))
}




#The following function splits a matrix into two groups
#and returns a list with matrix, ID, and previous clusters for complete data as the first group
#and a list of matrix, ID, and previous clusters for missing data as the 2nd group.
#This function takes in as parameters
#mat: A wide matrix
#id: a vector of IDs associated with each row of matrix
#column: The maximum column number to use in the provided matrix.
#cluster: the vector of cluster partition in the order of unique rows in the matrix
split_matrix3d <- function(mat, id, column, cluster) {

  check_missing_rows <- function(mat_list, col_index) {
    all_missing <- rep(FALSE, nrow(mat_list[[1]]))
    for (mat in mat_list) {
      all_missing <- all_missing | is.na(mat[, col_index])
    }
    return(all_missing)
  }

  # get the row indices where there are missing values in the provided column
  missing_rows <- which(check_missing_rows(mat, column))
  non_missing_rows <- which(!check_missing_rows(mat, column))
  # create two matrices: one with complete data and rows with intermittent missing values
  # and another with the remaining rows (i.e. rows missing the provided column)

  complete_missing <- function(x){
    complete_mat <- x[non_missing_rows, ]
    complete_id <- id[non_missing_rows]
    complete_clus <- cluster[non_missing_rows]
    missing_mat <- x[missing_rows, ]
    missing_id <- id[missing_rows]
    missing_clus <- cluster[missing_rows]

    list(complete_mat = complete_mat, complete_id = complete_id,complete_clus = complete_clus,
         missing_mat = missing_mat, missing_id = missing_id, missing_clus = missing_clus)
  }

  listcompletemiss <- lapply(X=mat, FUN=complete_missing)

  # return the two matrices and their corresponding ids as a list
  return(listcompletemiss)
}



#This function identifies if a matrix has intermittent missing values
#It takes in as parameters
#mat: the Matrix of interest
#column: up to what column to check for intermittent missing values
#note that the last column will not be considered intermittent if missing
has_intermittent_missing <- function(mat, column) {
  apply(mat[,1:column], 1, function(row) {
    any(is.na(row)) && !is.na(row[column])
  })
}

#The following function imputes intermittent missing data from a provided list
#It takes in as parameters
#matrix_list:  a list matrices of the same dimensions containing the missing values needing imputation
#imp: the imputation methods ("copyMean.locf", and other methods detailed in kML package document)
impute_intermittent_missing3d <- function(matrix_list, imp="copyMean.locf") {
  if (length(matrix_list) == 0) {
    stop("Input list is empty.")
  }

  # Get the number of rows and columns in the matrices
  num_rows <- nrow(matrix_list[[1]])
  num_cols <- ncol(matrix_list[[1]])

  # Initialize a boolean matrix to all FALSE
  result_vector<- rep(FALSE,length(matrix_list))

  # Iterate through each matrix in the list
  for (i in 1:length(result_vector)) {
    # Check for intermittent missing values (NA) in each column
    result_vector[i]  <- any(has_intermittent_missing(matrix_list[[i]], num_cols))
  }

  if (sum(result_vector) >= 1) {
    if(num_rows>3){
      for (i in 1:length(result_vector)) {
        # Check for intermittent missing values (NA) in each column
        matrix_list[[i]] <- imputation(matrix_list[[i]], method=imp, lowerBound="globalMin",upperBound="globalMax")
      }
      #submat[, 1:column] <- imputation(submat[, 1:column], method=imp, lowerBound="globalMin",upperBound="globalMax")
    }
    else{
      cat("Less than 10 Matrix numcol=", num_rows, "\n")
      for (i in 1:length(result_vector)) {
        row_means <- apply(matrix_list[[i]], 1, mean, na.rm = TRUE)

        imps <- matrix(row_means[rep(1:num_rows, each = num_cols)], nrow=num_rows, ncol=num_cols)

        matrix_list[[i]][is.na(matrix_list[[i]]) | is.nan(matrix_list[[i]])] <- imps[is.na(matrix_list[[i]]) | is.nan(matrix_list[[i]])]
        #submat[, 1:column] <- imputation(submat[, 1:column], method='linearInterpol.locf', lowerBound=0,upperBound="globalMax")
      }
    }

  }

  return(matrix_list)
}


#The following function imputes intermittent missing data from a provided matrix
#It takes in as parameters
#mat:  a matrix containing the missing values needing imputation
#id: a vector of IDs associated with each row of the matrix
#group: a vector identifying which group each row belongs to
#column: The Column of the matrix up to which to impute.
#imp: the imputation methods ("copyMean.locf", and other methods detailed in kML package document)
impute_column_by_group3d <- function(mat, id, group, column, nlimit, imp="copyMean.locf"){
  makemat2 <- function(X, numcol=ncol(mat)){
    return(matrix(X, ncol=numcol))
  }

  # Split the matrix into groups based on the grouping variable
  split_mat <-lapply(split(mat, group), makemat2)
  ids <- lapply(split(id, group), matrix)

  # Define a function to impute the column for each group
  impute_column <- function(submat) {
    test <- matrix(submat[, 1:column], ncol=column)
    num_rows <- nrow(test)
    num_cols <- ncol(test)
    # Check if the column contains any missing values
    if (any(is.na(submat[, 1:column]))) {
      if(nrow(test)>nlimit){
        #cat("Greater then 10 Matrix \n")
        #print(nrow(test)>3)
        # If the column contains missing values, impute using the impute function
        submat[, 1:column] <- imputation(submat[, 1:column], method=imp, lowerBound="globalMin",upperBound="globalMax")
        #print(submat[, 1:column] )
      }
      else{
        cat("Less than 10 Matrix numcol=", nrow(test), "\n")
        #row_means <- apply(test, 1, mean, na.rm = TRUE)
        #test[is.na(test) | is.nan(test)] <- row_means[rep(1:nrow(test), each = ncol(test))]
        row_means <- apply(test, 1, mean, na.rm = TRUE)
        imps <- matrix(row_means[rep(1:num_rows, each = num_cols)], nrow=num_rows, ncol=num_cols)
        test[is.na(test) | is.nan(test)] <- imps[is.na(test) | is.nan(test)]

        #submat[, 1:column] <- imputation(submat[, 1:column], method='linearInterpol.locf', lowerBound=0,upperBound="globalMax")
        submat[, 1:column] <- test
      }

    }

    # Return the updated submatrix
    return(submat)
  }

  # Apply the impute_column function to each group
  updated_groups <- lapply(split_mat, impute_column)

  # Combine the updated groups back into a single matrix
  updated_mat <- do.call(rbind, updated_groups)
  ids2 <- do.call(rbind, ids)


  # Sort the matrix back into the original order of the ids
  newframe <- data.frame(updated_mat, ids2)
  newframe <- arrange(newframe, ids2)
  #print(newframe)
  # Return the updated matrix
  #cat("End of imputation \n")
  return(as.matrix(newframe[1:ncol(mat)]))
}







# The following function finds the optimal number of clusters for a single array of variables of the same dimensions
#It takes in as parameters:

#' Clustering by Spectral decomposition
#'
#' Performs clustering using spectral decomposition for complete data
#' @param arr An array with the continuous variables of interest
#' @param id A vector of IDs associated with each individual in the order of the matrix
#' @param timevector The visits or time points associated with each column for each variable. Example: timevector = C(0,1,2,3)
#' @param B the number of uniform samples to be drawn for the GAP statistics
#' @param cor1 A square matrix with the correlations between the outcomes
#' @param var_a a boolean indicator of whether to account for the variances in the diagonal matrix
#' @param cor_0 a boolean indicator of whether to account for the correlations between the outcomes (obtained from cor1)
#' @param cor_e a boolean indicator of whether to account for the temporal correlations for the repeated measures
#' @param method the version of the Tibshirani Gap statistic to use "original" or "scaledPCA"
#' @param algo what algorithm to use to partition clusters: "kmeans", "kmedians", "kmeans2", "kmedians2", "pam", etc.... Note that kmeans2 and kmedians2 performs kmeans and kmedians 20 times with different starting points and selects the cluster partition with the best objective function.
#' @param dis the distance to use ("euclidean","manhattan" , "mahalanobis","minkowski", "maximum", "canberra", "binary","mahalanobis-ZCA", "mahalanobis-ZCA-cor","mahalanobis-PCA", "mahalanobis-PCA-cor", "mahalanobis-Cholesky")
#' @param gapdis the distance to use in the gap statistic function. Same options as "dis". Keep as euclidean if "dis" is already non-euclidean.
#' @param gapmethod which method to use to obtain optimal clusters ("Tibs2001SEmax", "globalSEmax", "fixed", "RamdomForrest")
#' @param covmethod Correlation method. "pearson" or "spearman".
#' @param standard Should each variable be globally standardized recursively before running algorithm at each iteration
#' @param fixedk The number of clusters to produce if you would like to fix the number of clusters. Otherwise, keep as FALSE
#' #'
#' @return returns a list of data frames containing ID, repeated variable 1, cluster partitions at each time point.
#' @export
#'
#' @examples find_optimal_clusters_single3d(arr=myarray, id=idvector, B=500, cor1=cor(myarray[[1]][,1],myarray[[2]][,1]) method="scaledPCA", algo="kmeans", dis="Cholesky", gapdis="euclidean", gapmethod="Tibs2001SEmax", covmethod="pearson", var_a = TRUE, cor_o=TRUE, cor_t=TRUE, novarcor=FALSE,standard=FALSE, fixmat=FALSE)
find_optimal_clusters_single3d <- function(arr, id, timevector, B, cor1, var_a = TRUE, cor_o=TRUE, cor_t=TRUE,
                                           method="scaledPCA", algo="kmeans",dis="Cholesky",gapdis="euclidean",
                                           gapmethod="Tibs2001SEmax", covmethod="pearson",standard=TRUE, fixedk=FALSE) {
  #browser()

  ind <- 1:length(timevector)

  #Subset matrix to timevector
  arr1 <- lapply(arr, function(x) x[, ind])

  if(standard==TRUE){
    for(r in 1:length(arr1)){
      arr1[[r]] <- (arr1[[r]] - mean(arr1[[r]],na.rm=TRUE)) / sd(arr1[[r]],na.rm=TRUE)
    }
  }

  if(length(timevector)==1){
    #mat1 <- t(t(mat1))
    arr1 <- lapply(arr1, function(x) t(t(x)))
  }

  mat1 <- listtomat(arr1,"wide")

  if(dis!="euclidean"){
    if(length(timevector)==1){
      mat1 <- whiten_cov(X=mat1, S=FALSE, center = FALSE, method = dis)
    }
    if(length(timevector)>1){
      cov2 <- varcov(arr1, cor1=cor1, var_a = var_a, cor_o=cor_o, cor_t=cor_t, covmethod=covmethod)
      mat1 <- whiten_cov(X=mat1, S=cov2, center = FALSE, method = dis)
    }
  }

  #if(max(timevector)>=4){print(mat1)}
  obtainclusvector2 <- function(x,k){
    obtainclusvectornew(mat=x, k=k, id=id, timevector=min(timevector):max(timevector), algo=algo)
  }
  funck <- function(num){
    getgap <- gapfunction(x=mat1, FUNcluster=obtainclusvector2, id=id, K=num, timevector=timevector, spaceH0=method, B=B, gapdis=gapdis)
    return(getgap)
  }

  if(gapmethod =="Tibs2001SEmax"){
    # Set an initial value for K
    K <- 1

    gapa <- funck(num=K)
    gapb <- funck(num=K+1)
    #print(gapa)
    #print(gapb)
    # # Loop until the condition is met
    while (gapa$gap < (gapb$gap - gapb$sdk)) {
      cat("K=", K, " Gap(K)=", gapa$gap, " Gap(K+1)=", gapb$gap, " sdk=", gapb$sdk,"\n")
      K <- K + 1

      gapa <- gapb
      gapb <- funck(K+1)

    }
    cat("K=", K, " Gap(K)=", gapa$gap, " Gap(K+1)=", gapb$gap, " sdk=", gapb$sdk,"\n")
    return(gapa)
  }
  else if(gapmethod =="globalSEmax"){
    # Set an initial value for K
    k_values <- 1:5

    gaps_sdk <- lapply(k_values, function(k) {funck(num=k)})
    gaps <- sapply(k_values, function(k) {gaps_sdk[[k]]$gap})
    sdks <- sapply(k_values, function(k) {gaps_sdk[[k]]$sdk})

    globalSEmax = {
      nc <- which.max(gaps)
      if (any(mp <- gaps[seq_len(nc - 1)] >= gaps[nc] - sdks[nc])) which(mp)[1] else nc
    }
    cat("K=", gaps_sdk[[globalSEmax]]$k, " Gap(K)=", gaps_sdk[[globalSEmax]]$gap, "\n")
    return(gaps_sdk[[globalSEmax]])

  }
  else if(gapmethod =="fixed"){
    # Set an initial value for K
    gaps_sdk <- funck(num=fixedk)
    cat("K=", gaps_sdk$k, " Gap(K)=", gaps_sdk$gap, "\n")
    return(gaps_sdk)
  }

  else if(gapmethod =="randomforrest"){

    kmeans_bootstrap <- function(data, indices, k) {

      sample_data <- data[indices, ]
      id2 <- id[indices]
      #arr2 <- list(sample_data[,1:6],sample_data[,7:12])


      obtainclusvector2 <- function(x,k){
        obtainclusvectornew(mat=x, k=k, id=id2, timevector=timevector, algo=algo)
      }
      funck <- function(num){
        getgap <- gapfunction(x=sample_data, FUNcluster=obtainclusvector2, id=id2, K=num, timevector=timevector, spaceH0=method, B=B, gapdis=gapdis)
        return(getgap)
      }

      result <- funck(num=k)
      return(result$gap)
    }

    # Set an initial value for K
    calculate_mean_wss <- function(data, k, num_bootstrap_samples) {
      # Set seed for reproducibility
      boot_results <- boot(data, kmeans_bootstrap, R = num_bootstrap_samples, k = k)
      return(mean(boot_results$t))
    }

    # Choose a range of k values
    k_values <- 1:6  # You can adjust this range based on your problem

    # Calculate mean WSS for each k
    mean_wss_values <- sapply(k_values, function(k) {
      calculate_mean_wss(mat1, k, B)
    })

    elbow_point <- k_values[which.min(diff(mean_wss_values))]

    cat("K=", elbow_point, " Gap(K)=", mean_wss_values[elbow_point], "\n")
    return(funck(num=elbow_point))

  }
  else if(gapmethod =="randomforrest2"){

    kmeans_bootstrap <- function(data, indices, k) {

      sample_data <- data[indices, ]
      id2 <- id[indices]
      #arr2 <- list(sample_data[,1:6],sample_data[,7:12])
      n <- nrow(sample_data)
      ii <- seq_len(n)
      obtainclusvector2 <- function(x,k){
        obtainclusvectornew(mat=x, k=k, id=id2, timevector=timevector, algo=algo)
      }
      funck <- function(num){
        W.k <- function(X, num) {
          clus <- if (num > 1)
            obtainclusvector2(X, num)$cluster
          else rep.int(1L, nrow(sample_data))
          list("w.k"=0.5 * sum(vapply(split(ii, clus), function(I) {
            xs <- X[I, , drop = FALSE]
            sum(distance(xs,gapdis)^2/nrow(xs))
          }, 0)),"clus"=clus)
        }
        w.k=W.k(sample_data, num)
        return(log(w.k[[1]]))
      }

      result <- funck(num=k)
      return(result)
    }


    calculate_mean_wss <- function(data, k, num_bootstrap_samples) {
      # Set seed for reproducibility
      boot_results <- boot(data, kmeans_bootstrap, R = num_bootstrap_samples, k = k)

      mean1 <- mean(boot_results$t)
      SE.sim <- sqrt((1 + 1/num_bootstrap_samples) * var(boot_results$t))
      return(list("logw"=mean1,"sdk"=SE.sim))
    }

    # Choose a range of k values
    k_values <- 1:6  # You can adjust this range based on your problem


    # Calculate mean WSS for each k
    mean_wss_values <- lapply(k_values, function(k) {
      calculate_mean_wss(mat1, k, B)
    })
    logws <- sapply(k_values, function(k) {mean_wss_values[[k]]$logw})
    gaps <- c(0,diff(logws))
    sdks <- sapply(k_values, function(k) {mean_wss_values[[k]]$sdk})
    #print(gaps)
    #print(sdks)
    globalSEmax = {
      nc <- which.min(gaps)
      if (any(mp <- gaps[seq_len(nc - 1)] < gaps[nc] - sdks[nc])) which(mp)[1] else nc
    }

    return(funck(num=globalSEmax))

  }
}










#The following function takes in an list of matrices, an ID vector, the associated time vector, performs R-CHRM and returns a list of data frames
# containing ID, repeated outcome, cluster partitions at each time point.
# It takes in as parameters:

#' Recursive Clustering by Heterogeneity of repeated Measures R-CHRM
#'
#' Performs clustering at each visit using the R-CHRM approach
#' @param arr a list of matrices of the same dimensions. Longitudinal variable in a wide format.
#' @param id a vector of IDs associated with each individual in the order of the matrix. 'id' and matrices should be sorted according to 'id'.
#' @param timevector The visits or time points associated with each column for each variable. Example: timevector = C(0,1,2,3)
#' @param init_cols The number of visits (columns) should the algorithm start with. 1,2,3
#' @param B the number of uniform samples to be drawn for the GAP statistics
#' @param method what version of the Tibshirani Gap statistic to use "original" or "scaledPCA"
#' @param nlimit number limit to restrict further clustering if number of sample is equal to or less than the specified.
#' @param imp Method of imputation ("copyMean.locf") for intermittent missing, if present. See 'kml' package for additional methods.
#' @param algo what algorithm to use to partition clusters: "kmeans", "kmedians", "kmeans2", "kmedians2", "pam", etc.... Note that kmeans2 and kmedians2 performs kmeans and kmedians 20 times with different starting points and selects the cluster partition with the best objective function.
#' @param dis the distance to use ("euclidean","manhattan" , "mahalanobis","minkowski", "maximum", "canberra", "binary","mahalanobis-ZCA", "mahalanobis-ZCA-cor","mahalanobis-PCA", "mahalanobis-PCA-cor", "mahalanobis-Cholesky")
#' @param gapdis the distance to use in the gap statistic function. Same options as "dis". Keep as euclidean if "dis" is already non-euclidean.
#' @param gapmethod which method to use to obtain optimal clusters ("Tibs2001SEmax", "globalSEmax", "fixed", "RamdomForrest")
#' @param covmethod Correlation method. "pearson" or "spearman".
#' @param var_a a boolean indicator of whether to account for the variances in the diagonal matrix
#' @param cor_o a boolean indicator of whether to account for the correlations between the outcomes (obtained from cor1)
#' @param cor_t a boolean indicator of whether to account for the temporal correlations for the repeated measures
#' @param novarcor A matrix with the between outcome correlation. Leave as FALSE if you want to algorithm to calculate the covariance matrix.
#' @param standard Should each variable be globally standardized recursively before running algorithm at each iteration
#' @param fixmat a matrix specifying where clusters deviate. Works only when the "fixed" option is specified in the gapmethod parameter.
#'
#' @return Returns a nested list with data frames containing ID, repeated outcome in wide format, cluster partitions at each visit.
#' @export
#'
#' @examples rchrmclus(arr=myarray, id=idvector, timevector=0:3, init_cols=1, B=500, method="scaledPCA", nlimit=20, imp="copyMean.locf", algo="kmeans", dis="Cholesky", gapdis="euclidean", gapmethod="globalSEmax", covmethod="pearson", var_a = TRUE, cor_o=TRUE, cor_t=TRUE, novarcor=FALSE,standard=FALSE, fixmat=FALSE)
#' Generate a MVN dataset (2 variable) using the createdata function
#' mymat <- createdata(p1_a=.5, E_a=c(.3,.4,.6,.4) ,s2_a=c(1,2,2,1,1,1,2,2), meanvalues=matrix(c(3,3,3,6,2,2,6,2,2,6,5,2,3,3,3,7,5,5,8,4,4,9,7,2), ncol=8), sizes=c(100,100,100))
#'
#' Create a list separating the two variables
#' mylist <- list(mymat[,1:4], mymat[,5:8])
#'
#' Create a vector of IDs associated with each row. Note that the vector should be sorted.
#' myid <- c(1001:1100,2001:2100, 3001:3100)
#'
#' Invoke the recursive algorithm
#' aaa <- rchrmclus(arr=mylist, id=myid, timevector=0:3, init_cols=1, B=200, method="scaledPCA", nlimit=20, imp="copyMean.locf", algo="kmeans", dis="Cholesky", gapdis="euclidean", gapmethod="Tibs2001SEmax", covmethod="pearson", var_a = TRUE, cor_o=TRUE, cor_t=TRUE, novarcor=FALSE,standard=TRUE, fixmat=FALSE)
#'
#' table(aaa[[1]][[1]]$cluster0)
#' table(aaa[[1]][[1]]$cluster1)
#' table(aaa[[1]][[1]]$cluster2)
#' table(aaa[[1]][[1]]$cluster3)
rchrmclus <- function(arr, id, timevector, init_cols=1, B=100, method="original", nlimit=20, imp="copyMean.locf", algo="kmeans", dis="euclidean", gapdis="euclidean", gapmethod="Tibs2001SEmax",
                                         covmethod="pearson", var_a = TRUE, cor_o=TRUE, cor_t=TRUE, novarcor=FALSE,standard=TRUE, fixmat=FALSE) {
  # Initialize the number of clusters and the previous clusters
  maxtime <- max(timevector)
  mintime <- min(timevector)

  if(dis !="euclidean"){
    first_columns <- sapply(arr, function(x){x[,1]})
    #print(first_columns)
    cor1 <- cor(first_columns, method=covmethod)

    if(novarcor == TRUE){
      cor1 <- diag(rep(1, length(diag(cor1))))
    }
    #print(cor1)
  }
  else{cor1=FALSE}

  K <- 1
  prev_clusters <- rep(1, nrow(arr[[1]]))


  makemat <- function(X, numcol=length(timevector)){
    return(matrix(X, ncol=numcol))
  }

  #Find intermittent missing values
  has_intermittent_missing <- function(mat, column) {
    apply(mat[,1:column], 1, function(row) {
      any(is.na(row)) && !is.na(row[column])
    })
  }


  createlistframe <- function(x){
    thedata <- data.frame(id, x, prev_clusters)
    return(thedata)
  }
  #print(mat)
  thedata <- lapply(X=arr, FUN=createlistframe)

  time <- mintime

  if(init_cols>1){
    time <- timevector[which(1:length(timevector)==init_cols)]
    for(v in 1:length(arr)){
      thedata[[v]]$cluster <- prev_clusters
    }
  }
  #time <- 1
  # Loop until no more time points or no more clusters can be found
  while (time <= maxtime) {
    # If there are previous clusters, break down the matrix by clusters
    cat("Time ", time, "; column ",which(timevector==time),"\n")
    if (!is.null(prev_clusters)) {

      #Split Matrix by non-missing column vs missing column.
      newmats <- split_matrix3d(mat=arr, id=id, column=which(timevector==time), cluster=prev_clusters)


      #Get IDs and clusters
      complete_id <- newmats[[1]]$complete_id
      complete_clus <- newmats[[1]]$complete_clus
      id_list <- lapply(split(complete_id, complete_clus), matrix)

      #Initialize length of lists
      complete_mat <- vector(mode = "list", length = length(newmats))
      complete_mat2 <- vector(mode = "list", length = length(newmats))
      mat_list <- vector(mode = "list", length = length(newmats))

      for(v in 1: length(newmats)){
        complete_mat[[v]] <- newmats[[v]]$complete_mat
        complete_mat2[[v]] <- newmats[[v]]$complete_mat

        if(which(timevector==time)>1){
          complete_mat[[v]] <- impute_column_by_group3d(mat=complete_mat[[v]], id=complete_id, group=complete_clus, column=which(timevector==time), nlimit=nlimit, imp="copyMean.locf")
        }

        mat_list[[v]] <-lapply(split(complete_mat[[v]], complete_clus), makemat)

      }

      #mat_list[[2]]
      allclusterlist2 <- list(NULL)
      clusteronly2 <- list(NULL)
      framelist <- list(NULL)


      #Combine list of data matrix into one list

      # Loop through the rest of the matrices in the list and find the optimal number of clusters
      for (i in 1:length(mat_list[[1]])) {
        #set.seed(123)

        obtain_k_fixed_cluster <- function(fixmat){
          aa <- apply(fixmat, 2, max) - apply(fixmat, 2, min)
          #diff(aa)+1
          #i = 4
          #timevector = 0:3
          #time=3
          #gapmethod="fixed"

          if(time==mintime){
            fixedk <- max(fixmat[,1])
          }
          else if(time>mintime){
            fixmat1 <- fixmat[,1:length(mintime:time)]
            fixmat1 <- unique(fixmat1)
            whichclus <- fixmat1[,which(timevector==time)-1]==i
            #(fixedk <-  min(sum(whichclus), min(fixmat1[whichclus,which(timevector==time)])))
            fixedk <-  min(sum(whichclus), min((diff(aa)+1)[which(timevector==time)-1]))
          }

          return(fixedk)
        }

        if(is.matrix(fixmat)){
          fixedk <- obtain_k_fixed_cluster(fixmat)
        }
        else{fixedk=FALSE}

        obtainclusvector2 <- function(x,k){
          obtainclusvectornew(mat=x, k=k, id=id, timevector=min(timevector):time, algo=algo)
        }

        mat_listb <- vector(mode = "list", length = length(arr))
        for(v in 1:length(arr)){
          mat_listb[[v]] <-  mat_list[[v]][[i]]
        }

        #check for sample size limit
        #if(length(id_list[[i]]) == 1) {
        #  allclusterlist2[[i]] <-  list("gap"=NULL, "sdk"=NULL, "k"=1, "logW"= NULL, "ElogW"=NULL, "clus"=1)
        #}
        if(length(id_list[[i]]) >= 1 & length(id_list[[i]]) < nlimit) {
          cat("Right before clustering for low counts \n")

          #mat_listc <- listtomat(mat_listb, "wide")
          #print(mat_listc)
          #allclusterlist2[[i]] <- gapfunction(mat=mat_listc, obtainclusvector2, id=id_list[[i]], K=1, timevector=mintime:time, method=method, gapdis=gapdis, B=B)
          allclusterlist2[[i]] <-  list("gap"=NULL, "sdk"=NULL, "k"=1, "logW"= NULL, "ElogW"=NULL, "clus"=1)
          cat("Right after clustering for low counts \n")
        }
        else if(length(id_list[[i]]) >= nlimit){
          cat("------------------------------------ \n")
          # print(id_list[[i]])
          allclusterlist2[[i]] <- find_optimal_clusters_single3d(arr=mat_listb, id=id_list[[i]], timevector=mintime:time, B=B,
                                                                 cor1=cor1, var_a=var_a, cor_o=cor_o, cor_t=cor_t, method=method, algo=algo, dis=dis, gapdis=gapdis,
                                                                 gapmethod=gapmethod,covmethod=covmethod,standard=standard, fixedk=fixedk)
        }

        clusteronly2[[i]] <- matrix(allclusterlist2[[i]]$clus+(100*i)+(1000*time),ncol=1)
        framelist[[i]] <- data.frame("id"=id_list[[i]], "cluster"=clusteronly2[[i]])
        #print(framelist[[i]])
      }


      #################Impute negative clusters for missing values in column.
      missing_mat <- vector(mode = "list", length = length(newmats))
      for(v in 1: length(newmats)){
        missing_mat[[v]] <- newmats[[v]]$missing_mat
      }
      #missing_mat <- newmats$missing_mat

      missing_id <- newmats[[1]]$missing_id
      missing_clus <- newmats[[1]]$missing_clus
      #framelist[[i+1]] <- data.frame("id"=missing_id, "cluster"=-1*((1000*time)+abs(missing_clus)))
      #framelist[[i+1]] <- data.frame("id"=missing_id, "cluster"=-1*((1000*time)+as.numeric(factor(missing_clus))))
      framelist[[i+1]] <- data.frame("id"=missing_id, "cluster"=-1*abs(missing_clus))


      # Combine the two matrices by rows
      mat_combined <- vector(mode = "list", length = length(arr))
      for(v in 1:length(arr)){
        mat_combined[[v]] <- rbind(complete_mat[[v]], missing_mat[[v]])
      }

      #mat_combined <- rbind(complete_mat, missing_mat)
      #cat("We done")
      # Combine the two vectors by concatenation
      id_combined <- c(complete_id, missing_id)
      # Get the order of the sorted combined id
      id_order <- order(id_combined)
      # Sort the combined matrix based on the sorted order of the combined id
      for(v in 1:length(arr)){
        arr[[v]] <- mat_combined[[v]][id_order,]
      }

      #mat <- mat_combined[id_order,]
      # print(nrow(mat))
      # Sort the combined id vector based on the sorted order of the combined id
      id <- id_combined[id_order]

      # Identify rows with intermittent missing values
      #intermittent_rows <- has_intermittent_missing(complete_mat2, which(timevector==time))
      has_intermittent_missing2 <- function(x){
        any(has_intermittent_missing(x, which(timevector==time)))
      }
      #cat("We done")
      #print(complete_mat2)

      if(time>mintime){
        intermittent_rows <- lapply(complete_mat2, has_intermittent_missing2)

        if (any(unlist(intermittent_rows))) {
          cat("Redo, go back to time ",time - 1)
          eval(parse(text=paste("thedata$cluster",time-1,"<- NULL", sep="")))
          time <- 1
          prev_clusters <- rep(1, length(prev_clusters))
          next
        }
      }

    }

    # Update the previous clusters and continue to the next iteration
    newmatframe <- data.frame(do.call(rbind, framelist))
    newids <- arrange(newmatframe, id)

    newids$prev_clusters <- ifelse(newids$cluster < 0, newids$cluster, as.numeric(factor(newids$cluster)))
    #eval(parse(text=paste("newids$prev_clusters",time,"<-newids$prev_clusters", sep="")))
    eval(parse(text=paste("newids$cluster",time,"<-newids$cluster", sep="")))

    prevprev_clusters <- prev_clusters

    for(v in 1:length(arr)){
      thedata[[v]] <- dplyr::right_join(thedata[[v]], newids, by="id")
    }
    #print(thedata[[1]])

    prev_clusters <- newids$prev_clusters
    id <- newids$id
    #mat <- as.matrix(newmatframe[,1:length(timevector)])

    time <- time + 1

  }
  return(list(thedata, arr))
}



#' Create multivariate normal data jointly using spectral decomposition.
#'
#' Data are generated using AR(1) correlations structure.
#' @param p1_a the between outcome correlation value
#' @param E_a the within outcome correlation value
#' @param s2_a a vector of variances for each visit for both outcome combined.
#' @param meanvalues Matrix containing the means for each visit cols for each cluster rows. this should have the same number of cols as the length of s2_a
#' @param sizes a vector with the sizes for each cluster. This should have the same number of rows as the matrix meanvalues
#'
#' @return a matrix with data geneterated using
#' @export
#'
#' @examples createdata(p1_a=.5, E_a=c(.3,.4,.6,.4) ,s2_a=c(1,2,2,1,1,1,2,2), meanvalues=matrix(c(3,3,3,6,2,2,6,2,2,6,5,2,3,3,3,7,5,5,8,4,4,9,7,2), ncol=8), sizes=c(100,100,100))
createdata <- function(p1_a, E_a ,s2_a, meanvalues, sizes){
  #source("./functions3d Multivariate v4.R")

  P1<- diag(c(1-p1_a,1-p1_a))+p1_a
  E <- ar1_cor(4,E_a)
  p2 <- P1 %x% E
  A <- diag(sqrt(s2_a))
  s1 <- A^(1/2)%*%p2%*%A^(1/2)

  M <- matrix(NA, nrow=sum(sizes), ncol=ncol(meanvalues))
  ind <- 0
  for(i in 1:nrow(meanvalues)){
    if(i==1){
      #cat("Here")
      M[1:sizes[i],] <- mvrnorm(n= sizes[i], mu= meanvalues[i,], Sigma=s1)
    }
    else{
      #cat("Here2 ")
      M[(ind+1):(ind+sizes[i]),] <- mvrnorm(n= sizes[i], mu= meanvalues[i,], Sigma=s1)
    }
    ind <- ind +sizes[i]
  }
  return(M)

}



mean_n <- function(frame, var1, clusvar){
  meanLED0 <- frame %>%
    group_by(eval(parse(text=paste(clusvar, sep="")))) %>%
    summarise_at(vars(var1), list(mean=mean, median=median), na.rm=T)
  meanLED0$count <- (frame %>%
                       group_by(eval(parse(text=paste(clusvar, sep="")))) %>%
                       summarise(count=n(), na.rm=T))$count
  meanLED0$meanb <- meanLED0$mean
  meanLED0$medianb <- meanLED0$median
  eval(parse(text=paste("meanLED0$",clusvar,"=meanLED0[[1]]", sep="")))
  eval(parse(text=paste("meanLED0[(meanLED0$",clusvar,"<0),]$mean<- NA", sep="")))
  eval(parse(text=paste("meanLED0[(meanLED0$",clusvar,"<0),]$median<- NA", sep="")))
  return(meanLED0)

}


meanclusterplot3d <- function(b, varnum, clus=c("cluster.x", "cluster1.y", "cluster2.y", "cluster3"),
                              colorclus="cluster3", imp=FALSE, title, cols2){
  c <- b[[1]][[varnum]]
  #c <- b[[1]]
  #d <- as.data.frame(b[[2]])
  d <- as.data.frame(b[[2]][[varnum]])
  names(d) <- cols2
  #print(names(d))

  if(imp==TRUE) {c[,2:5] <- d}

  names(c)[1:5] <- c("id", cols2)
  #print(names(c[,2:5]))

  meanLED0 <- mean_n(frame=c, var1=cols2[1], clusvar=clus[1])
  meanLED1 <- mean_n(frame=c, var1=cols2[2], clusvar=clus[2])
  meanLED2 <- mean_n(frame=c, var1=cols2[3], clusvar=clus[3])
  meanLED3 <- mean_n(frame=c, var1=cols2[4], clusvar=clus[4])

  #print(meanLED3)
  meanLED0a <- pivot_longer(as.data.frame(meanLED0[!is.nan(meanLED0$mean),]), cols= c('mean'), names_to = "year2", values_to = "mean0")
  meanLED0a$year <- 0
  meanLED1a <- pivot_longer(as.data.frame(meanLED1[!is.nan(meanLED1$mean),]), cols= c('mean'), names_to = "year2", values_to = "mean1")
  meanLED1a$year <- 1
  meanLED2a <- pivot_longer(as.data.frame(meanLED2[!is.nan(meanLED2$mean),]), cols= c('mean'), names_to = "year2", values_to = "mean2")
  meanLED2a$year <- 2
  meanLED3a <- pivot_longer(as.data.frame(meanLED3[!is.nan(meanLED3$mean),]), cols= c('mean'), names_to = "year2", values_to = "mean3")
  meanLED3a$year <- 3



  #finaldata3 <- pivot_longer(as.data.frame(finaldata2), cols= c('0','1', '2', '3', '4', '5'), names_to = "year", values_to = "LED")
  finaldata3 <- pivot_longer(as.data.frame(c), cols= cols2, names_to = "year", values_to = "Var")
  finaldata3$year <- as.numeric(gsub("X", "",finaldata3$year))

  finaldata3 <- left_join(x=finaldata3, y=meanLED0a, by=c('year', clus[1]))
  finaldata3 <- left_join(x=finaldata3, y=meanLED1a, by=c('year', clus[2]))
  finaldata3 <- left_join(x=finaldata3, y=meanLED2a, by=c('year', clus[3]))
  finaldata3 <- left_join(x=finaldata3, y=meanLED3a, by=c('year', clus[4]))


  ledmean <- apply(data.frame(finaldata3$mean0, finaldata3$mean1, finaldata3$mean2,finaldata3$mean3),1,min, na.rm=TRUE)
  ledmean[is.infinite(ledmean)] <- NA

  # print(ledmean)
  finaldata4 <- mutate(finaldata3, mean=ledmean)

  print(ggplot(data = finaldata4, aes(x = year, y = Var, group = id)) +
          geom_line() +
          theme_bw(base_size = 17) +
          geom_line(data = finaldata4, aes(x = year, y = mean, group = eval(parse(text=paste(clus[4]))), colour=factor(eval(parse(text=paste(colorclus)))), linewidth=2)) +
          theme(legend.position = "none") +
          labs(title = paste(title,"Over Time", sep=" "), subtitle = "Spaghetti Plot",
               x = "Yrs.", y = title, colour="Cluster",
               caption = paste0("nsub = ", length(unique(finaldata4$id)),
                                ", nobs = ", nrow(finaldata4))))

  return((finaldata4))
}



meanclusterplot3d5 <- function(b, varnum, clus=c("cluster.x", "cluster1.y", "cluster2.y", "cluster3"),
                               colorclus="cluster3", imp=TRUE, title, cols2){
  c <- b[[1]][[varnum]]
  #c <- b[[1]]
  #d <- as.data.frame(b[[2]])
  d <- as.data.frame(b[[2]][[varnum]])
  #names(d) <- cols2
  #print(names(d))

  if(imp==TRUE) {c[,2:7] <- d}

  names(c)[1:7] <- c("id", cols2)
  #print(names(c[,2:7]))

  #print(head(c))

  meanLED0 <- mean_n(frame=c, var1=cols2[1], clusvar=clus[1])
  meanLED1 <- mean_n(frame=c, var1=cols2[2], clusvar=clus[2])
  meanLED2 <- mean_n(frame=c, var1=cols2[3], clusvar=clus[3])
  meanLED3 <- mean_n(frame=c, var1=cols2[4], clusvar=clus[4])
  meanLED4 <- mean_n(frame=c, var1=cols2[5], clusvar=clus[5])
  meanLED5 <- mean_n(frame=c, var1=cols2[6], clusvar=clus[6])

  #print(meanLED3)
  meanLED0a <- pivot_longer(as.data.frame(meanLED0[!is.nan(meanLED0$mean) | !is.na(meanLED0$mean),]), cols= c('mean'), names_to = "year2", values_to = "mean0")
  meanLED0a$year <- 0
  meanLED1a <- pivot_longer(as.data.frame(meanLED1[!is.nan(meanLED1$mean) | !is.na(meanLED1$mean),]), cols= c('mean'), names_to = "year2", values_to = "mean1")
  meanLED1a$year <- 1
  meanLED2a <- pivot_longer(as.data.frame(meanLED2[!is.nan(meanLED2$mean) | !is.na(meanLED2$mean),]), cols= c('mean'), names_to = "year2", values_to = "mean2")
  meanLED2a$year <- 2
  meanLED3a <- pivot_longer(as.data.frame(meanLED3[!is.nan(meanLED3$mean) | !is.na(meanLED3$mean),]), cols= c('mean'), names_to = "year2", values_to = "mean3")
  meanLED3a$year <- 3
  meanLED4a <- pivot_longer(as.data.frame(meanLED4[!is.nan(meanLED4$mean) | !is.na(meanLED4$mean),]), cols= c('mean'), names_to = "year2", values_to = "mean4")
  meanLED4a$year <- 4
  meanLED5a <- pivot_longer(as.data.frame(meanLED5[!is.nan(meanLED5$mean) | !is.na(meanLED5$mean),]), cols= c('mean'), names_to = "year2", values_to = "mean5")
  meanLED5a$year <- 5


  #finaldata3 <- pivot_longer(as.data.frame(finaldata2), cols= c('0','1', '2', '3', '4', '5'), names_to = "year", values_to = "LED")
  finaldata3 <- pivot_longer(as.data.frame(c), cols= cols2, names_to = "year", values_to = "Var")
  finaldata3$year <- as.numeric(gsub("X", "",finaldata3$year))

  finaldata3 <- left_join(x=finaldata3, y=meanLED0a, by=c('year', clus[1]))
  finaldata3 <- left_join(x=finaldata3, y=meanLED1a, by=c('year', clus[2]))
  finaldata3 <- left_join(x=finaldata3, y=meanLED2a, by=c('year', clus[3]))
  finaldata3 <- left_join(x=finaldata3, y=meanLED3a, by=c('year', clus[4]))
  finaldata3 <- left_join(x=finaldata3, y=meanLED4a, by=c('year', clus[5]))
  finaldata3 <- left_join(x=finaldata3, y=meanLED5a, by=c('year', clus[6]))

  suppressWarnings({
    ledmean <- apply(data.frame(finaldata3$mean0, finaldata3$mean1, finaldata3$mean2,finaldata3$mean3,finaldata3$mean4,finaldata3$mean5),1,min, na.rm=TRUE)
    ledmean[is.infinite(ledmean)] <- NA
  })
  # print(ledmean)
  finaldata4 <- mutate(finaldata3, mean=ledmean)
  #finaldata4 <- subset(finaldata4, !is.na(eval(parse(text=paste(clus)))))
  # print(ggplot(data = finaldata4, aes(x = year, y = Var, group = id)) +
  #         geom_line(na.rm=TRUE) +
  #         theme_bw(base_size = 17) +
  #         geom_line(data = finaldata4, linewidth=2, aes(x = year, y = mean, group = eval(parse(text=paste(clus[6]))), colour=factor(eval(parse(text=paste(colorclus)))), linewidth=2)) +
  #         theme(legend.position = "none") +
  #         labs(title = paste(title,"Over Time", sep=" "), subtitle = "Spaghetti Plot",
  #              x = "Yrs.", y = title, colour="Cluster",
  #              caption = paste0("nsub = ", length(unique(finaldata4$id)),
  #                               ", nobs = ", nrow(finaldata4))))

  #finaldata4 <- subset(finaldata4, !is.na(Var))
  finaldata5 <- finaldata4[,c("id","year","cluster5","mean")]

  finaldata5 <- subset(finaldata5, eval(parse(text=paste(clus[6]))) > 0)


  print(ggplot(data = finaldata4, aes(x = year, y = Var, group = id)) +
          geom_line() +
          theme_bw(base_size = 17) +
          geom_line(data = finaldata5, aes(x = year, y = mean, group = eval(parse(text=paste(clus[6]))), colour=factor(eval(parse(text=paste(colorclus)))), linewidth=2)) +
          theme(legend.position = "none") +
          labs(title = paste(title,"Over Time", sep=" "), subtitle = "Spaghetti Plot",
               x = "Yrs.", y = title, colour="Cluster",
               caption = paste0("nsub = ", length(unique(finaldata4$id)),
                                ", nobs = ", nrow(finaldata4))))
  return((finaldata4))
}
################################################################################################

#' k_mean_clustering_different_ks, function that return plot of withinss for different k.
#' @params, data, matrix, dataset that should be clustered
#' @param , k integers, vector of different k values (integers) to cluster for 

k_mean_clustering_different_ks <- function(data, k = 1:30){
  if(is.data.frame(data))
    data <- as.matrix(data)
  if(!is.matrix(data))
    stop("Wrong input data, should be a data frame of matrix!")
  
  prop_within <- rep(NA, length(k))
  for(i in 1:length(k)){
    k_i <- k[i]
    clust_k <- kmeans(data, centers = k_i)
    prop_within[i] <- clust_k$tot.withinss
    print(paste(k_i, round(prop_within[i]), sep = ": "))
  }                  
  
  prop_within1 <- kmeans(data, centers = 1)$tot.withinss
  d <- data.frame(k = k, prop_within = prop_within/prop_within1)
  gg <- ggplot2::ggplot(d, ggplot2::aes(x = k, y = prop_within)) + 
    ggplot2::geom_point() + 
    ggplot2::ylim(0,1)
  
  return(gg)                  
}

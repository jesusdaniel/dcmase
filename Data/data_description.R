#' US airport data
#'
#' The US airport data are a collection of networks representing 
#' the number of monthly flights between some of the main US airports 
#' during the period January 2016 to September 2021. 
#' The nodes correspond to 343 airports in the continental US, 
#' and each edge encodes the total number of flights between two airports 
#' during a given month. The data are publicly available in the Bureau of 
#' Transportation Statistics website.
#' 
#' @docType data
#' 
#' @format \code{Adj_list} is a list with 69 adjacency matrices. Each adjacency matrix
#' encodes the total number of flights between two airports on a given month.
#' 
#' @format \code{airport_coordinates} is an array of size 343x2 containing the latitute and longitude 
#' of each airport. This is used for visualization purposes
#' 
#' @format \code{airport_names} is a vector with the airport codes, in the same order as they appear
#' in the adjacency matrices
#' Calculating genetic gain and rate of inbreeding
#'
#' This function can be used to calculate the genetic gain per generation for each trait and the rate of inbreeding from the output of the breeding program simulations (pedigree needs to be called ped).
#' Rate of inbreeding is calculated from the first simulated generation till the last generation, unless a different end generation is specified.
#' Genetic gain is an average genetic gain per generation, calculated over multiple generations.
#' @param Ntraits The number of simulated traits. Does not need to be specified if Ntraits is in the list 'BPdata'.
#' @param startgen The first generation from which the (average) genetic gain per generation needs to be calculated. By default, this is generation 2 for each trait.
#' @param endgen The last generation that needs to be included in the calculation of the (average) genetic gain per generation. By default, this is the last simulated generation for each trait.
#' @param endgenF The last generation that needs to be included in the calculateion of the rate of inbreeding. This is by default the last simulated generation for each trait.
#' @return A data frame with the rate of inbreeding and the genetic gain for each trait.
#' @export
#' @examples
#' \dontrun{
#' deltaG_F()
#'}


deltaG_F <- function(Ntraits = BPdata$Ntraits, startgen = c(rep(2, times=Ntraits)), endgen = c(rep(max(ped$generation), times=Ntraits)), endgenF = max(ped$generation)){
  mean_ped <- stats::aggregate(ped[, c(7, 12:(11+Ntraits))], list(ped$generation), mean)

  output <- data.frame(deltaF = 1-(1-mean_ped[endgenF+1,2])^(1/(endgenF-1)))

  for(trait in 1: Ntraits){
    output <- cbind(output, (mean_ped[(endgen[trait]+1),(2+trait)] - mean_ped[startgen[trait]+1,2+trait])/(endgen[trait]-startgen[trait]))
  }

  colnames <- c("deltaF")

  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("deltaG",traits), collapse=""))
  }

  names(output) <- colnames

  return (output)
}

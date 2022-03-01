#' Importing genetic parameters from excel file
#'
#' This function can be used to import data (heritabilities, variances, correlations) from a excel file with a specific format, and it produces a list. If this list is called 'BPdata', then the genetic parameters do not need to be specified anymore in the other functions.
#' @param nameexcelfile Name or path to the excel file.
#' @return A list with parameters that can be used in the other functions of this package.
#' @export
#' @examples
#' \dontrun{
#' BPdata <- gen_param("example.xlsx")
#'}


gen_param <- function(nameexcelfile){

  Nm <- as.numeric(readxl::read_excel(path=nameexcelfile, col_names=F, range="BP parameters!B3"))
  Nf <- as.numeric(readxl::read_excel(path=nameexcelfile, col_names=F, range="BP parameters!B2"))
  Nm2<- as.numeric(readxl::read_excel(path=nameexcelfile, col_names=F, range="BP parameters!B5"))
  Nf2<- as.numeric(readxl::read_excel(path=nameexcelfile, col_names=F, range="BP parameters!B4"))
  Ntraits<- as.numeric(readxl::read_excel(path=nameexcelfile, col_names=F, range="BP parameters!B8"))
  prob_male<- as.numeric(readxl::read_excel(path=nameexcelfile, col_names=F, range="BP parameters!B6"))
  prob_female<- as.numeric(readxl::read_excel(path=nameexcelfile, col_names=F, range="BP parameters!B7"))

  h2 <- as.matrix(readxl::read_excel(path=nameexcelfile, col_names=F, range=paste("Genetic parameters!C2:C", (1+Ntraits), sep="")))
  c2 <- as.matrix(readxl::read_excel(path=nameexcelfile, col_names=F, range=paste("Genetic parameters!D2:D", (1+Ntraits), sep="")))
  mean <- as.matrix(readxl::read_excel(path=nameexcelfile, col_names=F, range=paste("Genetic parameters!B2:B", (1+Ntraits), sep="")))
  phen_var <- as.matrix(readxl::read_excel(path=nameexcelfile, col_names=F, range=paste("Genetic parameters!E2:E", (1+Ntraits), sep="")))
  gen_var <- h2*phen_var
  com_var <- c2*phen_var
  res_var <- (1-h2-c2)*phen_var

  letters <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","AA")

  Rgen <- as.matrix(readxl::read_excel(path=nameexcelfile, col_names=F, range=paste("Genetic correlations!A1:", letters[Ntraits],(Ntraits), sep="")))
  Rcom <- as.matrix(readxl::read_excel(path=nameexcelfile, col_names=F, range=paste("Com correlations!A1:", letters[Ntraits], (Ntraits), sep="")))
  Rres <- as.matrix(readxl::read_excel(path=nameexcelfile, col_names=F, range=paste("Residual correlations!A1:", letters[Ntraits],(Ntraits), sep="")))

  BPdata <- list(Nm, Nf, Nm2, Nf2,Ntraits, prob_male, prob_female, h2, c2, mean, phen_var, gen_var, com_var, res_var, Rgen, Rcom, Rres)
  names(BPdata) <-c("Nm", "Nf", "Nm2", "Nf2", "Ntraits", "prob_male", "prob_female", "h2", "c2", "mean", "phen_var", "a_var", "c_var","e_var", "Rgen", "Rcom", "Rres")

  return(BPdata)
}

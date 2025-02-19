% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/founderpop_functions.R
\name{founderpopgroup}
\alias{founderpopgroup}
\title{Founder population group mating design}
\usage{
founderpopgroup(
  Nm = BPdata$Nm,
  Nf = BPdata$Nf,
  Nm2 = BPdata$Nm2,
  Nf2 = BPdata$Nf2,
  Nbatch = NA,
  Nbatch2 = NA,
  batch = NA,
  batch2 = NA,
  Ntraits = BPdata$Ntraits,
  TraitsIndex = c(1:Ntraits),
  Rgen = BPdata$Rgen,
  Rres = BPdata$Rres,
  Rcom = BPdata$Rcom,
  mean = BPdata$mean,
  a_var = BPdata$a_var,
  c_var = BPdata$c_var,
  e_var = BPdata$e_var,
  est_EBV = FALSE,
  EBV,
  accuracy,
  indexweights = c(1)
)
}
\arguments{
\item{Nm}{Number of males that will be used to breed the first generation of offspring. Males are listed first in the output data frame and coded with 1 for sex.}

\item{Nf}{Number of females that will be used to breed the first generation of offspring. Females are listed second in the output data frame and coded with 2 for sex.}

\item{Nm2}{Additional males that will not be used to breed the first generation of offspring, but can be used as additional selection candidates in futher generations. Default is zero.}

\item{Nf2}{Additional females that will not be used to breed the first generation of offspring, but can be used as additional selection candidates in futher generations. Default is zero.}

\item{Nbatch}{Number of batches over which the founder animals are divided. The number of founder animals that need to be simulated should be a multiple of the number of batches. If Nbatch is not specified, the parameter batch can be used to divide the founder animals over batches. If both Nbatch and batch are not specified, then all founder animals will be assigned to batch 0.}

\item{Nbatch2}{Number of batches over which the additional founder animals (Nm2 + Nf2) are divided. The number of additional founder animals that need to be simulated should be a multiple of the number of batches. If Nbatch2 is not specified, the parameter batch2 can be used to divide the additional founder animals over batches. If both Nbatch2 and batch2 are both not specified, then the additional founder animals are divided over batches in the same way as the other founder animals (Nm+Nf).}

\item{batch}{A vector with names of the batches. The number of founder animals that need to be simulated should be the same or a multiple of the length of the vector. If batch is not specified, the parameter Nbatch can be used to divide the founder animals over batches. If both Nbatch and batch are not specified, then all founder animals will be assigned to batch 0.}

\item{batch2}{A vector with names of the batches for the additional founder animals (Nm2 + Nf2). The number of additional founder animals that are simulated should be the same or a multiple of the length of the vector. If batch2 is not specified, the parameter Nbatch2 can be used to divide the additional founder animals over batches. If both Nbatch2 and batch2 are not specified, then the additional founder animals are divided over batches in the same way as the other founder animals (Nm+Nf).}

\item{Ntraits}{Number of traits to be simulated. Does not need to be specified if Ntraits is in the list 'BPdata'.}

\item{TraitsIndex}{Vector of traits that are in the index. By default, all traits are in the selection index.}

\item{Rgen}{Matrix of all genetic correlations between all Ntraits. Only needs to be specified if there is no matrix of genetic correlations named Rgen in the list called 'BPdata'.}

\item{Rres}{Matrix of all residual correlations between all Ntraits. Only needs to be specified if there is no matrix of residual correlations named Rres in the list called 'BPdata'.}

\item{Rcom}{Matrix of all common environmental correlations between all Ntraits. Only needs to be specified if there is no matrix of common environmental correlations named Rcom in the list called 'BPdata'.}

\item{mean}{Vector of means of all traits. Only needs to be specified if there is no vector of means named mean in the list 'BPdata'.}

\item{a_var}{Vector of genetic variances of all traits. Only needs to be specified if there is no vector of genetic variances named a_var in the list 'BPdata'.}

\item{c_var}{Vector of common environmental variances of all traits. Only needs to be specified if there is no vector of common environmental variances named c_var in the list 'BPdata'. If there is no common environmental effect, provide a vector of zero's.}

\item{e_var}{Vector of residual variances of all traits. Only needs to be specified if there is no vector of residual variances named e_var in the list 'BPdata'.}

\item{est_EBV}{TRUE or FALSe for estimating breeding values for all founder animals. The default is FALSE.}

\item{EBV}{If est_EBV is TRUE, a vector of methods for simulating EBVs for each trait in the selection index need to be specified. There are four options: 0 for giving each animal a breeding values of 0 for the specific trait; "mean_pop" for simulating EBVs for all founder populations that are equal to the mean of the population, the mean of the population is provided in the vector 'mean'; "pheno" for simulating EBVs that are equal to the phenotype of the animal for the specific trait and "EBV" for simulating EBVs as a value correlated to the true breeding value of the animal. This correlation equals the accuracy. Accuracies need to be specified in the parameter 'accuracy'.}

\item{accuracy}{If for one of the traits the EBVs need to be calculated with the method "EBV" specified in the parameter 'EBV', then accuracies need to be provided. A value need to be added for each trait, however this value can be zero for traits from which the EBVs are not simulated using method "EBV".}

\item{indexweights}{If traits need to be combined in an index, then desired gain indices need to be specified for each trait in the index.}
}
\value{
A data frame is returned, which should be called 'ped' in order to be able to use it in the other functions. The 'ped' data frame consist of all the simulated base animals, including their sex, generation, batch number, inbreeding level, phenotypes and genetic, common environmental and residual effects of each trait.
}
\description{
This function can be used to create a founder population of unrelated animals for a breeding program with group mating.
}
\examples{
ped <- founderpopgroup(Nm=60, Nm2=0,
                     Nf=60, Nf2 = 0,
                     batch = c(0,1,2),
                     Ntraits=2,
                     TraitsIndex = 2,
                     Rgen= matrix(c(1.00   , 0.48,
                                    0.48   , 1.00),
                                 nrow = 2),
                     Rcom= matrix(c(1.00  , 0,
                                    0   , 1.00),
                                 nrow = 2),
                     Rres= matrix(c(1.00   , 0.32,
                                    0.32   , 1.00),
                                 nrow = 2),
                     mean=c(50,500),
                     a_var=c(200,8000),
                     c_var=c(0,0),
                     e_var= c(250,12000))

              
ped <- founderpopgroup(Nm=60,
                     Nf=60,
                     Nm2=120,
                     Nf2=120,
                     Nbatch = 4,
                     batch2 = c(-3,-2,-1,0),
                     Ntraits=2,
                     TraitsIndex = c(1,2),
                     Rgen= matrix(c(1.00   , 0.48,
                                    0.48   , 1.00),
                                 nrow = 2),
                     Rcom= matrix(c(1.00  , 0,
                                    0   , 1.00),
                                 nrow = 2),
                     Rres= matrix(c(1.00   , 0.32,
                                    0.32   , 1.00),
                                 nrow = 2),
                     mean=c(50,500),
                     a_var=c(200,8000),
                     c_var=c(0,0),
                     e_var= c(250,12000),
                     est_EBV = TRUE,
                     EBV= c("pheno", "EBV"),
                     accuracy= c(NA,0.78),
                     indexweight= c(1,5))

}

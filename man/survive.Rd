% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection_functions.R
\name{survive}
\alias{survive}
\title{Survive}
\usage{
survive(gen, batch = 0, presel, surv = NA, fish_per_FSfam = NA)
}
\arguments{
\item{gen}{The generation of the fish}

\item{batch}{The batch of the fish. Default is 0. It is possible to provide a vector with multiple batches.}

\item{presel}{Identifies the environment. Fish first need to be pre-selected.}

\item{surv}{Proportion of fish that is assumed to survive. Either surv or fish_per_FSfam need to be provided.}

\item{fish_per_FSfam}{The number of fish that survives per full sib family. Fish are randomly selected within a full sib family if fish_per_FSfam is specified. Either surv or fish_per_FSfam need to be provided.}
}
\value{
This function will make changes to the data frame 'ped'. Pre-selected fish that do not survive will get a zero in their column 'preselected'.
}
\description{
Function to randomly select which fish from an environment 'survive': fish that don't are not pre-selected anymore for that environment. This is done by assigning 0 to the 'presel' variable of the fish.
}
\examples{
\donttest{ped <- founderpopfam(Nm = 60,
                     Nf = 60,
                     Nm2 = 0,
                     Nf2 = 0,
                     Ntraits = 2,
                     TraitsIndex = 2,
                     Rgen = matrix(c(1.00   , 0.48,
                                    0.48   , 1.00),
                                 nrow = 2),
                     Rcom = matrix(c(1.00  , 0.5,
                                    0.5   , 1.00),
                                 nrow = 2),
                     Rres = matrix(c(1.00   , 0.32,
                                    0.32   , 1.00),
                                 nrow = 2),
                     mean = c(50,500),
                     a_var = c(200,8000),
                     c_var = c(150,1000),
                     e_var = c(250,12000))
                     
Mating <- randommating(gen = 0,
                       Nfam_FS = 120)

for(fam in 1: nrow(Mating)){
ped <- offspringFSfam(gen = 1,
                      No = 100,
                      probmale = 0.5,
                      sire = Mating$Sire[fam],
                      dam = Mating$Dam[fam],
                      Ntraits = 2,
                      TraitsIndex = 2,
                      Rgen = matrix(c(1.00   , 0.48,
                                    0.48   , 1.00),
                                 nrow = 2),
                      Rcom = matrix(c(1.00  , 0.5,
                                    0.5   , 1.00),
                                 nrow = 2),
                      Rres = matrix(c(1.00   , 0.32,
                                    0.32   , 1.00),
                                 nrow = 2),
                      a_var = c(200,8000),
                      c_var = c(150,1000),
                      e_var = c(250,12000))
}

ped <- preselphen(gen = 1,
                  Nenv = 2,
                  Npresel = c(25,15),
                  trait = 1,
                  Ntraits = 2)
                  
ped <- survive(gen = 1,
               presel = 2,
               surv = 0.8)
}
}

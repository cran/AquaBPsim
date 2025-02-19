% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/offspring_functions.R
\name{groupmating}
\alias{groupmating}
\title{Group mating}
\usage{
groupmating(
  gen,
  batch = 0,
  batch_m = NA,
  batch_f = NA,
  No,
  contr_m,
  contr_f,
  distribution = "Gamma",
  shape = 0.75,
  scale = 0.11,
  selected = 1
)
}
\arguments{
\item{gen}{The generations of the sires and dams}

\item{batch}{The batches of the sires and dams. Default is 0.}

\item{batch_m}{The batch of the sires. Default is NA. If batch_m is specified, batch_f also needs to be specified and the parameter batch is not used.}

\item{batch_f}{The batch of the dams. Default is NA. If batch_f is specified, batch_m also needs to be specified and the parameter batch is not used.}

\item{No}{The total number of offspring of all matings.}

\item{contr_m}{Proportion of sires that contribute to the offspring.}

\item{contr_f}{Proportion of dams that contribute to the offspring.}

\item{distribution}{The distribution from which the contributions are drawn. Options are: "Gamma" (default) and "Uniform"}

\item{shape}{The shape of the gamma distribution. Default is 0.75.}

\item{scale}{The scale of the gamma distribution. Default is 0.11.}

\item{selected}{The value in ped$selected of the selected sires and dams. Default is 1.}
}
\value{
The output is a data frame with for each full sib family the sire and dam and the number of offspring per full sib family.
}
\description{
Function to determine which full sib families are produced in a group mating design.
}
\details{
By default, the contribution of the sires and dams that do reproduce come from a gamma distribution. The default shape and scale of the gamma distribution are 0.75 and 0.11, respectively. A uniform distribution can also be specified for the contributions of the sires and dams.
If not all sires and dams should contribute to the offspring, then the sires and dams that are going to reproduce are randomly chosen.
The output is a dataframe called Mating with the sire, dam and the size of each full sib family.

The sires and dams can come from multiple batches or generations. In that case, a vector of batches or generations need to be provided.
}
\examples{

{ped <- founderpopgroup(Nm=60,
                     Nf=60,
                     Nm2=120,
                     Nf2=120,
                     Nbatch = 4,
                     batch2 = c(-3,-2,-1,0),
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
                     
 Mating <- groupmating(gen = 0,
                       batch=-3,
                       No=1000,
                       contr_m = 0.5,
                       contr_f = 0.5)
}
}

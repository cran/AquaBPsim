## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  library(AquaBPsim)
#  
#  BPdata <- list(Ntraits = 3,
#               h2 = c(0.3,0.25,0.15),
#               c2 = c(0.1,0.05,0.02),
#               p_var = c(100,6400,4),
#               a_var = c(0.3,0.25,0.15)*c(100,6400,4),
#               c_var = c(0.1,0.05,0.02)*c(100,6400,4),
#               e_var = c(100,6400,4)*(1-c(0.1,0.05,0.02)-c(0.3,0.25,0.15)),
#               mean = c(50, 400, 10),
#               Rgen = matrix(c(1,   0.55,  0.1,
#                               0.55,   1,  0.3,
#                               0.1 ,  0.3,   1), nrow = 3),
#               Rres = matrix(c(1,   0.3,  0,
#                               0.3,   1,  0,
#                               0 ,  0, 1), nrow = 3),
#               Rcom = matrix(c(1,   0,  0,
#                               0,   1,  0,
#                               0 ,  0,  1), nrow = 3))
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  ped <- founderpopfam(Nm=200,
#                       Nf=200,
#                       batch=c(0,-1,-2,-3),
#                       Ntraits=3,
#                       TraitsIndex=c(2,3),
#                       a_var = c(0.3,0.25,0.15)*c(100,6400,4),
#                       c_var = c(0.1,0.05,0.02)*c(100,6400,4),
#                       e_var = c(100,6400,4)*(1-c(0.1,0.05,0.02)-c(0.3,0.25,0.15)),
#                       mean = c(50, 400, 10),
#                       Rgen = matrix(c(1,   0.55,  0.1,
#                                       0.55,   1,  0.3,
#                                       0.1 ,  0.3,   1), nrow = 3),
#                       Rres = matrix(c(1,   0.3,  0,
#                                       0.3,   1,  0,
#                                       0 ,  0, 1), nrow = 3),
#                       Rcom = matrix(c(1,   0,  0,
#                                       0,   1,  0,
#                                       0 ,  0,  1), nrow = 3))
#  
#  

## ----eval=FALSE---------------------------------------------------------------
#  Mating <- randommating(gen = 0,
#                         batch = -3,
#                         Nfam_FS = 100)
#  

## ----eval=FALSE---------------------------------------------------------------
#  for(mating in 1:nrow(Mating)){
#    ped <- offspringFSfam(gen=1,
#                          No=50,
#                          sire=Mating$Sire[mating],
#                          dam=Mating$Dam[mating],
#                          batch = 1,
#                          probmale = 0.5,
#                          TraitsIndex = c(2,3))
#  }
#  

## ----eval=FALSE---------------------------------------------------------------
#  ped <- preselphen(gen = 1,
#                    batch=1,
#                    Nenv = 2,
#                    Npresel = c(10,5),
#                    trait= 1)
#  

## ----eval=FALSE---------------------------------------------------------------
#  ped <- preselrandom(gen = 1,
#                      batch=1,
#                      Nenv = 2,
#                      Npresel = c(10,5))
#  

## ----eval=FALSE---------------------------------------------------------------
#  ped <- preselselcand(gen = 1,
#                    batch = 1,
#                    Nm =100,
#                    Nf = 100,
#                    max_FSfam = 15)

## ----eval=FALSE---------------------------------------------------------------
#  ped <- avail_selection(gen = 1,
#                         batch = 1,
#                         presel = 1,
#                         surv = 0.9)
#  

## ----eval=FALSE---------------------------------------------------------------
#  ped <- survive(gen = 1,
#                 batch = 1,
#                 presel = 2,
#                 surv = 0.9)
#  

## ----eval=FALSE---------------------------------------------------------------
#  EBV = r^2 * TBV + X * r * sqrt(1 - r^2)

## ----eval=FALSE---------------------------------------------------------------
#  ped <- breeding_values(gen=1,
#                         batch = 1,
#                         TraitsIndex=c(2,3),
#                         EBV=c("GEBV","GEBV"),
#                         GenomLength = 11.3,
#                         Ne = 100,
#                         SizeTraining = c(2000, 2000, 1000),
#                         indexweights=c(2,1))
#  

## ----eval=FALSE---------------------------------------------------------------
#   ped <- select(gen=1,
#                    batch = 1,
#                    Nm = 50,
#                    Nf = 50,
#                    mature_m = 0.5,
#                    mature_f = 0.5)
#  

## ----eval=FALSE---------------------------------------------------------------
#  library(AquaBPsim)
#  
#  BPdata <- list(Ntraits = 3,
#               h2 = c(0.3,0.25,0.15),
#               c2 = c(0.1,0.05,0.02),
#               p_var = c(100,6400,4),
#               a_var = c(0.3,0.25,0.15)*c(100,6400,4),
#               c_var = c(0.1,0.05,0.02)*c(100,6400,4),
#               e_var = c(100,6400,4)*(1-c(0.1,0.05,0.02)-c(0.3,0.25,0.15)),
#               mean = c(50, 400, 10),
#               Rgen = matrix(c(1,   0.55,  0.1,
#                               0.55,   1,  0.3,
#                               0.1 ,  0.3,   1), nrow = 3),
#               Rres = matrix(c(1,   0.3,  0,
#                               0.3,   1,  0,
#                               0 ,  0, 1), nrow = 3),
#               Rcom = matrix(c(1,   0,  0,
#                               0,   1,  0,
#                               0 ,  0,  1), nrow = 3))
#  
#  ped <- founderpopfam(Nm = 100,
#                       Nf = 100,
#                       TraitsIndex = c(2,3))
#  
#  # Selected males and females are randomly allocated to each other with the function randommating.
#  # 200 combinations are made in order to simulate 200 full sib families.
#  Mating <- randommating(gen = 0,
#                         Nfam_FS = 200)
#  
#  # Offspring are simualted for each full sib family separately.
#  for(mating in 1:nrow(Mating)){
#    ped <- offspringFSfam(gen = 1,
#                          No = 50,
#                          sire = Mating$Sire[mating],
#                          dam = Mating$Dam[mating],
#                          probmale = 0.5,
#                          TraitsIndex = c(2,3))
#  }
#  
#  # Pre-selection based on the phenotype of trait 1. Fish are pre-selected within a fullsib family: 10 for the nucleus and 10 for the production environment.
#  ped <- preselphen(gen = 1,
#                    Nenv = 2,
#                    Npresel = c(10,10),
#                    trait = 1)
#  
#  # 90% of the fish that are pre-selected for the nucleus will be available for selection. These 90% are randomly chosen.
#  ped <- avail_selection(gen = 1, presel = 1, surv = 0.9)
#  
#  # Simulating breeding values for all available selection candidates.
#  ped <- breeding_values(gen = 1,
#             TraitsIndex = c(2,3),
#             EBV = c("PEBV","sib_pheno"),
#             indexweights = c(2,1))
#  
#  # Selecting 100 males and 100 females based on their Index. It was assumed that 80% of the males and females are mature and therefore available for selection at the moment of selection
#  ped<- select(Nm = 100,
#               Nf = 100,
#               gen = 1,
#               mature_m = 0.8,
#               mature_f = 0.8)
#  
#  for(generation in 2:10){
#  Mating <- randommating(gen = generation-1,
#                         Nfam_FS = 200)
#  
#  for(mating in 1:nrow(Mating)){
#    ped <- offspringFSfam(gen = generation,
#                          No = 50,
#                          sire = Mating$Sire[mating],
#                          dam = Mating$Dam[mating],
#                          probmale = 0.5,
#                          TraitsIndex = c(2,3))
#  }
#  
#  ped <- preselphen(gen = generation,
#                    Nenv = 2,
#                    Npresel = c(10,10),
#                    trait = 1)
#  
#  ped <- avail_selection(gen = generation, presel = 1, surv = 0.9)
#  
#  ped <- breeding_values(gen = generation,
#                         TraitsIndex = c(2,3),
#                         EBV = c("PEBV","sib_pheno"),
#                         indexweights = c(2,1))
#  
#  ped<- select(Nm = 100,
#               Nf = 100,
#               gen = generation,
#               mature_m = 0.8,
#               mature_f = 0.8)
#  
#  }
#  
#  
#  # calculating genetic gain and rate of inbreeding
#  deltaG_F()
#  

## ----eval=FALSE---------------------------------------------------------------
#  library(AquaBPsim)
#  
#  BPdata <- list(Ntraits = 3,
#                 h2 = c(0.3,0.25,0.15),
#                 c2 = c(0.1,0.05,0.02),
#                 p_var = c(100,6400,4),
#                 a_var = c(0.3,0.25,0.15)*c(100,6400,4),
#                 c_var = c(0.1,0.05,0.02)*c(100,6400,4),
#                 e_var = c(100,6400,4)*(1-c(0.1,0.05,0.02)-c(0.3,0.25,0.15)),
#                 mean = c(50, 400, 10),
#                 Rgen = matrix(c(1,   0.55,  0.1,
#                                 0.55,   1,  0.3,
#                                 0.1 ,  0.3,   1), nrow = 3, byrow = TRUE),
#                 Rres = matrix(c(1,   0.3,  0,
#                                 0.3,   1,  0,
#                                 0 ,  0, 1), nrow = 3, byrow = TRUE),
#                 Rcom = matrix(c(1,   0,  0,
#                                 0,   1,  0,
#                                 0 ,  0,  1), nrow = 3, byrow = TRUE))
#  
#  next_batch <- c(2:6,1)
#  ped <- founderpopgroup(Nm = 120,
#                         Nf = 60,
#                         Nbatch = 6,
#                         TraitsIndex = c(2,3))
#  
#  # A loop is created in order to simulate each batch separately.
#  for(batch in 1:6){
#  
#    # founder animal from one batch are mated with each other. 50% of the males and females contributes to the offspring, with contribution drawn from a gamma distribution with shape 0.75 and scale 0.11. The total number of offspring per batch is approximately 1000.
#    Mating <- groupmating(gen = 0,
#                          batch = batch,
#                          No = 1000,
#                          contr_m = 0.5,
#                          contr_f =0.5)
#  
#    # Offspring are simualted for each full sib family.
#    for(mating in 1:nrow(Mating)){
#      ped <- offspringFSgroup(gen = 1,
#                              batch = batch,
#                              No = Mating$No[mating],
#                              sire = Mating$Sire[mating],
#                              dam = Mating$Dam[mating],
#                              probmale = 0.5,
#                              TraitsIndex = c(2,3))
#    }
#  
#    # Preselection based on the phenotype of trait 1. Fish are preselected within a batch: 300 for the nucleus and 100 for the production environment.
#    ped <- preselphen(gen = 1,
#                      batch = batch,
#                      Nenv = 2,
#                      Npresel = c(300,100),
#                      trait = 1,
#                      withinfam = F)
#  
#    # 90% of the fish that are preselected for the nucleus will be available for selection. These 90% are randomly chosen.
#    ped <- avail_selection(gen = 1, batch = batch, presel = 1, surv = 0.9)
#  
#    # Simulating genomic breeding values for all available selection candidates.
#    ped <- breeding_values(gen = 1,
#                           batch = batch,
#                           TraitsIndex = c(2,3),
#                           EBV = c("GEBV","GEBV"),
#                           accuracy = c(0.85,0.78),
#                           indexweights = c(2,1))
#  
#    # Selecting 20 males and 10 females based on their Index.
#    ped<- select(gen = 1,
#                 batch = batch,
#                 Nm = 20,
#                 Nf = 10)
#  }
#  
#  for(generation in 2:10){
#    for(batch in 1:6){
#  
#      Mating <- groupmating(gen = generation - 1,
#                            batch_m = batch,
#                            batch_f = next_batch[batch],
#                            No = 1000,
#                            contr_m = 0.5,
#                            contr_f = 0.5)
#  
#      for(mating in 1:nrow(Mating)){
#        ped <- offspringFSgroup(gen = generation,
#                                batch = batch,
#                                No = Mating$No[mating],
#                                sire = Mating$Sire[mating],
#                                dam = Mating$Dam[mating],
#                                probmale = 0.5,
#                                TraitsIndex = c(2,3))
#      }
#  
#      ped <- preselphen(gen = generation,
#                        batch = batch,
#                        Nenv = 2,
#                        Npresel = c(300,100),
#                        trait = 1,
#                        withinfam = F)
#  
#      ped <- avail_selection(gen = generation, batch = batch, presel = 1, surv = 0.9)
#  
#      ped <- breeding_values(gen = generation,
#                             batch = batch,
#                             TraitsIndex = c(2,3),
#                             EBV = c("GEBV","GEBV"),
#                             accuracy = c(0.85,0.78),
#                             indexweights = c(2,1))
#  
#      ped<- select(gen = generation,
#                   batch = batch,
#                   Nm = 20,
#                   Nf = 10)
#    }
#  }
#  
#  
#  # calculating genetic gain and rate of inbreeding
#  deltaG_F()
#  


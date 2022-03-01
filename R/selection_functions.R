# Functions:
#        - preselphen
#        - preselrandom
#        - avail_selection
#        - preselselcand
#        - select



#' Preselecting offspring based on phenotype
#'
#' This function can be used to preselect offspring based on the phenotype of one trait.
#'
#' @param gen The generation of the offspring.
#' @param batch The batch of the offspring. Default is 0.
#' @param withinfam Preselection within a full sib family or not, default is TRUE.
#' @param Nenv The number of environments the fish need to be preselected for, for example a certain number of fish need to be preselected for a nucleus and a certain number of fish need to be preselected for a production environment in order to measure sib traits. Preselection is based on the same trait.
#' @param Npresel A vector of the number of fish that needs to be preselected for each environment.
#' @param trait The trait on which the phenotypic preselection is based, for example trait = 2.
#' @param Ntraits The total number of simulated traits. Needs to be specified if Ntraits is not specified in a list called 'BPdata'.
#' @return This function will make changes to the data frame called 'ped'. Fish that are pre-selected will get a number assigned to their column 'preselected'.
#' @export
#' @examples
#' 
#' \donttest{ped <- founderpopfam(Nm = 60,
#'                      Nf = 60,
#'                      Nm2 = 0,
#'                      Nf2 = 0,
#'                      Ntraits = 2,
#'                      TraitsIndex = 2,
#'                      Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                      Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                      Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                      mean = c(50,500),
#'                      a_var = c(200,8000),
#'                      c_var = c(150,1000),
#'                      e_var = c(250,12000))
#'                      
#' Mating <- randommating(gen = 0,
#'                        Nfam_FS = 120)
#'
#' for(fam in 1: nrow(Mating)){
#' ped <- offspringFSfam(gen = 1,
#'                       No = 100,
#'                       probmale = 0.5,
#'                       sire = Mating$Sire[fam],
#'                       dam = Mating$Dam[fam],
#'                       Ntraits = 2,
#'                       TraitsIndex = 2,
#'                       Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                       Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                       Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                       a_var = c(200,8000),
#'                       c_var = c(150,1000),
#'                       e_var = c(250,12000))
#' }
#' 
#' ped <- preselphen(gen = 1,
#'                   Nenv = 2,
#'                   Npresel = c(25,15),
#'                   trait = 1,
#'                   Ntraits = 2)
#' }                  
#'                   
#'                   
#' \donttest{ped <- founderpopgroup(Nm = 60,
#'                      Nf = 60,
#'                      Nm2 = 0,
#'                      Nf2 = 0,
#'                      batch = c(-3,-2,-1,0),
#'                      Ntraits = 2,
#'                      TraitsIndex = 2,
#'                      Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                      Rcom = matrix(c(1.00  , 0,
#'                                     0   , 1.00),
#'                                  nrow = 2),
#'                      Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                      mean = c(50,500),
#'                      a_var = c(200,8000),
#'                      c_var = c(0,0),
#'                      e_var = c(250,12000))
#'  
#' Mating <- groupmating(gen = 0,
#'                        batch =-3,
#'                        No = 1000,
#'                        contr_m = 0.5,
#'                        contr_f = 0.5)                    
#'                      
#' for(fam in 1: nrow(Mating)){
#' ped <- offspringFSgroup(gen = 1,
#'                       No = Mating$No[fam],
#'                       probmale = 0.5,
#'                       sire = Mating$Sire[fam],
#'                       dam = Mating$Dam[fam],
#'                       batch = 1,
#'                       Ntraits = 2,
#'                       TraitsIndex = 2,
#'                       Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                       Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                       Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                       a_var = c(200,8000),
#'                       c_var = c(150,1000),
#'                       e_var = c(250,12000))
#' }
#'
#' ped <- preselphen(gen = 1,
#'                   batch = 1,
#'                   withinfam = FALSE,
#'                   Nenv = 2,
#'                   Npresel = c(400,150),
#'                   trait = 1,
#'                   Ntraits = 2)
#'                   }



preselphen <- function(gen, batch=0, withinfam=TRUE, Nenv, Npresel, trait, Ntraits=BPdata$Ntraits){
  if(withinfam==T){
    preselcand <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$preselected==0,]
    for(fam in min(preselcand$FSfam):max(preselcand$FSfam)){
      for(N in 1:Nenv){
        famavail <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$preselected==0 & ped$FSfam==fam,]
        famavail <- famavail[order(famavail[c(11+3*Ntraits+trait)], decreasing=T),]
        famavail <- famavail[1:Npresel[N],]
        ped$preselected[ped$id %in% famavail$id] <- N
      }
    }
  }else{
    for(N in 1:Nenv){
      preselcand <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$preselected==0,]
      preselcand <- preselcand[order(preselcand[c(11+3*Ntraits+trait)],decreasing=T),]
      preselcand <- preselcand[1:Npresel[N],1]
      ped$preselected[ped$id %in% preselcand] <- N
    }
  }
  return(ped)
}


#' Randomly preselecting offspring
#'
#' This function can be used to randomly preselect offspring or randomly allocate offspring to for example the nucleus and prodction environment.
#'
#' @param gen The generation of the offspring.
#' @param batch The batch of the offspring that need to be randomly preselected. Default is 0.
#' @param withinfam Preselection within a full sib family or not, default is TRUE.
#' @param Nenv The number of environments the fish need to be preselected for, for example a certain number of fish need to be preselected for a nucleus and a certain number of fish need to be preselected for a production environment in order to measure sib traits.
#' @param Npresel The number of fish that needs to be preselected for each environment.
#' @return This function will make changes to the data frame called 'ped'. Fish that are pre-selected will get a number assigned to their column 'preselected'.
#' @export
#' @examples
#' \donttest{ped <- founderpopfam(Nm = 60,
#'                      Nf = 60,
#'                      Nm2 = 0,
#'                      Nf2 = 0,
#'                      Ntraits = 2,
#'                      TraitsIndex = 2,
#'                      Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                      Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                      Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                      mean = c(50,500),
#'                      a_var = c(200,8000),
#'                      c_var = c(150,1000),
#'                      e_var = c(250,12000))
#'                      
#' Mating <- randommating(gen = 0,
#'                        Nfam_FS = 120)
#'
#' for(fam in 1: nrow(Mating)){
#' ped <- offspringFSfam(gen = 1,
#'                       No = 100,
#'                       probmale = 0.5,
#'                       sire = Mating$Sire[fam],
#'                       dam = Mating$Dam[fam],
#'                       Ntraits = 2,
#'                       TraitsIndex = 2,
#'                       Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                       Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                       Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                       a_var = c(200,8000),
#'                       c_var = c(150,1000),
#'                       e_var = c(250,12000))
#' }
#' 
#' ped <- preselrandom(gen = 1,
#'                     Nenv = 2,
#'                     Npresel = c(25,15))
#'}


preselrandom <- function(gen, batch=0, withinfam=TRUE, Nenv, Npresel){
  preselcand <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$preselected==0,]
   if(withinfam==TRUE){
    for(fam in min(preselcand$FSfam):max(preselcand$FSfam)){
      for(N in 1:Nenv){
        famavail <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$preselected==0 & ped$FSfam==fam,]
        famavail <- sample(famavail$id, Npresel[N])
        ped$preselected[ped$id %in% famavail] <- N
      }
    }
  }else{
    for(N in 1:Nenv){
      avail <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$preselected==0,]
      avail <- sample(avail$id, Npresel[N])
      ped$preselected[ped$id %in% avail] <- N
    }
  }
  return(ped)
}


#' Available as selection candidates
#'
#' Function to determine which fish are available for selection. Fish are randomly chosen.
#'
#' @param gen The generation of the fish.
#' @param batch The batch of the fish. Default is 0. It is possible to provide a vector with multiple batches.
#' @param presel Identifies which preselected fish are available for selection. If no pre-selection took place, then presel should be 0 (= default).
#' @param surv Proportion of fish that is assumed to survive till the moment of selection. Either surv or fish_per_FSfam need to be provided.
#' @param fish_per_FSfam The number of fish available for selection per full sib family. Fish are randomly selected within a full sib family if fish_per_FSfam is specified. Either surv or fish_per_FSfam need to be provided.
#' @return This function will make changes to the data frame called 'ped'. Fish that become available as selection candidates will be assigned a 1 to their column 'selcand'.
#' @export
#' @examples
#' \donttest{ped <- founderpopfam(Nm = 60,
#'                      Nf = 60,
#'                      Nm2 = 0,
#'                      Nf2 = 0,
#'                      Ntraits = 2,
#'                      TraitsIndex = 2,
#'                      Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                      Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                      Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                      mean = c(50,500),
#'                      a_var = c(200,8000),
#'                      c_var = c(150,1000),
#'                      e_var = c(250,12000))
#'                      
#' Mating <- randommating(gen = 0,
#'                        Nfam_FS = 120)
#'
#' for(fam in 1: nrow(Mating)){
#' ped <- offspringFSfam(gen = 1,
#'                       No = 100,
#'                       probmale = 0.5,
#'                       sire = Mating$Sire[fam],
#'                       dam = Mating$Dam[fam],
#'                       Ntraits = 2,
#'                       TraitsIndex = 2,
#'                       Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                       Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                       Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                       a_var = c(200,8000),
#'                       c_var = c(150,1000),
#'                       e_var = c(250,12000))
#' }
#' 
#' ped <- preselphen(gen = 1,
#'                   Nenv = 2,
#'                   Npresel = c(25,15),
#'                   trait = 1,
#'                   Ntraits = 2)
#'                   
#' ped <- avail_selection(gen = 1,
#'                        presel = 1,
#'                        surv = 0.9)
#'}

avail_selection <- function(gen, batch=0, presel = 0, surv=NA, fish_per_FSfam=NA){
  if(!is.na(surv)){
    avail <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$preselected %in% presel,]
    avail <- sample(avail$id, (nrow(avail)*surv))
    ped$selcand[ped$id %in% avail] <- 1
  }

  if(!is.na(fish_per_FSfam)){
    avail <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$preselected %in% presel,]
    for(fam in (min(avail$FSfam):max(avail$FSfam))){
      availfam <- avail[avail$FSfam==fam,]
      availfam <- sample(availfam$id, fish_per_FSfam)
      ped$selcand[ped$id %in% availfam] <- 1
    }
  }
  return (ped)
}


#' Survive
#'
#' Function to randomly select which fish from an environment 'survive': fish that don't are not pre-selected anymore for that environment. This is done by assigning 0 to the 'presel' variable of the fish.
#'
#' @param gen The generation of the fish
#' @param batch The batch of the fish. Default is 0. It is possible to provide a vector with multiple batches.
#' @param presel Identifies the environment. Fish first need to be pre-selected.
#' @param surv Proportion of fish that is assumed to survive. Either surv or fish_per_FSfam need to be provided.
#' @param fish_per_FSfam The number of fish that survives per full sib family. Fish are randomly selected within a full sib family if fish_per_FSfam is specified. Either surv or fish_per_FSfam need to be provided.
#' @return This function will make changes to the data frame 'ped'. Pre-selected fish that do not survive will get a zero in their column 'preselected'. 
#' @export
#' @examples
#' \donttest{ped <- founderpopfam(Nm = 60,
#'                      Nf = 60,
#'                      Nm2 = 0,
#'                      Nf2 = 0,
#'                      Ntraits = 2,
#'                      TraitsIndex = 2,
#'                      Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                      Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                      Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                      mean = c(50,500),
#'                      a_var = c(200,8000),
#'                      c_var = c(150,1000),
#'                      e_var = c(250,12000))
#'                      
#' Mating <- randommating(gen = 0,
#'                        Nfam_FS = 120)
#'
#' for(fam in 1: nrow(Mating)){
#' ped <- offspringFSfam(gen = 1,
#'                       No = 100,
#'                       probmale = 0.5,
#'                       sire = Mating$Sire[fam],
#'                       dam = Mating$Dam[fam],
#'                       Ntraits = 2,
#'                       TraitsIndex = 2,
#'                       Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                       Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                       Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                       a_var = c(200,8000),
#'                       c_var = c(150,1000),
#'                       e_var = c(250,12000))
#' }
#' 
#' ped <- preselphen(gen = 1,
#'                   Nenv = 2,
#'                   Npresel = c(25,15),
#'                   trait = 1,
#'                   Ntraits = 2)
#'                   
#' ped <- survive(gen = 1,
#'                presel = 2,
#'                surv = 0.8)
#'}

survive <- function(gen, batch=0, presel, surv=NA, fish_per_FSfam=NA){
  if(!is.na(surv)){
    avail <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$preselected %in% presel,]
    avail <- sample(avail$id, (nrow(avail)*(1-surv)))
    ped$preselected[ped$id %in% avail] <- 0
  }

  if(!is.na(fish_per_FSfam)){
    avail <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$preselected %in% presel,]
    for(fam in (min(avail$FSfam):max(avail$FSfam))){
      availfam <- avail[avail$FSfam==fam,]
      availfam <- sample(availfam$id, (nrow(availfam) - fish_per_FSfam))
      ped$preselected[ped$id %in% availfam] <- 0
    }
  }
  return (ped)
}

#' Preselection of selection candidates
#'
#' Function to preselect selection candidates based on their Index/EBV/phenotype. Fish are preselected from the fish that are available for selection (ped$selcand==1). Fish that are not preselected will be assigned a value of 2 for ped$selcand. These fish will not be used in the function select.
#'
#' @param gen The generation of the selection candidates. A vector of generations can be provided.
#' @param batch The batch of the selection candidates. Default is 0. It is possible to provide a vector with multiple batches.
#' @param select_on Options: "Index", "EBV" or "Phenotype". For "EBV" and "Phenotype", also the trait need to be specified in 'trait'. Default is Index.
#' @param trait Which trait or EBV the selection is based on when option "EBV" or "Phenotype" is choosen in 'select on'.
#' @param Ntraits Number of simulated traits. Does not need to be specified if Ntraits is specified in a list called 'BPdata'.
#' @param Nm Number of males to preselect (In total or per full sib family, depending on within_FSfam).
#' @param Nf Number of females to preselect (In total or per full sib family, depending on within_FSfam).
#' @param N Total number of pre-selected fish. Does not need to be specified if Nm and Nf are specified.
#' @param within_FSfam If True, pre-selection takes place within a full sib family. Default is False. Only use within family selection in a family design. Each full sib family must have at least Nm male sibs and Nf female sibs (or N sibs).
#' @param max_FSfam Maximum number of sibs that can be selected per full sib familie, in case selection does not take place within a full sib family. Default is 'all'.
#' @return This function will change the data frame called 'ped'. Fish that are not preselected will be assigned a value of 2 for ped$selcand.
#' @export
#' @examples
#' \donttest{ped <- founderpopfam(Nm = 60,
#'                      Nf = 60,
#'                      Nm2 = 0,
#'                      Nf2 = 0,
#'                      Ntraits = 2,
#'                      TraitsIndex = 2,
#'                      Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                      Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                      Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                      mean = c(50,500),
#'                      a_var = c(200,8000),
#'                      c_var = c(150,1000),
#'                      e_var = c(250,12000))
#'                      
#' Mating <- randommating(gen = 0,
#'                        Nfam_FS = 120)
#'
#' for(fam in 1: nrow(Mating)){
#' ped <- offspringFSfam(gen = 1,
#'                       No = 100,
#'                       probmale = 0.5,
#'                       sire = Mating$Sire[fam],
#'                       dam = Mating$Dam[fam],
#'                       Ntraits = 2,
#'                       TraitsIndex = 2,
#'                       Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                       Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                       Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                       a_var = c(200,8000),
#'                       c_var = c(150,1000),
#'                       e_var = c(250,12000))
#' }
#' 
#' ped <- preselphen(gen = 1,
#'                   Nenv = 2,
#'                   Npresel = c(25,15),
#'                   trait = 1,
#'                   Ntraits = 2)
#'                   
#' ped <- avail_selection(gen = 1,
#'                        presel = 1,
#'                        surv = 0.9)
#'                        
#' ped <- breeding_values(gen = 1,
#'                        TraitsIndex = 2,
#'                        EBV = "GEBV",
#'                        GenomLength = 11.3,
#'                        Ne = 100,
#'                        SizeTraining = nrow(ped[ped$preselected ==2,]),
#'                        Ntraits = 2,
#'                        a_var = c(200,8000),
#'                        h2 = c(0.33,0.38))
#'                        
#' ped <- preselselcand(gen = 1,
#'                  Nm = 300,
#'                  Nf = 300,
#'                  max_FSfam = 15,
#'                  Ntraits = 2)
#'}

preselselcand <- function(gen, batch=0, select_on="Index", trait, Ntraits=BPdata$Ntraits, Nm, Nf, N, within_FSfam = FALSE, max_FSfam="all" ){

# Between family selection
  if(within_FSfam == F){

  ## Between family selection, Nm and Nf specifies
    if(!(is.na(Nm))){
  selcand_m <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$sex==1 & ped$selcand==1,]
  selcand_f <- ped[ped$generation %in% gen &ped$batch %in% batch & ped$sex==2 & ped$selcand==1,]

  if(max_FSfam=="all"){
    if(select_on=="Index"){
      selcand_m <- selcand_m[order(selcand_m$Index, decreasing=T),]
      selcand_f <- selcand_f[order(selcand_f$Index, decreasing=T),]
    }
    if(select_on=="EBV"){
      selcand_m <- selcand_m[order(selcand_m[,c(11+4*Ntraits+trait)], decreasing=T),]
      selcand_f <- selcand_f[order(selcand_f[,c(11+4*Ntraits+trait)], decreasing=T),]
    }
    if(select_on=="Phenotype"){
      selcand_m <- selcand_m[order(selcand_m[,c(11+3*Ntraits+trait)], decreasing=T),]
      selcand_f <- selcand_f[order(selcand_f[,c(11+3*Ntraits+trait)], decreasing=T),]
    }

    selm <- selcand_m[1:Nm,]
    self <- selcand_f[1:Nf,]
    ped$selcand[!(ped$id %in% selm$id) & ped$generation %in% gen & ped$batch %in% batch & ped$sex==1 & ped$selcand==1] <- 2
    ped$selcand[!(ped$id %in% self$id) & ped$generation %in% gen & ped$batch %in% batch & ped$sex==2 & ped$selcand==1] <- 2
  }else{
    selcand <- rbind(selcand_m, selcand_f)
    if(select_on=="Index"){selcand <- selcand[order(selcand$Index, decreasing=T),]}
    if(select_on=="EBV"){selcand <- selcand[order(selcand[,c(11+4*Ntraits+trait)], decreasing=T),]}
    if(select_on=="Phenotype"){selcand <- selcand[order(selcand[,c(11+3*Ntraits+trait)], decreasing=T),] }

    for(i in 1:(Nm + Nf)){
      ped$selcand[ped$id %in% selcand$id[i]] <- 99
      selectedFS <- ped[ped$selcand==99 & ped$generation %in% gen & ped$batch %in% batch & ped$FSfam==selcand$FSfam[i],]
      selected_m <- ped[ped$selcand==99 & ped$generation %in% gen & ped$batch %in% batch & ped$sex==1,]
      selected_f <- ped[ped$selcand==99 & ped$generation %in% gen & ped$batch %in% batch & ped$sex==2,]

      if(nrow(selectedFS) == max_FSfam){
        sel_candidates1 <- selcand[1:i,]
        sel_candidates2 <- selcand[(i+1):nrow(selcand),]
        sel_candidates2 <- sel_candidates2[sel_candidates2$FSfam != selcand$FSfam[i],]
        selcand <-  rbind(sel_candidates1,sel_candidates2)
      }

      if(nrow(selected_m) == Nm){
        sel_candidates1 <- selcand[1:i,]
        sel_candidates2 <- selcand[(i+1):nrow(selcand),]
        sel_candidates2 <- sel_candidates2[sel_candidates2$sex == 2,]
        selcand <-  rbind(sel_candidates1,sel_candidates2)
      }
      if(nrow(selected_f) == Nf){
        sel_candidates1 <- selcand[1:i,]
        sel_candidates2 <- selcand[(i+1):nrow(selcand),]
        sel_candidates2 <- sel_candidates2[sel_candidates2$sex == 1,]
        selcand <-  rbind(sel_candidates1,sel_candidates2)
      }

    }
    ped$selcand[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1] <- 2
    ped$selcand[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==99] <- 1
  }
    }else{

      ## Between family selection, no Nm and Nf specifies
      selcand <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,]
      if(select_on=="Index"){selcand <- selcand[order(selcand$Index, decreasing=T),]}
      if(select_on=="EBV"){selcand <- selcand[order(selcand[,c(11+4*Ntraits+trait)], decreasing=T),]}
      if(select_on=="Phenotype"){selcand <- selcand[order(selcand[,c(11+3*Ntraits+trait)], decreasing=T),] }

      if(max_FSfam=="all"){
        sel <- selcand[1:N,]
        ped$selcand[!(ped$id %in% sel$id) & ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1] <- 2
      }else{

        for(i in 1:N){
          ped$selcand[ped$id %in% selcand$id[i]] <- 99
          selectedFS <- ped[ped$selcand==99 & ped$generation %in% gen & ped$batch %in% batch & ped$FSfam==selcand$FSfam[i],]

          if(nrow(selectedFS) == max_FSfam){
            sel_candidates1 <- selcand[1:i,]
            sel_candidates2 <- selcand[(i+1):nrow(selcand),]
            sel_candidates2 <- sel_candidates2[sel_candidates2$FSfam != selcand$FSfam[i],]
            selcand <-  rbind(sel_candidates1,sel_candidates2)
          }

        }
        ped$selcand[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1] <- 2
        ped$selcand[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==99] <- 1
      }
    }
  }

  if(within_FSfam == T){

    ## Within family selection, Nm and Nf specifies
    if(!(is.na(Nm))){
      allselcand_m <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$sex==1 & ped$selcand==1,]
      allselcand_f <- ped[ped$generation %in% gen &ped$batch %in% batch & ped$sex==2 & ped$selcand==1,]

      for(fam in (min(allselcand_m$FSfam):max(allselcand_m$FSfam))){

        selcand_m <- allselcand_m[allselcand_m$FSfam == fam]
        selcand_f <- allselcand_f[allselcand_f$FSfam == fam]

        if(select_on=="Index"){
          selcand_m <- selcand_m[order(selcand_m$Index, decreasing=T),]
          selcand_f <- selcand_f[order(selcand_f$Index, decreasing=T),]
        }
        if(select_on=="EBV"){
          selcand_m <- selcand_m[order(selcand_m[,c(11+4*Ntraits+trait)], decreasing=T),]
          selcand_f <- selcand_f[order(selcand_f[,c(11+4*Ntraits+trait)], decreasing=T),]
        }
        if(select_on=="Phenotype"){
          selcand_m <- selcand_m[order(selcand_m[,c(11+3*Ntraits+trait)], decreasing=T),]
          selcand_f <- selcand_f[order(selcand_f[,c(11+3*Ntraits+trait)], decreasing=T),]
        }

        selm <- selcand_m[1:Nm,]
        self <- selcand_f[1:Nf,]
        ped$selcand[!(ped$id %in% selm$id) & ped$generation %in% gen & ped$batch %in% batch & ped$sex==1 & ped$selcand==1 & ped$FSfam == fam] <- 2
        ped$selcand[!(ped$id %in% self$id) & ped$generation %in% gen & ped$batch %in% batch & ped$sex==2 & ped$selcand==1 & ped$FSfam == fam] <- 2
    }

    }else{

      ## Within family selection, no Nm and Nf specifies
      allselcand <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,]

      for(fam in (min(allselcand$FSfam):max(allselcand$FSfam))){

        selcand <- allselcand[allselcand$FSfam == fam]

      if(select_on=="Index"){selcand <- selcand[order(selcand$Index, decreasing=T),]}
      if(select_on=="EBV"){selcand <- selcand[order(selcand[,c(11+4*Ntraits+trait)], decreasing=T),]}
      if(select_on=="Phenotype"){selcand <- selcand[order(selcand[,c(11+3*Ntraits+trait)], decreasing=T),] }

        sel <- selcand[1:N,]
        ped$selcand[!(ped$id %in% sel$id) & ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1 & ped$FSfam == fam] <- 2
      }
    }
  }

   return(ped)
}



#' Selection
#'
#' Function to select selection candidates based on their Index, EBV or phenotype. ped$selcand should be 1 for the selection candidates.
#'
#' @param gen The generation of the selection candidates. A vector with multiple generations can be provided.
#' @param batch The batch of the selection candidates. Default is 0. It is possible to provide a vector with multiple batches.
#' @param select_on Options: "Index", "EBV" or "Phenotype". For "EBV" and "Phenotype", also the trait need to be specified in 'trait'. Default is Index.
#' @param trait Which trait or EBV the selection is base on when option "EBV" or "Phenotype" is choosen in 'select_on'.
#' @param Ntraits Number of simulated traits. Does not need to be specified if Ntraits is specified in a list called 'BPdata'.
#' @param Nm Number of males to select.
#' @param Nf Number of females to select.
#' @param max_FSfam Maximum number of sibs that can be selected per full sib familie. Default is 'all'.
#' @param mature_m Proportion of male selection candidates that is assumed to be mature and available at the moment of selection. Default is 1
#' @param mature_f Proportion of female selection candidates that is assumed to be mature and available at the moment of selection. Default is 1
#' @param selected The value assigned to ped$selected for the selected animals, default is 1.
#' @return This function will change the data frame called 'ped'. Fish that are selected will be assigned a value to their column 'selected'.
#' @export
#' @examples
#' \donttest{ped <- founderpopfam(Nm = 60,
#'                      Nf = 60,
#'                      Nm2 = 0,
#'                      Nf2 = 0,
#'                      Ntraits = 2,
#'                      TraitsIndex = 2,
#'                      Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                      Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                      Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                      mean = c(50,500),
#'                      a_var = c(200,8000),
#'                      c_var = c(150,1000),
#'                      e_var = c(250,12000))
#'                      
#' Mating <- randommating(gen = 0,
#'                        Nfam_FS = 120)
#'
#' for(fam in 1: nrow(Mating)){
#' ped <- offspringFSfam(gen = 1,
#'                       No = 100,
#'                       probmale = 0.5,
#'                       sire = Mating$Sire[fam],
#'                       dam = Mating$Dam[fam],
#'                       Ntraits = 2,
#'                       TraitsIndex = 2,
#'                       Rgen = matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                       Rcom = matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                       Rres = matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                       a_var = c(200,8000),
#'                       c_var = c(150,1000),
#'                       e_var = c(250,12000))
#' }
#' 
#' ped <- preselphen(gen = 1,
#'                   Nenv = 2,
#'                   Npresel = c(25,15),
#'                   trait = 1,
#'                   Ntraits = 2)
#'                   
#' ped <- avail_selection(gen = 1,
#'                        presel = 1,
#'                        surv = 0.9)
#'                        
#' ped <- breeding_values(gen = 1,
#'                        TraitsIndex = 2,
#'                        EBV = "GEBV",
#'                        GenomLength = 11.3,
#'                        Ne = 100,
#'                        SizeTraining = nrow(ped[ped$preselected ==2,]),
#'                        Ntraits = 2,
#'                        a_var = c(200,8000),
#'                        h2 = c(0.33,0.38))
#'                        
#' ped <- select(gen=1,
#'                  Nm = 60,
#'                  Nf = 60,
#'                  mature_m = 0.5,
#'                  mature_f = 0.4,
#'                  Ntraits = 2)
#'}

select <- function(gen, batch=0, select_on="Index", trait, Ntraits=BPdata$Ntraits, Nm=BPdata$Nm, Nf=BPdata$Nf, max_FSfam="all", mature_m=1, mature_f=1, selected=1){
  avail_selcand_m <- ped[ped$generation %in% gen &ped$batch %in% batch & ped$sex==1 & ped$selcand==1 & ped$selected==0,]
  avail_selcand_m <- sample(avail_selcand_m$id, mature_m*nrow(avail_selcand_m))
  avail_selcand_f <- ped[ped$generation %in% gen &ped$batch %in% batch & ped$sex==2 & ped$selcand==1 & ped$selected==0,]
  avail_selcand_f <- sample(avail_selcand_f$id, mature_f*nrow(avail_selcand_f))

  selcand_m <- ped[ped$id %in% avail_selcand_m,]
  selcand_f <- ped[ped$id %in% avail_selcand_f,]

    if(max_FSfam=="all"){
      if(select_on=="Index"){
        selcand_m <- selcand_m[order(selcand_m$Index, decreasing=T),]
        selcand_f <- selcand_f[order(selcand_f$Index, decreasing=T),]
      }
      if(select_on=="EBV"){
        selcand_m <- selcand_m[order(selcand_m[,c(11+4*Ntraits+trait)], decreasing=T),]
        selcand_f <- selcand_f[order(selcand_f[,c(11+4*Ntraits+trait)], decreasing=T),]
      }
      if(select_on=="Phenotype"){
        selcand_m <- selcand_m[order(selcand_m[,c(11+3*Ntraits+trait)], decreasing=T),]
        selcand_f <- selcand_f[order(selcand_f[,c(11+3*Ntraits+trait)], decreasing=T),]
      }

    selm <- selcand_m[1:Nm,]
    self <- selcand_f[1:Nf,]
    ped$selected[ped$id %in% selm$id] <- selected
    ped$selected[ped$id %in% self$id] <- selected
  }else{
    selcand <- rbind(selcand_m, selcand_f)
    if(select_on=="Index"){selcand <- selcand[order(selcand$Index, decreasing=T),]}
    if(select_on=="EBV"){selcand <- selcand[order(selcand[,c(11+4*Ntraits+trait)], decreasing=T),]}
    if(select_on=="Phenotype"){selcand <- selcand[order(selcand[,c(11+3*Ntraits+trait)], decreasing=T),] }

    for(i in 1:(Nm + Nf)){
      ped$selected[ped$id %in% selcand$id[i]] <- selected
      selectedFS <- ped[ped$selected==selected & ped$generation %in% gen & ped$batch %in% batch & ped$FSfam==selcand$FSfam[i],]
      selected_m <- ped[ped$selected==selected & ped$generation %in% gen & ped$batch %in% batch & ped$sex==1,]
      selected_f <- ped[ped$selected==selected & ped$generation %in% gen & ped$batch %in% batch & ped$sex==2,]

      if(nrow(selectedFS) == max_FSfam){
        sel_candidates1 <- selcand[1:i,]
        sel_candidates2 <- selcand[(i+1):nrow(selcand),]
        sel_candidates2 <- sel_candidates2[sel_candidates2$FSfam != selcand$FSfam[i],]
        selcand <-  rbind(sel_candidates1,sel_candidates2)
      }

      if(nrow(selected_m) == Nm){
        sel_candidates1 <- selcand[1:i,]
        sel_candidates2 <- selcand[(i+1):nrow(selcand),]
        sel_candidates2 <- sel_candidates2[sel_candidates2$sex == 2,]
        selcand <-  rbind(sel_candidates1,sel_candidates2)
      }
      if(nrow(selected_f) == Nf){
        sel_candidates1 <- selcand[1:i,]
        sel_candidates2 <- selcand[(i+1):nrow(selcand),]
        sel_candidates2 <- sel_candidates2[sel_candidates2$sex == 1,]
        selcand <-  rbind(sel_candidates1,sel_candidates2)
      }

      }
   }

  return(ped)
}


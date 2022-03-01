### Functies: 1) offspringFSfam -> maken van 1 full sib familie, dus for(1:nrow(Mating)){offspringFSfam.....}. Voor family design
###           1b)offspringFSgroup -> zelfde maar dan voor group mating design (enige verschil FSfam -> contribution)
###           2) randommating -> voor family design allocating sires to dams. Output tabel Mating
###           3) matinggroup -> Met group mating, welke full sib families maken, met contribution etc



#' Creating offspring for family design
#'
#' This function can be used to create the offspring of one full sib family in a family design, with genetic, common environmental and residual effects for each offspring. Offspring are added to the ped file.
#' @param gen The generation of the offspring.
#' @param No The number of offspring per fullsib family.
#' @param sire The sire of the full sib family. Sire should also be in the ped file with information on the inbreeding level and true breeding values for each trait.
#' @param dam The dam of the full sib family. Dam should also be in the ped file with information on the inbreeding level and the true breeding values for each trait.
#' @param batch The batch of the offspring. Default is 0.
#' @param probmale The probability that a offspring is male. The probability that the offspring is female is calculated as 1 - probmale. Probmale does not need to be specified if it is in the list called 'BPdata'. The default is 0.5.
#' @param Ntraits Number of traits to be simulated. Does not need to be specified if Ntraits is in the list 'BPdata'.
#' @param TraitsIndex Vector of traits that are in the selection index. By default, all traits are in the index.
#' @param Rgen Matrix of all genetic correlations between all Ntraits. Only needs to be specified if there is no matrix of genetic correlations named Rgen in the list called 'BPdata'.
#' @param Rcom Matrix of all common environmental correlations between all Ntraits. Only needs to be specified if there is no matrix of common environmental correlations named Rcom in the list called 'BPdata'.
#' @param Rres Matrix of all residual correlations between all Ntraits. Only needs to be specified if there is no matrix of residual correlations named Rres in the list called 'BPdata'.
#' @param a_var Vector of genetic variances of all traits. Only needs to be specified if there is no vector of genetic variances named a_var in the list called 'BPdata'.
#' @param c_var Vector of common environmental variances of all traits. Only needs to be specified if there is no vector of common environmental variances named c_var in the list called 'BPdata'.
#' @param e_var Vector of residual variances of all traits. Only needs to be specified if there is no vector of residual variances named e_var in the list called 'BPdata'.
#' @param inbreeding If TRUE (default), then the inbreeding level is calculated for each offspring.
#' @export
#' @return This function returns the 'ped' data frame: the new offspring are added to this data frame.
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
#' }


offspringFSfam <- function(gen, No, sire, dam, batch=0, probmale=BPdata$prob_male, Ntraits=BPdata$Ntraits,  TraitsIndex=c(1:Ntraits), Rgen=BPdata$Rgen,Rres=BPdata$Rres,Rcom=BPdata$Rcom, a_var=BPdata$a_var, c_var = BPdata$c_var, e_var=BPdata$e_var, inbreeding=TRUE){
  if(is.null(probmale)){probmale <- 0.5}

  if(No==1){
    offspring <- data.frame( id = (max(ped$id)+(1:(5*No))) ,
                           sire  = sire ,
                           dam = dam,
                           sex = sample(c(1,2), size = 5*No, replace = T, prob=c(probmale,(1-probmale))),
                           generation = gen,
                           batch = batch,
                           inbreeding = NA,
                           FSfam = max(ped$FSfam)+1,
                           selected = 0,
                           preselected = 0,
                           selcand = 0)
  }else{
    offspring <- data.frame( id = (max(ped$id)+(1:(No))) ,
                             sire  = sire ,
                             dam = dam,
                             sex = sample(c(1,2), size = No, replace = T, prob=c(probmale,(1-probmale))),
                             generation = gen,
                             batch = batch,
                             inbreeding = NA,
                             FSfam = max(ped$FSfam)+1,
                             selected = 0,
                             preselected = 0,
                             selcand = 0)
  }

  offspring <- cbind(offspring, MASS::mvrnorm(n = nrow(offspring), mu = rep(0, nrow(Rgen)), Sigma = Rgen))

  colnames <- c("id",
                "sire",
                "dam",
                "sex",
                "generation",
                "batch",
                "inbreeding",
                "FSfam",
                "selected",
                "preselected",
                "selcand")

  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("g_trait",traits), collapse=""))
  }

  names(offspring)<- colnames


  FSfam <- data.frame(FSfam=max(ped$FSfam)+1)
  FSfam <- cbind(FSfam, matrix(MASS::mvrnorm(n = nrow(FSfam), mu = rep(0, Ntraits), Sigma = Rcom),nrow=nrow(FSfam)))
  offspring<- merge(offspring, FSfam, by="FSfam", all = T)

  offspring <- (cbind(offspring, MASS::mvrnorm(n = nrow(offspring), mu = rep(0, Ntraits), Sigma = Rres)))


  offspring <- cbind(offspring,matrix(c(rep(0, times=nrow(offspring)*(Ntraits))),nrow=nrow(offspring)))
  offspring <- cbind(offspring,matrix(c(rep(0, times=nrow(offspring)*length(TraitsIndex))),nrow=nrow(offspring)))

  colnames <- c("FSfam",
                "id",
                "sire",
                "dam",
                "sex",
                "generation",
                "batch",
                "inbreeding",
                "selected",
                "preselected",
                "selcand")

  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("g_trait",traits), collapse=""))
  }

  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("c_trait",traits), collapse=""))
  }
  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("e_trait",traits), collapse=""))
  }

  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("p_trait",traits), collapse=""))
  }
  for(traits in 1:length(TraitsIndex)){
    colnames <- c(colnames, paste(c("EBV",traits), collapse=""))
  }

  names(offspring)<- colnames
  offspring$Index <- 0


  offspring <- (cbind(offspring, MASS::mvrnorm(n = nrow(offspring), mu = rep(0, Ntraits), Sigma = Rgen)))
  offspring <- cbind(offspring,matrix(c(rep(0, times=nrow(offspring)*(Ntraits))),nrow=nrow(offspring)))

  colnames <- c(colnames, "Index")
  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("MS",traits), collapse=""))
  }
  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("g1_trait",traits), collapse=""))
  }

  names(offspring)<- colnames

  for(trait in 1: Ntraits){
    offspring[,c(11+4*Ntraits+1+length(TraitsIndex)+trait)] <- offspring[,c(11+4*Ntraits+1+length(TraitsIndex)+trait)]* sqrt(a_var[trait]*0.5*(1-0.5*(ped$inbreeding[match(offspring$sire, ped$id)] + ped$inbreeding[match(offspring$dam, ped$id)])))
    offspring[,c(11+5*Ntraits+1+length(TraitsIndex)+trait)] <- (ped[match(offspring$sire, ped$id),c(11+trait) ]/2 +
                                                                  ped[match(offspring$dam, ped$id),c(11+trait)]/2)
  }

  for(trait in 1:Ntraits){
    #print(trait)
    offspring[,c(11+trait)] <- offspring[,c(11+4*Ntraits+1+length(TraitsIndex)+trait)] + offspring[,c(11+5*Ntraits+1+length(TraitsIndex)+trait)]
    offspring[,c(11+Ntraits+trait)] <- offspring[,c(11+Ntraits+trait)]*sqrt(c_var[trait])
    offspring[,c(11+2*Ntraits+trait)] <- offspring[,c(11+2*Ntraits+trait)]*sqrt(e_var[trait])
    offspring[,c(11+3*Ntraits+trait)] <- offspring[,c(11+trait)] + offspring[,c(11+Ntraits+trait)] + offspring[,c(11+2*Ntraits+trait)]
  }


  if(No==1){
    offspring <- offspring[1,c(1:(11+4*Ntraits+1+length(TraitsIndex)))]
  }else{
    offspring <- offspring[,c(1:(11+4*Ntraits+1+length(TraitsIndex)))]
  }

  ped <- rbind(ped, offspring)

  if(inbreeding == T){ped$inbreeding <- pedigree::calcInbreeding(ped)}

  return(ped)
}

#' Creating offspring for a group mating design
#'
#' This function can be used to create the offspring of one full sib family in a group mating design, with genetic, common environmental and residual effects for each offspring. Offspring are added to the ped file.
#' @param gen The generation of the offspring.
#' @param No The number of offspring per fullsib family.
#' @param sire The sire of the full sib family. Sire should also be in the ped file with information on the inbreeding level and true breeding values for each trait.
#' @param dam The dam of the full sib family. Dam should also be in the ped file with information on the inbreeding level and the true breeding values for each trait.
#' @param batch The batch of the offspring. Default is 0.
#' @param probmale The probability that a offspring is male. The probability that the offspring is female is calculated as 1 - probmale. Probmale does not need to be specified if it is in the list called 'BPdata'. The default is 0.5.
#' @param Ntraits Number of traits to be simulated. Does not need to be specified if Ntraits is in the list 'BPdata'.
#' @param TraitsIndex Vector of traits that are in the selection index. By default, all traits are in the index.
#' @param Rgen Matrix of all genetic correlations between all Ntraits. Only needs to be specified if there is no matrix of genetic correlations named Rgen in the list called 'BPdata'.
#' @param Rcom Matrix of all common environmental correlations between all Ntraits. Only needs to be specified if there is no matrix of common environmental correlations named Rcom in the list called 'BPdata'.
#' @param Rres Matrix of all residual correlations between all Ntraits. Only needs to be specified if there is no matrix of residual correlations named Rres in the list called 'BPdata'.
#' @param a_var Vector of genetic variances of all traits. Only needs to be specified if there is no vector of genetic variances named a_var in the list called 'BPdata'.
#' @param c_var Vector of common environmental variances of all traits. Only needs to be specified if there is no vector of common environmental variances named c_var in the list called 'BPdata'.
#' @param e_var Vector of residual variances of all traits. Only needs to be specified if there is no vector of residual variances named e_var in the list called 'BPdata'.
#' @param inbreeding If True (default), then the inbreeding level is calculated for each offspring.
#' @return This function returns the 'ped' data frame: the new offspring are added to this data frame.
#' @export
#' @examples
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
#'                        batch = -3,
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
#'                       Ntraits =2,
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
#' }


offspringFSgroup <- function(gen, No, sire, dam, batch=0, probmale=BPdata$prob_male, Ntraits=BPdata$Ntraits,  TraitsIndex=c(1:Ntraits), Rgen=BPdata$Rgen,Rres=BPdata$Rres,Rcom=BPdata$Rcom, a_var=BPdata$a_var, c_var = BPdata$c_var, e_var=BPdata$e_var, inbreeding=TRUE){
  if(is.null(probmale)){probmale <- 0.5}

  if(No==1){
    offspring <- data.frame( id = (max(ped$id)+(1:(5*No))) ,
                             sire  = sire ,
                             dam = dam,
                             sex = sample(c(1,2), size = 5*No, replace = T, prob=c(probmale,(1-probmale))),
                             generation = gen,
                             batch = batch,
                             inbreeding = NA,
                             Contribution = 0,
                             selected = 0,
                             preselected = 0,
                             selcand = 0)

    offspring <- cbind(offspring, MASS::mvrnorm(n = nrow(offspring), mu = rep(0, nrow(Rgen)), Sigma = Rgen))
    offspring <- cbind(offspring, MASS::mvrnorm(n = nrow(offspring), mu = rep(0, nrow(Rcom)), Sigma = Rcom))
    offspring <- cbind(offspring, MASS::mvrnorm(n = nrow(offspring), mu = rep(0, Ntraits), Sigma = Rres))
    offspring <- cbind(offspring,matrix(c(rep(0, times=nrow(offspring)*(Ntraits))),nrow=nrow(offspring)))
    offspring <- cbind(offspring,matrix(c(rep(0, times=nrow(offspring)*length(TraitsIndex))),nrow=nrow(offspring)))
  }else{
    offspring <- data.frame( id = (max(ped$id)+(1:(No))) ,
                             sire  = sire ,
                             dam = dam,
                             sex = sample(c(1,2), size = No, replace = T, prob=c(probmale,(1-probmale))),
                             generation = gen,
                             batch = batch,
                             inbreeding = NA,
                             Contribution = 0,
                             selected = 0,
                             preselected = 0,
                             selcand = 0)

    offspring <- cbind(offspring, MASS::mvrnorm(n = nrow(offspring), mu = rep(0, nrow(Rgen)), Sigma = Rgen))
    offspring <- cbind(offspring, MASS::mvrnorm(n = nrow(offspring), mu = rep(0, nrow(Rcom)), Sigma = Rcom))
    offspring <- cbind(offspring, MASS::mvrnorm(n = nrow(offspring), mu = rep(0, Ntraits), Sigma = Rres))
    offspring <- cbind(offspring,matrix(c(rep(0, times=nrow(offspring)*(Ntraits))),nrow=nrow(offspring)))
    offspring <- cbind(offspring,matrix(c(rep(0, times=nrow(offspring)*length(TraitsIndex))),nrow=nrow(offspring)))
    }


  colnames <- c("id",
                "sire",
                "dam",
                "sex",
                "generation",
                "batch",
                "inbreeding",
                "Contribution",
                "selected",
                "preselected",
                "selcand")

  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("g_trait",traits), collapse=""))
  }

  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("c_trait",traits), collapse=""))
  }
  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("e_trait",traits), collapse=""))
  }

  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("p_trait",traits), collapse=""))
  }
  for(traits in 1:length(TraitsIndex)){
    colnames <- c(colnames, paste(c("EBV",traits), collapse=""))
  }

  names(offspring)<- colnames
  offspring$Index <- 0


  offspring <- (cbind(offspring, MASS::mvrnorm(n = nrow(offspring), mu = rep(0, Ntraits), Sigma = Rgen)))
  offspring <- cbind(offspring,matrix(c(rep(0, times=nrow(offspring)*(Ntraits))),nrow=nrow(offspring)))

  colnames <- c(colnames, "Index")
  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("MS",traits), collapse=""))
  }
  for(traits in 1:Ntraits){
    colnames <- c(colnames, paste(c("g1_trait",traits), collapse=""))
  }

  names(offspring)<- colnames

  for(trait in 1: Ntraits){
    offspring[,c(11+4*Ntraits+1+length(TraitsIndex)+trait)] <- offspring[,c(11+4*Ntraits+1+length(TraitsIndex)+trait)]* sqrt(a_var[trait]*0.5*(1-0.5*(ped$inbreeding[match(offspring$sire, ped$id)] + ped$inbreeding[match(offspring$dam, ped$id)])))
    offspring[,c(11+5*Ntraits+1+length(TraitsIndex)+trait)] <- (ped[match(offspring$sire, ped$id),c(11+trait) ]/2 +
                                                                  ped[match(offspring$dam, ped$id),c(11+trait)]/2)
  }

  for(trait in 1:Ntraits){
    #print(trait)
    offspring[,c(11+trait)] <- offspring[,c(11+4*Ntraits+1+length(TraitsIndex)+trait)] + offspring[,c(11+5*Ntraits+1+length(TraitsIndex)+trait)]
    offspring[,c(11+Ntraits+trait)] <- offspring[,c(11+Ntraits+trait)]*sqrt(c_var[trait])
    offspring[,c(11+2*Ntraits+trait)] <- offspring[,c(11+2*Ntraits+trait)]*sqrt(e_var[trait])
    offspring[,c(11+3*Ntraits+trait)] <- offspring[,c(11+trait)] + offspring[,c(11+Ntraits+trait)] + offspring[,c(11+2*Ntraits+trait)]
  }


  if(No==1){
    offspring <- offspring[1,c(1:(11+4*Ntraits+1+length(TraitsIndex)))]
  }else{
    offspring <- offspring[,c(1:(11+4*Ntraits+1+length(TraitsIndex)))]
    }

  ped <- rbind(ped, offspring)
  if(inbreeding == T){ped$inbreeding <- pedigree::calcInbreeding(ped)}

  return(ped)
}



#' Random mating family design
#'
#' Function to randomly allocate sires to dams.
#' 
#' @details
#' A dataframe called ped needs to be present in the data. Ped needs to contain all the sires and dams that need to be allocated to each other and the columns sex, selected, generation, batch and id (first column).
#'
#' The sires and dams can come from multiple batches or generations. In that case, a vector of batches or generations need to be provided.
#'
#' Optionally, a column with the number of offspring per full sib family can be added to the dataframe. To do this, either the argument No_FSfam or No needs to be provided. In case No_FSfam is provided, each mating will have the same number of offspring (namely the value provided with No_FSfam). When No is provided, then the total number of offspring will be evenly divided among each mating if possible.
#' @param gen The generations of the sires and dams.
#' @param batch The batch of the sires and dams. Default is 0.
#' @param batch_m The batch of the sires. Default is NA. If batch_m is specified, batch_f also needs to be specified and the parameter batch is not used.
#' @param batch_f The batch of the dams. Default is NA. If batch_f is specified, batch_m also needs to be specified and the parameter batch is not used.
#' @param Nfam_FS The number of full sib families.
#' @param No_FSfam The number of offspring in each full sib family. Default is NA.
#' @param No The total number of offspring of all matings. Default is NA.
#' @param selected The value in ped$selected of the selected sires and dams. Default is 1.
#' @return The output is a data frame with for each full sib family the sire and dam and (optional) the number of offspring per full sib family.
#' @export
#' @examples
#'{ ped <- founderpopfam(Nm=60,
#'                      Nf=60,
#'                      Nm2=0,
#'                      Nf2=0,
#'                      Ntraits=2,
#'                      TraitsIndex = 2,
#'                      Rgen= matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                      Rcom= matrix(c(1.00  , 0.5,
#'                                     0.5   , 1.00),
#'                                  nrow = 2),
#'                      Rres= matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                      mean=c(50,500),
#'                      a_var=c(200,8000),
#'                      c_var=c(150,1000),
#'                      e_var= c(250,12000))
#'                      
#' Mating <- randommating(gen = 0,
#'                        Nfam_FS = 120)
#'}

randommating <- function(gen, batch = 0, batch_m = NA, batch_f = NA, Nfam_FS, No = NA, No_FSfam = NA, selected = 1){
  if(is.na(batch_m)){
    sires <- ped[ped$sex==1 &  ped$generation %in% gen & ped$batch %in% batch & ped$selected==selected, c(1,4)]
    dams <- ped[ped$sex==2 & ped$generation %in% gen & ped$batch %in% batch & ped$selected==selected, c(1,4)]
  }else{
    sires <- ped[ped$sex==1 & ped$generation %in% gen & ped$batch %in% batch_m & ped$selected==selected, c(1,4)]
    dams <- ped[ped$sex==2 & ped$generation %in% gen & ped$batch %in% batch_f & ped$selected==selected, c(1,4)]
  }

  Nmating_f <- Nfam_FS/nrow(dams)
  Nmating_m <- Nfam_FS/nrow(sires)

  dams$n <- c(sample(c(rep((1:(nrow(dams)/Nmating_m)), times=Nmating_m)), nrow(dams)))
  sires$n <- c(sample(c(rep((1:(nrow(sires)/Nmating_f)), times=Nmating_f)), nrow(sires)))

  curr_fam <- 1
  fams <- sires$id[sires$n==curr_fam]
  famd <- dams$id[dams$n==curr_fam]

  Mating <- expand.grid(fams, famd)

  for(curr_fam in 2:Nfam_FS){
    fams <- sires$id[sires$n==curr_fam]
    famd <- dams$id[dams$n==curr_fam]

    Matingnew <- expand.grid(fams, famd)
    Mating <- rbind(Mating,Matingnew)
  }

  names(Mating) <- c("Sire", "Dam")

  if(!is.na(No_FSfam)){Mating$No <- No_FSfam }

  if(!is.na(No)){
    Mating$No <- 0
    for(n in 1:nrow(Mating)){
      Mating$No[n] <- round((No-sum(Mating$No))/(nrow(Mating)-(n-1)))
    }
  }

  return(Mating)
}



#' Group mating
#'
#' Function to determine which full sib families are produced in a group mating design.
#' @details
#' By default, the contribution of the sires and dams that do reproduce come from a gamma distribution. The default shape and scale of the gamma distribution are 0.75 and 0.11, respectively. A uniform distribution can also be specified for the contributions of the sires and dams.
#' If not all sires and dams should contribute to the offspring, then the sires and dams that are going to reproduce are randomly chosen.
#' The output is a dataframe called Mating with the sire, dam and the size of each full sib family.
#'
#' The sires and dams can come from multiple batches or generations. In that case, a vector of batches or generations need to be provided.
#' @param gen The generations of the sires and dams
#' @param batch The batches of the sires and dams. Default is 0.
#' @param batch_m The batch of the sires. Default is NA. If batch_m is specified, batch_f also needs to be specified and the parameter batch is not used.
#' @param batch_f The batch of the dams. Default is NA. If batch_f is specified, batch_m also needs to be specified and the parameter batch is not used.
#' @param No The total number of offspring of all matings.
#' @param contr_m Proportion of sires that contribute to the offspring.
#' @param contr_f Proportion of dams that contribute to the offspring.
#' @param distribution The distribution from which the contributions are drawn. Options are: "Gamma" (default) and "Uniform"
#' @param shape The shape of the gamma distribution. Default is 0.75.
#' @param scale The scale of the gamma distribution. Default is 0.11.
#' @param selected The value in ped$selected of the selected sires and dams. Default is 1.
#' @return The output is a data frame with for each full sib family the sire and dam and the number of offspring per full sib family.
#' @export
#' @examples
#' 
#' {ped <- founderpopgroup(Nm=60,
#'                      Nf=60,
#'                      Nm2=120,
#'                      Nf2=120,
#'                      Nbatch = 4,
#'                      batch2 = c(-3,-2,-1,0),
#'                      Ntraits=2,
#'                      TraitsIndex = 2,
#'                      Rgen= matrix(c(1.00   , 0.48,
#'                                     0.48   , 1.00),
#'                                  nrow = 2),
#'                      Rcom= matrix(c(1.00  , 0,
#'                                     0   , 1.00),
#'                                  nrow = 2),
#'                      Rres= matrix(c(1.00   , 0.32,
#'                                     0.32   , 1.00),
#'                                  nrow = 2),
#'                      mean=c(50,500),
#'                      a_var=c(200,8000),
#'                      c_var=c(0,0),
#'                      e_var= c(250,12000))
#'                      
#'  Mating <- groupmating(gen = 0,
#'                        batch=-3,
#'                        No=1000,
#'                        contr_m = 0.5,
#'                        contr_f = 0.5)
#'}


groupmating <- function(gen, batch=0, batch_m=NA, batch_f=NA, No, contr_m, contr_f, distribution = "Gamma", shape=0.75, scale=0.11, selected = 1){
  if(is.na(batch_m)){
    sires <- ped[ped$sex==1 &  ped$generation %in% gen & ped$batch %in% batch & ped$selected==selected, c(1,4,8)]
    dams <- ped[ped$sex==2 & ped$generation %in% gen & ped$batch %in% batch & ped$selected==selected, c(1,4,8)]
  }else{
    sires <- ped[ped$sex==1 & ped$generation %in% gen & ped$batch %in% batch_m & ped$selected==selected, c(1,4,8)]
    dams <- ped[ped$sex==2 & ped$generation %in% gen & ped$batch %in% batch_f & ped$selected==selected, c(1,4,8)]
  }

  if(distribution == "Gamma"){
    sires$Contribution <- stats::rgamma(n=nrow(sires) ,shape = shape, scale = scale)
    non_mature_m <- sample(sires$id, (1-contr_m)*nrow(sires))
    sires$Contribution[sires$id %in% non_mature_m] <- 0
    sires$Contribution <- sires$Contribution/sum(sires$Contribution)

    dams$Contribution <- stats::rgamma(n=nrow(dams) ,shape = shape, scale = scale)
    non_mature_f <- sample(dams$id, (1-contr_f)*nrow(dams))
    dams$Contribution[dams$id %in% non_mature_f] <- 0
    dams$Contribution <- dams$Contribution/sum(dams$Contribution)
  }

  if(distribution == "Uniform"){
    sires$Contribution <- 1
    non_mature_m <- sample(sires$id, (1-contr_m)*nrow(sires))
    sires$Contribution[sires$id %in% non_mature_m] <- 0
    sires$Contribution <- sires$Contribution/sum(sires$Contribution)

    dams$Contribution <- 1
    non_mature_f <- sample(dams$id, (1-contr_f)*nrow(dams))
    dams$Contribution[dams$id %in% non_mature_f] <- 0
    dams$Contribution <- dams$Contribution/sum(dams$Contribution)
  }

  Mating <- expand.grid(sires$id, dams$id)
  Mating$offspring <- kronecker(sires$Contribution, dams$Contribution)
  Mating$offspring <- round((Mating$offspring/sum(Mating$offspring))*No)
  Mating <- Mating[Mating$offspring >0,]
  names(Mating) <- c("Sire", "Dam", "No")

  return(Mating)
}




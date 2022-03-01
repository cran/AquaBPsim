#' Founder population family design
#'
#' This function can be used to create a founder population of unrelated animals for a breeding program with a family design.
#' @param Nm Number of males that will be used to breed the first generation of offspring. Males are listed first in the output data frame and coded with 1 for sex.
#' @param Nf Number of females that will be used to breed the first generation of offspring. Females are listed second in the output data frame and coded with 2 for sex.
#' @param Nm2 Additional males that will not be used to breed the first generation of offspring, but can be used as additional selection candidates in further generations. Default is zero.
#' @param Nf2 Additional females that will not be used to breed the first generation of offspring, but can be used as additional selection candidates in further generations. Default is zero.
#' @param Nbatch Number of batches over which the founder animals are divided. The number of founder animals that need to be simulated should be a multiple of the number of batches. If Nbatch is not specified, the parameter batch can be used to divide the founder animals over batches. If both Nbatch and batch are not specified, then all founder animals will be assigned to batch 0.
#' @param Nbatch2 Number of batches over which the additional founder animals (Nm2 + Nf2) are divided. The number of additional founder animals that need to be simulated should be a multiple of the number of batches. If Nbatch2 is not specified, the parameter batch2 can be used to divide the additional founder animals over batches. If both Nbatch2 and batch2 are both not specified, then the additional founder animals are divided over batches in the same way as the other founder animals (Nm+Nf).
#' @param batch A vector with names of the batches. The number of founder animals that need to be simulated should be the same or a multiple of the length of the vector. If batch is not specified, the parameter Nbatch can be used to divide the founder animals over batches. If both Nbatch and batch are not specified, then all founder animals will be assigned to batch 0.
#' @param batch2 A vector with names of the batches for the additional founder animals (Nm2 + Nf2). The number of additional founder animals that are simulated should be the same or a multiple of the length of the vector. If batch2 is not specified, the parameter Nbatch2 can be used to divide the additional founder animals over batches. If both Nbatch2 and batch2 are not specified, then the additional founder animals are divided over batches in the same way as the other founder animals (Nm+Nf).
#' @param Ntraits Number of traits to be simulated. Does not need to be specified if Ntraits is in the list 'BPdata'.
#' @param TraitsIndex Vector of traits that are in the index. By default, all traits are in the selection index.
#' @param Rgen Matrix of all genetic correlations between all Ntraits. Only needs to be specified if there is no matrix of genetic correlations named Rgen in the list called 'BPdata'.
#' @param Rcom Matrix of all common environmental correlations between all Ntraits. Only needs to be specified if there is no matrix of common environmental correlations named Rcom in the list called 'BPdata'.
#' @param Rres Matrix of all residual correlations between all Ntraits. Only needs to be specified if there is no matrix of residual correlations named Rres in the list called 'BPdata'.
#' @param mean Vector of means of all traits. Only needs to be specified if there is no vector of means named mean in the list 'BPdata'.
#' @param a_var Vector of genetic variances of all traits. Only needs to be specified if there is no vector of genetic variances named a_var in the list 'BPdata'.
#' @param c_var Vector of common environmental variances of all traits. Only needs to be specified if there is no vector of common environmental variances named c_var in the list 'BPdata'. If there is no common environmental effect, provide a vector of zero's.
#' @param e_var Vector of residual variances of all traits. Only needs to be specified if there is no vector of residual variances named e_var in the list 'BPdata'.
#' @param est_EBV TRUE or FALSe for estimating breeding values for all founder animals. The default is FALSE.
#' @param EBV If est_EBV is TRUE, a vector of methods for simulating EBVs for each trait in the selection index need to be specified. There are four options: 0 for giving each animal a breeding values of 0 for the specific trait; "mean_pop" for simulating EBVs for all founder populations that are equal to the mean of the population, the mean of the population is provided in the vector 'mean'; "pheno" for simulating EBVs that are equal to the phenotype of the animal for the specific trait and "EBV" for simulating EBVs as a value correlated to the true breeding value of the animal. This correlation equals the accuracy. Accuracies need to be specified in the parameter 'accuracy'.
#' @param accuracy If for one of the traits the EBVs need to be calculated with the method "EBV", specified in the parameter 'EBV', then accuracies need to be provided. A value need to be added for each trait, however this value can be NA for traits from which the EBVs are not simulated using method "EBV".
#' @param indexweights If traits need to be combined in an index, then desired gain indices need to be specified for each trait in the index.
#' @return A data frame is returned, which should be called 'ped' in order to be able to use it in the other functions. The 'ped' data frame consist of all the simulated base animals, including their sex, generation, batch number, inbreeding level, phenotypes and genetic, common environmental and residual effects of each trait.
#' @export
#' @examples
#'
#' ped <- founderpopfam(Nm=60, Nm2=0,
#'                      Nf=60, Nf2=0,
#'                      batch = c(0,1,2),
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
#'               
#' ped <- founderpopfam(Nm=60, 
#'                      Nf=60, 
#'                      Nm2=120,
#'                      Nf2=120,
#'                      Nbatch = 4,
#'                      batch2 = c(-3,-2,-1,0),
#'                      Ntraits=2,
#'                      TraitsIndex = c(1,2),
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
#'                      e_var= c(250,12000),
#'                      est_EBV = TRUE,
#'                      EBV= c("pheno", "EBV"),
#'                      accuracy= c(NA,0.78),
#'                      indexweight= c(1,5))
#'


founderpopfam <- function(Nm=BPdata$Nm, Nf=BPdata$Nf, Nm2=BPdata$Nm2, Nf2=BPdata$Nf2, Nbatch=NA, Nbatch2=NA, batch=NA, batch2=NA, Ntraits=BPdata$Ntraits, TraitsIndex=c(1:Ntraits), Rgen=BPdata$Rgen,Rres=BPdata$Rres,Rcom=BPdata$Rcom, mean=BPdata$mean, a_var=BPdata$a_var, c_var = BPdata$c_var, e_var=BPdata$e_var, est_EBV=FALSE, EBV, accuracy, indexweights=c(rep(1, Ntraits))){

  if(is.null(Nm2)){Nm2 <- 0}
  if(is.null(Nf2)){Nf2 <- 0}

  # Simulating the founders
  ped <- data.frame(id = 1:(Nm+Nf+Nm2+Nf2),
                    sire = 0,
                    dam = 0,
                    sex = c(rep(1, times = Nm), rep(2, times = Nf),rep(1, times = Nm2), rep(2, times = Nf2)),
                    generation = 0,
                    batch = 0,
                    inbreeding = 0,
                    FSfam = 0,
                    selected = c(rep(1,times=(Nm+Nf)), rep(0,times=Nm2+Nf2)),
                    preselected = c(rep(0,times=(Nm+Nf)), rep(1,times=Nm2+Nf2)),
                    selcand = c(rep(0,times=(Nm+Nf)), rep(1,times=Nm2+Nf2)))

  # Adding the right batch numbers/names
  if(is.na(Nbatch) & is.na(Nbatch2) & sum(is.na(batch))==0 & sum(is.na(batch2))>0){
    ped$batch <- rep(x = batch, times = (Nm+Nf+Nm2+Nf2)/length(batch))}

  if(is.na(Nbatch) & is.na(Nbatch2) & sum(is.na(batch))==0 & sum(is.na(batch2))==0){
    ped$batch <- c(rep(x = batch, times = (Nm+Nf)/length(batch)),rep(x = batch2, times = (Nm2+Nf2)/length(batch2)))}

  if(is.na(Nbatch) & !is.na(Nbatch2) & sum(is.na(batch))==0 & sum(is.na(batch2))>0){
    ped$batch <- c(rep(x = batch, times = (Nm+Nf)/length(batch)),rep(x = 1:Nbatch2, times = (Nm2+Nf2)/Nbatch2))}

  if(!is.na(Nbatch) & is.na(Nbatch2) & sum(is.na(batch))>0 & sum(is.na(batch2))>0){
    ped$batch <- rep(x = 1:Nbatch, times = (Nm+Nf+Nm2+Nf2)/Nbatch)}

  if(!is.na(Nbatch) & is.na(Nbatch2) & sum(is.na(batch))>0 & sum(is.na(batch2))==0){
    ped$batch <- c(rep(x = 1:Nbatch, times = (Nm+Nf)/Nbatch),rep(x = batch2, times = (Nm2+Nf2)/length(batch2)))}

  if(!is.na(Nbatch) & !is.na(Nbatch2) & sum(is.na(batch))>0 & sum(is.na(batch2))>0){
    ped$batch <- c(rep(x = 1:Nbatch, times = (Nm+Nf)/Nbatch),rep(x = 1:Nbatch2, times = (Nm2+Nf2)/Nbatch2))}

  # Adding normally distributed correlated values for the traits and colums for the phenotypes and EBVs for traits in the index
  ped <- cbind(ped, MASS::mvrnorm(n = nrow(ped), mu = rep(0, nrow(Rgen)), Sigma = Rgen))
  ped <- cbind(ped, MASS::mvrnorm(n = nrow(ped), mu = rep(0, nrow(Rcom)), Sigma = Rcom))
  ped <- cbind(ped, MASS::mvrnorm(n = nrow(ped), mu = rep(0, nrow(Rres)), Sigma = Rres))
  ped <- cbind(ped,matrix(c(rep(0, times=nrow(ped)*(Ntraits))),nrow=nrow(ped)))
  ped <- cbind(ped,matrix(c(rep(0, times=nrow(ped)*(length(TraitsIndex)))),nrow=nrow(ped)))

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

  names(ped) <- colnames

  ped$Index <- 0

  # # Adding rights variances and mean for genetic, residual and common environmental effects + calculating phenotype
  for(trait in 1:Ntraits){
    #print(trait)
    ped[,c(11+trait)] <- ped[,c(11+trait)]*sqrt(a_var[trait]) + mean[trait]
    ped[,c(11+Ntraits+trait)] <- ped[,c(11+Ntraits+trait)]*sqrt(c_var[trait])
    ped[,c(11+2*Ntraits+trait)] <- ped[,c(11+2*Ntraits+trait)]*sqrt(e_var[trait])
    ped[,c(11+3*Ntraits+trait)] <- ped[,c(11+trait)] + ped[,c(11+Ntraits+trait)] + ped[,c(11+2*Ntraits+trait)]
  }

  if(est_EBV==T){
  # If needed, breeding values are simulated, 4 options for that
    for(trait in 1:length(TraitsIndex)){
      if(EBV[trait]==0){ ped[,c(11+4*Ntraits+trait)]<-0 }
      if(EBV[trait]=="mean_pop"){ ped[,c(11+4*Ntraits+trait)]<- mean[TraitsIndex[trait]] }
      if(EBV[trait]=="pheno"){ ped[,c(11+4*Ntraits+trait)]<- ped[,c(11+3*Ntraits+TraitsIndex[trait])] }
      if(EBV[trait]=="EBV"){
        pe <- ped[,c(1,2,3)]
        pe$pe <- stats::rnorm(nrow(pe)) * sqrt(accuracy[trait]) * sqrt(1 - (accuracy[trait])) * sqrt(a_var[TraitsIndex[trait]])
        ped[,c(11+4*Ntraits+trait)] <- ((ped[,c(11+trait)]-mean[TraitsIndex[trait]])*accuracy[trait]) + pe$pe + mean[TraitsIndex[trait]]
      }
    }

    # If EBVs are simulated, then these EBVs are combined in an Index
    weights <- indexweights[1]/sqrt(a_var[TraitsIndex[1]])

    if(length(TraitsIndex) > 1){
      for(i in 2:length(TraitsIndex)){
        weights <- c(weights, indexweights[i]/sqrt(a_var[TraitsIndex[i]]))
      }
    }

      weights <- weights/weights[1]

      for(i in 1:length(TraitsIndex)){
        ped$Index <- ped$Index + weights[i]* ped[,c(11+4*Ntraits+i)]
      }
  }

  return(ped)
}






#' Founder population group mating design
#'
#' This function can be used to create a founder population of unrelated animals for a breeding program with group mating.
#' @param Nm Number of males that will be used to breed the first generation of offspring. Males are listed first in the output data frame and coded with 1 for sex.
#' @param Nf Number of females that will be used to breed the first generation of offspring. Females are listed second in the output data frame and coded with 2 for sex.
#' @param Nm2 Additional males that will not be used to breed the first generation of offspring, but can be used as additional selection candidates in futher generations. Default is zero.
#' @param Nf2 Additional females that will not be used to breed the first generation of offspring, but can be used as additional selection candidates in futher generations. Default is zero.
#' @param Nbatch Number of batches over which the founder animals are divided. The number of founder animals that need to be simulated should be a multiple of the number of batches. If Nbatch is not specified, the parameter batch can be used to divide the founder animals over batches. If both Nbatch and batch are not specified, then all founder animals will be assigned to batch 0.
#' @param Nbatch2 Number of batches over which the additional founder animals (Nm2 + Nf2) are divided. The number of additional founder animals that need to be simulated should be a multiple of the number of batches. If Nbatch2 is not specified, the parameter batch2 can be used to divide the additional founder animals over batches. If both Nbatch2 and batch2 are both not specified, then the additional founder animals are divided over batches in the same way as the other founder animals (Nm+Nf).
#' @param batch A vector with names of the batches. The number of founder animals that need to be simulated should be the same or a multiple of the length of the vector. If batch is not specified, the parameter Nbatch can be used to divide the founder animals over batches. If both Nbatch and batch are not specified, then all founder animals will be assigned to batch 0.
#' @param batch2 A vector with names of the batches for the additional founder animals (Nm2 + Nf2). The number of additional founder animals that are simulated should be the same or a multiple of the length of the vector. If batch2 is not specified, the parameter Nbatch2 can be used to divide the additional founder animals over batches. If both Nbatch2 and batch2 are not specified, then the additional founder animals are divided over batches in the same way as the other founder animals (Nm+Nf).
#' @param Ntraits Number of traits to be simulated. Does not need to be specified if Ntraits is in the list 'BPdata'.
#' @param TraitsIndex Vector of traits that are in the index. By default, all traits are in the selection index.
#' @param Rgen Matrix of all genetic correlations between all Ntraits. Only needs to be specified if there is no matrix of genetic correlations named Rgen in the list called 'BPdata'.
#' @param Rcom Matrix of all common environmental correlations between all Ntraits. Only needs to be specified if there is no matrix of common environmental correlations named Rcom in the list called 'BPdata'.
#' @param Rres Matrix of all residual correlations between all Ntraits. Only needs to be specified if there is no matrix of residual correlations named Rres in the list called 'BPdata'.
#' @param mean Vector of means of all traits. Only needs to be specified if there is no vector of means named mean in the list 'BPdata'.
#' @param a_var Vector of genetic variances of all traits. Only needs to be specified if there is no vector of genetic variances named a_var in the list 'BPdata'.
#' @param c_var Vector of common environmental variances of all traits. Only needs to be specified if there is no vector of common environmental variances named c_var in the list 'BPdata'. If there is no common environmental effect, provide a vector of zero's.
#' @param e_var Vector of residual variances of all traits. Only needs to be specified if there is no vector of residual variances named e_var in the list 'BPdata'.
#' @param est_EBV TRUE or FALSe for estimating breeding values for all founder animals. The default is FALSE.
#' @param EBV If est_EBV is TRUE, a vector of methods for simulating EBVs for each trait in the selection index need to be specified. There are four options: 0 for giving each animal a breeding values of 0 for the specific trait; "mean_pop" for simulating EBVs for all founder populations that are equal to the mean of the population, the mean of the population is provided in the vector 'mean'; "pheno" for simulating EBVs that are equal to the phenotype of the animal for the specific trait and "EBV" for simulating EBVs as a value correlated to the true breeding value of the animal. This correlation equals the accuracy. Accuracies need to be specified in the parameter 'accuracy'.
#' @param accuracy If for one of the traits the EBVs need to be calculated with the method "EBV" specified in the parameter 'EBV', then accuracies need to be provided. A value need to be added for each trait, however this value can be zero for traits from which the EBVs are not simulated using method "EBV".
#' @param indexweights If traits need to be combined in an index, then desired gain indices need to be specified for each trait in the index.
#' @return A data frame is returned, which should be called 'ped' in order to be able to use it in the other functions. The 'ped' data frame consist of all the simulated base animals, including their sex, generation, batch number, inbreeding level, phenotypes and genetic, common environmental and residual effects of each trait.
#' @export
#' @examples
#' ped <- founderpopgroup(Nm=60, Nm2=0,
#'                      Nf=60, Nf2 = 0,
#'                      batch = c(0,1,2),
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
#'               
#' ped <- founderpopgroup(Nm=60,
#'                      Nf=60,
#'                      Nm2=120,
#'                      Nf2=120,
#'                      Nbatch = 4,
#'                      batch2 = c(-3,-2,-1,0),
#'                      Ntraits=2,
#'                      TraitsIndex = c(1,2),
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
#'                      e_var= c(250,12000),
#'                      est_EBV = TRUE,
#'                      EBV= c("pheno", "EBV"),
#'                      accuracy= c(NA,0.78),
#'                      indexweight= c(1,5))
#'

founderpopgroup <- function(Nm=BPdata$Nm, Nf=BPdata$Nf, Nm2=BPdata$Nm2, Nf2=BPdata$Nf2, Nbatch=NA, Nbatch2=NA, batch=NA, batch2=NA, Ntraits=BPdata$Ntraits, TraitsIndex=c(1:Ntraits), Rgen=BPdata$Rgen,Rres=BPdata$Rres,Rcom=BPdata$Rcom, mean=BPdata$mean, a_var=BPdata$a_var, c_var = BPdata$c_var, e_var=BPdata$e_var, est_EBV=FALSE, EBV, accuracy, indexweights=c(1)){
# Same as founderpopgroup, only the column FSfam is replaced with Contribution

  if(is.null(Nm2)){Nm2 <- 0}
  if(is.null(Nf2)){Nf2 <- 0}

   ped <- data.frame(id = 1:(Nm+Nf+Nm2+Nf2),
                    sire = 0,
                    dam = 0,
                    sex = c(rep(1, times = Nm), rep(2, times = Nf),rep(1, times = Nm2), rep(2, times = Nf2)),
                    generation = 0,
                    batch = 0,
                    inbreeding = 0,
                    Contribution = 0,
                    selected = c(rep(1,times=(Nm+Nf)), rep(0,times=Nm2+Nf2)),
                    preselected = c(rep(0,times=(Nm+Nf)), rep(1,times=Nm2+Nf2)),
                    selcand = c(rep(0,times=(Nm+Nf)), rep(1,times=Nm2+Nf2)))


  if(is.na(Nbatch) & is.na(Nbatch2) & sum(is.na(batch))==0 & sum(is.na(batch2))>0){
    ped$batch <- rep(x = batch, times = (Nm+Nf+Nm2+Nf2)/length(batch))}

  if(is.na(Nbatch) & is.na(Nbatch2) & sum(is.na(batch))==0 & sum(is.na(batch2))==0){
    ped$batch <- c(rep(x = batch, times = (Nm+Nf)/length(batch)),rep(x = batch2, times = (Nm2+Nf2)/length(batch2)))}

  if(is.na(Nbatch) & !is.na(Nbatch2) & sum(is.na(batch))==0 & sum(is.na(batch2))>0){
    ped$batch <- c(rep(x = batch, times = (Nm+Nf)/length(batch)),rep(x = 1:Nbatch2, times = (Nm2+Nf2)/Nbatch2))}

  if(!is.na(Nbatch) & is.na(Nbatch2) & sum(is.na(batch))>0 & sum(is.na(batch2))>0){
    ped$batch <- rep(x = 1:Nbatch, times = (Nm+Nf+Nm2+Nf2)/Nbatch)}

  if(!is.na(Nbatch) & is.na(Nbatch2) & sum(is.na(batch))>0 & sum(is.na(batch2))==0){
    ped$batch <- c(rep(x = 1:Nbatch, times = (Nm+Nf)/Nbatch),rep(x = batch2, times = (Nm2+Nf2)/length(batch2)))}

  if(!is.na(Nbatch) & !is.na(Nbatch2) & sum(is.na(batch))>0 & sum(is.na(batch2))>0){
    ped$batch <- c(rep(x = 1:Nbatch, times = (Nm+Nf)/Nbatch),rep(x = 1:Nbatch2, times = (Nm2+Nf2)/Nbatch2))}


    ped <- cbind(ped, MASS::mvrnorm(n = nrow(ped), mu = rep(0, nrow(Rgen)), Sigma = Rgen))
    ped <- cbind(ped, MASS::mvrnorm(n = nrow(ped), mu = rep(0, nrow(Rcom)), Sigma = Rcom))
    ped <- cbind(ped, MASS::mvrnorm(n = nrow(ped), mu = rep(0, nrow(Rres)), Sigma = Rres))
    ped <- cbind(ped,matrix(c(rep(0, times=nrow(ped)*(Ntraits))),nrow=nrow(ped)))
    ped <- cbind(ped,matrix(c(rep(0, times=nrow(ped)*(length(TraitsIndex)))),nrow=nrow(ped)))

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

  names(ped) <- colnames

  ped$Index <- 0

  for(trait in 1:Ntraits){
    #print(trait)
    ped[,c(11+trait)] <- ped[,c(11+trait)]*sqrt(a_var[trait]) + mean[trait]
    ped[,c(11+Ntraits+trait)] <- ped[,c(11+Ntraits+trait)]*sqrt(c_var[trait])
    ped[,c(11+2*Ntraits+trait)] <- ped[,c(11+2*Ntraits+trait)]*sqrt(e_var[trait])
    ped[,c(11+3*Ntraits+trait)] <- ped[,c(11+trait)] + ped[,c(11+Ntraits+trait)] + ped[,c(11+2*Ntraits+trait)]
  }

  if(est_EBV==T){

    for(trait in 1:length(TraitsIndex)){
      if(EBV[trait]==0){ ped[,c(11+4*Ntraits+trait)]<-0 }
      if(EBV[trait]=="mean_pop"){ ped[,c(11+4*Ntraits+trait)]<- mean[TraitsIndex[trait]] }
      if(EBV[trait]=="pheno"){ ped[,c(11+4*Ntraits+trait)]<- ped[,c(11+3*Ntraits+TraitsIndex[trait])] }
      if(EBV[trait]=="EBV"){
        pe <- ped[,c(1,2,3)]
        pe$pe <- stats::rnorm(nrow(pe)) * sqrt(accuracy[trait]) * sqrt(1 - (accuracy[trait])) * sqrt(a_var[TraitsIndex[trait]])
        ped[,c(11+4*Ntraits+trait)] <- ((ped[,c(11+trait)]-mean[TraitsIndex[trait]])*accuracy[trait]) + pe$pe + mean[TraitsIndex[trait]]
      }
    }

    weights <- indexweights[1]/sqrt(a_var[TraitsIndex[1]])

    if(length(TraitsIndex) > 1){
      for(i in 2:length(TraitsIndex)){
        weights <- c(weights, indexweights[i]/sqrt(a_var[TraitsIndex[i]]))
      }
    }

    weights <- weights/weights[1]

    for(i in 1:length(TraitsIndex)){
      ped$Index <- ped$Index + weights[i]* ped[,c(11+4*Ntraits+i)]
    }
  }

  return(ped)
}

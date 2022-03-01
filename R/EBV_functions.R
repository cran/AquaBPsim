#' simulating correlated variable
#'
#' @param x Vector of values to which the other vector needs to be correlated
#' @param cor correlation
#' @return Returns a vector with values that are correlated to vector x.
#' @export
#' @examples
#'
#'cor_var(c(2,4,2,2,6,7,5,6,6,7,9,4,5), 0.5)



cor_var <- function(x, cor){
  x1 <- stats::rnorm(length(x))
  res <- stats::residuals(stats::lm(x1 ~ x))
  val <- cor * stats::sd(res) * x + res * stats::sd(x) * sqrt(1 - cor^2)
  val <- (val / stats::sd(val)) * sqrt(stats::var(x) * cor^2)
  return(val)}


#' Simulating estimated breeding values
#'
#' This function can be used to simulate estimated breeding values (EBV) for the selection candidates.
#'
#' @details
#' Breeding values will be simulated only for fish with ped$selcand == 1. EBVs are not estimated but simulated as a value correlated to the true breeding value. The correlation is equal to the accuracy, which can be calculated or provided by the user.
#'
#' There are three options for simulating the EBVs, namely:
#'
#'     - "pheno": EBV equal to phenotype.
#'
#'     - "PEBV": Pedigree estimated breeding values: either an accuracy needs to be provided or the accuracy will be calculated using information about the number of full sibs and half sibs present in the ped file. When the accuracy is calculated, it is assumed that the selection candidate has an own performance of the trait and that no common environmental effects are present (for example in group mating design). Only sibs that are also selection candidates are used in the calculation of the accuracy. When an accuracy is provided, prediction errors are correlated to the common environmental effects.
#'
#'     - "GEBV": Genomically estimated breeding values. Either an accuracy needs to be provided or the genome length, effective population size and size of the training population need to be provided in order to calculate the accuracy using the formula of Deatwyler et al. (2010).
#'
#'     - "sib_pheno": EBVs are calculated from the phenotypes of the full sibs and half sibs. If this option is choosen, then the parameter presel_sibs needs to be used to specify which sibs are going to be used to calculate the breeding values. For example, if presel_sibs = 2 is specified, then only the sibs that are preselected for 'environment 2' are used.
#'
#'
#' If EBVs are simulated for more than one trait, then the EBVs are combined in an index. Combining the traits in an index using the desired gains can be done in two ways:
#'
#'     - Method 1 (default): For each trait, the desired gains are divided by the genetic standard deviations of the trait and then multipplied with the EBV. The index is a summation of these values of each trait.
#'
#'     - Method 2: The EBVs are first standardized to a mean of 100 and a standard deviation of 10 and then multiplied with the desired gains. The index is a summation of these values of each trait.
#'
#'
#' @param gen The generation of the fish for which breeding values need to be simulated.
#' @param batch The batch of the fish for which breeding values need to be simulated, default is 0.
#' @param TraitsIndex A vector of traits that are in the index. For these traits, breeding values will be simulated.
#' @param Ntraits Total number of simulated traits. Does not need to be specified if a list called 'BPdata' including the variable Ntraits is available.
#' @param EBV A vector of methods for simulating EBVs for each trait in the selection index needs to be specified. Options: "pheno" for simulating an EBV that is equal to the phenotype of the trait; "GEBV" for simulating GEBVs as a value correlated to the true breeding value of the animal, this correlation equals the accuracy; "PEBV" for simulating EBVs as a value correlated to the true breeding value of the animal (correlation = accuracy) and the prediction errors are correlated to the common environmental effects; "sib_pheno" for simualting EBVs using only information of the sibs of the selection candidate. Accuracies can be specified in the parameter 'accuracy', or in the case of GEBVs, accuracies can be calculated with the formula of Deatwyler et al. (2010) when the effective population size, genome length and size of the trainingspopulation for each trait are specified.
#' @param accuracy Vector of accuracies with the length equal to the number of traits in the Index. For traits for which the accuracy is not needed or traits for which the accuracy needs to be calculated, NA needs to be specified in the vector.
#' @param presel_sibs Vector to indicate which preselected fish can be used to calculate the EBV when EBV = "sibs_pheno". Needs to be specified for each trait in index (if not applicable, then NA can be specified).
#' @param GenomLength Genome length in Morgan.
#' @param Ne Effective population size.
#' @param SizeTraining Vector with the size of the training population for each trait in the index.
#' @param h2 Vector with heritabilities. Does not need to be specified when the heritabilities are already provided in the list 'BPdata'.
#' @param c2 Vector with common environmental effects. Does not need to be specified when the common environmental effects are already provided in the list 'BPdata'.
#' @param indexweights If traits need to be combined in an index, then desired gain indices need to be specified for each trait in the index.
#' @param a_var Vector of genetic variances of all traits, does not need to be specified if already provided in the list 'BPdata'.
#' @param method_indexweights Either 1 or 2, default is 1, see Details.
#' @return This function will change the data frame called 'ped'. Simulated breeding values will be added to the EBV columns and values will be added to the column 'Index'. 
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
#'                       e_var= c(250,12000))
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
#'}


breeding_values <- function(gen, batch=0, TraitsIndex=c(1:Ntraits), Ntraits=BPdata$Ntraits, EBV, accuracy=NA, GenomLength, Ne, SizeTraining, h2 = BPdata$h2, presel_sibs = c(rep(1, times=length(TraitsIndex))), c2=BPdata$c2, indexweights=c(1), a_var=BPdata$a_var, method_indexweights=1){
  #estimation EBVs, loop for each trait in index seperately
  for(trait in 1:length(TraitsIndex)){
    # If EBV equals teh phenotype/own performance fish
    if(EBV[trait]=="pheno"){ ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+4*Ntraits+trait)]<- ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+3*Ntraits+TraitsIndex[trait])] }

    #Simulating genomic estimated breeding values, either with accuracies that are already provided or calculating accuracy (genome length, ne, size trainingspopulation need to be given)
    if(EBV[trait]=="GEBV"){
      if(is.na(accuracy[trait])){
        Me <- (2*Ne*GenomLength)/(log(4*Ne*GenomLength, base=10))
        acc <- sqrt(SizeTraining[trait]*h2[TraitsIndex[trait]]/(SizeTraining[trait]* h2[TraitsIndex[trait]] + Me))

        pe <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(1,2,3)]
        pe$pe <- stats::rnorm(nrow(pe)) * sqrt(acc) * sqrt(1 - (acc)) * sqrt(a_var[TraitsIndex[trait]])
        mean <- mean(ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+TraitsIndex[trait])])
        ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+4*Ntraits+trait)] <- ((ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+TraitsIndex[trait])]-mean)*acc) + pe$pe + mean

      }else{
        pe <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(1,2,3)]
        pe$pe <- stats::rnorm(nrow(pe)) * sqrt(accuracy[trait]) * sqrt(1 - (accuracy[trait])) * sqrt(a_var[TraitsIndex[trait]])
        mean <- mean(ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+TraitsIndex[trait])])
        ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+4*Ntraits+trait)] <- ((ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+TraitsIndex[trait])]-mean)*accuracy[trait]) + pe$pe + mean
      }
    }

    # Simulating pedigree breeding values, either with accuracies that are already provided or the accuracy is calculated according to the number of full sibs
    # and half sibs present (assumed that fish also has own performance)
    if(EBV[trait]=="PEBV"){
      if(is.na(accuracy[trait])){
        preselected <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,]
        families <- preselected[,2:3]
        families <- unique(families)

        for(family in 1:nrow(families)){
          FS_family <- preselected[preselected$sire==families$sire[family] & preselected$dam==families$dam[family],]
          Mfs <- nrow(FS_family) - 1
          dam <- preselected[preselected$sire!=families$sire[family] & preselected$dam==families$dam[family],]
          sire <- preselected[preselected$sire==families$sire[family] & preselected$dam!=families$dam[family],]
          Mhs_dam <- nrow(dam)
          Mhs_sire <- nrow(sire)

          # Met fullsibs en half sibs van sire en dam
          if(Mfs>0 & Mhs_dam >0 & Mhs_sire>0){
            P <- matrix(c(1,           (0.5*h2[TraitsIndex[trait]]),                (0.25*h2[TraitsIndex[trait]]),                       (0.25*h2[TraitsIndex[trait]]),
                          (0.5*h2[TraitsIndex[trait]]) , ((1+(Mfs-1)*0.5*h2[TraitsIndex[trait]])/Mfs),(0.25*h2[TraitsIndex[trait]]),                       (0.25*h2[TraitsIndex[trait]]),
                          (0.25*h2[TraitsIndex[trait]]), (0.25*h2[TraitsIndex[trait]]),               ((1+(Mhs_dam-1)*0.5*h2[TraitsIndex[trait]])/Mhs_dam), 0,
                          (0.25*h2[TraitsIndex[trait]]), (0.25*h2[TraitsIndex[trait]]),               0,                                 ((1+(Mhs_sire-1)*0.5*h2[TraitsIndex[trait]])/Mhs_sire)),
                        nrow = 4, byrow = TRUE)

            G <- matrix(c(h2[TraitsIndex[trait]], 0.5*h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]]), nrow=4, byrow=TRUE)

            P_inv <- matlib::inv(P)
            b <- P_inv %*% G
            b_1 <- matrix(b,nrow=1)
            acc <- sqrt((b_1 %*% G)/h2[TraitsIndex[trait]])}

          # Met alleen half sibs van sire en dam
          if(Mfs==0 & Mhs_dam >0 & Mhs_sire>0){
            P <- matrix(c(1,           (0.25*h2[TraitsIndex[trait]]),                       (0.25*h2[TraitsIndex[trait]]),
                          (0.25*h2[TraitsIndex[trait]]), ((1+(Mhs_dam-1)*0.5*h2[TraitsIndex[trait]])/Mhs_dam), 0,
                          (0.25*h2[TraitsIndex[trait]]), 0,                                 ((1+(Mhs_sire-1)*0.5*h2[TraitsIndex[trait]])/Mhs_sire)),
                        nrow = 3, byrow = TRUE)

            G <- matrix(c(h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]]), nrow=3, byrow=TRUE)

            P_inv <- matlib::inv(P)
            b <- P_inv %*% G
            b_1 <- matrix(b,nrow=1)
            acc <- sqrt((b_1 %*% G)/h2[TraitsIndex[trait]])}

          # Met  half sibs van sire
          if(Mfs==0 & Mhs_dam ==0 & Mhs_sire>0){
            P <- matrix(c(1,            (0.25*h2[TraitsIndex[trait]]),
                          (0.25*h2[TraitsIndex[trait]]),  ((1+(Mhs_sire-1)*0.5*h2[TraitsIndex[trait]])/Mhs_sire)),
                        nrow = 2, byrow = TRUE)

            G <- matrix(c(h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]]), nrow=2, byrow=TRUE)

            P_inv <- matlib::inv(P)
            b <- P_inv %*% G
            b_1 <- matrix(b,nrow=1)
            acc <- sqrt((b_1 %*% G)/h2[TraitsIndex[trait]])}

          # Met alleen own record
          if(Mfs==0 & Mhs_dam ==0 & Mhs_sire==0){
            acc<- matrix(c(sqrt(h2[TraitsIndex[trait]])),nrow=1)}

          # Met fullsibs en half sibs van sire
          if(Mfs>0 & Mhs_dam ==0 & Mhs_sire>0){
            P <- matrix(c(1,           (0.5*h2[TraitsIndex[trait]]),                (0.25*h2[TraitsIndex[trait]]),
                          (0.5*h2[TraitsIndex[trait]]) , ((1+(Mfs-1)*0.5*h2[TraitsIndex[trait]])/Mfs),(0.25*h2[TraitsIndex[trait]]),
                          (0.25*h2[TraitsIndex[trait]]), (0.25*h2[TraitsIndex[trait]]),               ((1+(Mhs_sire-1)*0.5*h2[TraitsIndex[trait]])/Mhs_sire)),
                        nrow = 3, byrow = TRUE)

            G <- matrix(c(h2[TraitsIndex[trait]], 0.5*h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]]), nrow=3, byrow=TRUE)

            P_inv <- matlib::inv(P)
            b <- P_inv %*% G
            b_1 <- matrix(b,nrow=1)
            acc <- sqrt((b_1 %*% G)/h2[TraitsIndex[trait]])}

          # Met fullsibs
          if(Mfs>0 & Mhs_dam ==0 & Mhs_sire==0){
            P <- matrix(c(1,           (0.5*h2[TraitsIndex[trait]]),
                          (0.5*h2[TraitsIndex[trait]]) , ((1+(Mfs-1)*0.5*h2[TraitsIndex[trait]])/Mfs)),
                        nrow = 2, byrow = TRUE)

            G <- matrix(c(h2[TraitsIndex[trait]], 0.5*h2[TraitsIndex[trait]]), nrow=2, byrow=TRUE)

            P_inv <- matlib::inv(P)
            b <- P_inv %*% G
            b_1 <- matrix(b,nrow=1)
            acc <- sqrt((b_1 %*% G)/h2[TraitsIndex[trait]])}

          # Met fullsibs en half sibs dam
          if(Mfs>0 & Mhs_dam >0 & Mhs_sire==0){
            P <- matrix(c(1,           (0.5*h2[TraitsIndex[trait]]),                (0.25*h2[TraitsIndex[trait]]),
                          (0.5*h2[TraitsIndex[trait]]) , ((1+(Mfs-1)*0.5*h2[TraitsIndex[trait]])/Mfs),(0.25*h2[TraitsIndex[trait]]),
                          (0.25*h2[TraitsIndex[trait]]), (0.25*h2[TraitsIndex[trait]]),               ((1+(Mhs_dam-1)*0.5*h2[TraitsIndex[trait]])/Mhs_dam)),
                        nrow = 3, byrow = TRUE)

            G <- matrix(c(h2[TraitsIndex[trait]], 0.5*h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]]), nrow=3, byrow=TRUE)

            P_inv <- matlib::inv(P)
            b <- P_inv %*% G
            b_1 <- matrix(b,nrow=1)
            acc <- sqrt((b_1 %*% G)/h2[TraitsIndex[trait]])}

          # Met half sibs  dam
          if(Mfs==0 & Mhs_dam >0 & Mhs_sire==0){
            P <- matrix(c(1,                        (0.25*h2[TraitsIndex[trait]]),
                          (0.25*h2[TraitsIndex[trait]]),               ((1+(Mhs_dam-1)*0.5*h2[TraitsIndex[trait]])/Mhs_dam)),
                        nrow = 2, byrow = TRUE)

            G <- matrix(c(h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]]), nrow=2, byrow=TRUE)

            P_inv <- matlib::inv(P)
            b <- P_inv %*% G
            b_1 <- matrix(b,nrow=1)
            acc <- sqrt((b_1 %*% G)/h2[TraitsIndex[trait]])}

          preselected$acc[preselected$id %in% FS_family$id] <- acc[1,1]
        }

        preselected$pe <- stats::rnorm(nrow(preselected)) * sqrt(preselected$acc) * sqrt(1 - (preselected$acc)) * sqrt(a_var[TraitsIndex[trait]])
        mean <- mean(ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+TraitsIndex[trait])])
        ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+4*Ntraits+trait)] <- ((ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+TraitsIndex[trait])]-mean)*preselected$acc) + preselected$pe + mean

      }else{
        cor_com <- (c2[TraitsIndex[trait]]*h2[TraitsIndex[trait]]*(1-(accuracy[trait])^2))
        pe <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,]
        pe$pe <- cor_var(pe[,c(11+Ntraits+TraitsIndex[trait])], cor_com)
        pe$pe <- (pe$pe/stats::sd(pe$pe)) * accuracy[trait] * sqrt(1 - accuracy[trait]^2) * sqrt(a_var[TraitsIndex[trait]])
        mean <- mean(ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+TraitsIndex[trait])])
        ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+4*Ntraits+trait)] <- ((ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(11+TraitsIndex[trait])]-mean)*accuracy[trait]) + pe$pe + mean
      }
    }

    if(EBV[trait]=="sib_pheno"){
      sibtrait <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$preselected==presel_sibs[trait],c(1,2,3,11+3*Ntraits+TraitsIndex[trait])]
      selcand <- ped[ped$generation %in% gen & ped$batch %in% batch & ped$selcand==1,c(1,2,3)]
      selcand$EBV <- 0

      for(fish in 1:nrow(selcand)){
        full_sib <- sibtrait[sibtrait$sire==selcand$sire[fish] & sibtrait$dam==selcand$dam[fish],]
        half_dam <- sibtrait[sibtrait$sire==selcand$sire[fish] & sibtrait$dam !=selcand$dam[fish],]
        half_sire <- sibtrait[sibtrait$sire !=selcand$sire[fish] & sibtrait$dam==selcand$dam[fish],]

        mean_FS <- mean(full_sib[,c(4)])
        mean_HS_dam <- mean(half_dam[,c(4)])
        mean_HS_sire <- mean(half_sire[,c(4)])

        Mfs <- nrow(full_sib)
        Mhs_dam <- nrow(half_dam)
        Mhs_sire <- nrow(half_sire)

        # Met fullsibs en half sibs van sire en dam
        if(Mfs>0 & Mhs_dam >0 & Mhs_sire>0){
          P <- matrix(c(((1+(Mfs-1)*0.5*h2[TraitsIndex[trait]])/Mfs),(0.25*h2[TraitsIndex[trait]]),                       (0.25*h2[TraitsIndex[trait]]),
                        (0.25*h2[TraitsIndex[trait]]),               ((1+(Mhs_dam-1)*0.5*h2[TraitsIndex[trait]])/Mhs_dam), 0,
                        (0.25*h2[TraitsIndex[trait]]),               0,                                 ((1+(Mhs_sire-1)*0.5*h2[TraitsIndex[trait]])/Mhs_sire)),
                      nrow = 3, byrow = TRUE)

          G <- matrix(c(0.5*h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]]), nrow=3, byrow=TRUE)

          P_inv <- matlib::inv(P)
          b <- P_inv %*% G

          b_sum <- b[1] + b[2] + b[3]
          b1 <- b[1]/b_sum
          b2 <- b[2]/b_sum
          b3 <- b[3]/b_sum

          selcand$EBV[fish] <- b1*mean_FS + b2*mean_HS_dam + b3* mean_HS_sire
        }

        # Met alleen half sibs van sire en dam
        if(Mfs==0 & Mhs_dam >0 & Mhs_sire>0){
          P <- matrix(c(((1+(Mhs_dam-1)*0.5*h2[TraitsIndex[trait]])/Mhs_dam), 0,
                        0,                                 ((1+(Mhs_sire-1)*0.5*h2[TraitsIndex[trait]])/Mhs_sire)),
                      nrow = 2, byrow = TRUE)

          G <- matrix(c(0.25*h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]]), nrow=2, byrow=TRUE)

          P_inv <- matlib::inv(P)
          b <- P_inv %*% G

          b_sum <- b[1] + b[2]
          b1 <- b[1]/b_sum
          b2 <- b[2]/b_sum

          selcand$EBV[fish] <- b1*mean_HS_dam + b2* mean_HS_sire
        }

        # Met  half sibs van sire
        if(Mfs==0 & Mhs_dam ==0 & Mhs_sire>0){
          selcand$EBV[fish] <- mean_HS_sire
        }

        # Met fullsibs en half sibs van sire
        if(Mfs>0 & Mhs_dam ==0 & Mhs_sire>0){
          P <- matrix(c(((1+(Mfs-1)*0.5*h2[TraitsIndex[trait]])/Mfs),(0.25*h2[TraitsIndex[trait]]),
                        (0.25*h2[TraitsIndex[trait]]),               ((1+(Mhs_sire-1)*0.5*h2[TraitsIndex[trait]])/Mhs_sire)),
                      nrow = 2, byrow = TRUE)

          G <- matrix(c(0.5*h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]]), nrow=2, byrow=TRUE)

          P_inv <- matlib::inv(P)
          b <- P_inv %*% G

          b_sum <- b[1] + b[2]
          b1 <- b[1]/b_sum
          b2 <- b[2]/b_sum

          selcand$EBV[fish] <- b1*mean_FS + b2* mean_HS_sire
        }

        # Met fullsibs
        if(Mfs>0 & Mhs_dam ==0 & Mhs_sire==0){
          selcand$EBV[fish] <- mean_FS
        }

        # Met fullsibs en half sibs dam
        if(Mfs>0 & Mhs_dam >0 & Mhs_sire==0){
          P <- matrix(c(((1+(Mfs-1)*0.5*h2[TraitsIndex[trait]])/Mfs),(0.25*h2[TraitsIndex[trait]]),
                        (0.25*h2[TraitsIndex[trait]]),               ((1+(Mhs_dam-1)*0.5*h2[TraitsIndex[trait]])/Mhs_dam)),
                      nrow = 2, byrow = TRUE)

          G <- matrix(c(0.5*h2[TraitsIndex[trait]], 0.25*h2[TraitsIndex[trait]]), nrow=2, byrow=TRUE)

          P_inv <- matlib::inv(P)
          b <- P_inv %*% G

          b_sum <- b[1] + b[2]
          b1 <- b[1]/b_sum
          b2 <- b[2]/b_sum

          selcand$EBV[fish] <- b1*mean_FS + b2* mean_HS_dam
        }

        # Met half sibs  dam
        if(Mfs==0 & Mhs_dam >0 & Mhs_sire==0){
          selcand$EBV[fish] <- mean_HS_dam
        }

        if(Mfs==0 & Mhs_dam ==0 & Mhs_sire==0){
          selcand$EBV[fish] <- (ped[ped$id==selcand$sire, c(11+4*Ntraits+trait)] + ped[ped$id==selcand$dam,c(11+4*Ntraits+trait)])/2
        }

        ped[ped$id==selcand$id[fish],c(11+4*Ntraits+trait)] <- selcand$EBV[fish]
      }
    }
  }


  #indexweight
  if(method_indexweights==1){
    weights <- indexweights[1]/sqrt(a_var[TraitsIndex[1]])

    if(length(TraitsIndex) > 1){
      for(i in 2:length(TraitsIndex)){
        weights <- c(weights, indexweights[i]/sqrt(a_var[TraitsIndex[i]]))
      }
    }

    weights <- weights/weights[1]

    for(i in 1:length(TraitsIndex)){
      ped$Index[ped$generation==gen & ped$batch==batch & ped$selcand==1] <- ped$Index[ped$generation==gen & ped$batch==batch & ped$selcand==1] + weights[i]* ped[ped$generation==gen & ped$batch==batch & ped$selcand==1,c(11+4*Ntraits+i)]
    }
  }


  if(method_indexweights==2){
    for(EBV in 1:length(TraitsIndex)){
      ped[ped$generation==gen& ped$batch==batch& ped$selcand==1, c(11+4*Ntraits+EBV)] <- ped[ped$generation==gen& ped$batch==batch& ped$selcand==1, c(11+4*Ntraits+EBV)] /sqrt(stats::var(ped[ped$generation==gen& ped$batch==batch& ped$selcand==1, c(11+4*Ntraits+EBV)] )) *10
      ped[ped$generation==gen& ped$batch==batch& ped$selcand==1, c(11+4*Ntraits+EBV)] <- ped[ped$generation==gen& ped$batch==batch& ped$selcand==1, c(11+4*Ntraits+EBV)] - mean(ped[ped$generation==gen& ped$batch==batch& ped$selcand==1, c(11+4*Ntraits+EBV)]) + 100
    }

    ped$Index[ped$generation==gen & ped$batch==batch & ped$selcand==1] <- ped[ped$generation==gen& ped$batch==batch& ped$selcand==1, c(11+4*Ntraits+1)]*indexweights[1]

    if(length(TraitsIndex)>1){
      for(i in 2:length(TraitsIndex)){
        ped$Index[ped$generation==gen & ped$batch==batch & ped$selcand==1] <-ped$Index[ped$generation==gen & ped$batch==batch & ped$selcand==1] + ped[ped$generation==gen& ped$batch==batch& ped$selcand==1, c(11+4*Ntraits+i)]*indexweights[i]
      }
    }
  }


  return(ped)
}


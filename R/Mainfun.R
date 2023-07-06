#' MF_single
#'
#' \code{MF_single}:\cr
#' Multi-functionality measures in a single ecosystem. This is a function in the framework of Hill-Chao numbers
#'
#' @param func_data the ecosystem function data can be input as a vector of functions (for a single assemblage), matrix/data.frame (assemblages by functions).\cr
#' The data input must be normalized between 0 and 1 already and must contain only the ecosystem function columns.\cr
#' For \code{species_data} is not \code{NULL}, the rownames of func_data should be names of plotID.
#' @param species_data the species abundance data must include three columns: 'plotID', 'species' and 'abundance'. Default is \code{NULL}.
#' @param q a numerical vector specifying the diversity orders. Default is 0, 1 and 2.

#' @import devtools
#' @import ggplot2
#' @import dplyr
#' @import iNEXT.3D
#' @import tidyverse
#' @import tidyr
#' @importFrom stats cor

#' @return  a data.frame with columns 'Type', 'Order.q' and 'qMF'. For \code{species_data} is not \code{NULL}, the data.frame will contain 'Species.diversity'.
#' @export


MF_single <- function(func_data, species_data = NULL, q = c(0,1,2)){

  if ((length(func_data) == 1) && (inherits(func_data, c("numeric", "integer"))))
    stop("Error: Your data does not have enough information.")

  if (FALSE %in% (0 <= func_data && func_data <= 1))
    stop("Error: Your data does not be normalized, please transform the data between 0 and 1 first.")
  # if (is.vector(data) & (FALSE %in% (0 <= data & data <= 1)))
  #   stop("Error: Your data does not have enough information to be normalized automatically, please transform the data between 0 and 1 first.")

  if(!is.null(species_data)){
    if(!all(names(species_data) %in% c("plotID","species","abundance")))
      stop("Error: The species_data should include the following three columns: 'plotID', 'species' and 'abundance'.")
    else if(!all(rownames(func_data) %in% unique(species_data$plotID)))
      stop("Error: The column 'plotID' of species_data should include all of the assemblages (i.e, row names) in func_data.")
    else if(sum(is.na(species_data))!=0)
      stop("Error: There exists NA values in species_data.")
  }

  if (is.vector(func_data)) func_data <- matrix(func_data,nrow=1)


  # if (nrow(fun_data)<2 & (FALSE %in% (0 <= fun_data & fun_data <= 1)))
  #   stop("Error: Your data does not have enough information to be normalized automatically, please transform the data between 0 and 1 first.")
  # if ((FALSE %in% (0 <= fun_data & fun_data <= 1)) && is_normalized)
  #   stop("Error: Your data does not be normalized between 0 and 1, please change the argument is_normalized to be FALSE.")


  func_data[is.na(func_data)] <- 0
  id_data <- data.frame(plotID = rownames(func_data))

  qMF_output <- sapply(q,function(i) apply(func_data,1,function(x) qMF(v=x,q=i))) %>% as.data.frame()
  # if(nrow(data)==1) qMF_output <- t(qMF_output)
  names(qMF_output) <- paste0("qMF_Uncorrelated_",q)

  if(ncol(func_data)>1){
    tau <- seq(0,1,0.01)
    transform_D <- sqrt(1-abs(cor(func_data)))
    qMF_tau_output <- sapply(q,function(i) {
      out <- lapply(tau,function(t){
        apply(func_data,1,function(x) qMF(v=x,q=i,tau=t,di=transform_D))
      }) %>% do.call(cbind,.)
      AUC <- apply(out,1,function(x){
        AUC_L <- sum(x[seq_along(x[-1])]*diff(tau))
        AUC_R <- sum(x[-1]*diff(tau))
        (AUC_L+AUC_R)/2
      })
      AUC
    }) %>% as.data.frame()
    # if(nrow(norm_data)==1) qMF_tau_output <- t(qMF_tau_output)
    names(qMF_tau_output) <- paste0("qMF_Correlated_",q)

    output <- cbind(id_data,qMF_output,qMF_tau_output) %>%
      pivot_longer(cols = starts_with("qMF"),
                   names_to = c("Type", "Order.q"),
                   names_pattern = "qMF_?(.*)_(.*)",
                   values_to = "qMF") %>%
      mutate(Order.q = factor(paste0("q = ",Order.q)))
  }
  else{
    output <- cbind(id_data,qMF_output) %>%
      pivot_longer(cols = starts_with("qMF"),
                   names_to = c("Type", "Order.q"),
                   names_pattern = "qMF_?(.*)_(.*)",
                   values_to = "qMF") %>%
      mutate(Order.q = factor(paste0("q = ",Order.q)))
  }



  if(!is.null(species_data)){
    species_div <- lapply(id_data$plotID,function(p){
      data_s <- species_data %>% filter(plotID == p) %>% group_by(species) %>%
        summarise(abundance = sum(abundance))
      div = iNEXT.3D::AO3D(data_s$abundance,q=q,nboot = 0,method = "Observed")
      return(data.frame("Species diversity"=rep(div$qD,ifelse(ncol(func_data)>1,2,1))))
    }) %>% do.call(rbind,.)

    output <- cbind(output,species_div)
  }

  return(output)
}


#' \code{MF_multiple}:\cr
#' Multi-functionality measures in multiple ecosystems. This is a function in the framework of Hill-Chao numbers
#'
#' @param func_data the ecosystem function data can be input as a vector of functions (for a single assemblage), matrix/data.frame (assemblages by functions).\cr
#' The data input must be normalized between 0 and 1 already.\cr
#' For \code{by_group = NULL} the func_data must contain only the ecosystem function columns. Otherwise, you must add the \code{by_group} column in data.
#' For \code{species_data} is not \code{NULL}, the rownames of func_data should be names of plotID.
#' @param species_data the species abundance data must include three columns: 'plotID', 'species' and 'abundance'. Default is \code{NULL}.
#' @param q a numerical vector specifying the diversity orders. Default is 0, 1 and 2.
#' @param by_group name of the column to be paired by group. Default is \code{NULL}.

#' @return a data.frame with columns 'Order.q' and 'MF_gamma', 'MF_alpha' and 'MF_beta' for uncorrelated output, 'MF_gamma_cor', 'MF_alpha_cor' and 'MF_beta_cor' for correlated output.\cr
#' For \code{species_data} is not \code{NULL}, the data.frame will contain 'Species_gamma', 'Species_alpha' and 'Species_beta'.
#' @export

MF_multiple <- function(func_data, species_data = NULL, q = c(0,1,2), by_group = NULL){

  if ((length(func_data) == 1) && (inherits(func_data, c("numeric", "integer"))))
    stop("Error: Your data does not have enough information.")

  if (FALSE %in% (0 <= (func_data %>% dplyr::select(-by_group)) && (func_data %>% dplyr::select(-by_group)) <= 1))
    stop("Error: Your data does not be normalized, please transform the data between 0 and 1 first.")

  if(!is.null(species_data)){
    if(!all(names(species_data) %in% c("plotID","species","abundance")))
      stop("Error: The species_data should include the following three columns: 'plotID', 'species' and 'abundance'.")
    else if(!all(rownames(func_data) %in% unique(species_data$plotID)))
      stop("Error: The column 'plotID' of species_data should include all of the assemblages (i.e, row names) in func_data.")
    else if(sum(is.na(species_data))!=0)
      stop("Error: There exists NA values in species_data.")
  }

  if (is.vector(func_data)) func_data <- matrix(func_data,nrow=1)


  if(is.null(by_group)){
    func_data <- ungroup(func_data)

    spe_data <- func_data
    two_idx <- combn(1:nrow(spe_data),2)

    fun_data <- func_data
    if(ncol(fun_data)>1){
      tau <- seq(0,1,0.01)
      fun_data[is.na(fun_data)] <- 0
      transform_D <- sqrt(1-abs(cor(fun_data)))
    }

    output <- lapply(1:ncol(two_idx),function(i){
      d_x <- cbind(x1=t(spe_data[two_idx[1,i],]),
                   x2=t(spe_data[two_idx[2,i],]))
      rF <- sapply(q,function(y) qMF_diversity(d_x,y,diversity = "gamma"))
      aF <- sapply(q,function(y) qMF_diversity(d_x,y,diversity = "alpha"))
      bF <- rF/aF

      two_plot <- rownames(spe_data)[c(two_idx[1,i],two_idx[2,i])]
      result <- data.frame("plotID"=rep(paste(two_plot[1],two_plot[2],sep = " v.s. "),length(q)),
                           "Order.q"=paste0("q = ",q),
                           "Uncorrelated_Gamma"=rF,
                           "Uncorrelated_Alpha"=aF,
                           "Uncorrelated_Beta"=bF)


      if(ncol(fun_data)>1){
        rF_tau <- lapply(tau,function(t){
          out <- data.frame("Order.q"=paste0("q = ",q),
                            "value"=sapply(q,function(y) qMF_diversity(d_x,y,t,transform_D,diversity = "gamma")),
                            "Tau"=t)
          out
        }) %>% do.call(rbind,.)

        aF_tau <- parallel::mclapply(tau,function(t){
          out <- data.frame("Order.q"=paste0("q = ",q),
                            "value"=sapply(q,function(y) qMF_diversity(d_x,y,t,transform_D,diversity = "alpha")),
                            "Tau"=t)
        }) %>% do.call(rbind,.)
        bF_tau <- rF_tau$value/aF_tau$value

        result_tau <- data.frame("Order.q"=rF_tau$Order.q,
                                 "Tau"=rF_tau$Tau,
                                 "MF_g"=rF_tau$value,
                                 "MF_a"=aF_tau$value,
                                 "MF_b"=bF_tau) %>% group_by(Order.q) %>%
          summarise(L_g = sum(MF_g[seq_along(MF_g[-1])]*diff(Tau)),
                    R_g = sum(MF_g[-1]*diff(Tau)),
                    L_a = sum(MF_a[seq_along(MF_a[-1])]*diff(Tau)),
                    R_a = sum(MF_a[-1]*diff(Tau)),
                    L_b = sum(MF_b[seq_along(MF_b[-1])]*diff(Tau)),
                    R_b = sum(MF_b[-1]*diff(Tau))) %>% ungroup %>%
          mutate(Correlated_Gamma = (L_g+R_g)/2,
                 Correlated_Alpha = (L_a+R_a)/2,
                 Correlated_Beta = (L_b+R_b)/2) %>%
          dplyr::select(c(Correlated_Gamma:Correlated_Beta))

        result <- cbind(result,result_tau)
      }

      if(!is.null(species_data)){

        Obs_div <- function(data,q){
          data_gamma = rowSums(data)
          data_gamma = data_gamma[data_gamma > 0]
          data_alpha = as.matrix(data) %>% as.vector
          data_alpha = data_alpha[data_alpha > 0]

          gamma = iNEXT.3D::AO3D(data_gamma,q=q,nboot = 0,method = "Observed")
          g = gamma$qD
          alpha = iNEXT.3D::AO3D(data_alpha,q=q,nboot = 0,method = "Observed")
          a = alpha$qD/2
          b = g/a
          return(list(gamma=g,alpha=a,beta=b))
        }

        spec_data1 <- (species_data %>% filter(plotID==two_plot[1])) %>%
          group_by(species) %>% summarise(x1 = sum(abundance)) %>% select(c(species,x1))
        spec_data2 <- (species_data %>% filter(plotID==two_plot[2])) %>%
          group_by(species) %>% summarise(x2 = sum(abundance)) %>% select(c(species,x2))
        spec_data <- dplyr::full_join(spec_data1,spec_data2,by="species") %>% select(-species)
        spec_data[is.na(spec_data)] <- 0
        div_out <- Obs_div(spec_data,q=q)
        s_r <- div_out$gamma
        s_a <- div_out$alpha
        s_b <- div_out$beta

        result <- cbind(result,data.frame("Species_Gamma"=s_r,
                                          "Species_Alpha"=s_a,
                                          "Species_Beta"=s_b))
      }

      result
    }) %>% do.call(rbind,.)

  }
  else if(length(by_group)!=1) stop("Error: The number of the group variable should not be more than 1.")
  else if(!(by_group %in% colnames(func_data))) stop("Error: The group variable is not included in the given data.")
  else{
    func_data <- ungroup(func_data)
    output <- lapply(unique(unlist(func_data %>% dplyr::select(by_group)) %>% as.character()),function(group){
      spe_data <- func_data %>% filter(across(by_group) == group) %>% select(-by_group)
      two_idx <- combn(1:nrow(spe_data),2)

      fun_data <- func_data %>% select(-by_group)
      if(ncol(fun_data)>1){
        tau <- seq(0,1,0.01)
        fun_data[is.na(fun_data)] <- 0
        transform_D <- sqrt(1-abs(cor(fun_data)))
      }

      out <- lapply(1:ncol(two_idx),function(i){
        d_x <- cbind(x1=t(spe_data[two_idx[1,i],]),
                     x2=t(spe_data[two_idx[2,i],]))
        rF <- sapply(q,function(y) qMF_diversity(d_x,y,diversity = "gamma"))
        aF <- sapply(q,function(y) qMF_diversity(d_x,y,diversity = "alpha"))
        bF <- rF/aF

        two_plot <- rownames(spe_data)[c(two_idx[1,i],two_idx[2,i])]
        result <- data.frame("plotID"=rep(paste(two_plot[1],two_plot[2],sep = " v.s. "),length(q)),
                             "group"=rep(group,length(q)),
                             "Order.q"=paste0("q = ",q),
                             "Uncorrelated_Gamma"=rF,
                             "Uncorrelated_Alpha"=aF,
                             "Uncorrelated_Beta"=bF)
        colnames(result)[2] <- by_group


        if(ncol(fun_data)>1){
          rF_tau <- lapply(tau,function(t){
            out <- data.frame("Order.q"=paste0("q = ",q),
                              "value"=sapply(q,function(y) qMF_diversity(d_x,y,t,transform_D,diversity = "gamma")),
                              "Tau"=t)
            out
          }) %>% do.call(rbind,.)

          aF_tau <- parallel::mclapply(tau,function(t){
            out <- data.frame("Order.q"=paste0("q = ",q),
                              "value"=sapply(q,function(y) qMF_diversity(d_x,y,t,transform_D,diversity = "alpha")),
                              "Tau"=t)
          }) %>% do.call(rbind,.)
          bF_tau <- rF_tau$value/aF_tau$value

          result_tau <- data.frame("Order.q"=rF_tau$Order.q,
                                   "Tau"=rF_tau$Tau,
                                   "MF_g"=rF_tau$value,
                                   "MF_a"=aF_tau$value,
                                   "MF_b"=bF_tau) %>% group_by(Order.q) %>%
            summarise(L_g = sum(MF_g[seq_along(MF_g[-1])]*diff(Tau)),
                      R_g = sum(MF_g[-1]*diff(Tau)),
                      L_a = sum(MF_a[seq_along(MF_a[-1])]*diff(Tau)),
                      R_a = sum(MF_a[-1]*diff(Tau)),
                      L_b = sum(MF_b[seq_along(MF_b[-1])]*diff(Tau)),
                      R_b = sum(MF_b[-1]*diff(Tau))) %>% ungroup %>%
            mutate(Correlated_Gamma = (L_g+R_g)/2,
                   Correlated_Alpha = (L_a+R_a)/2,
                   Correlated_Beta = (L_b+R_b)/2) %>%
            dplyr::select(c(Correlated_Gamma:Correlated_Beta))

          result <- cbind(result,result_tau)
        }

        if(!is.null(species_data)){

          Obs_div <- function(data,q){
            data_gamma = rowSums(data)
            data_gamma = data_gamma[data_gamma > 0]
            data_alpha = as.matrix(data) %>% as.vector
            data_alpha = data_alpha[data_alpha > 0]

            gamma = iNEXT.3D::AO3D(data_gamma,q=q,nboot = 0,method = "Observed")
            g = gamma$qD
            alpha = iNEXT.3D::AO3D(data_alpha,q=q,nboot = 0,method = "Observed")
            a = alpha$qD/2
            b = g/a
            return(list(gamma=g,alpha=a,beta=b))
          }

          spec_data1 <- (species_data %>% filter(plotID==two_plot[1])) %>%
            group_by(species) %>% summarise(x1 = sum(abundance)) %>% select(c(species,x1))
          spec_data2 <- (species_data %>% filter(plotID==two_plot[2])) %>%
            group_by(species) %>% summarise(x2 = sum(abundance)) %>% select(c(species,x2))
          spec_data <- dplyr::full_join(spec_data1,spec_data2,by="species") %>% select(-species)
          spec_data[is.na(spec_data)] <- 0
          div_out <- Obs_div(spec_data,q=q)
          s_r <- div_out$gamma
          s_a <- div_out$alpha
          s_b <- div_out$beta

          result <- cbind(result,data.frame("Species_Gamma"=s_r,
                                            "Species_Alpha"=s_a,
                                            "Species_Beta"=s_b))
        }
        # other_data <- data[,-fun_cols] %>% filter(Country == country) %>% dplyr::select(!Country)
        # other <- paste(other_data[two_idx[1,i],],other_data[two_idx[2,i],],sep = "-") %>%
        #   rep(each=length(q)) %>% matrix(nrow = length(q)) %>% as.data.frame()
        # names(other) <- names(other_data)

        result
      }) %>% do.call(rbind,.)
    }) %>% do.call(rbind,.)
  }

  output <- output %>% tidyr::pivot_longer(cols = starts_with(c("Uncorrelated","Correlated")),
                                           names_to = c("Type","Scale"),
                                           names_pattern = "(.*)_(.*)",
                                           values_to = "qMF") %>%
    dplyr::mutate(Species.diversity=ifelse(Scale=="Gamma",Species_Gamma,
                                           ifelse(Scale=="Alpha",Species_Alpha,
                                                  Species_Beta))) %>%
    dplyr::select(-c(Species_Gamma:Species_Beta))

  return(output)
}



#' For functions normalized
#'
#' \code{func_normalized}:\cr
#' To get the normalized functions between 0 and 1.
#'
#' @param data data can be input as a matrix/data.frame (multiple plots by functions).
#' @param fun_cols the order number of the columns which be used as the ecosystem function.
#' @param negative name of columns to be normalized negatively.
#' @param by_group name of the column to normalized by group. Default is \code{NULL}.

#' @return a data.frame with normalized functions.
#' @export


func_normalized <- function(data, fun_cols = 1:ncol(data), negative = NULL, by_group = NULL){
  data <- dplyr::ungroup(data)
  fun_data <- data[,fun_cols]
  if(is.null(negative)) neg_cols <- 0
  else if (!all(negative %in% colnames(fun_data))) stop("Error: Negative columns must be included in function columns.")
  else neg_cols <- which(colnames(fun_data) %in% negative)

  trans_fun <- function(x, positive=TRUE){
    mi <- min(x,na.rm = T)
    ma <- max(x,na.rm = T)
    if(positive) y <- (x-mi)/(ma-mi)
    else y <- (ma-x)/(ma-mi)
    y[x==0] <- NA
    return(y)
  }

  if(is.null(by_group)){
    trans_out <- lapply(1:ncol(fun_data),function(i){
      if(i %in% neg_cols) trans_fun(fun_data[,i],F)
      else trans_fun(fun_data[,i],T)
    }) %>% do.call(cbind,.)
    colnames(trans_out) <- colnames(fun_data)
  }
  else{
    if(length(by_group)!=1) stop("Error: The number of the group variable for normalization should not be more than 1.")
    else if(!(by_group %in% colnames(data))) stop("Error: The group variable is not included in the given data.")

    data <- data %>% dplyr::arrange(across(by_group))
    gr_data <- cbind(data %>% dplyr::select(all_of(by_group)),data[,fun_cols])
    trans_out <- lapply(unique(gr_data[,1]),function(g){
      sub_gr <- gr_data[gr_data[,1]==g,] %>% .[,-1] %>% as.data.frame()
      out <- lapply(1:ncol(sub_gr),function(i){
        if(i %in% neg_cols) trans_fun(sub_gr[,i],F)
        else trans_fun(sub_gr[,i],T)
      }) %>% do.call(cbind,.)
      out
    }) %>% do.call(rbind,.)
  }
  data[,fun_cols] <- trans_out
  return(data)
}


# Calculate multi-functionality measures in a single ecosystem.
#
# \code{qMF} Calculate uncorrelated or correlated MF measures using special case of Hill-Chao numbers.
#
# @param w a vector of the weighted.
# @param v a vector of the ecosystem with several functions.
# @param q a numerical vector specifying the diversity orders.
# @param tau a single value used in the correlated MF measures.
# @param di a distance matrix used in the correlated MF measures.
# @return a value of correlated/uncorrelated MF measure in a single ecosystem.

qMF <- function(w=1,v,q,tau=0.5,di=NULL){
  if(is.null(di) | tau == 0){
    ai_tau <- v[v>0]
    V <- w*ai_tau
  }
  else{
    d_tau <- ifelse(di<tau,di,tau)
    ai <- (1-d_tau/tau)%*%v
    ai_tau <- ai[ai>0]
    v <- v[ai>0]
    V <- w*v*v/ai_tau
  }
  if(q==1) exp(-sum(V*(ai_tau/sum(V*ai_tau))*log(ai_tau/sum(V*ai_tau))))
  else sum(V*(ai_tau/sum(V*ai_tau))^q)^(1/(1-q))
}


# Calculate multi-functionality measures in multiple ecosystems.
#
# \code{qMF_diversity} Calculate uncorrelated or correlated MF measures using special case of Hill-Chao numbers.
#
# @param v a matrix contains two columns of the ecosystem with several functions.
# @param q a numerical vector specifying the diversity orders.
# @param tau a single value used in the correlated MF measures.
# @param di a distance matrix used in the correlated MF measures.
# @param diversity decide to compute whether 'gamma' diversity or 'alpha' diversity.
# @return a value of correlated/uncorrelated MF measure in multiple ecosystems.

qMF_diversity <- function(v,q,tau=0.5,di=NULL,diversity="gamma"){
  F_alpha <- function(w=1,v,q,R=c(0.5,0.5),N=2){
    v <- as.matrix(v)
    v1 <- v
    v[is.na(v)] <- 0
    V <- w*(v%*%R)
    v_plus <- c()
    for (i in 1:N){
      v_plus <- cbind(v_plus,v[,i]*R[i])
    }
    if (q==1) {
      Vv <- ifelse(v_plus==0,0,(v_plus/sum(V*(v%*%R)))*log(v_plus/sum(V*(v%*%R))))
      1/N*exp(-sum(t(V)%*%Vv))
    }
    else {
      Vv <- ifelse(is.na(v1),0,(v_plus/sum(V*(v%*%R)))^q)
      1/N*(sum(V*rowSums(Vv)))^(1/(1-q))
    }
  }
  F_gamma <- function(w=1,v,q,R=c(0.5,0.5),N=2){
    v <- as.matrix(v)
    v[is.na(v)] <- 0
    V <- w*(v%*%R)
    if(q==1) {
      Vv <- ifelse(v%*%R==0,0,((v%*%R)/sum(V*(v%*%R)))*log((v%*%R)/sum(V*(v%*%R))))
      exp(-sum(V*Vv))
    }
    else (sum(V*((v%*%R)/sum(V*(v%*%R)))^q))^(1/(1-q))
  }
  F_alpha_tau <- function(w=1,v,q,tau,di,R=c(0.5,0.5),N=2){
    v <- as.matrix(v)
    di <- di[rowSums(is.na(v))!=N,rowSums(is.na(v))!=N]
    v <- v[rowSums(is.na(v))!=N,]
    v1 <- v
    v[is.na(v)] <- 0
    d_tau <- ifelse(di<tau,di,tau)
    ai <- (1-d_tau/tau)%*%v
    f_bar <- rowSums(v)/N
    ai_bar <- (1-d_tau/tau)%*%f_bar
    v_tau <- ifelse(f_bar==0,0,f_bar*f_bar/ai_bar)
    V <- w*v_tau

    if (q==1) {
      Vv <- ifelse(ai==0,0,(ai/sum(V%*%ai))*log(ai/sum(V%*%ai)))
      1/N*exp(-sum(t(V)%*%Vv))
    }
    else {
      Vv <- ifelse(is.na(v1) & ai==0,0,(ai/sum(V%*%ai))^q)
      1/N*(sum(V*rowSums(Vv)))^(1/(1-q))
    }
  }
  F_gamma_tau <- function(w=1,v,q,tau,di,R=c(0.5,0.5),N=2){
    v <- as.matrix(v)
    v[is.na(v)] <- 0
    d_tau <- ifelse(di<tau,di,tau)
    ai <- (1-d_tau/tau)%*%v
    f_bar <- rowSums(v)/N
    ai_bar <- (1-d_tau/tau)%*%f_bar
    v_tau <- ifelse(f_bar==0,0,f_bar*f_bar/ai_bar)
    V <- w*v_tau

    if(q==1) {
      Vv <- ifelse(ai_bar==0,0,((ai_bar)/sum(V*ai_bar))*log((ai_bar)/sum(V*ai_bar)))
      exp(-sum(V*Vv))
    }
    else (sum(V*((ai_bar)/sum(V*ai_bar))^q))^(1/(1-q))
  }

  if(is.null(di) | tau == 0){
    if(diversity=="gamma") out = F_gamma(v=v,q=q)
    else out = F_alpha(v=v,q=q)
  }
  else{
    if(diversity=="gamma") out = F_gamma_tau(v=v,q=q,tau = tau,di=di)
    else out = F_alpha_tau(v=v,q=q,tau = tau,di=di)
  }
  return(out)
}





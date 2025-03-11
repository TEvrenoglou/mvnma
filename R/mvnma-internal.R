
catch <- function(argname, matchcall, data, encl)
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)


multi_arm <- function(dat){
  
  u <- unique(dat$studlab)    
  
  r <- list()
  
  t <- list()
  
  E <- list()
  
  for(i in 1:length(u)){
    
    r[[i]] <- dat %>% 
      filter(studlab == u[i])
    
    if(nrow(r[[i]])>2){
      
      t[[i]] <- as.data.frame(table(r[[i]]$treat2))
      
      E[[i]] <- t[[i]]$Var1[which(t[[i]]$Freq==max(t[[i]]$Freq))]
      
      r[[i]] <- r[[i]] %>% 
        filter(treat2 == E[[i]])
    }
    
  }
  #### df harmonized in terms of treat2 
  
  df <- list_rbind(r)
  
  return(df)
}


create_T <- function(data1,max_arms){
  
  u <- unique(data1$studlab)  
  
  mat <- matrix(NA,nrow = length(u),ncol=max_arms)
  
  r <- list()
  
  treats <- list()
  
  for (i in 1:length(u)){
    
    r[[i]] <- data1 %>% 
      filter(studlab==u[i])
    
    
    treats[[i]] <- c(unique(r[[i]]$label2),unique(r[[i]]$label1))
    
    mat[i,][1:length(treats[[i]])] <- treats[[i]]
    
    # replace NA's with 0
    #T[i,] <- ifelse(is.na(T[i,]),0,T[i,])
    
    
  }  
  
  return(mat)
  
}


'%!in%' <- function(x,y)!('%in%'(x,y))

###### Helpers for new mvdata

create_data <- function(p,...){
  
  if(any(class(p)=="list")){
    
    dat1 <- list()
    
    dat2 <- list()
    
    studies <- list()
    
    n_outcomes <- length(p)
    
    for(i in 1:n_outcomes){
      
      p[[i]] <- p[[i]] %>% 
        dplyr::select(studlab,TE,seTE,treat1,treat2) 
      
      p[[i]] <- multi_arm(p[[i]])
      
      p[[i]]$outcome= i
      
      p[[i]] <- add_arms(p[[i]])
    }
    
    ## combine all p's
    
    comb_p <- list_rbind(p)
    
    for(i in 1:n_outcomes){
      
      dat1[[i]] <- comb_p %>% 
        filter(outcome==i)
      
      dat2[[i]] <- comb_p %>% 
        filter(outcome!=i) 
      
      dat2[[i]]$TE = dat2[[i]]$seTE = NA
      
      dat2[[i]]$new_outcome = i
      
      studies[[i]] <- unique(which(dat2[[i]]$studlab %!in% dat1[[i]]$studlab))
      
      if(length(studies[[i]])>0){
        
        dat2[[i]] <- dat2[[i]][studies[[i]],]
        
        dat2[[i]]$outcome <- dat2[[i]]$new_outcome
        
        dat2[[i]]$new_outcome <- NULL
        
        row.names(dat2[[i]]) <- NULL
        
      }else{
        
        dat2 <- list()
      }
      
    }
    
    if(length(dat2)!=0){
      
      dat2 <- list_rbind(dat2)
      
      comb_p <- rbind.data.frame(comb_p,dat2)
      
      comb_p <- comb_p %>% 
        arrange(outcome)
      
      row.names(comb_p) = NULL
      
    }else{
      
      comb_p <- comb_p %>% 
        arrange(outcome)
      
    }
    
    data_final<- comb_p %>%
      group_by(studlab) %>%
      arrange(desc(treat1),.by_group = T) %>% 
      arrange(n_arms)
    
    #data_final$n_arms <- NULL
    
    data_final <- add_labels(data_final)
    
  }else{
    
    stop("Argument 'p' must be a list of pairwise objects.",
         call. = FALSE)
    
  }
  
  return(data_final)
  
}

add_arms <- function(dat,...){
  
  u <- unique(dat$studlab)
  
  dat1 <- list()
  
  arms <- list()
  
  for(i in 1:length(u)){
    
    dat1[[i]] <- dat %>% 
      filter(studlab == u[i])
    
    arms[[i]] <- length(unique(c(dat1[[i]]$treat1,dat1[[i]]$treat2)))
    
    
    dat1[[i]]$n_arms <- arms[[i]]
    
  }
  
  dat1 <- list_rbind(dat1)  
  
  return(dat1)
}

add_labels <- function(data){
  
  all_treats <- unique(c(data$treat1,data$treat2))  
  
  levels_treats <- as.data.frame(levels(as.factor(all_treats)))
  
  names(levels_treats) <- c("treat")
  
  levels_treats$level <- 1:nrow(levels_treats)  
  
  data$label1 <- NA
  
  data$label2 <- NA
  
  for(i in 1:nrow(data)){
    
    for(j in 1:nrow(levels_treats)){
      
      if(data$treat1[i]==levels_treats$treat[j]){
        
        data$label1[i] = levels_treats$level[j]
        
      }
      
      if(data$treat2[i]==levels_treats$treat[j]){
        
        data$label2[i] = levels_treats$level[j]
        
      } 
      
    }
    
  }  
  
  return(data)
  
}

make_jags_data <- function(dat){
  
  ## number of studies  
  
  Ns <- length(unique(dat$studlab))    
  
  ## labtreat 
  
  labtreat1 <- cbind.data.frame(dat$treat1,dat$label1)
  
  labtreat2 <- cbind.data.frame(dat$treat2,dat$label2)
  
  names(labtreat1) <- names(labtreat2) <- c("treat","label")
  
  labtreat <- rbind.data.frame(labtreat1,labtreat2)
  
  labtreat <- labtreat %>% 
    distinct() %>% 
    arrange(label)
  
  #NT 
  
  NT <- length(labtreat$treat)
  
  ## number of outcomes
  
  n_outcomes <- length(unique(dat$outcome))
  
  ## arms per study
  arm_data <- dat[!duplicated(dat$studlab), ]
  
  max_arms <- max(arm_data$n_arms)
  
  ## number of 2 arms studies
  
  two_arm <- length(which(arm_data$n_arms==2))
  
  treat_data <- create_T(dat,max_arms=max_arms)
  
  
  ## extract vector with treatment effects
  
  y <- dat$TE
  
  y <- ifelse(is.na(y),0,y)
  
  dat_out <- list()
  
  var <- list()
  
  treat_out <- list()
  
  names_vec <- c()
  
  for(i in 1:n_outcomes){
    
    dat_out[[i]] <- dat[dat$outcome==i,]
    
    var[[i]] <- dat_out[[i]]$seTE^2
    
    var[[i]] <- ifelse(is.na(var[[i]]),10000,var[[i]])
    
    names_vec[i] <- paste("var",i,sep = "")
    
    dat_out[[i]] <- dat_out[[i]][complete.cases(dat_out[[i]]$TE),]
    
    treat_out[[i]] <- unique(c(dat_out[[i]]$treat1,dat_out[[i]]$treat2))
    
  }
  
  var_f <- list.cbind(var) 
  
  var_f <- as.data.frame(var_f)
  
  names(var_f) <- names_vec
  
  
  dat_f <- list("y"=y,
                "var" = var_f,
                "T" = treat_data,
                "Ns" = Ns,
                "N2h" = two_arm,
                #"na" = arms,
                "labtreat" = labtreat,
                "treat_out" = treat_out,
                "NT" = NT
  )
  
  return(dat_f)
  
}

is.list.pairwise <- function(p,...){
  
  all_class <- sapply(p,class)   
  
  check <- c()
  
  for(i in 1:ncol(all_class)){
    
    check[i] <- isTRUE("pairwise" %in% all_class[,i])  
    
  }
  
  pair <- ifelse(sum(check)>=2,"pairwise",NA)
  
  return(pair)
  
}



gather_results <- function(res,n.out,labtreat,ref,reference.group,
                           data,sims_bugs_out,outlab,
                           
                           ...){
  
  reference <- c()
  
  ds <- c()
  
  basic_comp <- list()
  
  dat_treat <- list()
  
  psi <- list()
  
  rho <- list()
  
  d <- list()
  
  for(i in 1:n.out){
    
    reference[i] <- paste("ref",i,sep = "")
    
    ds[i] <- paste("d",i,sep = "")
    
    basic_comp[[i]] <- res %>% 
      filter(grepl(reference[i], ind))
    
    row.names(basic_comp[[i]]) <- labtreat[-ref]
    
    dat_treat[[i]] <- data$treat_out[[i]][-which(data$treat_out[[i]]==reference.group)]
    
    basic_comp[[i]] <- basic_comp[[i]][which(row.names(basic_comp[[i]]) %in% dat_treat[[i]]),]
    
    basic_comp[[i]] <- basic_comp[[i]] %>% 
      dplyr::select(mean,sd,`2.5%`,`97.5%`,Rhat)
    
    names(basic_comp[[i]]) <- c("TE","sd","lb.ci","ub.ci","Rhat")
    
    ## psi's (similar to tau)
    
    psi[[i]] <- res %>%
      filter(grepl(c("psi"), ind)) %>% 
      dplyr::select(mean,sd,`2.5%`,`97.5%`,Rhat)
    
    names(psi[[i]]) <- c("psi","sd","lb.ci","ub.ci","Rhat")
    
    ## rho
    rho[[i]] <- res %>% 
      filter(grepl("rho",ind))
    
    d[[i]] <- sims_bugs_out[,c(grepl(ds[i], names(sims_bugs_out)))]
    
    names(d[[i]]) <- labtreat
    
    d[[i]] <- d[[i]] %>% 
      dplyr::select(all_of(data$treat_out[[i]]))
    
    
  }
  
  ## prepare output
  
  outcome_correlation <- rho[[1]]
  
  psi <- psi[[1]]
  
  row.names(psi) <- paste(outlab)
  
  outcome_correlation <- outcome_correlation %>% 
    dplyr::select(mean,sd,`2.5%`,`97.5%`,Rhat)
  
  names(outcome_correlation) <- c("corr_coef","sd","lb.ci","ub.ci","Rhat")
  
  ## create row.names for outcome_correlation
  r1 <- seq(1:length(outlab))
  
  r1 <- t(combn(r1,2))
  
  r.names <- matrix(ncol=1,nrow=nrow(outcome_correlation))
  
  for(i in 1:nrow(r.names)){
    
    for(j in 1:ncol(r1)){
      
      r.names[i,] <- paste(outlab[r1[i,1]],outlab[r1[i,2]],sep = "/")    
    }
    
  }
  
  row.names(outcome_correlation) <- r.names
  
  if(n.out==2){
    
    out1 <- list("basic_estimates"=basic_comp[[1]],
                 "heterogeneity"=psi[1,],
                 "samples"=d[[1]]
    )
    
    out2 <- list("basic_estimates"=basic_comp[[2]],
                 "heterogeneity"=psi[2,],
                 "samples"=d[[2]]
    )
    
    res.final <- list(out1,out2,outcome_correlation)
    
    names(res.final) <- c(paste(outlab[1]),paste(outlab[2]),"outcome_correlation")
    
    
    
  }
  
  if(n.out==3){
    
    out1 <- list("basic_estimates"=basic_comp[[1]],
                 "heterogeneity"=psi[1,],
                 "samples"=d[[1]]
    )
    
    out2 <- list("basic_estimates"=basic_comp[[2]],
                 "heterogeneity"=psi[2,],
                 "samples"=d[[2]]
    )
    
    out3 <- list("basic_estimates"=basic_comp[[3]],
                 "heterogeneity"=psi[3,],
                 "samples"=d[[3]]
    )
    
    res.final <- list(out1,out2,out3,outcome_correlation)
    
    names(res.final) <- c(paste(outlab[1]),paste(outlab[2]),paste(outlab[3]),"outcome_correlation")
    
  }
  
  if(n.out==4){
    
    out1 <- list("basic_estimates"=basic_comp[[1]],
                 "heterogeneity"=psi[1,],
                 "samples"=d[[1]]
    )
    
    out2 <- list("basic_estimates"=basic_comp[[2]],
                 "heterogeneity"=psi[2,],
                 "samples"=d[[2]]
    )
    
    out3 <- list("basic_estimates"=basic_comp[[3]],
                 "heterogeneity"=psi[3,],
                 "samples"=d[[3]]
    )
    
    out4 <- list("basic_estimates"=basic_comp[[4]],
                 "heterogeneity"=psi[4,],
                 "samples"=d[[4]]
    )
    
    res.final <- list(out1,out2,out3,out4,outcome_correlation)
    
    names(res.final) <- c(paste(outlab[1]),paste(outlab[2]),
                          paste(outlab[3]),paste(outlab[4]),
                          "outcome_correlation")
    
  }
  
  if(n.out==5){
    
    out1 <- list("basic_estimates"=basic_comp[[1]],
                 "heterogeneity"=psi[1,],
                 "samples"=d[[1]]
    )
    
    out2 <- list("basic_estimates"=basic_comp[[2]],
                 "heterogeneity"=psi[2,],
                 "samples"=d[[2]]
    )
    
    out3 <- list("basic_estimates"=basic_comp[[3]],
                 "heterogeneity"=psi[3,],
                 "samples"=d[[3]]
    )
    
    out4 <- list("basic_estimates"=basic_comp[[4]],
                 "heterogeneity"=psi[4,],
                 "samples"=d[[4]]
    )
    
    out5 <- list("basic_estimates"=basic_comp[[5]],
                 "heterogeneity"=psi[5,],
                 "samples"=d[[5]]
    )
    
    res.final <- list(out1,out2,out3,out4,out5,outcome_correlation)
    
    names(res.final) <- c(paste(outlab[1]),paste(outlab[2]),
                          paste(outlab[3]),paste(outlab[4]),
                          paste(outlab[5]),
                          "outcome_correlation")
    
  }
  
  return(res.final)
}


get_all_estimates <- function(x){
  
  x <- x[names(x) != "outcome_correlation"]  
  
  all <- list()
  
  all_TE <- list()
  
  all_sd <- list()
  
  for(i in 1:length(x)){
    
    all[[i]] <- get_diff(x[[i]]$samples)
    
    all_TE[[i]] <- all[[i]]$TE.random
    
    all_sd[[i]] <- all[[i]]$sd.random
    
  }
  
  res <- list("TE.random"=all_TE,
              "sd.random" = all_sd
  )
  
  return(res)
}


get_diff <- function(samples){
  
  diff <- list()
  
  E <- names(samples)
  
  ests <- list()
  
  means <- list()
  
  sds <- list()
  
  for(i in 1:ncol(samples)){
    
    diff[[i]] <- samples %>% 
      dplyr::select(any_of(E[i]))
    
    ests[[i]] <- samples[,E[i]]-samples
    
    ests[[i]][E[[i]]] <- 0
    
    
  }
  
  for(k in 1:length(ests)){
    
    means[[k]] <- colMeans(ests[[k]])
    
    sds[[k]] <- colSds(as.matrix(ests[[k]]))
    
  }
  
  TE <- matrix(NA,ncol = ncol(samples),nrow=ncol(samples))
  
  sd <- matrix(NA,ncol = ncol(samples),nrow=ncol(samples))
  
  row.names(TE) <- colnames(TE) <- names(TE)
  
  row.names(sd) <- colnames(sd) <- names(sd)
  
  for(j in 1:nrow(TE)){
    
    TE[j,] <- means[[j]]
    
    sd[j,] <- sds[[j]]
    
  }
  
  TE <- as.data.frame(TE)
  
  sd <- as.data.frame(sd)
  
  row.names(TE) <- row.names(sd) <- names(TE) <- names(sd) <- E
  
  
  res_f <- list("TE.random"=TE,
                "sd.random"=sd
  )
  
  
  return(res_f) 
  
}   


# functions:
#` Computes beta parameters from mean and variance
#`
beta.params <- function(m, v){
  a <- ((1 - m) * (m^2))/v - m
  a[a < 0] <- NA
  b <- a * (1 - m)/m
  ab <- list(alpha = a, beta = b)
  return(ab)
}

summarize.counts <- function(data.1){
  
  # First half
  data.1 %>%
    filter(hr < 7) %>%
    group_by(sector, year, date) %>%
    select(sector, year, date, survey_duration_hr, counts, trench_length_m, date_num) %>% #-> tmp.1
    summarize(across(.cols = counts, 
                     .fns = sum,  
                     na.rm = T, 
                     .names = "counts_1"),
              length_m = first(trench_length_m),
              date_num = first(date_num),
              duration_hr = first(survey_duration_hr)) -> sum.0.6.hrs.12hr
  
  # Second half
  data.1 %>%
    filter(hr > 6) %>%
    group_by(sector, year, date) %>%
    summarize(across(.cols = counts, 
                     .fns = sum, 
                     na.rm = T, 
                     .names = "counts_2")) -> sum.7.12.hrs.12hr
  
  # put them together and convert counts_2 = NA when only 6-hr survey
  sum.0.6.hrs.12hr %>% 
    left_join(sum.7.12.hrs.12hr, by = c("sector", "year", "date")) %>%
    mutate(counts_2a = ifelse(duration_hr == 6, NA, counts_2)) -> counts.1st.2nd
  
  return(counts.1st.2nd)  
}

save.fig.fcn <- function(fig, file.name, replace = F, dpi = 600, device = "png",
                         height = 4, 
                         width = 6, 
                         units = "in",
                         bg = "white"){
  if (isTRUE(replace) | (!isTRUE(replace) & !file.exists(file.name)))
    ggsave(fig, 
           filename = file.name, 
           dpi = dpi, device = device,
           height = height, width = width, units = units,
           bg = bg)
  
}

compute.LOOIC <- function(loglik.mat, data.vector, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  #loglik.vec <- as.vector(loglik)
  
  # each column corresponds to a data point and rows are MCMC samples
  #loglik.mat <- matrix(loglik.vec, nrow = n.per.chain * MCMC.params$n.chains)
  # take out the columns that correspond to missing data points
  loglik.mat <- loglik.mat[, !is.na(data.vector)]
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = 4)
  
  loo.out <- rstanarm::loo(loglik.mat, 
                           r_eff = Reff, 
                           cores = 4, k_threshold = 0.7)
  
  out.list <- list(Reff = Reff,
                   loo.out = loo.out)
  
  return(out.list)  
}


# a function to extract posterior samples from jags output
extract.samples <- function(varname, zm){
  dev <- unlist(lapply(zm, FUN = function(x) x[, varname]))
  return(dev)
}

cross.validataion.fcn <- function(model.base.name, model.ver, jags.data, jags.params, MCMC.params, original.data){
  
  counts.1st.2nd.jags.data <- original.data
  model.file <- paste0(model.base.name, model.ver, ".txt")
  
  #cross.validation.out.file.name <- paste0("RData/cross-validation_out_", model.ver, ".rds")
  cross.validation.out.file.name <- paste0("RData/cross-validation_out_log_counts_", model.ver, ".rds")
  
  if (!file.exists(cross.validation.out.file.name)){
    start.time <- now()
    # Prepare the output dataframe
    counts.X.val <- counts.1st.2nd.jags.data %>%
      mutate(estim.m = NA,
             q2.5.m = NA,
             q50.m = NA,
             q97.5.m = NA,
             estim.lambda = NA,
             q2.5.lambda = NA,
             q50.lambda = NA,
             q97.5.lambda = NA,
             estim.p = NA,
             q2.5.p = NA,
             q50.p = NA,
             q97.5.p = NA)
    
    jags.params <- c(jags.params, "m")
    
    # take one m out at a time and run jags on the modified dataset:
    k <- 1
    for (k in 1:nrow(counts.1st.2nd.jags.data)){
      print(paste( k, " of ", nrow(counts.1st.2nd.jags.data)))
      
      x.val.jags.data <- counts.1st.2nd.jags.data
      # Turn one to NA
      x.val.jags.data$counts_2[k] <- NA
      # Year index
      Y <- x.val.jags.data$year_num[k]
      # Sector index
      S <- x.val.jags.data$sector_num[k]
      # day index
      D <- x.val.jags.data$seq_date[k]
      
      # Poisson data
      # jags.data$m <- x.val.jags.data$counts_2
      
      # log Normal data
      jags.data$m <- log(x.val.jags.data$counts_2 + 1)
      
      jm <- jags(jags.data,
                 inits = NULL,
                 parameters.to.save= jags.params,
                 model.file = model.file,
                 n.chains = MCMC.params$n.chains,
                 n.burnin = MCMC.params$n.burnin,
                 n.thin = MCMC.params$n.thin,
                 n.iter = MCMC.params$n.samples,
                 DIC = T, parallel=T)
      
      # calculate m[k] using the posterior samples 
      if (str_split(model.ver, "-", simplify = TRUE)[1] == "v5"){
        p.samples <-  extract.samples("p",
                                      jm$samples) 
      } else {
        p.samples <-  extract.samples(paste0("p[", jags.data$year[k], "]"),
                                      jm$samples)
        
      }
      
      # find m value from observed n and posterior samples of p
      m.k <- log(counts.X.val$counts_1[k] + 1) + p.samples
      
      qtiles.m <- quantile(m.k, c(0.025, 0.5, 0.975))
      
      # Stats from calculated m values:
      counts.X.val$estim.m[k] <- mean(m.k)
      counts.X.val$q2.5.m[k] <- qtiles.m[1]
      counts.X.val$q50.m[k] <- qtiles.m[2]
      counts.X.val$q97.5.m[k] <- qtiles.m[3]
      
      # find appropriate lambda in mean, q2.5 and q97.5. 
      if (str_split(model.ver, "-", simplify = TRUE)[2] == "6"){
        counts.X.val$estim.lambda[k] <- jm$mean$lambda[Y,S]
        counts.X.val$q2.5.lambda[k] <- jm$q2.5$lambda[Y,S]
        counts.X.val$q50.lambda[k] <- jm$q50$lambda[Y,S]
        counts.X.val$q97.5.lambda[k] <- jm$q97.5$lambda[Y,S]
        
        if (str_split(model.ver, "-", simplify = TRUE)[1] == "v5"){
          counts.X.val$estim.p[k] <- jm$mean$p
          counts.X.val$q2.5.p[k] <- jm$q2.5$p
          counts.X.val$q50.p[k] <- jm$q50$p
          counts.X.val$q97.5.p[k] <- jm$q97.5$p
        } else {
          counts.X.val$estim.p[k] <- jm$mean$p[Y,S]
          counts.X.val$q2.5.p[k] <- jm$q2.5$p[Y,S]
          counts.X.val$q50.p[k] <- jm$q50$p[Y,S]
          counts.X.val$q97.5.p[k] <- jm$q97.5$p[Y,S]
        }
      } else if (str_split(model.ver, "-", simplify = TRUE)[2] == "1"){
        counts.X.val$estim.lambda[k] <- jm$mean$lambda[Y,S,D]
        counts.X.val$q2.5.lambda[k] <- jm$q2.5$lambda[Y,S,D]
        counts.X.val$q50.lambda[k] <- jm$q50$lambda[Y,S,D]
        counts.X.val$q97.5.lambda[k] <- jm$q97.5$lambda[Y,S,D]
        
        if (str_split(model.ver, "-", simplify = TRUE)[1] == "v5"){
          counts.X.val$estim.p[k] <- jm$mean$p
          counts.X.val$q2.5.p[k] <- jm$q2.5$p
          counts.X.val$q50.p[k] <- jm$q50$p
          counts.X.val$q97.5.p[k] <- jm$q97.5$p
        } else {
          counts.X.val$estim.p[k] <- jm$mean$p[Y]
          counts.X.val$q2.5.p[k] <- jm$q2.5$p[Y]
          counts.X.val$q50.p[k] <- jm$q50$p[Y]
          counts.X.val$q97.5.p[k] <- jm$q97.5$p[Y]
        }      
      } else {
        counts.X.val$estim.lambda[k] <- jm$mean$lambda[Y,S,D]
        counts.X.val$q2.5.lambda[k] <- jm$q2.5$lambda[Y,S,D]
        counts.X.val$q50.lambda[k] <- jm$q50$lambda[Y,S,D]
        counts.X.val$q97.5.lambda[k] <- jm$q97.5$lambda[Y,S,D]
        
        if (str_split(model.ver, "-", simplify = TRUE)[1] == "v5"){
          counts.X.val$estim.p[k] <- jm$mean$p
          counts.X.val$q2.5.p[k] <- jm$q2.5$p
          counts.X.val$q50.p[k] <- jm$q50$p
          counts.X.val$q97.5.p[k] <- jm$q97.5$p
        } else {
          counts.X.val$estim.p[k] <- jm$mean$p[Y,S,D]
          counts.X.val$q2.5.p[k] <- jm$q2.5$p[Y,S,D]
          counts.X.val$q50.p[k] <- jm$q50$p[Y,S,D]
          counts.X.val$q97.5.p[k] <- jm$q97.5$p[Y,S,D]
        }
        
      }
    }
    end.time <- now()
    cross.validation.out <- list(output = counts.X.val,
                                 input.jags.data = counts.1st.2nd.jags.data,
                                 run.date = Sys.Date(),
                                 system = Sys.info(),
                                 run.time = end.time - start.time)
    saveRDS(cross.validation.out,
            file = cross.validation.out.file.name)  
  } else {
    
    cross.validation.out <- readRDS(file = cross.validation.out.file.name)
  }
  
  return(cross.validation.out)
}


calculate.ratio <- function(data.1){
  
  # only include 12-hr data:
  data.1 %>% 
    filter(survey_duration_hr == 12) -> data.1.12hr
  
  # Create a sequential ID number for each survey
  # date_num_df <- data.frame(date_num = unique(data.1.12hr$date_num),
  #                           date_seq = 1:length(unique(data.1.12hr$date_num)))
  # 
  # data.1.12hr %>%
  #   left_join(date_num_df, by = "date_num") -> data.1.12hr
  # 
  # 
  # First half
  data.1.12hr %>%
    filter(hr < 7) %>%
    group_by(sector, year, date) %>%
    select(sector, year, date, counts, trench_length_m, date_num) %>% #-> tmp.1
    summarize(across(.cols = counts, 
                     .fns = sum, 
                     na.rm = T, 
                     .names = "counts_1"),
              length = first(trench_length_m),
              date_num = first(date_num)) -> sum.0.6.hrs.12hr
  
  data.1.12hr %>%
    filter(hr < 7) %>%
    group_by(sector, year, date) %>%
    summarize(across(.cols = counts_m, 
                     .fns = sum, 
                     na.rm = T, 
                     .names = "counts_m_1")) -> counts.m.1
  
  sum.0.6.hrs.12hr$counts_m_1 <- counts.m.1$counts_m_1
  
  # Second half
  data.1.12hr %>%
    filter(hr > 6) %>%
    group_by(sector, year, date) %>%
    summarize(across(.cols = counts, 
                     .fns = sum, 
                     na.rm = T, 
                     .names = "counts_2")) -> sum.7.12.hrs.12hr
  
  data.1.12hr %>%
    filter(hr > 6) %>%
    group_by(sector, year, date) %>%
    summarize(across(.cols = counts_m, 
                     .fns = sum, 
                     na.rm = T, .names = "counts_m_2")) -> counts.m.2
  
  sum.7.12.hrs.12hr$counts_m_2 <- counts.m.2$counts_m_2
  
  # put them together and compute the ratio
  sum.0.6.hrs.12hr %>% 
    left_join(sum.7.12.hrs.12hr, by = c("sector", "year", "date")) %>% 
    mutate(ratio = counts_m_1/(counts_m_1 + counts_m_2))-> counts.1st.2nd
  
  return(counts.1st.2nd)  
}

define.df_q50 <- function(x) { 
  y <- data.frame(t(x))
  if (ncol(y) == 5){
    colnames(y) <- c("2014R_q50", "2017R_q50", 
                     "North_q50", "South_q50", "West_q50")
  } else if (ncol(y) == 6) {
    colnames(y) <- c("2014R_q50", "2017R_q50", "2019R_q50",
                     "North_q50", "South_q50", "West_q50")
    
  }
  y$DOSeason <- seq(1, nrow(y))
  return(y)}

define.df_q2.5 <- function(x) { 
  y <- data.frame(t(x))
  if (ncol(y) == 5){
    colnames(y) <- c("2014R_q2.5", "2017R_q2.5", 
                     "North_q2.5", "South_q2.5", "West_q2.5")
  } else if (ncol(y) == 6) {
    colnames(y) <- c("2014R_q2.5", "2017R_q2.5", "2019R_q2.5",
                     "North_q2.5", "South_q2.5", "West_q2.5")
    
  }
  y$DOSeason <- seq(1, nrow(y))
  return(y)}

define.df_q97.5 <- function(x) { 
  y <- data.frame(t(x))
  if (ncol(y) == 5){
    colnames(y) <- c("2014R_q97.5", "2017R_q97.5", 
                     "North_q97.5", "South_q97.5", "West_q97.5")
  } else if (ncol(y) == 6) {
    colnames(y) <- c("2014R_q97.5", "2017R_q97.5", "2019R_q97.5",
                     "North_q97.5", "South_q97.5", "West_q97.5")
    
  }
  y$DOSeason <- seq(1, nrow(y))
  return(y)}


extract.statistics <- function(jm){
  ys_q50 <- apply(jm$q50$y, MARGIN = 1, FUN = define.df_q50)
  ys_q2.5 <- apply(jm$q2.5$y, MARGIN = 1, FUN = define.df_q2.5)
  ys_q97.5 <- apply(jm$q97.5$y, MARGIN = 1, FUN = define.df_q97.5)
  
  # ys_q2.5.long <- melt(ys_q2.5, value.name = "y_q2.5")
  # ys_q50.long <- melt(ys_q50, value.name = "y_q50")
  # ys_q97.5.long <- melt(ys_q97.5, value.name = "y_q97.5")
  
  Xs_q50 <- apply(jm$q50$X, MARGIN = 1, 
                  FUN = define.df_q50)
  
  Xs_q2.5 <- apply(jm$q2.5$X, MARGIN = 1, FUN = define.df_q2.5)
  
  Xs_q97.5 <- apply(jm$q97.5$X, MARGIN = 1, FUN = define.df_q97.5)
  
  # Xs_q2.5.long <- melt(Xs_q2.5, value.name = "X_q2.5")
  # 
  # Xs_q50.long <- melt(Xs_q50, value.name = "X_q50")
  # 
  # Xs_q97.5.long <- melt(Xs_q97.5, value.name = "y_q97.5")
  
  seasons <- c(2014, 2015, 2016, 2017)
  for (s in 1:length(seasons)){
    Xs_q50[[s]]$Season <- seasons[s]
    Xs_q2.5[[s]]$Season  <- seasons[s]
    Xs_q97.5[[s]]$Season  <- seasons[s]
    
    ys_q50[[s]]$Season <- seasons[s]
    ys_q2.5[[s]]$Season  <- seasons[s]
    ys_q97.5[[s]]$Season  <- seasons[s]
  }
  
  Xs_q2.5_df <- do.call(rbind, Xs_q2.5)
  Xs_q50_df <- do.call(rbind, Xs_q50)
  Xs_q97.5_df <- do.call(rbind, Xs_q97.5)
  
  Xs_2014R <- cbind(Xs_q2.5_df$`2014R_q2.5`, Xs_q50_df$`2014R_q50`,
                    Xs_q97.5_df$`2014R_q97.5`, Xs_q2.5_df[, c("DOSeason", "Season")])
  colnames(Xs_2014R) <- c("Xs_q2.5", "Xs_q50", "Xs_q97.5", "DOSeason", "Season")
  Xs_2014R$Sector <- "2014R"
  
  Xs_2017R <- cbind(Xs_q2.5_df$`2017R_q2.5`, Xs_q50_df$`2017R_q50`,
                    Xs_q97.5_df$`2017R_q97.5`, Xs_q2.5_df[, c("DOSeason", "Season")])
  colnames(Xs_2017R) <- c("Xs_q2.5", "Xs_q50", "Xs_q97.5", "DOSeason", "Season")
  Xs_2017R$Sector <- "2017R"
  
  Xs_North <- cbind(Xs_q2.5_df$North_q2.5, Xs_q50_df$North_q50,
                    Xs_q97.5_df$North_q97.5, Xs_q2.5_df[, c("DOSeason", "Season")])
  colnames(Xs_North) <- c("Xs_q2.5", "Xs_q50", "Xs_q97.5", "DOSeason", "Season")
  Xs_North$Sector <- "North"
  
  Xs_South <- cbind(Xs_q2.5_df$South_q2.5, Xs_q50_df$South_q50,
                    Xs_q97.5_df$South_q97.5, Xs_q2.5_df[, c("DOSeason", "Season")])
  colnames(Xs_South) <- c("Xs_q2.5", "Xs_q50", "Xs_q97.5", "DOSeason", "Season")
  Xs_South$Sector <- "South"
  
  Xs_West <- cbind(Xs_q2.5_df$West_q2.5, Xs_q50_df$West_q50,
                   Xs_q97.5_df$West_q97.5, Xs_q2.5_df[, c("DOSeason", "Season")])
  colnames(Xs_West) <- c("Xs_q2.5", "Xs_q50", "Xs_q97.5", "DOSeason", "Season")
  Xs_West$Sector <- "West"
  
  return(out = list(Xs_2014R = Xs_2014R,
                    Xs_2017R = Xs_2017R,
                    Xs_North = Xs_North,
                    Xs_South = Xs_South,
                    Xs_West = Xs_West))
  
}

extract.statistics.2021 <- function(jm){
  ys_q50 <- apply(jm$q50$y, MARGIN = 1, FUN = define.df_q50)
  ys_q2.5 <- apply(jm$q2.5$y, MARGIN = 1, FUN = define.df_q2.5)
  ys_q97.5 <- apply(jm$q97.5$y, MARGIN = 1, FUN = define.df_q97.5)
  
  Xs_q50 <- apply(jm$q50$X, MARGIN = 1, 
                  FUN = define.df_q50)
  
  Xs_q2.5 <- apply(jm$q2.5$X, MARGIN = 1, FUN = define.df_q2.5)
  
  Xs_q97.5 <- apply(jm$q97.5$X, MARGIN = 1, FUN = define.df_q97.5)
  
  # Xs_q2.5.long <- melt(Xs_q2.5, value.name = "X_q2.5")
  # 
  # Xs_q50.long <- melt(Xs_q50, value.name = "X_q50")
  # 
  # Xs_q97.5.long <- melt(Xs_q97.5, value.name = "y_q97.5")
  
  seasons <- c(2014, 2015, 2016, 2017, 2018, 2019, 2020)
  for (s in 1:length(seasons)){
    Xs_q50[[s]]$Season <- seasons[s]
    Xs_q2.5[[s]]$Season  <- seasons[s]
    Xs_q97.5[[s]]$Season  <- seasons[s]
    
    ys_q50[[s]]$Season <- seasons[s]
    ys_q2.5[[s]]$Season  <- seasons[s]
    ys_q97.5[[s]]$Season  <- seasons[s]
  }
  
  Xs_q2.5_df <- do.call(rbind, Xs_q2.5)
  Xs_q50_df <- do.call(rbind, Xs_q50)
  Xs_q97.5_df <- do.call(rbind, Xs_q97.5)
  
  Xs_2014R <- cbind(Xs_q2.5_df$`2014R_q2.5`, Xs_q50_df$`2014R_q50`,
                    Xs_q97.5_df$`2014R_q97.5`, Xs_q2.5_df[, c("DOSeason", "Season")])
  colnames(Xs_2014R) <- c("Xs_q2.5", "Xs_q50", "Xs_q97.5", "DOSeason", "Season")
  Xs_2014R$Sector <- "2014R"
  
  Xs_2017R <- cbind(Xs_q2.5_df$`2017R_q2.5`, Xs_q50_df$`2017R_q50`,
                    Xs_q97.5_df$`2017R_q97.5`, Xs_q2.5_df[, c("DOSeason", "Season")])
  colnames(Xs_2017R) <- c("Xs_q2.5", "Xs_q50", "Xs_q97.5", "DOSeason", "Season")
  Xs_2017R$Sector <- "2017R"
  
  Xs_2019R <- cbind(Xs_q2.5_df$`2019R_q2.5`, Xs_q50_df$`2019R_q50`,
                    Xs_q97.5_df$`2019R_q97.5`, Xs_q2.5_df[, c("DOSeason", "Season")])
  colnames(Xs_2019R) <- c("Xs_q2.5", "Xs_q50", "Xs_q97.5", "DOSeason", "Season")
  Xs_2019R$Sector <- "2019R"

  Xs_North <- cbind(Xs_q2.5_df$North_q2.5, Xs_q50_df$North_q50,
                    Xs_q97.5_df$North_q97.5, Xs_q2.5_df[, c("DOSeason", "Season")])
  colnames(Xs_North) <- c("Xs_q2.5", "Xs_q50", "Xs_q97.5", "DOSeason", "Season")
  Xs_North$Sector <- "North"
  
  Xs_South <- cbind(Xs_q2.5_df$South_q2.5, Xs_q50_df$South_q50,
                    Xs_q97.5_df$South_q97.5, Xs_q2.5_df[, c("DOSeason", "Season")])
  colnames(Xs_South) <- c("Xs_q2.5", "Xs_q50", "Xs_q97.5", "DOSeason", "Season")
  Xs_South$Sector <- "South"
  
  Xs_West <- cbind(Xs_q2.5_df$West_q2.5, Xs_q50_df$West_q50,
                   Xs_q97.5_df$West_q97.5, Xs_q2.5_df[, c("DOSeason", "Season")])
  colnames(Xs_West) <- c("Xs_q2.5", "Xs_q50", "Xs_q97.5", "DOSeason", "Season")
  Xs_West$Sector <- "West"
  
  return(out = list(Xs_2014R = Xs_2014R,
                    Xs_2017R = Xs_2017R,
                    Xs_2019R = Xs_2019R,
                    Xs_North = Xs_North,
                    Xs_South = Xs_South,
                    Xs_West = Xs_West))
  
}

# Beach length for 2019R needs to be corrected when I have the info. 2021-06-02
Beach.Length <- data.frame(name = c("2019R", "2017R", "2014R", "South", "North", "West"),
                           length = c(160, 250, 150, 350, 700, 350))

summary.stats <- function(out.stats){
  out.stats$Xs_2014R %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum(exp(Xs_q50)),
              q2.5 = sum(exp(Xs_q2.5)),
              q97.5 = sum(exp(Xs_q97.5)),
              total_q50 = sum(exp(Xs_q50)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100,
              total_q2.5 = sum(exp(Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100,
              total_q97.5 = sum(exp(Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100) %>%
    mutate(Sector = "2014R") -> Xs_2014R_summary
  
  colnames(Xs_2014R_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_North %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum(exp(Xs_q50)),
              q2.5 = sum(exp(Xs_q2.5)),
              q97.5 = sum(exp(Xs_q97.5)),
              total_q50 = sum(exp(Xs_q50)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100,
              total_q2.5 = sum(exp(Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100,
              total_q97.5 = sum(exp(Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100) %>%
    mutate(Sector = "North") -> Xs_North_summary
  colnames(Xs_North_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_South %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum(exp(Xs_q50)),
              q2.5 = sum(exp(Xs_q2.5)),
              q97.5 = sum(exp(Xs_q97.5)),
              total_q50 = sum(exp(Xs_q50)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100,
              total_q2.5 = sum(exp(Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100,
              total_q97.5 = sum(exp(Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100) %>%
    mutate(Sector = "South") -> Xs_South_summary
  colnames(Xs_South_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_West %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum(exp(Xs_q50)),
              q2.5 = sum(exp(Xs_q2.5)),
              q97.5 = sum(exp(Xs_q97.5)),
              total_q50 = sum(exp(Xs_q50)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100,
              total_q2.5 = sum(exp(Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100,
              total_q97.5 = sum(exp(Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100) %>%
    mutate(Sector = "West") -> Xs_West_summary
  colnames(Xs_West_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                 "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_2017R %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum(exp(Xs_q50)),
              q2.5 = sum(exp(Xs_q2.5)),
              q97.5 = sum(exp(Xs_q97.5)),
              total_q50 = sum(exp(Xs_q50)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100,
              total_q2.5 = sum(exp(Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100,
              total_q97.5 = sum(exp(Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100) %>%
    mutate(Sector = "2017R") -> Xs_2017R_summary
  colnames(Xs_2017R_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  all.summary <- rbind(Xs_2017R_summary, 
                       Xs_2014R_summary,
                       Xs_North_summary,
                       Xs_South_summary,
                       Xs_West_summary)
  
  return(all.summary)
}

summary.stats.2021 <- function(out.stats){
  out.stats$Xs_2014R %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum(exp(Xs_q50)),
              q2.5 = sum(exp(Xs_q2.5)),
              q97.5 = sum(exp(Xs_q97.5)),
              total_q50 = sum(exp(Xs_q50)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100,
              total_q2.5 = sum(exp(Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100,
              total_q97.5 = sum(exp(Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100) %>%
    mutate(Sector = "2014R") -> Xs_2014R_summary
  
  colnames(Xs_2014R_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_North %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum(exp(Xs_q50)),
              q2.5 = sum(exp(Xs_q2.5)),
              q97.5 = sum(exp(Xs_q97.5)),
              total_q50 = sum(exp(Xs_q50)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100,
              total_q2.5 = sum(exp(Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100,
              total_q97.5 = sum(exp(Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100) %>%
    mutate(Sector = "North") -> Xs_North_summary
  colnames(Xs_North_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_South %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum(exp(Xs_q50)),
              q2.5 = sum(exp(Xs_q2.5)),
              q97.5 = sum(exp(Xs_q97.5)),
              total_q50 = sum(exp(Xs_q50)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100,
              total_q2.5 = sum(exp(Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100,
              total_q97.5 = sum(exp(Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100) %>%
    mutate(Sector = "South") -> Xs_South_summary
  colnames(Xs_South_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_West %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum(exp(Xs_q50)),
              q2.5 = sum(exp(Xs_q2.5)),
              q97.5 = sum(exp(Xs_q97.5)),
              total_q50 = sum(exp(Xs_q50)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100,
              total_q2.5 = sum(exp(Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100,
              total_q97.5 = sum(exp(Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100) %>%
    mutate(Sector = "West") -> Xs_West_summary
  colnames(Xs_West_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                 "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_2017R %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum(exp(Xs_q50)),
              q2.5 = sum(exp(Xs_q2.5)),
              q97.5 = sum(exp(Xs_q97.5)),
              total_q50 = sum(exp(Xs_q50)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100,
              total_q2.5 = sum(exp(Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100,
              total_q97.5 = sum(exp(Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100) %>%
    mutate(Sector = "2017R") -> Xs_2017R_summary
  colnames(Xs_2017R_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_2019R %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum(exp(Xs_q50)),
              q2.5 = sum(exp(Xs_q2.5)),
              q97.5 = sum(exp(Xs_q97.5)),
              total_q50 = sum(exp(Xs_q50)) * Beach.Length[which(Beach.Length$name == "2019R"), "length"]/100,
              total_q2.5 = sum(exp(Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "2019R"), "length"]/100,
              total_q97.5 = sum(exp(Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "2019R"), "length"]/100) %>%
    mutate(Sector = "2019R") -> Xs_2019R_summary
  colnames(Xs_2019R_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  all.summary <- rbind(Xs_2014R_summary, 
                       Xs_2017R_summary,
                       Xs_2019R_summary,
                       Xs_North_summary,
                       Xs_South_summary,
                       Xs_West_summary)
  
  return(all.summary)
}

summary.stats.nolog <- function(out.stats){
  # estimated numbers are in per 100 m. So, the total number for the entire sector is computed
  # by multiplying the estimates with (total length)/100.  
  out.stats$Xs_2014R %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum((Xs_q50)),
              q2.5 = sum((Xs_q2.5)),
              q97.5 = sum((Xs_q97.5)),
              total_q50 = sum((Xs_q50)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100,
              total_q2.5 = sum((Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100,
              total_q97.5 = sum((Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100) %>%
    mutate(Sector = "2014R") -> Xs_2014R_summary
  
  colnames(Xs_2014R_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_North %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum((Xs_q50)),
              q2.5 = sum((Xs_q2.5)),
              q97.5 = sum((Xs_q97.5)),
              total_q50 = sum((Xs_q50)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100,
              total_q2.5 = sum((Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100,
              total_q97.5 = sum((Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100) %>%
    mutate(Sector = "North") -> Xs_North_summary
  colnames(Xs_North_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_South %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum((Xs_q50)),
              q2.5 = sum((Xs_q2.5)),
              q97.5 = sum((Xs_q97.5)),
              total_q50 = sum((Xs_q50)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100,
              total_q2.5 = sum((Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100,
              total_q97.5 = sum((Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100) %>%
    mutate(Sector = "South") -> Xs_South_summary
  colnames(Xs_South_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_West %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum((Xs_q50)),
              q2.5 = sum((Xs_q2.5)),
              q97.5 = sum((Xs_q97.5)),
              total_q50 = sum((Xs_q50)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100,
              total_q2.5 = sum((Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100,
              total_q97.5 = sum((Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100) %>%
    mutate(Sector = "West") -> Xs_West_summary
  colnames(Xs_West_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                 "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_2017R %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum((Xs_q50)),
              q2.5 = sum((Xs_q2.5)),
              q97.5 = sum((Xs_q97.5)),
              total_q50 = sum((Xs_q50)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100,
              total_q2.5 = sum((Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100,
              total_q97.5 = sum((Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100) %>%
    mutate(Sector = "2017R") -> Xs_2017R_summary
  colnames(Xs_2017R_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  all.summary <- rbind(Xs_2017R_summary,
                       Xs_2014R_summary,
                       Xs_North_summary,
                       Xs_South_summary,
                       Xs_West_summary)
  
  all.summary[all.summary < 0] <- 0
  
  return(all.summary)
}


summary.stats.nolog.2021 <- function(out.stats){
  # estimated numbers are in per 100 m. So, the total number for the entire sector is computed
  # by multiplying the estimates with (total length)/100.  
  out.stats$Xs_2014R %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum((Xs_q50)),
              q2.5 = sum((Xs_q2.5)),
              q97.5 = sum((Xs_q97.5)),
              total_q50 = sum((Xs_q50)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100,
              total_q2.5 = sum((Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100,
              total_q97.5 = sum((Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "2014R"), "length"]/100) %>%
    mutate(Sector = "2014R") -> Xs_2014R_summary
  
  colnames(Xs_2014R_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_North %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum((Xs_q50)),
              q2.5 = sum((Xs_q2.5)),
              q97.5 = sum((Xs_q97.5)),
              total_q50 = sum((Xs_q50)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100,
              total_q2.5 = sum((Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100,
              total_q97.5 = sum((Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "North"), "length"]/100) %>%
    mutate(Sector = "North") -> Xs_North_summary
  colnames(Xs_North_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_South %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum((Xs_q50)),
              q2.5 = sum((Xs_q2.5)),
              q97.5 = sum((Xs_q97.5)),
              total_q50 = sum((Xs_q50)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100,
              total_q2.5 = sum((Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100,
              total_q97.5 = sum((Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "South"), "length"]/100) %>%
    mutate(Sector = "South") -> Xs_South_summary
  colnames(Xs_South_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_West %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum((Xs_q50)),
              q2.5 = sum((Xs_q2.5)),
              q97.5 = sum((Xs_q97.5)),
              total_q50 = sum((Xs_q50)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100,
              total_q2.5 = sum((Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100,
              total_q97.5 = sum((Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "West"), "length"]/100) %>%
    mutate(Sector = "West") -> Xs_West_summary
  colnames(Xs_West_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                 "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_2017R %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum((Xs_q50)),
              q2.5 = sum((Xs_q2.5)),
              q97.5 = sum((Xs_q97.5)),
              total_q50 = sum((Xs_q50)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100,
              total_q2.5 = sum((Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100,
              total_q97.5 = sum((Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "2017R"), "length"]/100) %>%
    mutate(Sector = "2017R") -> Xs_2017R_summary
  colnames(Xs_2017R_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  
  out.stats$Xs_2019R %>% group_by(as.factor(Season)) %>%
    summarise(q50 = sum((Xs_q50)),
              q2.5 = sum((Xs_q2.5)),
              q97.5 = sum((Xs_q97.5)),
              total_q50 = sum((Xs_q50)) * Beach.Length[which(Beach.Length$name == "2019R"), "length"]/100,
              total_q2.5 = sum((Xs_q2.5)) * Beach.Length[which(Beach.Length$name == "2019R"), "length"]/100,
              total_q97.5 = sum((Xs_q97.5)) * Beach.Length[which(Beach.Length$name == "2019R"), "length"]/100) %>%
    mutate(Sector = "2019R") -> Xs_2019R_summary
  
  colnames(Xs_2019R_summary) <- c("Season", "q50", "q2.5", "q97.5", 
                                  "total_q50", "total_q2.5", "total_q97.5", "Sector")
  all.summary <- rbind(Xs_2014R_summary,
                       Xs_2017R_summary,
                       Xs_2019R_summary,
                       Xs_North_summary,
                       Xs_South_summary,
                       Xs_West_summary)
  
  all.summary[all.summary < 0] <- 0
  
  return(all.summary)
}

Girondot_fcn <- function(d, S, K, P, min, max){
  K <- abs(K)
  S <- abs(S)
  S1 <- -S
  M1 <- (1 + (2 * exp(K) - 1) * exp((1/S1) * (P - d))) ^ (-1/exp(K))
  M2 <- (1 + (2 * exp(K) - 1) * exp((1/S) * (P - d))) ^ (-1/exp(K))
  N <- min + (max - min) * (M1 * M2)
  return(N)
}

# Extracting posterior samples of deviance or any other variable from jags output:
extract.samples <- function(varname, zm){
  dev <- unlist(lapply(zm, FUN = function(x) x[, varname]))
  return(dev)
}



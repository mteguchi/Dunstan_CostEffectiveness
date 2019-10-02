

define.df_q50 <- function(x) { 
  y <- data.frame(t(x))
  colnames(y) <- c("2014R_q50", "2017R_q50", 
                   "North_q50", "South_q50", "West_q50")
  y$DOSeason <- seq(1, nrow(y))
  return(y)}

define.df_q2.5 <- function(x) { 
  y <- data.frame(t(x))
  colnames(y) <- c("2014R_q2.5", "2017R_q2.5", 
                   "North_q2.5", "South_q2.5", "West_q2.5")
  y$DOSeason <- seq(1, nrow(y))
  return(y)}

define.df_q97.5 <- function(x) { 
  y <- data.frame(t(x))
  colnames(y) <- c("2014R_q97.5", "2017R_q97.5", 
                   "North_q97.5", "South_q97.5", "West_q97.5")
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
  for (s in 1:4){
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

Beach.Length <- data.frame(name = c("2017R", "2014R", "South", "North", "West"),
                           length = c(250, 150, 350, 700, 350))

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
  
  all.summary <- rbind(Xs_2014R_summary,
                       Xs_North_summary,
                       Xs_South_summary,
                       Xs_West_summary)
  
  return(all.summary)
}


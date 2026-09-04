library(rstan)
library(ggplot2)
library(patchwork)


# Fit files for each species and model

fit_files <- list(
  Arabiensis = list(
    `Model 1` = "ArabiensisNewResults/Model1/Model1_fit.rds",
    `Model 2` = "ArabiensisNewResults/Model2/Model2_fit.rds",
    `Model 3` = "ArabiensisNewResults/Model3/Model3_fit.rds",
    `Model 4` = "ArabiensisNewResults/Model4/Model4_fit.rds"
  ),
  
  Culex = list(
    `Model 1` = "CulexNewResults/Model1/Model1_fit.rds",
    `Model 2` = "CulexNewResults/Model2/Model2_fit.rds",
    `Model 3` = "CulexNewResults/Model3/Model3_fit.rds",
    `Model 4` = "CulexNewResults/Model4/Model4_fit.rds"
  )
)


# Observed control and intervention data files

data_files <- list(
  Arabiensis = list(
    control = "Data/arabiensis_cleanEH_BIT046_data0.csv",
    intervention = "Data/arabiensis_cleanEH_BIT046_data1.csv"
  ),
  
  Culex = list(
    control = "Data/culex_cleanEH_BIT046_data0.csv",
    intervention = "Data/culex_cleanEH_BIT046_data1.csv"
  )
)


# Labels used in the final PPC figures

model_labels <- c(
  "Model 1" = "Model one",
  "Model 2" = "Model two",
  "Model 3" = "Model three",
  "Model 4" = "Model four"
)


# Number of posterior predictive replicates

n_rep <- 500
set.seed(12345)


# Extract posterior matrix columns matching a parameter name

extract_cols <- function(post, pattern) {
  
  post[
    ,
    grep(pattern, colnames(post)),
    drop = FALSE
  ]
}


# Summarise observed and posterior predictive count frequencies
# for each possible mosquito count

make_ppc_table <- function(observed, replicated) {
  
  max_count <- max(
    observed,
    replicated,
    na.rm = TRUE
  )
  
  vals <- 0:max_count
  
  
  # Frequency of each count in every posterior predictive replicate
  
  rep_freq <- sapply(
    vals,
    function(k) {
      rowSums(replicated == k)
    }
  )
  
  
  # Frequency of each count in the observed data
  
  observed_freq <- as.integer(
    table(
      factor(
        observed,
        levels = vals
      )
    )
  )
  
  
  # Posterior predictive mean and 95% interval for each frequency
  
  data.frame(
    count = vals,
    
    observed = observed_freq,
    
    predicted_mean = apply(
      rep_freq,
      2,
      mean
    ),
    
    predicted_lower = apply(
      rep_freq,
      2,
      quantile,
      probs = 0.025
    ),
    
    predicted_upper = apply(
      rep_freq,
      2,
      quantile,
      probs = 0.975
    )
  )
}


# Plot observed frequencies against posterior predictive frequencies

plot_ppc <- function(tab, title, x_label, x_limit) {
  
  ggplot(
    tab,
    aes(x = count)
  ) +
    
    # 95% posterior predictive interval
    
    geom_ribbon(
      aes(
        ymin = predicted_lower,
        ymax = predicted_upper
      ),
      fill = "orange",
      alpha = 0.30
    ) +
    
    # Observed frequencies
    
    geom_col(
      aes(y = observed),
      fill = "steelblue",
      width = 0.85,
      alpha = 0.85
    ) +
    
    # Posterior predictive mean
    
    geom_line(
      aes(y = predicted_mean),
      colour = "darkorange",
      linewidth = 0.8
    ) +
    
    geom_point(
      aes(y = predicted_mean),
      colour = "darkorange",
      size = 1.3
    ) +
    
    labs(
      title = title,
      x = x_label,
      y = "Number of observations"
    ) +
    
    coord_cartesian(
      xlim = c(0, x_limit)
    ) +
    
    theme_minimal(
      base_size = 10
    ) +
    
    theme(
      plot.title = element_text(
        size = 11,
        face = "bold"
      ),
      
      panel.grid.minor = element_blank()
    )
}


# Produce total-count and fed-count PPCs for one fitted model

make_model_ppcs <- function(
    fit_file,
    data0_file,
    data1_file,
    species,
    model_name,
    n_rep = 500
) {
  
  # Load observed control and intervention data
  
  Data0 <- read.csv(data0_file)
  Data1 <- read.csv(data1_file)
  
  
  # Observed total mosquito counts
  # UA mosquitoes are excluded from the fitted outcome
  
  total0_obs <-
    Data0$FA +
    Data0$FD +
    Data0$UD
  
  total1_obs <-
    Data1$FA +
    Data1$FD +
    Data1$UD
  
  
  # Observed fed mosquito counts
  
  fed0_obs <-
    Data0$FA +
    Data0$FD
  
  fed1_obs <-
    Data1$FA +
    Data1$FD
  
  
  # Combine control and intervention observations
  
  observed_total <- c(
    total0_obs,
    total1_obs
  )
  
  observed_fed <- c(
    fed0_obs,
    fed1_obs
  )
  
  
  # Load saved Stan fit and convert posterior draws to a matrix
  
  fit <- readRDS(fit_file)
  
  posterior <- as.matrix(fit)
  
  
  # Randomly select posterior draws for the PPC
  
  idx <- sample(
    seq_len(nrow(posterior)),
    size = n_rep,
    replace = FALSE
  )
  
  
  # Extract posterior quantities required to generate replicated data
  
  log_totalC <- extract_cols(
    posterior,
    "^log_totalC\\["
  )
  
  log_qC <- extract_cols(
    posterior,
    "^log_qC\\["
  )
  
  log_totalI <- extract_cols(
    posterior,
    "^log_totalI\\["
  )
  
  log_pI <- extract_cols(
    posterior,
    "^log_pI\\["
  )
  
  phiC <- posterior[
    ,
    "phi_valuesC",
    drop = FALSE
  ]
  
  phiI <- extract_cols(
    posterior,
    "^phiI\\["
  )
  
  
  # Empty matrices for replicated total mosquito counts
  
  totalC_rep <- matrix(
    NA_integer_,
    nrow = n_rep,
    ncol = length(total0_obs)
  )
  
  totalI_rep <- matrix(
    NA_integer_,
    nrow = n_rep,
    ncol = length(total1_obs)
  )
  
  
  # Generate posterior predictive total counts
  # using the negative binomial model
  
  for (s in seq_len(n_rep)) {
    
    d <- idx[s]
    
    totalC_rep[s, ] <- rnbinom(
      length(total0_obs),
      mu = exp(log_totalC[d, ]),
      size = phiC[d, 1]
    )
    
    totalI_rep[s, ] <- rnbinom(
      length(total1_obs),
      mu = exp(log_totalI[d, ]),
      size = phiI[d, ]
    )
  }
  
  
  # Combine control and intervention replicated totals
  
  total_rep <- cbind(
    totalC_rep,
    totalI_rep
  )
  
  
  # Empty matrices for replicated fed mosquito counts
  
  fedC_rep <- matrix(
    NA_integer_,
    nrow = n_rep,
    ncol = length(fed0_obs)
  )
  
  fedI_rep <- matrix(
    NA_integer_,
    nrow = n_rep,
    ncol = length(fed1_obs)
  )
  
  
  # Generate posterior predictive fed counts
  # conditional on the observed total mosquito counts
  
  for (s in seq_len(n_rep)) {
    
    d <- idx[s]
    
    fedC_rep[s, ] <- rbinom(
      length(fed0_obs),
      size = total0_obs,
      prob = pmin(
        pmax(
          exp(log_qC[d, ]),
          0
        ),
        1
      )
    )
    
    fedI_rep[s, ] <- rbinom(
      length(fed1_obs),
      size = total1_obs,
      prob = pmin(
        pmax(
          exp(log_pI[d, ]),
          0
        ),
        1
      )
    )
  }
  
  
  # Combine control and intervention replicated fed counts
  
  fed_rep <- cbind(
    fedC_rep,
    fedI_rep
  )
  
  
  # Summarise posterior predictive frequencies
  
  fed_tab <- make_ppc_table(
    observed_fed,
    fed_rep
  )
  
  total_tab <- make_ppc_table(
    observed_total,
    total_rep
  )
  
  
  # Create panel title
  
  plot_title <- paste(
    species,
    "—",
    model_labels[model_name]
  )
  
  
  # Fed-count PPC
  
  fed_plot <- plot_ppc(
    fed_tab,
    title = plot_title,
    x_label = "Fed count",
    x_limit = 30
  )
  
  
  # Total-count PPC
  
  total_plot <- plot_ppc(
    total_tab,
    title = plot_title,
    x_label = "Total count",
    x_limit = 80
  )
  
  
  list(
    fed = fed_plot,
    total = total_plot
  )
}


# Run the PPC procedure for all four models and both species

results <- list()

for (species in names(fit_files)) {
  
  results[[species]] <- list()
  
  for (model_name in names(fit_files[[species]])) {
    
    results[[species]][[model_name]] <-
      make_model_ppcs(
        fit_file = fit_files[[species]][[model_name]],
        data0_file = data_files[[species]]$control,
        data1_file = data_files[[species]]$intervention,
        species = species,
        model_name = model_name,
        n_rep = n_rep
      )
  }
}


# Combine fed-count PPCs
# Rows = models, columns = species

fed_plot <-
  results$Arabiensis$`Model 1`$fed +
  results$Culex$`Model 1`$fed +
  results$Arabiensis$`Model 2`$fed +
  results$Culex$`Model 2`$fed +
  results$Arabiensis$`Model 3`$fed +
  results$Culex$`Model 3`$fed +
  results$Arabiensis$`Model 4`$fed +
  results$Culex$`Model 4`$fed +
  
  plot_layout(
    ncol = 2,
    nrow = 4
  ) +
  
  plot_annotation(
    title = "Posterior predictive checks for fed mosquito counts",
    subtitle = paste(
      "Blue bars: observed frequencies; orange line: posterior predictive mean;",
      "orange band: 95% posterior predictive interval"
    ),
    
    theme = theme(
      plot.title = element_text(
        size = 14
      ),
      
      plot.subtitle = element_text(
        size = 10
      )
    )
  )


# Combine total-count PPCs
# Rows = models, columns = species

total_plot <-
  results$Arabiensis$`Model 1`$total +
  results$Culex$`Model 1`$total +
  results$Arabiensis$`Model 2`$total +
  results$Culex$`Model 2`$total +
  results$Arabiensis$`Model 3`$total +
  results$Culex$`Model 3`$total +
  results$Arabiensis$`Model 4`$total +
  results$Culex$`Model 4`$total +
  
  plot_layout(
    ncol = 2,
    nrow = 4
  ) +
  
  plot_annotation(
    title = "Posterior predictive checks for total mosquito counts",
    subtitle = paste(
      "Blue bars: observed frequencies; orange line: posterior predictive mean;",
      "orange band: 95% posterior predictive interval"
    ),
    
    theme = theme(
      plot.title = element_text(
        size = 14
      ),
      
      plot.subtitle = element_text(
        size = 10
      )
    )
  )


# Display final figures

fed_plot
total_plot



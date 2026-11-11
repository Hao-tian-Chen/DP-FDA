Sys.setenv(LANG = "en")

library(ggplot2)
library(patchwork)

###############################################################################
## Load Data and Preprocess

#### Data for Plot 1 ####
plot1_data_1 <- readRDS("data/Poly/Poly_sigma0_0_B1p1.rds")|>
  dplyr::filter(mu %in% c(1))
plot1_data_2 <- readRDS("data/Poly/Poly_sigma0_0_B1p1.rds")|>
  dplyr::filter(mu %in% c(0.1))

## Reference line y = 1/m for each basis dimension
ref_data_1_plot1 <- data.frame(
  m = sort(unique(plot1_data_1$m))
)
ref_data_1_plot1$ref_y <- 1 / ref_data_1_plot1$m
ref_data_1_plot1$x_min <- min(plot1_data_1$n)
ref_data_1_plot1$x_max <- max(plot1_data_1$n)

ref_data_2_plot1 <- data.frame(
  m = sort(unique(plot1_data_2$m))
)
ref_data_2_plot1$ref_y <- 1 / ref_data_2_plot1$m
ref_data_2_plot1$x_min <- min(plot1_data_2$n)
ref_data_2_plot1$x_max <- max(plot1_data_2$n)

#### Data for Plot 2 ####
plot2_data <- readRDS("data/Poly/Poly_Bprime.rds")
ref_data_plot2 <- data.frame(
  m = sort(unique(plot2_data$m))
)
ref_data_plot2$ref_y <- 1 / ref_data_plot2$m
ref_data_plot2$x_min <- min(plot2_data$n)
ref_data_plot2$x_max <- max(plot2_data$n)

#### Data for Plot 3 ####
plot3_data <- readRDS("data/Poly/Poly_noisy_sqrtlog.rds")
plot3_data_1 <- subset(plot3_data, plot3_data$sigma_obs == 0.1)
plot3_data_2 <- subset(plot3_data, plot3_data$sigma_obs == 1)
plot3_data_3 <- subset(plot3_data, plot3_data$sigma_obs == 10)

## Reference line y = 1/m for each basis dimension
ref_data_1_plot3 <- data.frame(
  m = sort(unique(plot3_data_1$m))
)
ref_data_1_plot3$ref_y <- 1 / ref_data_1_plot3$m
ref_data_1_plot3$x_min <- min(plot3_data_1$n)
ref_data_1_plot3$x_max <- max(plot3_data_1$n)

ref_data_2_plot3 <- data.frame(
  m = sort(unique(plot3_data_2$m))
)
ref_data_2_plot3$ref_y <- 1 / ref_data_2_plot3$m
ref_data_2_plot3$x_min <- min(plot3_data_2$n)
ref_data_2_plot3$x_max <- max(plot3_data_2$n)

ref_data_3_plot3 <- data.frame(
  m = sort(unique(plot3_data_3$m))
)
ref_data_3_plot3$ref_y <- 1 / ref_data_3_plot3$m
ref_data_3_plot3$x_min <- min(plot3_data_3$n)
ref_data_3_plot3$x_max <- max(plot3_data_3$n)

###############################################################################
## Figure 1: IMSE Ratio of FTP to PTF vs n,m (sigma_0 = 0, mu = 0.1, 1)
#### mu = 1 ####
plot1_ratio_1 <- ggplot(
  plot1_data_1,
  aes(
    x = n,
    y = Ratio_A3_A2,
    color = factor(m),
    shape = factor(m),
    group = factor(m)
  )
) +
  # Reference lines
  geom_segment(
    data = ref_data_1_plot1,
    aes(
      x = x_min,
      xend = x_max,
      y = ref_y,
      yend = ref_y,
      color = factor(m)
    ),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Curves and points
  geom_line(
    aes(color = factor(m)),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = factor(m)),
    size = 3,
    alpha = 0.8
  ) +
  
  # x-axis
  #scale_x_log10(breaks = c(5, 6, 7, 8, 9, 10, 15, 20, 30)) +
  scale_x_log10(
    breaks = 2^seq(2, 12, by = 2),
    labels = parse(text = paste0("2^", seq(2, 12, by = 2)))
  ) +
  scale_y_log10(
    limits = c(0.08, 200),
    breaks = c(0.1, 1, 10, 100),
    expand = c(0, 0)
  ) +
  labs(
    title = expression(mu == 1),
    x = expression("Sample size:" ~ n),
    y = expression(IMSE[FTP] / IMSE[PTF]),
    color = expression("Basis dimension:" ~ m)
  ) +
  scale_shape_manual(
    values = c(
      "2" = 16,  # Filled circle
      "4" = 17,  # Filled triangle
      "6" = 15,  # Filled square
      "8" = 18   # Filled diamond
    )
  ) +
  
  # Paper-style background
  theme_light(base_size = 18) +
  
  theme(
    # Panel background and boundary
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    # Gridlines
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # Axis appearance
    axis.text.x = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) + guides(
    color = guide_legend(
      nrow = 1, byrow = TRUE, 
      override.aes = list( 
        shape = c(16, 17, 15, 18), alpha = 1 )
      ),
    shape = "none"
  )

#### mu = 0.1 ####
plot1_ratio_2 <- ggplot(
  plot1_data_2,
  aes(
    x = n,
    y = Ratio_A3_A2,
    color = factor(m),
    shape = factor(m),
    group = factor(m)
  )
) +
  # Reference lines
  geom_segment(
    data = ref_data_2_plot1,
    aes(
      x = x_min,
      xend = x_max,
      y = ref_y,
      yend = ref_y,
      color = factor(m)
    ),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Curves and points
  geom_line(
    aes(color = factor(m)),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = factor(m)),
    size = 3,
    alpha = 0.8
  ) +
  
  # x-axis
  #scale_x_log10(breaks = c(5, 6, 7, 8, 9, 10, 15, 20, 30)) +
  scale_x_log10(
    breaks = 2^seq(2, 12, by = 2),
    labels = parse(text = paste0("2^", seq(2, 12, by = 2)))
  ) +
  scale_y_log10(
    limits = c(0.08, 200),
    breaks = c(0.1, 1, 10, 100),
    expand = c(0, 0)
  ) +
  labs(
    title = expression(mu == 0.1),
    x = expression("Sample size:" ~ n),
    y = expression(IMSE[FTP] / IMSE[PTF]),
    color = expression("Basis dimension:" ~ m)
  ) +
  scale_shape_manual(
    values = c(
      "2" = 16,  # Filled circle
      "4" = 17,  # Filled triangle
      "6" = 15,  # Filled square
      "8" = 18   # Filled diamond
    )
  ) +
  
  # Paper-style background
  theme_light(base_size = 18) +
  
  theme(
    # Panel background and boundary
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    # Gridlines
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # Axis appearance
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    axis.ticks = element_line(color = "black"),
    
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) + guides(
    color = guide_legend(
      nrow = 1, byrow = TRUE, 
      override.aes = list( 
        shape = c(16, 17, 15, 18), alpha = 1 )
    ),
    shape = "none"
  )

combined_plot1 <- plot1_ratio_2 + plot1_ratio_1 + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom",
        legend.box = "horizontal")

print(combined_plot1)

ggsave(filename = "figures/Poly/IMSE_Ratio_mu_e-1-1.pdf", plot = combined_plot1, 
       width = 12, height = 6, units = "in", device = "pdf")


###############################################################################
## Figure 2: IMSE of FTP and PTF vs n, Bprime (mu = 1, sigma_0 = 1)
#### FTP ####
plot2_IMSE_FTP <- ggplot(
  plot2_data,
  aes(
    x = n, y = n * IMSE_A3_FTP_constrained,
    color = factor(Bprime_choice),
    linetype = factor(Bprime_choice),
    shape = factor(Bprime_choice),
    group = factor(Bprime_choice)
  )
) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.8) +
  
  # Axes
  scale_x_log10(
    breaks = 2^seq(2, 12, by = 2),
    labels = parse(text = paste0("2^", seq(2, 12, by = 2)))
  ) +
  scale_y_log10(
    limits = c(0.15, 1500),
    breaks = c(1, 10, 100, 1000),
    expand = c(0, 0)
  ) +
  
  # Labels
  labs(
    #title = expression(n %*% IMSE[FTP]),
    x = expression("Sample size:" ~ n),
    y = expression(n %*% IMSE[FTP]),
    color = expression(B * "'" ~ "choices:")
  ) +
  
  # Manual scales
  scale_color_discrete(
    breaks = c("B", "sqrtlog", "log", "sqrtn", "n_2_3"),
    labels = parse(text = c(
      "B", "sqrt(log(n))", "log(n)", "sqrt(n)", "n^{2/3}"
    ))
  ) +
  scale_linetype_manual(values = c(
    B = "solid", sqrtlog = "dashed", log = "dotted",
    sqrtn = "dotdash", n_2_3 = "longdash"
  )) +
  scale_shape_manual(values = c(
    B = 16, sqrtlog = 17, log = 15, sqrtn = 18, n_2_3 = 8
  )) +
  
  # Theme
  theme_light(base_size = 18) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) +
  guides(
    color = guide_legend(
      nrow = 1, byrow = TRUE, position = "bottom",
      override.aes = list( 
        linetype = c("solid", "dashed", "dotted", "dotdash", "longdash"), 
        shape = c(16, 17, 15, 18, 8), alpha = 1 )),
    linetype = "none",
    shape = "none"
  )

#### PTF ####
plot2_IMSE_PTF <- ggplot(
  plot2_data,
  aes(
    x = n, y = n * IMSE_A2_PTF_constrained,
    color = factor(Bprime_choice),
    linetype = factor(Bprime_choice),
    shape = factor(Bprime_choice),
    group = factor(Bprime_choice)
  )
) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.8) +
  
  # Axes
  scale_x_log10(
    breaks = 2^seq(2, 12, by = 2),
    labels = parse(text = paste0("2^", seq(2, 12, by = 2)))            
  ) +
  scale_y_log10(
    limits = c(0.15, 1500),
    breaks = c(1, 10, 100, 1000),
    expand = c(0, 0)
  ) +
  
  # Labels
  labs(
    #title = expression(n %*% IMSE[PTF]),
    x = expression("Sample size:" ~ n),
    y = expression(n %*% IMSE[PTF]),
    color = expression(B * "'" ~ "choices:")
  ) +
  
  # Manual scales
  scale_color_discrete(
    breaks = c("B", "sqrtlog", "log", "sqrtn", "n_2_3"),
    labels = parse(text = c(
      "B", "sqrt(log(n))", "log(n)", "sqrt(n)", "n^{2/3}"
    ))
  ) +
  scale_linetype_manual(values = c(
    B = "solid", sqrtlog = "dashed", log = "dotted",
    sqrtn = "dotdash", n_2_3 = "longdash"
  )) +
  scale_shape_manual(values = c(
    B = 16, sqrtlog = 17, log = 15, sqrtn = 18, n_2_3 = 8
  )) +
  
  # Theme
  theme_light(base_size = 18) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) +
  guides(
    color = guide_legend(
      nrow = 1, byrow = TRUE, position = "bottom",
      override.aes = list( 
        linetype = c("solid", "dashed", "dotted", "dotdash", "longdash"), 
        shape = c(16, 17, 15, 18, 8), alpha = 1 )),
    linetype = "none",
    shape = "none"
  )

#### Ratio ####

# plot2_ratio <- ggplot(
#   plot2_data,
#   aes(
#     x = n,
#     y = Ratio_A3_A2,
#     color = factor(Bprime_choice),
#     linetype = factor(Bprime_choice),
#     shape = factor(Bprime_choice),
#     group = factor(Bprime_choice)
#   )
# ) +
#   # Reference line
#   geom_segment(
#     data = ref_data_plot2,
#     aes(x = x_min, xend = x_max, y = ref_y, yend = ref_y),
#     color = "black",
#     linetype = "dashed",
#     linewidth = 0.7,
#     show.legend = FALSE,
#     inherit.aes = FALSE
#   ) +
#   
#   # Curves and points
#   geom_line(linewidth = 0.7, alpha = 0.8) +
#   geom_point(size = 3, alpha = 0.8) +
#   
#   # Axes
#   scale_x_log10(breaks = 2^(2:11)) +
#   scale_y_log10(breaks = c(0.25, 0.5, 1), expand = c(0, 0)) +
#   
#   # Labels
#   labs(
#     #title = expression(IMSE[FTP] / IMSE[PTF]),
#     x = expression("Sample size:" ~ n),
#     y = expression(IMSE[FTP] / IMSE[PTF]),
#     color = expression(B * "'" ~ "choices:")
#   ) +
#   
#   # Same color, line type, and shape settings as plot_IMSE_FTP
#   scale_color_discrete(
#     breaks = c("B", "sqrtlog", "log", "sqrtn", "n_2_3"),
#     labels = parse(text = c(
#       "B", "sqrt(log(n))", "log(n)", "sqrt(n)", "n^{2/3}"
#     ))
#   ) +
#   scale_linetype_manual(values = c(
#     B = "solid", sqrtlog = "dashed", log = "dotted",
#     sqrtn = "dotdash", n_2_3 = "longdash"
#   )) +
#   scale_shape_manual(values = c(
#     B = 16, sqrtlog = 17, log = 15, sqrtn = 18, n_2_3 = 8
#   )) +
#   
#   # Theme
#   theme_light(base_size = 14) +
#   theme(
#     panel.background = element_rect(fill = "white", color = NA),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.text = element_text(color = "black", size = 12),
#     axis.title = element_text(color = "black", size = 16),
#     axis.ticks = element_line(color = "black"),
#     plot.title = element_text(face = "bold", hjust = 0.5)
#   ) +
#   guides(
#     color = guide_legend(
#       nrow = 1,
#       byrow = TRUE,
#       position = "bottom",
#       override.aes = list( 
#         linetype = c("solid", "dashed", "dotted", "dotdash", "longdash"), 
#         shape = c(16, 17, 15, 18, 8), alpha = 1 )),
#     linetype = "none",
#     shape = "none"
#   )

combined_plot2 <- plot2_IMSE_FTP + plot2_IMSE_PTF +
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom",
        legend.box = "horizontal")

print(combined_plot2)

ggsave(filename = "figures/Poly/Poly_Bprime_choices.pdf", plot = combined_plot2, 
       width = 12, height = 6, units = "in", device = "pdf")

###############################################################################
## Figure 3: IMSE Ratio of FTP to PTF vs n,m (mu = 1, sigma_0 = 0.01, 1, 10)
#### sigma_0 = 0.1 ####
plot3_ratio_1 <- ggplot(
  plot3_data_1,
  aes(
    x = n,
    y = Ratio_A3_A2,
    color = factor(m),
    shape = factor(m),
    group = factor(m)
  )
) +
  # Reference lines
  geom_segment(
    data = ref_data_1_plot3,
    aes(
      x = x_min,
      xend = x_max,
      y = ref_y,
      yend = ref_y,
      color = factor(m)
    ),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Curves and points
  geom_line(
    aes(color = factor(m)),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = factor(m)),
    size = 3,
    alpha = 0.8
  ) +
  
  # x-axis
  #scale_x_log10(breaks = c(5, 6, 7, 8, 9, 10, 15, 20, 30)) +
  scale_x_log10(
    breaks = 2^seq(2, 12, by = 2),
    labels = parse(text = paste0("2^", seq(2, 12, by = 2)))
  ) +
  scale_y_log10(
    limits = c(0.08, 2.8),
    breaks = c(0.125, 0.25, 0.5, 1, 2),
    expand = c(0, 0)
  ) +
  labs(
    title = expression(sigma[0] == 0.1),
    x = expression("Sample size:" ~ n),
    y = expression(IMSE[FTP] / IMSE[PTF]),
    color = expression("Basis dimension:" ~ m)
  ) + scale_shape_manual(
    values = c(
      "2" = 16,  # Filled circle
      "4" = 17,  # Filled triangle
      "6" = 15,  # Filled square
      "8" = 18   # Filled diamond
    )
  ) +
  
  
  # Paper-style background
  theme_light(base_size = 18) +
  
  theme(
    # Panel background and boundary
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    # Gridlines
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # Axis appearance
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    axis.ticks = element_line(color = "black"),
    
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) + guides(
    color = guide_legend(
      nrow = 1, byrow = TRUE, 
      override.aes = list( 
        shape = c(16, 17, 15, 18), alpha = 1 )
    ),
    shape = "none"
  )

#### sigma_0 = 1 ####
plot3_ratio_2 <- ggplot(
  plot3_data_2,
  aes(
    x = n,
    y = Ratio_A3_A2,
    color = factor(m),
    shape = factor(m),
    group = factor(m)
  )
) +
  # Reference lines
  geom_segment(
    data = ref_data_2_plot3,
    aes(
      x = x_min,
      xend = x_max,
      y = ref_y,
      yend = ref_y,
      color = factor(m)
    ),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Curves and points
  geom_line(
    aes(color = factor(m)),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = factor(m)),
    size = 3,
    alpha = 0.8
  ) +
  
  # x-axis
  #scale_x_log10(breaks = c(5, 6, 7, 8, 9, 10, 15, 20, 30)) +
  scale_x_log10(
    breaks = 2^seq(2, 12, by = 2),
    labels = parse(text = paste0("2^", seq(2, 12, by = 2)))
  ) +
  scale_y_log10(
    limits = c(0.08, 2.8),
    breaks = c(0.125, 0.25, 0.5, 1, 2),
    expand = c(0, 0)
  ) +
  labs(
    title = expression(sigma[0] == 1),
    x = expression("Sample size:" ~ n),
    y = expression(IMSE[FTP] / IMSE[PTF]),
    color = expression("Basis dimension:" ~ m)
  ) + scale_shape_manual(
    values = c(
      "2" = 16,  # Filled circle
      "4" = 17,  # Filled triangle
      "6" = 15,  # Filled square
      "8" = 18   # Filled diamond
    )
  ) +
  
  # Paper-style background
  theme_light(base_size = 18) +
  
  theme(
    # Panel background and boundary
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    # Gridlines
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # Axis appearance
    axis.text.x = element_text(size  = 18, color = "black"),
    axis.title.x = element_text(size  = 18, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    plot.title = element_text(
      size  = 18,
      face = "bold",
      hjust = 0.5
    ),
    
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) + guides(
    color = guide_legend(
      nrow = 1, byrow = TRUE, 
      override.aes = list( 
        shape = c(16, 17, 15, 18), alpha = 1 )
    ),
    shape = "none"
  )

#### sigma_0 = 10 ####
plot3_ratio_3 <- ggplot(
  plot3_data_3,
  aes(
    x = n,
    y = Ratio_A3_A2,
    color = factor(m),
    shape = factor(m),
    group = factor(m)
  )
) +
  # Reference lines
  geom_segment(
    data = ref_data_3_plot3,
    aes(
      x = x_min,
      xend = x_max,
      y = ref_y,
      yend = ref_y,
      color = factor(m)
    ),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Curves and points
  geom_line(
    aes(color = factor(m)),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = factor(m)),
    size = 3,
    alpha = 0.8
  ) +
  
  # x-axis
  #scale_x_log10(breaks = c(5, 6, 7, 8, 9, 10, 15, 20, 30)) +
  scale_x_log10(
    breaks = 2^seq(2, 12, by = 2),
    labels = parse(text = paste0("2^", seq(2, 12, by = 2)))
  ) +
  scale_y_log10(
    limits = c(0.08, 2.8),
    breaks = c(0.125, 0.25, 0.5, 1, 2),
    expand = c(0, 0)
  ) +
  labs(
    title = expression(sigma[0] == 10),
    x = expression("Sample size:" ~ n),
    y = expression(IMSE[FTP] / IMSE[PTF]),
    color = expression("Basis dimension:" ~ m)
  ) + scale_shape_manual(
    values = c(
      "2" = 16,  # Filled circle
      "4" = 17,  # Filled triangle
      "6" = 15,  # Filled square
      "8" = 18   # Filled diamond
    )
  ) +
  
  # Paper-style background
  theme_light(base_size = 18) +
  
  theme(
    # Panel background and boundary
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    # Gridlines
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # Axis appearance
    axis.text.x = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) + guides(
    color = guide_legend(
      nrow = 1, byrow = TRUE, 
      override.aes = list( 
        shape = c(16, 17, 15, 18), alpha = 1 )
    ),
    shape = "none"
  )

combined_plot3 <- plot3_ratio_1 + plot3_ratio_2 + plot3_ratio_3 + 
  plot_layout(ncol = 3, guides = "collect") & 
  theme(legend.position = "bottom",
        legend.box = "horizontal")

print(combined_plot3)

ggsave(filename = "figures/Poly/IMSE_Ratio_noisy.pdf", plot = combined_plot3, 
       width = 12, height = 6, units = "in", device = "pdf")

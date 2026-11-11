library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

##### RKHS Space

######## Part 1: Load Data and Preprocess

#### Data for Plot 1 ####
plot1_data_1 <- readRDS("data/RKHS/RKHS_sigma0_0_B1p1.rds")|>
  dplyr::filter(mu %in% c(1))
plot1_data_2 <- readRDS("data/RKHS/RKHS_sigma0_0_B1p1.rds")|>
  dplyr::filter(mu %in% c(0.1))

## Reference line y = 1/T for each basis dimension
ref_data_1_plot1 <- data.frame(
  T_pts = sort(unique(plot1_data_1$T_pts))
)
ref_data_1_plot1$ref_y <- 1 / ref_data_1_plot1$T_pts
ref_data_1_plot1$x_min <- min(plot1_data_1$n)
ref_data_1_plot1$x_max <- max(plot1_data_1$n)

ref_data_2_plot1 <- data.frame(
  T_pts = sort(unique(plot1_data_2$T_pts))
)
ref_data_2_plot1$ref_y <- 1 / ref_data_2_plot1$T_pts
ref_data_2_plot1$x_min <- min(plot1_data_2$n)
ref_data_2_plot1$x_max <- max(plot1_data_2$n)

#### Data for Plot 2 ####
plot2_data <- readRDS("data/RKHS/RKHS_Bprime.rds")
ref_data_plot2 <- data.frame(
  T_pts = sort(unique(plot2_data$T_pts))
)
ref_data_plot2$ref_y <- 1 / ref_data_plot2$T_pts
ref_data_plot2$x_min <- min(plot2_data$n)
ref_data_plot2$x_max <- max(plot2_data$n)

#### Data for Plot 3 ####
plot3_data <- readRDS("data/RKHS/RKHS_noisy.rds")
plot3_data_1 <- subset(plot3_data, plot3_data$sigma_0 == 0.1)
plot3_data_2 <- subset(plot3_data, plot3_data$sigma_0 == 1)
plot3_data_3 <- subset(plot3_data, plot3_data$sigma_0 == 10)

## Reference line y = 1/T for each basis dimension
ref_data_1_plot3 <- data.frame(
  T_pts = sort(unique(plot3_data_1$T_pts))
)
ref_data_1_plot3$ref_y <- 1 / ref_data_1_plot3$T_pts
ref_data_1_plot3$x_min <- min(plot3_data_1$n)
ref_data_1_plot3$x_max <- max(plot3_data_1$n)

ref_data_2_plot3 <- data.frame(
  T_pts = sort(unique(plot3_data_2$T_pts))
)
ref_data_2_plot3$ref_y <- 1 / ref_data_2_plot3$T_pts
ref_data_2_plot3$x_min <- min(plot3_data_2$n)
ref_data_2_plot3$x_max <- max(plot3_data_2$n)

ref_data_3_plot3 <- data.frame(
  T_pts = sort(unique(plot3_data_3$T_pts))
)
ref_data_3_plot3$ref_y <- 1 / ref_data_3_plot3$T_pts
ref_data_3_plot3$x_min <- min(plot3_data_3$n)
ref_data_3_plot3$x_max <- max(plot3_data_3$n)

###############################################################################

######## Figure 1: IMSE Ratio of FTP to PTF vs n,T (sigma_0 = 0, mu = 0.1, 1)
#### mu = 1 ####
plot1_ratio_1 <- ggplot(
  plot1_data_1,
  aes(
    x = n,
    y = Ratio_A3_A2,
    color = factor(T_pts),
    shape = factor(T_pts),
    group = factor(T_pts)
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
      color = factor(T_pts)
    ),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Curves and points
  geom_line(
    aes(color = factor(T_pts)),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = factor(T_pts)),
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
    limits = c(0.08, 100),
    breaks = c(0.1, 1, 10, 100),
    expand = c(0, 0)
  ) +
  labs(
    title = expression(mu == 1),
    x = expression("Sample size:" ~ n),
    y = expression(MSE[FTP] / MSE[PTF]),
    color = expression("Grid Size:" ~ T)
  ) + scale_shape_manual(
    values = c(
      "2" = 16,  # Filled circle
      "4" = 17,  # Filled triangle
      "6" = 15,  # Filled square
      "8" = 18   # Filled diamond
    )
  ) +
  theme_light(base_size = 18)+
  # Paper-style background
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
    color = factor(T_pts),
    shape = factor(T_pts),
    group = factor(T_pts)
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
      color = factor(T_pts)
    ),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Curves and points
  geom_line(
    aes(color = factor(T_pts)),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = factor(T_pts)),
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
    limits = c(0.08, 100),
    breaks = c(0.1, 1, 10, 100),
    expand = c(0, 0)
  ) +
  labs(
    title = expression(mu == 0.1),
    x = expression("Sample size:" ~ n),
    y = expression(MSE[FTP] / MSE[PTF]),
    color = expression("Grid Size:" ~ T)
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
  theme_light(base_size = 18)+
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

ggsave(filename = "figures/RKHS/MSE_Ratio_mu_e-1-1.pdf", plot = combined_plot1, 
       width = 12, height = 6, units = "in", device = "pdf")


###############################################################################
#### Table display ####
raw <- readRDS("data/RKHS/RKHS_sigma0_0_B1p1.rds")

table_data <- raw |>
  filter(
    n %in% 2^seq(2, 12, by = 2),
    mu %in% c(0.1, 1)
  ) |>
  transmute(
    T = T_pts,
    n,
    mu,
    MSE_FTP = IMSE_A3_FTP_constrained,
    MSE_PTF = IMSE_A2_PTF_constrained,
    Ratio = Ratio_A3_A2,
    Prob_PTF = A2_prop_constrained
  )

# Each (T, n, mu) should identify exactly one result.
stopifnot(
  !anyDuplicated(table_data[c("T", "n", "mu")])
)

table_wide <- table_data |>
  pivot_wider(
    id_cols = c(T, n),
    names_from = mu,
    values_from = c(MSE_FTP, MSE_PTF, Ratio, Prob_PTF),
    names_glue = "{.value}_mu_{mu}",
    names_sort = TRUE,
    names_vary = "fastest"
  ) |>
  arrange(T, n)

# Full-precision numeric values
View(table_wide)
print(table_wide, n = Inf, width = Inf)

table_display <- table_wide |>
  mutate(
    n = paste0("2^", as.integer(log2(n))),
    across(starts_with("MSE_"), ~ sprintf("%.2e", .x)),
    across(
      starts_with("Ratio_") | starts_with("Prob_PTF_"),
      ~ sprintf("%.3f", .x)
    )
  )

View(table_display)
print(table_display, n = Inf, width = Inf)

######## Figure 2: IMSE of FTP and PTF vs n, Bprime (mu = 1, sigma_0 = 1)
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
    limits = c(0.8, 350),
    breaks = c(1, 10, 100),
    expand = c(0, 0)
  ) +
  
  # Labels
  labs(
    #title = expression(n %*% MSE[FTP]),
    x = expression("Sample size:" ~ n),
    y = expression(n %*% MSE[FTP]),
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
    limits = c(0.8, 350),
    breaks = c(1, 10, 100),
    expand = c(0, 0)
  ) +
  
  # Labels
  labs(
    #title = expression(n %*% MSE[PTF]),
    x = expression("Sample size:" ~ n),
    y = expression(n %*% MSE[PTF]),
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

combined_plot2 <- plot2_IMSE_FTP + plot2_IMSE_PTF +
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom",
        legend.box = "horizontal")

print(combined_plot2)

ggsave(filename = "figures/RKHS/RKHS_Bprime_choices.pdf", plot = combined_plot2, 
       width = 12, height = 6, units = "in", device = "pdf")

###############################################################################

######## Figure 3: IMSE Ratio of FTP to PTF vs n,m (mu = 1, sigma_0 = 0.01, 1, 10)
#### sigma_0 = 0.1 ####
plot3_ratio_1 <- ggplot(
  plot3_data_1,
  aes(
    x = n,
    y = Ratio_A3_A2,
    color = factor(T_pts),
    shape = factor(T_pts),
    group = factor(T_pts)
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
      color = factor(T_pts)
    ),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Curves and points
  geom_line(
    aes(color = factor(T_pts)),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = factor(T_pts)),
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
    limits = c(0.1, 2.5),
    breaks = c(0.125, 0.25, 0.5, 1, 2),
    expand = c(0, 0)
  ) +
  labs(
    title = expression(sigma[0] == 0.1),
    x = expression("Sample size:" ~ n),
    y = expression(MSE[FTP] / MSE[PTF]),
    color = expression("Grid Size:" ~ T)
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
  theme_light(base_size = 18)+
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
    color = factor(T_pts),
    shape = factor(T_pts),
    group = factor(T_pts)
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
      color = factor(T_pts)
    ),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Curves and points
  geom_line(
    aes(color = factor(T_pts)),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = factor(T_pts)),
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
    limits = c(0.1, 2.5),
    breaks = c(0.125, 0.25, 0.5, 1, 2),
    expand = c(0, 0)
  ) +
  labs(
    title = expression(sigma[0] == 1),
    x = expression("Sample size:" ~ n),
    y = expression(MSE[FTP] / MSE[PTF]),
    color = expression("Grid Size:" ~ T)
  ) +
  scale_shape_manual(
    values = c(
      "2" = 16,  # Filled circle
      "4" = 17,  # Filled triangle
      "6" = 15,  # Filled square
      "8" = 18   # Filled diamond
    )
  ) +
  theme_light(base_size = 18)+
  # Paper-style background
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
    color = factor(T_pts),
    shape = factor(T_pts),
    group = factor(T_pts)
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
      color = factor(T_pts)
    ),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Curves and points
  geom_line(
    aes(color = factor(T_pts)),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = factor(T_pts)),
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
    limits = c(0.1, 2.5),
    breaks = c(0.125, 0.25, 0.5, 1, 2),
    expand = c(0, 0)
  ) +
  labs(
    title = expression(sigma[0] == 10),
    x = expression("Sample size:" ~ n),
    y = expression(MSE[FTP] / MSE[PTF]),
    color = expression("Grid Size:" ~ T)
  ) +
  scale_shape_manual(
    values = c(
      "2" = 16,  # Filled circle
      "4" = 17,  # Filled triangle
      "6" = 15,  # Filled square
      "8" = 18   # Filled diamond
    )
  ) +
  theme_light(base_size = 18)+
  # Paper-style background
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

ggsave(filename = "figures/RKHS/MSE_Ratio_noisy.pdf", plot = combined_plot3, 
       width = 12, height = 6, units = "in", device = "pdf")

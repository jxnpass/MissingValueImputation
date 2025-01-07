source("sim_functions.R")
full_results <- readRDS("simulation_results.RDS")

add_theme <- theme_light() + 
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = 'bold'),
    axis.title.y = element_text(size = 14, face = 'bold'),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

### See individual results
# Mean graphic

ggplot(data = full_results$means %>% 
         pivot_longer(
           cols = c(-Method, -type, -rho),
           names_to = "Statistic",
           values_to = "Mean"
         ), mapping = aes(x = fct_inorder(Method), y = Mean, fill = Method)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_violin(show.legend = FALSE) +
  facet_grid(
    cols = vars(fct_inorder(type)), 
    rows = vars(rho), 
    labeller = labeller(
      rho = function(x) paste("\u03c1 =", x)
    )
  ) +
  labs(
    title = "Estimation of Means",
    subtitle = paste("\u03bc = 0"),
    caption = "Each density curve holds 1000 x 3 estimations",
    x = "Method", 
    y = "Difference"
  ) +
  add_theme +
  scale_fill_manual(
    values = c("LWD" = "red3", 
               "EM" = "green3",
               "CS" = "skyblue3")
  )

# variance
ggplot(data = full_results$cov %>% 
         pivot_longer(
           cols = c(-Method, -type, -rho),
           names_to = "Statistic",
           values_to = "Mean"
         ) %>% 
         filter(Statistic %in% c("X1.X1", "X2.X2", "X3.X3")), 
       mapping = aes(x = fct_inorder(Method), y = Mean, fill = Method)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_violin(show.legend = FALSE) +
  facet_grid(
    cols = vars(fct_inorder(type)), 
    rows = vars(rho), 
    labeller = labeller(
      rho = function(x) paste("\u03c1 =", x)
    )
  ) +
  labs(
    title = "Estimation of Variances",
    subtitle = paste("\u03c3_ii = 1"),
    caption = "Each density curve holds 1000 x 3 estimations",
    x = "Method", 
    y = "Difference"
  ) +
  add_theme +
  scale_fill_manual(
    values = c("LWD" = "red3", 
               "EM" = "green3",
               "CS" = "skyblue3")
  )

# covariance
ggplot(data = full_results$cov %>% 
         pivot_longer(
           cols = c(-Method, -type, -rho),
           names_to = "Statistic",
           values_to = "Mean"
         ) %>% 
         filter(!(Statistic %in% c("X1.X1", "X2.X2", "X3.X3"))), 
       mapping = aes(x = fct_inorder(Method), y = Mean, fill = Method)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_violin(show.legend = FALSE) +
  facet_grid(
    cols = vars(fct_inorder(type)), 
    rows = vars(rho), 
    labeller = labeller(
      rho = function(x) paste("\u03c1 =", x)
    )
  ) +
  labs(
    title = "Estimation of Covariances",
    subtitle = paste("\u03c3_ij = \u03c3_ii\u03c1"),
    caption = "Each density curve holds 1000 x 3 estimations",
    x = "Method", 
    y = "Difference"
  ) +
  add_theme +
  scale_fill_manual(
    values = c("LWD" = "red3", 
               "EM" = "green3",
               "CS" = "skyblue3")
  )

# corr
ggplot(data = full_results$cor %>% 
         pivot_longer(
           cols = c(-Method, -type, -rho),
           names_to = "Statistic",
           values_to = "Mean"
         ), 
       mapping = aes(x = fct_inorder(Method), y = Mean, fill = Method)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_violin(show.legend = FALSE) +
  facet_grid(
    cols = vars(fct_inorder(type)), 
    rows = vars(rho), 
    labeller = labeller(
      rho = function(x) paste("\u03c1 =", x)
    )
  ) +
  labs(
    title = "Estimation of Correlations",
    subtitle = paste("\u03c1 = 0, 0.3, 0.7"),
    caption = "Each density curve holds 1000 x 3 estimations",
    x = "Method", 
    y = "Difference"
  ) +
  add_theme +
  scale_fill_manual(
    values = c("LWD" = "red3", 
               "EM" = "green3",
               "CS" = "skyblue3")
  )

# MCAR test distribution
ggplot(data = full_results$tests, 
       mapping = aes(x = MCARstat)) +
  geom_density() +
  facet_grid(scales = "free_x",
    cols = vars(fct_inorder(type)), 
    rows = vars(rho), 
    labeller = labeller(
      rho = function(x) paste("\u03c1 =", x)
    )
  ) +
  labs(
    title = "Little's MCAR Test Statistics",
    subtitle = bquote("Deriving from " * chi^2 * " Distribution"),
    caption = "A Test Statistic Greater Than 16.83 Indicates Significance at p < .05",
    x = "Test Statistic", 
    y = "Density"
  ) +
  add_theme +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())


# Table 1

calculate_statistics(full_results$means) %>% 
  mutate(Method = factor(Method, c("LWD", "EM", "CS")),
         type = factor(type, c("MCAR", "MAR", "MNAR"))) %>% 
  arrange(rho, type, Method) %>% 
  mutate_if(is.numeric, ~ round(.x, digits = 4)) %>% View()

# Table 2
full_results$cov %>% 
  pivot_longer(
    cols = c(-Method, -type, -rho),
    names_to = "Statistic",
    values_to = "Mean"
  ) %>% 
  mutate(VarStat = ifelse(Statistic %in% c("X1.X1", "X2.X2", "X3.X3"), "Var", "Cov")) %>% 
  group_by(rho, type, Method, VarStat) %>% 
  summarise(Bias = mean(Mean),
            Var = var(Mean)) %>% 
  arrange(rho, type, Method) %>% 
  mutate_if(is.numeric, ~ round(.x, digits = 4)) %>% 
  View()

full_results$cor %>% 
  pivot_longer(
    cols = c(-Method, -type, -rho),
    names_to = "Statistic",
    values_to = "Mean"
  ) %>% 
  group_by(rho, type, Method) %>% 
  summarise(Bias = mean(Mean),
            Var = var(Mean)) %>% 
  arrange(rho, type, Method) %>% 
  mutate_if(is.numeric, ~ round(.x, digits = 4)) %>% 
  View()






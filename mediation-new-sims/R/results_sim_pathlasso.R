library(tidyverse)
library(qs)
full_results <- qread("results/pathlasso_results.qs")
method_colors <- c("magenta3", "steelblue1")
method_lvls <- c("Proposed", "Pathway Lasso")
names(method_colors) <- method_lvls

sel_results <- full_results %>% unnest(data) %>%
  select(-ests, -method_summary) %>%
  unnest(sel_result) %>%
  mutate(method = if_else(method=="Ours", "Proposed", method),
         method=factor(method, method_lvls),
         sel_type = factor(sel_idx, levels=c(1:7),
                           labels=c(
                             "Mediator 1",
                             "Mediator 2",
                             "Outcome Predictor",
                             "Effected by Treatment",
                             "Outcome Predictor",
                             "Outcome Predictor",
                             "Noise"
                           )),
         # n = factor(n, labels=sprintf("n=%i", c(2,4)*100))
         )

p1 <- sel_results %>%
  ggplot(aes(x=sel_type, y=Shat)) +
  geom_col(aes(fill=method), position=position_dodge()) +
  theme_bw() +
  scale_y_continuous(limits=c(0,1), labels = ~sprintf("%i%%", .x*100)) +
  scale_fill_manual(values=method_colors) +
  # facet_grid(cols=vars(n)) +
  labs(y="Selection Probability",
       fill="\u03ba",
       x="Variable Type") +
  theme(axis.text.x = element_text(hjust = 1, angle=15),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank())
p1

bias_table <-
  full_results %>% unnest(data) %>%
  unnest(ests) %>%
  mutate(NDE_star = 2, bias=abs(NDE_hat_mean - NDE_star)) %>%
  mutate(lam_type = if_else(method != "Ours", NA, lam_type)) %>%
  select(n:rho, method, lam_type,
         NDE_star, NDE_hat_mean, bias, SE=NDE_hat_sd) %>%
  mutate(lower95 = bias-qnorm(0.975, sd=SE),
         upper95=bias+qnorm(0.975, sd=SE))

bias_table1 <- bias_table


bias_table <-
  full_results %>% unnest(data) %>%
  unnest(ests) %>%
  mutate(NIE_star = 4+(n**(alpha_rate+beta_rate)), bias=abs(NIE_hat_mean - NIE_star)) %>%
  mutate(lam_type = if_else(method != "Ours", NA, lam_type)) %>%
  select(n:rho, method, lam_type,
         NIE_star, NIE_hat_mean, bias, SE=NIE_hat_sd) %>%
  mutate(lower95 = bias-qnorm(0.975, sd=SE),
         upper95=bias+qnorm(0.975, sd=SE))
bias_table2 <- bias_table

bias_table1$metric <- "NDE"
bias_table2$metric <- "NIE"

p2 <- rbind(bias_table1 %>% select(method, metric, bias, lower95, upper95),
      bias_table2 %>% select(method, metric, bias, lower95, upper95)) %>%
  mutate(method = if_else(method=="Ours", "Proposed", method),
         method = factor(method, method_lvls)) %>%
  ggplot(aes(x=method, y=bias)) +
  geom_col(aes(fill=method, group=method)) +
  theme_bw() +
  scale_fill_manual(values=method_colors) +
  facet_wrap(vars(metric)) +
  labs(y="Bias", x="") +
  theme(axis.text.x = element_text(hjust = 1, angle=15),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

cowplot::plot_grid(p1, p2, nrow=2, labels=LETTERS[1:2])

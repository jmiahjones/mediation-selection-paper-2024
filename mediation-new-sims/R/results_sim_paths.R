library(qs)
library(tidyverse)
full_results <- qread("results/full_results.qs")

plot_sel <- full_results %>% unnest(data) %>%
  select(-ests) %>%
  mutate(sel_plot =
           map(sel_result, \(x)
               ggplot(x, aes(x=sel_idx, y=Shat)) + geom_point() +
                 theme_bw() +
                 scale_y_continuous(limits=c(0,1))
           )
  ) %>%
  arrange(n, p, m1_size, use_lm)

library(cowplot)
plot_sel %>%
  mutate(name=paste(n, p, m1_size, use_lm, sep="-")) %>%
  pull(name) %>%
  cowplot::plot_grid(plotlist=plot_sel$sel_plot, labels = .,
                     label_y = 1.02,
                     # label_x = -0.5,
                     label_size=10,
                     nrow = 3, ncol=4
                     )

library(kableExtra)
full_results %>% unnest(data) %>%
  select(-ests) %>%
  unnest(sel_result) %>%
  group_by(n,p,m1_size,sel_idx) %>%
  summarize(Shat = mean(Shat)) %>% ungroup %>%
  mutate(sel_idx=paste0("M", sel_idx)) %>%
  tidyr::pivot_wider(id_cols=c("n","p","m1_size"),
                     names_from="sel_idx", values_from="Shat") %>%
  kbl(format="latex", digits = 3)

full_results %>% unnest(data) %>%
  select(-sel_result) %>% unnest(ests) %>%
  mutate(NIE_star = m1_size*2.8,
         bias=abs(NIE_hat_mean - NIE_star),
         rel_bias = bias/NIE_star,
         p=factor(p)) %>%
  select(n:m1_size, NIE_star, NIE_hat_mean, bias, rel_bias) %>%
  ggplot(aes(x=m1_size, y=rel_bias)) +
  geom_line(aes(color=p, group=p)) +
  theme_bw() +
  labs(y="", x="\u03b1")

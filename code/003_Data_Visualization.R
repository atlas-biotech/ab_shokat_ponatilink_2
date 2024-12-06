# 003 Validation and Visualization #
library(dplyr)
library(ggplot2)
library(ggpubr)

df = read.csv('processed_data/pl2_master_datafile.cvs')
print(paste0('Baseline coverage: ',sum(!is.na(df$freq.bl)),' out of 8405, ',100*round(sum(!is.na(df$freq.bl))/8405,3),'%'))

df.view = df %>%
  select(mut,freq.bl,l2fc.a.1000,l2fc.pl2.100,l2fc.p.30)

df.view.strict = df.view %>% 
  filter(freq.bl >= 1/10000)

df.view.mid = df.view %>% 
  filter(freq.bl >= 1/100000)


resmuts = c('T315I','T315L','T315Q','T315M','T315E','T315A','F311I','F317V',
            'F317I','F317C','V299L','Y253H','E255K','E255V','G250E','Q252H',
            'F359C','F359I','F359V','H396R','A344P','A337T','E355K','E459K',
            'P465S','V468F','M351T','M351V','I502M','I502L')

df.resmuts = df.view[df.view$mut %in% resmuts,]

figure_aes = theme_classic2() +
  theme(title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

asc.resmuts = df.view.strict %>% 
  filter(l2fc.a.1000 >= 1)

# Asciminib Resmuts #

ggplot(asc.resmuts, aes(l2fc.a.1000, l2fc.pl2.100)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  labs(title = paste0("Asc vs PL2 L2FC"),
       x = "Asciminib 1000 nM",
       y = "PonatiLink-2 100 nM") +
  ylim(c(-8,5)) +
  xlim(c(-8,5)) +
  figure_aes

ggplot(asc.resmuts, aes(l2fc.a.1000, l2fc.p.30)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  labs(title = paste0("Asc vs Pon L2FC"),
       x = "Asciminib 1000 nM",
       y = "Ponatinib 30 nM") +
  ylim(c(-8,5)) +
  xlim(c(-8,5)) +
  figure_aes

# Ponatinib Resmuts #



# Other plots #

ggplot(df, aes(l2fc.a.1000.a, l2fc.a.1000.b)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  labs(title = paste0("Asciminib 1000 nM L2FC"),
       x = "Rep A",
       y = "Rep B") +
  figure_aes

ggplot(df, aes(l2fc.pl2.100.a, l2fc.pl2.100.b)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  stat_cor(method = "pearson", label.x = 3, label.y = -8) +
  labs(title = paste0("PonatiLink-2 100 nM L2FC"),
       x = "Rep A",
       y = "Rep B") +
  figure_aes

ggplot(df, aes(l2fc.p.30.a, l2fc.p.30.b)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  stat_cor(method = "pearson", label.x = 3, label.y = -8) +
  labs(title = paste0("Ponatinib 30 nM L2FC"),
       x = "Rep A",
       y = "Rep B") +
  figure_aes

ggplot(df.shokat, aes(freq.bl.a, freq.bl.b)) +
  geom_point(fill = 'grey60', alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  stat_cor(method = "pearson") +
  # # xlim(c(1e-06,1e-02)) +
  # ylim(c(1e-06,1e-02)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = paste0("ABL Library Baseline Frequency"),
       x = "Rep A",
       y = "Rep B") +
  figure_aes

# BASELINE CORRELATION

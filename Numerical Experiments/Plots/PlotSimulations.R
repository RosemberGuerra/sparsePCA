## Plotting the results 2 ##
# Rosember Guerra #
# update: 10-03-2020

rm(list = ls(all.names = TRUE))

library(data.table)
library(ggplot2)
library(ggpubr)

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

load("JoinDataSet2.RData")
class(JoinDataSet2)

# Dimension for each plot: width=1200, hight=858 #
Width = 12.50
Height = 8.94

### Squared relative error - loadings/weights and scores ###

# Experiment I:

plot_I <- ggplot(data = JoinDataSet2[p_sparse != "Sparsity 0%" & n_components == "K = 2" &
                                       VAFx =="vafx = 80%"],
               aes(x = Methodology, y = ExpI))
plot_I + geom_boxplot(aes(fill = n_variables), outlier.size=0.3, lwd=0.2) +
  facet_grid(rows = vars(Measurements), cols = vars(s_size,p_sparse)) +
  # geom_hline(yintercept = 0.85, col = "black", linetype = 2) +
  geom_vline(xintercept = 3.5, col = "black", linetype = 2) +
  guides(fill=guide_legend(title="N. variables")) +
  xlab("") + 
  ylim(0,1)+
  scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
  ylab("Condition type I") +
  coord_fixed(ratio = 6/1) +
  theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=8),
                              axis.text.y = element_text(size=10),
                              strip.text = element_text(size=10),
                              legend.text = element_text(size=10),
                              legend.title = element_text(size=10),
                              legend.key.size = unit(0.5, "cm"))

ggsave("ExpI-Measurements.eps",path = "plots",width = Width, height = Height,
       limitsize = FALSE)

# experiment II #
plot_II <- ggplot(data = JoinDataSet2[ n_components == "K = 2"&
                                         VAFx =="vafx = 80%"],
               aes(x = Methodology, y = ExpII))
plot_II + geom_boxplot(aes(fill = n_variables), outlier.size=0.3, lwd=0.2) +
  facet_grid(rows = vars(Measurements), cols = vars(s_size,p_sparse2)) +
  # geom_hline(yintercept = 0.85, col = "black", linetype = 2) +
  geom_vline(xintercept = 3.5, col = "black", linetype = 2) +
  guides(fill=guide_legend(title="N. variables")) +
  xlab("") + 
  ylim(0,1)+
  scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
  ylab("Condition type II") +
  coord_fixed(ratio = 8/1) +
  theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=8),
                              axis.text.y = element_text(size=10),
                              strip.text = element_text(size=10),
                              legend.text = element_text(size=10),
                              legend.title = element_text(size=10),
                              legend.key.size = unit(0.5, "cm"))

ggsave("ExpII-Measurements.eps",path = "plots",width = Width, height = Height,
       limitsize = FALSE)

# Experiment III #
plot_III <- ggplot(data = JoinDataSet2[p_sparse != "Sparsity 0%" & n_components == "K = 2" &
                                       VAFx =="vafx = 80%"],
                 aes(x = Methodology, y = ExpIII))
plot_III + geom_boxplot(aes(fill = n_variables), outlier.size=0.3, lwd=0.2) +
  facet_grid(rows = vars(Measurements2), cols = vars(s_size,p_sparse)) +
  # geom_hline(yintercept = 0.85, col = "black", linetype = 2) +
  geom_vline(xintercept = 3.5, col = "black", linetype = 2) +
  guides(fill=guide_legend(title="N. variables")) +
  xlab("") + 
  ylim(0,1)+
  scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
  ylab("Condition type III") +
  coord_fixed(ratio = 6/1) +
  theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=8),
                              axis.text.y = element_text(size=10),
                              strip.text = element_text(size=10),
                              legend.text = element_text(size=10),
                              legend.title = element_text(size=10),
                              legend.key.size = unit(0.5, "cm"))

ggsave("ExpIII-Measurements.eps",path = "plots",width = Width, height = Height,
       limitsize = FALSE)

### --- Appendix plots ---  ###
# Experiment I, K=3 #

plot_I <- ggplot(data = JoinDataSet2[p_sparse != "Sparsity 0%" & n_components == "K = 3" &
                                       VAFx =="vafx = 80%"],
                 aes(x = Methodology, y = ExpI))
plot_I + geom_boxplot(aes(fill = n_variables), outlier.size=0.3, lwd=0.2) +
  facet_grid(rows = vars(Measurements), cols = vars(s_size,p_sparse)) +
  # geom_hline(yintercept = 0.85, col = "black", linetype = 2) +
  geom_vline(xintercept = 3.5, col = "black", linetype = 2) +
  guides(fill=guide_legend(title="N. variables")) +
  xlab("") + 
  ylim(0,1)+
  scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
  ylab("Condition type I") +
  coord_fixed(ratio = 6/1) +
  theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=8),
                              axis.text.y = element_text(size=10),
                              strip.text = element_text(size=10),
                              legend.text = element_text(size=10),
                              legend.title = element_text(size=10),
                              legend.key.size = unit(0.5, "cm"))

ggsave("AppxExpI-Measurements.eps",path = "plots",width = Width, height = Height,
       limitsize = FALSE)

# experiment II, K=3 #
plot_II <- ggplot(data = JoinDataSet2[ n_components == "K = 3" &
                                         VAFx =="vafx = 80%"],
                  aes(x = Methodology, y = ExpII))
plot_II + geom_boxplot(aes(fill = n_variables), outlier.size=0.3, lwd=0.2) +
  facet_grid(rows = vars(Measurements), cols = vars(s_size,p_sparse2)) +
  # geom_hline(yintercept = 0.85, col = "black", linetype = 2) +
  geom_vline(xintercept = 3.5, col = "black", linetype = 2) +
  guides(fill=guide_legend(title="N. variables")) +
  xlab("") + 
  ylim(0,1)+
  scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
  ylab("Condition type II") +
  coord_fixed(ratio = 8/1) +
  theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=8),
                              axis.text.y = element_text(size=10),
                              strip.text = element_text(size=10),
                              legend.text = element_text(size=10),
                              legend.title = element_text(size=10),
                              legend.key.size = unit(0.5, "cm"))

ggsave("AppxExpII-Measurements.eps",path = "plots",width = Width, height = Height,
       limitsize = FALSE)

# Experiment III, K=3 #
plot_III <- ggplot(data = JoinDataSet2[p_sparse != "Sparsity 0%" & n_components == "K = 3" &
                                         VAFx =="vafx = 80%"],
                   aes(x = Methodology, y = ExpIII))
plot_III + geom_boxplot(aes(fill = n_variables), outlier.size=0.3, lwd=0.2) +
  facet_grid(rows = vars(Measurements2), cols = vars(s_size,p_sparse)) +
  # geom_hline(yintercept = 0.85, col = "black", linetype = 2) +
  geom_vline(xintercept = 3.5, col = "black", linetype = 2) +
  guides(fill=guide_legend(title="N. variables")) +
  xlab("") + 
  ylim(0,1)+
  scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
  ylab("Condition type III") +
  coord_fixed(ratio = 6/1) +
  theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=8),
                              axis.text.y = element_text(size=10),
                              strip.text = element_text(size=10),
                              legend.text = element_text(size=10),
                              legend.title = element_text(size=10),
                              legend.key.size = unit(0.5, "cm"))

ggsave("AppxExpIII-Measurements.eps",path = "plots",width = Width, height = Height,
       limitsize = FALSE)

# Experiment I and III, with PS =  0% #
plot <- ggplot(data = JoinDataSet2[p_sparse == "Sparsity 0%" &  VAFx =="vafx = 80%" &
                                     Measurements != "MR"],
               aes(x = Methodology, y = ExpI))
plot + geom_boxplot(aes(fill = n_variables), outlier.size=0.3, lwd=0.2) +
  facet_grid(rows = vars(Measurements), cols = vars(s_size,n_components)) +
  # geom_hline(yintercept = 0.85, col = "black", linetype = 2) +
  geom_vline(xintercept = 3.5, col = "black", linetype = 2) +
  guides(fill=guide_legend(title="N. variables")) +
  xlab("") + 
  ylim(0,1)+
  scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
  ylab("SRE - Loadings & Weights") +
  coord_fixed(ratio = 6/1) +
  theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=8),
                              axis.text.y = element_text(size=10),
                              strip.text = element_text(size=10),
                              legend.text = element_text(size=10),
                              legend.title = element_text(size=10),
                              legend.key.size = unit(0.5, "cm"))
ggsave("AppxExpI0-Measurements.eps",path = "plots",width = Width, height = Height,
       limitsize = FALSE)

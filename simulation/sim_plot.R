# sim plot
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggpubr)

### bootstrap parametric integral dr
## tau=4
# case1
load("survival/scenario 1/true_survival_1.RData")
load("survival/scenario 1/sim1-n4000.RData")
f.sim.case1 <- f.sim.result$tau4
f.sim.case1 <- f.sim.case1[grep("direct.mc.ha", rownames(f.sim.case1), value = TRUE), ]
rm(f.sim.result)

true.value.case1.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("survival/scenario 1/sim1-baseline-n4000.RData")
# f.sim.result <- f.sim.result[!(rownames(f.sim.result) %in% c("dr.integral.driptw", "dr.integral.driptw.trunc", "dr.integral.drow", "dr.integral.drmw", "dr.integral.drenw")), ]
f.sim.case1 <- rbind(f.sim.case1, f.sim.result)
rm(f.sim.result)

# case2
load("survival/scenario 2/true_survival_2.RData")
load("survival/scenario 2/sim2-n4000.RData")
f.sim.case2 <- f.sim.result$tau4
f.sim.case2 <- f.sim.case2[grep("direct.mc.ha", rownames(f.sim.case2), value = TRUE), ]
rm(f.sim.result)

true.value.case2.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("survival/scenario 2/sim2-baseline-n4000.RData")
# f.sim.result <- f.sim.result[!(rownames(f.sim.result) %in% c("dr.integral.driptw", "dr.integral.driptw.trunc", "dr.integral.drow", "dr.integral.drmw", "dr.integral.drenw")), ]
f.sim.case2 <- rbind(f.sim.case2, f.sim.result)
rm(f.sim.result)

# case3
load("competing/scenario 1/true_competing_7.RData")
load("competing/scenario 1/sim7-n4000.RData")
f.sim.case3 <- f.sim.result$tau4
f.sim.case3 <- f.sim.case3[grep("direct.mc.ha", rownames(f.sim.case3), value = TRUE), ]
rm(f.sim.result)

true.value.case3.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("competing/scenario 1/sim7-baseline-n4000.RData")
# f.sim.result <- f.sim.result[!(rownames(f.sim.result) %in% c("dr.integral.driptw", "dr.integral.driptw.trunc", "dr.integral.drow", "dr.integral.drmw", "dr.integral.drenw")), ]
f.sim.case3 <- rbind(f.sim.case3, f.sim.result)
rm(f.sim.result)

# case4
load("competing/scenario 2/true_competing_8.RData")
load("competing/scenario 2/sim8-n4000.RData")
f.sim.case4 <- f.sim.result$tau4
f.sim.case4 <- f.sim.case4[grep("direct.mc.ha", rownames(f.sim.case4), value = TRUE), ]
rm(f.sim.result)

true.value.case4.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("competing/scenario 2/sim8-baseline-n4000.RData")
# f.sim.result <- f.sim.result[!(rownames(f.sim.result) %in% c("dr.integral.driptw", "dr.integral.driptw.trunc", "dr.integral.drow", "dr.integral.drmw", "dr.integral.drenw")), ]
f.sim.case4 <- rbind(f.sim.case4, f.sim.result)
rm(f.sim.result)

#
nsim <- 1000
psiLj.case1.tau4 <- psiLj.case2.tau4 <- psiLj.case3.tau4 <- psiLj.case4.tau4 <- data.frame(matrix(NA, nrow=19, ncol=5))
colnames(psiLj.case1.tau4) <- colnames(psiLj.case2.tau4) <- colnames(psiLj.case3.tau4) <- colnames(psiLj.case4.tau4) <- c("true","effect","ese","ase","cp")
rownames(psiLj.case1.tau4) <- rownames(psiLj.case2.tau4) <-
  c("RMST.dml", "RMST.dml.t", "RMST.dml.ow", "RMST.dml.mw", "RMST.dml.enw",
    "S.or", "S.or.ow", "S.or.mw", "S.or.enw", "S.ipcw", "S.ipcw.t", "S.ipcw.ow", "S.ipcw.mw", "S.ipcw.enw", "S.dr", "S.dr.t", "S.dr.ow", "S.dr.mw", "S.dr.enw")
# rownames(psiLj.case1.tau4) <- rownames(psiLj.case2.tau4) <-
#   c("RMST.driptw", "RMST.driptw.t", "RMST.drow", "RMST.drmw", "RMST.drenw",
#     "S.or", "S.or.ow", "S.or.mw", "S.or.enw", "S.ipcw.iptw", "S.ipcw.iptw.t", "S.ipcw.ow", "S.ipcw.mw", "S.ipcw.enw", "S.driptw", "S.driptw.t", "S.drow", "S.drmw", "S.drenw")

rownames(psiLj.case3.tau4) <- rownames(psiLj.case4.tau4) <-
  c("RMTLj.dml", "RMTLj.dml.t", "RMTLj.dml.ow", "RMTLj.dml.mw", "RMTLj.dml.enw",
    "Fj.or", "Fj.or.ow", "Fj.or.mw", "Fj.or.enw", "Fj.ipcw", "Fj.ipcw.t", "Fj.ipcw.ow", "Fj.ipcw.mw", "Fj.ipcw.enw", "Fj.dr", "Fj.dr.t", "Fj.dr.ow", "Fj.dr.mw", "Fj.dr.enw")
# rownames(psiLj.case3.tau4) <- rownames(psiLj.case4.tau4) <-
#   c("RMTLj.driptw", "RMTLj.driptw.t", "RMTLj.drow", "RMTLj.drmw", "RMTLj.drenw",
#     "Fj.or", "Fj.or.ow", "Fj.or.mw", "Fj.or.enw", "Fj.ipcw.iptw", "Fj.ipcw.iptw.t", "Fj.ipcw.ow", "Fj.ipcw.mw", "Fj.ipcw.enw", "Fj.driptw", "Fj.driptw.t", "Fj.drow", "Fj.drmw", "Fj.drenw")

group.survival <- c("S.or", "S.ipcw", "S.dr", "RMST.dml",
                    "S.ipcw.t", "S.dr.t", "RMST.dml.t",
                    "S.or.ow", "S.ipcw.ow", "S.dr.ow", "RMST.dml.ow",
                    "S.or.mw", "S.ipcw.mw", "S.dr.mw", "RMST.dml.mw",
                    "S.or.enw", "S.ipcw.enw", "S.dr.enw", "RMST.dml.enw")

group.competing <- c("Fj.or", "Fj.ipcw", "Fj.dr", "RMTLj.dml",
                     "Fj.ipcw.t", "Fj.dr.t", "RMTLj.dml.t",
                     "Fj.or.ow", "Fj.ipcw.ow", "Fj.dr.ow", "RMTLj.dml.ow",
                     "Fj.or.mw", "Fj.ipcw.mw", "Fj.dr.mw", "RMTLj.dml.mw",
                     "Fj.or.enw", "Fj.ipcw.enw", "Fj.dr.enw", "RMTLj.dml.enw")
## ate
# case1
psiLj.case1.tau4$true <- as.numeric(true.value.case1.tau4[, c(c("rmst.diff", "rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff"), c("rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff"), rep(c("rmst.diff", "rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff"), 2))])
psiLj.case1.tau4$effect <- apply(f.sim.case1[,grep("estimate.diff.v",colnames(f.sim.case1))], 1, mean) # - psiLj.case1.tau4$true
psiLj.case1.tau4$ese <- apply(f.sim.case1[,grep("estimate.diff.v",colnames(f.sim.case1))], 1, sd)
psiLj.case1.tau4$ase <- apply(f.sim.case1[,grep("estimate.diff.se",colnames(f.sim.case1))], 1, mean)
psiLj.case1.tau4$cp <- apply(f.sim.case1[,grep("estimate.diff.v",colnames(f.sim.case1))]+qnorm(0.05/2)*f.sim.case1[,grep("estimate.diff.se",colnames(f.sim.case1))] < psiLj.case1.tau4$true
                             & f.sim.case1[,grep("estimate.diff.v",colnames(f.sim.case1))]+qnorm(1-0.05/2)*f.sim.case1[,grep("estimate.diff.se",colnames(f.sim.case1))] > psiLj.case1.tau4$true, 1,sum)/nsim
psiLj.case1.tau4 <- psiLj.case1.tau4[group.survival,]

# case2
psiLj.case2.tau4$true <- as.numeric(true.value.case2.tau4[, c(c("rmst.diff", "rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff"), c("rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff"), rep(c("rmst.diff", "rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff"), 2))])
psiLj.case2.tau4$effect <- apply(f.sim.case2[,grep("estimate.diff.v",colnames(f.sim.case2))], 1, mean) # - psiLj.case2.tau4$true
psiLj.case2.tau4$ese <- apply(f.sim.case2[,grep("estimate.diff.v",colnames(f.sim.case2))], 1, sd)
psiLj.case2.tau4$ase <- apply(f.sim.case2[,grep("estimate.diff.se",colnames(f.sim.case2))], 1, mean)
psiLj.case2.tau4$cp <- apply(f.sim.case2[,grep("estimate.diff.v",colnames(f.sim.case2))]+qnorm(0.05/2)*f.sim.case2[,grep("estimate.diff.se",colnames(f.sim.case2))] < psiLj.case2.tau4$true
                             & f.sim.case2[,grep("estimate.diff.v",colnames(f.sim.case2))]+qnorm(1-0.05/2)*f.sim.case2[,grep("estimate.diff.se",colnames(f.sim.case2))] > psiLj.case2.tau4$true, 1,sum)/nsim
psiLj.case2.tau4 <- psiLj.case2.tau4[group.survival,]

# case3
psiLj.case3.tau4$true <- as.numeric(true.value.case3.tau4[, c(c("rmtlj.diff", "rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff"), c("rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff"), rep(c("rmtlj.diff", "rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff"), 2))])
psiLj.case3.tau4$effect <- apply(f.sim.case3[,grep("estimate.diff.v",colnames(f.sim.case3))], 1, mean) # - psiLj.case3.tau4$true
psiLj.case3.tau4$ese <- apply(f.sim.case3[,grep("estimate.diff.v",colnames(f.sim.case3))], 1, sd)
psiLj.case3.tau4$ase <- apply(f.sim.case3[,grep("estimate.diff.se",colnames(f.sim.case3))], 1, mean)
psiLj.case3.tau4$cp <- apply(f.sim.case3[,grep("estimate.diff.v",colnames(f.sim.case3))]+qnorm(0.05/2)*f.sim.case3[,grep("estimate.diff.se",colnames(f.sim.case3))] < psiLj.case3.tau4$true
                             & f.sim.case3[,grep("estimate.diff.v",colnames(f.sim.case3))]+qnorm(1-0.05/2)*f.sim.case3[,grep("estimate.diff.se",colnames(f.sim.case3))] > psiLj.case3.tau4$true, 1,sum)/nsim
psiLj.case3.tau4 <- psiLj.case3.tau4[group.competing,]

# case4
psiLj.case4.tau4$true <- as.numeric(true.value.case4.tau4[, c(c("rmtlj.diff", "rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff"), c("rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff"), rep(c("rmtlj.diff", "rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff"), 2))])
psiLj.case4.tau4$effect <- apply(f.sim.case4[,grep("estimate.diff.v",colnames(f.sim.case4))], 1, mean) # - psiLj.case4.tau4$true
psiLj.case4.tau4$ese <- apply(f.sim.case4[,grep("estimate.diff.v",colnames(f.sim.case4))], 1, sd)
psiLj.case4.tau4$ase <- apply(f.sim.case4[,grep("estimate.diff.se",colnames(f.sim.case4))], 1, mean)
psiLj.case4.tau4$cp <- apply(f.sim.case4[,grep("estimate.diff.v",colnames(f.sim.case4))]+qnorm(0.05/2)*f.sim.case4[,grep("estimate.diff.se",colnames(f.sim.case4))] < psiLj.case4.tau4$true
                             & f.sim.case4[,grep("estimate.diff.v",colnames(f.sim.case4))]+qnorm(1-0.05/2)*f.sim.case4[,grep("estimate.diff.se",colnames(f.sim.case4))] > psiLj.case4.tau4$true, 1,sum)/nsim
psiLj.case4.tau4 <- psiLj.case4.tau4[group.competing,]

## plot
# effect
psiLj.effect.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=effect))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  annotate("text", x=4, y=-0.05, label= "ATE")+
  annotate("text", x=9.5, y=-0.05, label= "ATO")+
  annotate("text", x=13.5, y=-0.05, label= "ATM")+
  annotate("text", x=17.5, y=-0.05, label= "ATEN")+
  scale_y_continuous(limits=c(-0.8, 0), breaks=c(-0.8, -0.6, -0.4, -0.2, 0))

psiLj.effect.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=effect))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  annotate("text", x=4, y=-0.05, label= "ATE")+
  annotate("text", x=9.5, y=-0.05, label= "ATO")+
  annotate("text", x=13.5, y=-0.05, label= "ATM")+
  annotate("text", x=17.5, y=-0.05, label= "ATEN")+
  scale_y_continuous(limits=c(-0.8, 0), breaks=c(-0.8, -0.6, -0.4, -0.2, 0))

psiLj.effect.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=effect))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  annotate("text", x=4, y=0.8, label= "ATE")+
  annotate("text", x=9.5, y=0.8, label= "ATO")+
  annotate("text", x=13.5, y=0.8, label= "ATM")+
  annotate("text", x=17.5, y=0.8, label= "ATEN")+
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8))

psiLj.effect.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=effect))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  annotate("text", x=4, y=0.8, label= "ATE")+
  annotate("text", x=9.5, y=0.8, label= "ATO")+
  annotate("text", x=13.5, y=0.8, label= "ATM")+
  annotate("text", x=17.5, y=0.8, label= "ATEN")+
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8))

# ese
psiLj.ese.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=ese))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  annotate("text", x=4, y=0.4, label= "ATE")+
  annotate("text", x=9.5, y=0.4, label= "ATO")+
  annotate("text", x=13.5, y=0.4, label= "ATM")+
  annotate("text", x=17.5, y=0.4, label= "ATEN")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=ese))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  annotate("text", x=4, y=0.4, label= "ATE")+
  annotate("text", x=9.5, y=0.4, label= "ATO")+
  annotate("text", x=13.5, y=0.4, label= "ATM")+
  annotate("text", x=17.5, y=0.4, label= "ATEN")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=ese))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  annotate("text", x=4, y=0.4, label= "ATE")+
  annotate("text", x=9.5, y=0.4, label= "ATO")+
  annotate("text", x=13.5, y=0.4, label= "ATM")+
  annotate("text", x=17.5, y=0.4, label= "ATEN")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=ese))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  annotate("text", x=4, y=0.4, label= "ATE")+
  annotate("text", x=9.5, y=0.4, label= "ATO")+
  annotate("text", x=13.5, y=0.4, label= "ATM")+
  annotate("text", x=17.5, y=0.4, label= "ATEN")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# ase
psiLj.ase.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=ase))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  annotate("text", x=4, y=0.4, label= "ATE")+
  annotate("text", x=9.5, y=0.4, label= "ATO")+
  annotate("text", x=13.5, y=0.4, label= "ATM")+
  annotate("text", x=17.5, y=0.4, label= "ATEN")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=ase))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  annotate("text", x=4, y=0.4, label= "ATE")+
  annotate("text", x=9.5, y=0.4, label= "ATO")+
  annotate("text", x=13.5, y=0.4, label= "ATM")+
  annotate("text", x=17.5, y=0.4, label= "ATEN")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=ase))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  annotate("text", x=4, y=0.4, label= "ATE")+
  annotate("text", x=9.5, y=0.4, label= "ATO")+
  annotate("text", x=13.5, y=0.4, label= "ATM")+
  annotate("text", x=17.5, y=0.4, label= "ATEN")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=ase))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  annotate("text", x=4, y=0.4, label= "ATE")+
  annotate("text", x=9.5, y=0.4, label= "ATO")+
  annotate("text", x=13.5, y=0.4, label= "ATM")+
  annotate("text", x=17.5, y=0.4, label= "ATEN")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# cp
psiLj.cp.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=cp))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  annotate("text", x=4, y=1.05, label= "ATE")+
  annotate("text", x=9.5, y=1.05, label= "ATO")+
  annotate("text", x=13.5, y=1.05, label= "ATM")+
  annotate("text", x=17.5, y=1.05, label= "ATEN")+
  scale_y_continuous(labels = scales::percent, limits=c(0.3, 1.05),
                     breaks=c(0.30, 0.50, 0.70, 0.90, 0.95, 1))

psiLj.cp.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=cp))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  annotate("text", x=4, y=1.05, label= "ATE")+
  annotate("text", x=9.5, y=1.05, label= "ATO")+
  annotate("text", x=13.5, y=1.05, label= "ATM")+
  annotate("text", x=17.5, y=1.05, label= "ATEN")+
  scale_y_continuous(labels = scales::percent, limits=c(0.3, 1.05),
                     breaks=c(0.30, 0.50, 0.70, 0.90, 0.95, 1))

psiLj.cp.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=cp))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  annotate("text", x=4, y=1.05, label= "ATE")+
  annotate("text", x=9.5, y=1.05, label= "ATO")+
  annotate("text", x=13.5, y=1.05, label= "ATM")+
  annotate("text", x=17.5, y=1.05, label= "ATEN")+
  scale_y_continuous(labels = scales::percent, limits=c(0.3, 1.05),
                     breaks=c(0.30, 0.50, 0.70, 0.90, 0.95, 1))

psiLj.cp.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=cp))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 16, 17), 3)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  annotate("text", x=4, y=1.05, label= "ATE")+
  annotate("text", x=9.5, y=1.05, label= "ATO")+
  annotate("text", x=13.5, y=1.05, label= "ATM")+
  annotate("text", x=17.5, y=1.05, label= "ATEN")+
  scale_y_continuous(labels = scales::percent, limits=c(0.3, 1.05),
                     breaks=c(0.30, 0.50, 0.70, 0.90, 0.95, 1))
## save
# survival
ggarrange(psiLj.effect.case1.tau4, psiLj.effect.case2.tau4,
          psiLj.ese.case1.tau4, psiLj.ese.case2.tau4,
          psiLj.cp.case1.tau4, psiLj.cp.case2.tau4, nrow=3, ncol=2,
          labels=c("(a) Setting 1: good overlap, survival",
                   "(b) Setting 2: poor overlap, survival",
                   "(c)", "(d)", "(e)", "(f)"),
          hjust=-0.1, font.label=list(size = 12, face = "plain"))
ggsave("sim_survival_main.pdf", width=8, height=12)

ggarrange(psiLj.ase.case1.tau4, psiLj.ase.case2.tau4, nrow=1, ncol=2,
          labels=c("(a) Setting 1: good overlap, survival",
                   "(b) Setting 2: poor overlap, survival"),
          hjust=-0.1, font.label=list(size = 12, face = "plain"))
ggsave("sim_survival_ase.pdf", width=8, height=4)

# defense slides
ggarrange(psiLj.effect.case1.tau4, psiLj.effect.case2.tau4,
          nrow=1, ncol=2,
          labels=c("(a) Setting 1: good overlap, survival",
                   "(b) Setting 2: poor overlap, survival"),
          hjust=-0.1, font.label=list(size = 12, face = "plain"))
ggsave("sim_survival_effect.pdf", width=8, height=4)

ggarrange(psiLj.ese.case1.tau4, psiLj.ese.case2.tau4,
          nrow=1, ncol=2,
          labels=c("(c) Setting 1: good overlap, survival",
                   "(d) Setting 2: poor overlap, survival"),
          hjust=-0.1, font.label=list(size = 12, face = "plain"))
ggsave("sim_survival_se.pdf", width=8, height=4)

ggarrange(psiLj.cp.case1.tau4, psiLj.cp.case2.tau4, nrow=1, ncol=2,
          labels=c("(e) Setting 1: good overlap, survival",
                   "(f) Setting 2: poor overlap, survival"),
          hjust=-0.1, font.label=list(size = 12, face = "plain"))
ggsave("sim_survival_cp.pdf", width=8, height=4)

# competing
ggarrange(psiLj.effect.case3.tau4, psiLj.effect.case4.tau4,
          psiLj.ese.case3.tau4, psiLj.ese.case4.tau4,
          psiLj.cp.case3.tau4, psiLj.cp.case4.tau4, nrow=3, ncol=2,
          labels=c("(a) Setting 1: good overlap, competing",
                   "(b) Setting 2: poor overlap, competing",
                   "(c)", "(d)", "(e)", "(f)"),
          hjust=-0.1, font.label=list(size = 12, face = "plain"))
ggsave("sim_competing_main.pdf", width=8, height=12)

ggarrange(psiLj.ase.case3.tau4, psiLj.ase.case4.tau4, nrow=1, ncol=2,
          labels=c("(a) Setting 1: good overlap, competing",
                   "(b) Setting 2: poor overlap, competing"),
          hjust=-0.1, font.label=list(size = 12, face = "plain"))
ggsave("sim_competing_ase.pdf", width=8, height=4)






# # effect
# ggarrange(psiLj.effect.case1.tau4, psiLj.effect.case2.tau4, nrow=1, ncol=2,
#           labels=c("(a)             Setting 1: good overlap, survival",
#                    "(b)             Setting 2: poor overlap, survival"),
#           hjust=-0.1, font.label=list(size = 12, face = "plain"))
# ggsave("sim_survival_effect.pdf", width=8, height=4)
# 
# ggarrange(psiLj.effect.case3.tau4, psiLj.effect.case4.tau4, nrow=1, ncol=2,
#           labels=c("(a)           Setting 3: good overlap, competing",
#                    "(b)           Setting 4: poor overlap, competing"),
#           hjust=-0.1, font.label=list(size = 12, face = "plain"))
# ggsave("sim_competing_effect.pdf", width=8, height=4)
# 
# # se
# ggarrange(psiLj.ese.case1.tau4, psiLj.ese.case2.tau4,
#           psiLj.ase.case1.tau4, psiLj.ase.case2.tau4, nrow=2, ncol=2,
#           labels=c("(a)             Setting 1: good overlap, survival",
#                    "(b)             Setting 2: poor overlap, survival",
#                    "(c)             Setting 1: good overlap, survival",
#                    "(d)             Setting 2: poor overlap, survival"),
#           hjust=-0.1, font.label=list(size = 12, face = "plain"))
# ggsave("sim_survival_se.pdf", width=8, height=8)
# 
# ggarrange(psiLj.ese.case3.tau4, psiLj.ese.case4.tau4,
#           psiLj.ase.case3.tau4, psiLj.ase.case4.tau4, nrow=2, ncol=2,
#           labels=c("(a)           Setting 3: good overlap, competing",
#                    "(b)           Setting 4: poor overlap, competing",
#                    "(c)           Setting 3: good overlap, competing",
#                    "(d)           Setting 4: poor overlap, competing"),
#           hjust=-0.1, font.label=list(size = 12, face = "plain"))
# ggsave("sim_competing_se.pdf", width=8, height=8)
# 
# #cp
# ggarrange(psiLj.cp.case1.tau4, psiLj.cp.case2.tau4, nrow=1, ncol=2,
#           labels=c("(a)             Setting 1: good overlap, survival",
#                    "(b)             Setting 2: poor overlap, survival"),
#           hjust=-0.1, font.label=list(size = 12, face = "plain"))
# ggsave("sim_survival_cp.pdf", width=8, height=4)
# 
# ggarrange(psiLj.cp.case3.tau4, psiLj.cp.case4.tau4, nrow=1, ncol=2,
#           labels=c("(a)           Setting 3: good overlap, competing",
#                    "(b)           Setting 4: poor overlap, competing"),
#           hjust=-0.1, font.label=list(size = 12, face = "plain"))
# ggsave("sim_competing_cp.pdf", width=8, height=4)


### dml integral dr
## tau=4
# case1
load("survival/scenario 1/true_survival_1.RData")
load("survival/scenario 1/sim1.RData")
f.sim.case1 <- f.sim.result$tau4
f.sim.case1 <- f.sim.case1[grep("mc.ha", rownames(f.sim.case1), value = TRUE), ]
rm(f.sim.result)

true.value.case1.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("survival/scenario 1/sim1-baseline.RData")
f.sim.result <- f.sim.result[!(rownames(f.sim.result) %in% c("dr.integral.driptw", "dr.integral.driptw.trunc", "dr.integral.drow", "dr.integral.drmw", "dr.integral.drenw")), ]
f.sim.case1 <- rbind(f.sim.case1, f.sim.result)
rm(f.sim.result)

# case2
load("survival/scenario 2/true_survival_2.RData")
load("survival/scenario 2/sim2.RData")
f.sim.case2 <- f.sim.result$tau4
f.sim.case2 <- f.sim.case2[grep("mc.ha", rownames(f.sim.case2), value = TRUE), ]
rm(f.sim.result)

true.value.case2.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("survival/scenario 2/sim2-baseline.RData")
f.sim.result <- f.sim.result[!(rownames(f.sim.result) %in% c("dr.integral.driptw", "dr.integral.driptw.trunc", "dr.integral.drow", "dr.integral.drmw", "dr.integral.drenw")), ]
f.sim.case2 <- rbind(f.sim.case2, f.sim.result)
rm(f.sim.result)

# case3
load("competing/scenario 1/true_competing_7.RData")
load("competing/scenario 1/sim7.RData")
f.sim.case3 <- f.sim.result$tau4
f.sim.case3 <- f.sim.case3[grep("mc.ha", rownames(f.sim.case3), value = TRUE), ]
rm(f.sim.result)

true.value.case3.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("competing/scenario 1/sim7-baseline.RData")
f.sim.result <- f.sim.result[!(rownames(f.sim.result) %in% c("dr.integral.driptw", "dr.integral.driptw.trunc", "dr.integral.drow", "dr.integral.drmw", "dr.integral.drenw")), ]
f.sim.case3 <- rbind(f.sim.case3, f.sim.result)
rm(f.sim.result)

# case4
load("competing/scenario 2/true_competing_8.RData")
load("competing/scenario 2/sim8.RData")
f.sim.case4 <- f.sim.result$tau4
f.sim.case4 <- f.sim.case4[grep("mc.ha", rownames(f.sim.case4), value = TRUE), ]
rm(f.sim.result)

true.value.case4.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("competing/scenario 2/sim8-baseline.RData")
f.sim.result <- f.sim.result[!(rownames(f.sim.result) %in% c("dr.integral.driptw", "dr.integral.driptw.trunc", "dr.integral.drow", "dr.integral.drmw", "dr.integral.drenw")), ]
f.sim.case4 <- rbind(f.sim.case4, f.sim.result)
rm(f.sim.result)

## updated tau 4
nsim <- 1000
psiLj.case1.tau4 <- psiLj.case2.tau4 <- psiLj.case3.tau4 <- psiLj.case4.tau4 <- data.frame(matrix(NA, nrow=19, ncol=5))
colnames(psiLj.case1.tau4) <- colnames(psiLj.case2.tau4) <- colnames(psiLj.case3.tau4) <- colnames(psiLj.case4.tau4) <- c("true","effect","ese","ase","cp")
rownames(psiLj.case1.tau4) <- rownames(psiLj.case2.tau4) <-
  c("S.driptw", "S.driptw.t", "S.drow", "S.drmw", "S.drenw", "RMST.driptw", "RMST.driptw.t", "RMST.drow", "RMST.drmw", "RMST.drenw",
                                                              "S.or", "S.or.ow", "S.or.mw", "S.or.enw", "S.ipcw.iptw", "S.ipcw.iptw.t", "S.ipcw.ow", "S.ipcw.mw", "S.ipcw.enw")
rownames(psiLj.case3.tau4) <- rownames(psiLj.case4.tau4) <-
  c("Fj.driptw", "Fj.driptw.t", "Fj.drow", "Fj.drmw", "Fj.drenw", "RMTLj.driptw", "RMTLj.driptw.t", "RMTLj.drow", "RMTLj.drmw", "RMTLj.drenw",
    "Fj.or", "Fj.or.ow", "Fj.or.mw", "Fj.or.enw", "Fj.ipcw.iptw", "Fj.ipcw.iptw.t", "Fj.ipcw.ow", "Fj.ipcw.mw", "Fj.ipcw.enw")

group.survival <- c("S.or", "S.ipcw.iptw", "S.driptw", "RMST.driptw",
                    "S.ipcw.iptw.t", "S.driptw.t", "RMST.driptw.t",
                    "S.or.ow", "S.ipcw.ow", "S.drow", "RMST.drow",
                    "S.or.mw", "S.ipcw.mw", "S.drmw", "RMST.drmw",
                    "S.or.enw", "S.ipcw.enw", "S.drenw", "RMST.drenw")

group.competing <- c("Fj.or", "Fj.ipcw.iptw", "Fj.driptw", "RMTLj.driptw",
                     "Fj.ipcw.iptw.t", "Fj.driptw.t", "RMTLj.driptw.t",
                     "Fj.or.ow", "Fj.ipcw.ow", "Fj.drow", "RMTLj.drow",
                     "Fj.or.mw", "Fj.ipcw.mw", "Fj.drmw", "RMTLj.drmw",
                     "Fj.or.enw", "Fj.ipcw.enw", "Fj.drenw", "RMTLj.drenw")
## ate
# case1
psiLj.case1.tau4$true <- as.numeric(true.value.case1.tau4[, c(rep(c("rmst.diff", "rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff"), 2), c("rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff", c("rmst.diff", "rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff")))])
psiLj.case1.tau4$effect <- apply(f.sim.case1[,grep("estimate.diff.v",colnames(f.sim.case1))], 1, mean) # - psiLj.case1.tau4$true
psiLj.case1.tau4$ese <- apply(f.sim.case1[,grep("estimate.diff.v",colnames(f.sim.case1))], 1, sd)
psiLj.case1.tau4$ase <- apply(f.sim.case1[,grep("estimate.diff.se",colnames(f.sim.case1))], 1, mean)
psiLj.case1.tau4$cp <- apply(f.sim.case1[,grep("estimate.diff.v",colnames(f.sim.case1))]+qnorm(0.05/2)*f.sim.case1[,grep("estimate.diff.se",colnames(f.sim.case1))] < psiLj.case1.tau4$true
                             & f.sim.case1[,grep("estimate.diff.v",colnames(f.sim.case1))]+qnorm(1-0.05/2)*f.sim.case1[,grep("estimate.diff.se",colnames(f.sim.case1))] > psiLj.case1.tau4$true, 1,sum)/nsim
psiLj.case1.tau4 <- psiLj.case1.tau4[group.survival,]

# case2
psiLj.case2.tau4$true <- as.numeric(true.value.case2.tau4[, c(rep(c("rmst.diff", "rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff"), 2), c("rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff", c("rmst.diff", "rmst.diff", "rmst.ow.diff", "rmst.mw.diff", "rmst.enw.diff")))])
psiLj.case2.tau4$effect <- apply(f.sim.case2[,grep("estimate.diff.v",colnames(f.sim.case2))], 1, mean) # - psiLj.case2.tau4$true
psiLj.case2.tau4$ese <- apply(f.sim.case2[,grep("estimate.diff.v",colnames(f.sim.case2))], 1, sd)
psiLj.case2.tau4$ase <- apply(f.sim.case2[,grep("estimate.diff.se",colnames(f.sim.case2))], 1, mean)
psiLj.case2.tau4$cp <- apply(f.sim.case2[,grep("estimate.diff.v",colnames(f.sim.case2))]+qnorm(0.05/2)*f.sim.case2[,grep("estimate.diff.se",colnames(f.sim.case2))] < psiLj.case2.tau4$true
                             & f.sim.case2[,grep("estimate.diff.v",colnames(f.sim.case2))]+qnorm(1-0.05/2)*f.sim.case2[,grep("estimate.diff.se",colnames(f.sim.case2))] > psiLj.case2.tau4$true, 1,sum)/nsim
psiLj.case2.tau4 <- psiLj.case2.tau4[group.survival,]

# case3
psiLj.case3.tau4$true <- as.numeric(true.value.case3.tau4[, c(rep(c("rmtlj.diff", "rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff"), 2), c("rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff", c("rmtlj.diff", "rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff")))])
psiLj.case3.tau4$effect <- apply(f.sim.case3[,grep("estimate.diff.v",colnames(f.sim.case3))], 1, mean) # - psiLj.case3.tau4$true
psiLj.case3.tau4$ese <- apply(f.sim.case3[,grep("estimate.diff.v",colnames(f.sim.case3))], 1, sd)
psiLj.case3.tau4$ase <- apply(f.sim.case3[,grep("estimate.diff.se",colnames(f.sim.case3))], 1, mean)
psiLj.case3.tau4$cp <- apply(f.sim.case3[,grep("estimate.diff.v",colnames(f.sim.case3))]+qnorm(0.05/2)*f.sim.case3[,grep("estimate.diff.se",colnames(f.sim.case3))] < psiLj.case3.tau4$true
                             & f.sim.case3[,grep("estimate.diff.v",colnames(f.sim.case3))]+qnorm(1-0.05/2)*f.sim.case3[,grep("estimate.diff.se",colnames(f.sim.case3))] > psiLj.case3.tau4$true, 1,sum)/nsim
psiLj.case3.tau4 <- psiLj.case3.tau4[group.competing,]

# case4
psiLj.case4.tau4$true <- as.numeric(true.value.case4.tau4[, c(rep(c("rmtlj.diff", "rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff"), 2), c("rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff", c("rmtlj.diff", "rmtlj.diff", "rmtlj.ow.diff", "rmtlj.mw.diff", "rmtlj.enw.diff")))])
psiLj.case4.tau4$effect <- apply(f.sim.case4[,grep("estimate.diff.v",colnames(f.sim.case4))], 1, mean) # - psiLj.case4.tau4$true
psiLj.case4.tau4$ese <- apply(f.sim.case4[,grep("estimate.diff.v",colnames(f.sim.case4))], 1, sd)
psiLj.case4.tau4$ase <- apply(f.sim.case4[,grep("estimate.diff.se",colnames(f.sim.case4))], 1, mean)
psiLj.case4.tau4$cp <- apply(f.sim.case4[,grep("estimate.diff.v",colnames(f.sim.case4))]+qnorm(0.05/2)*f.sim.case4[,grep("estimate.diff.se",colnames(f.sim.case4))] < psiLj.case4.tau4$true
                             & f.sim.case4[,grep("estimate.diff.v",colnames(f.sim.case4))]+qnorm(1-0.05/2)*f.sim.case4[,grep("estimate.diff.se",colnames(f.sim.case4))] > psiLj.case4.tau4$true, 1,sum)/nsim
psiLj.case4.tau4 <- psiLj.case4.tau4[group.competing,]

## plot
# effect
psiLj.effect.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=effect))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  scale_y_continuous(limits=c(-0.8, 0.8), breaks=c(-0.8, -0.4, 0, 0.4, 0.8))

psiLj.effect.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=effect))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  scale_y_continuous(limits=c(-0.8, 0.8), breaks=c(-0.8, -0.4, 0, 0.4, 0.8))

psiLj.effect.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=effect))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  scale_y_continuous(limits=c(-0.8, 0.8), breaks=c(-0.8, -0.4, 0, 0.4, 0.8))

psiLj.effect.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=effect))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  scale_y_continuous(limits=c(-0.8, 0.8), breaks=c(-0.8, -0.4, 0, 0.4, 0.8))

# ese
psiLj.ese.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=ese))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=ese))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=ese))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=ese))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# ase
psiLj.ase.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=ase))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=ase))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=ase))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=ase))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# cp
psiLj.cp.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=cp))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.3,1),
                     breaks=c(0.30, 0.50, 0.70, 0.90, 0.95, 1))

psiLj.cp.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=cp))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.3,1),
                     breaks=c(0.30, 0.50, 0.70, 0.90, 0.95, 1))

psiLj.cp.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=cp))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.3,1),
                     breaks=c(0.30, 0.50, 0.70, 0.90, 0.95, 1))

psiLj.cp.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=cp))+
  geom_point(size=2, shape=c(16, 16, 16, 17, 16, 16, 17, rep(c(16, 16, 17, 17), 3)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  geom_vline(xintercept=c(7.5, 11.5, 15.5), linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.3,1),
                     breaks=c(0.30, 0.50, 0.70, 0.90, 0.95, 1))

# ggsave("sim_tau4a.pdf", width = 20, height = 20,
#        plot_grid(psiLj.effect.case1.tau4, psiLj.ese.case1.tau4, psiLj.ase.case1.tau4, psiLj.cp.case1.tau4,
#                  psiLj.effect.case2.tau4, psiLj.ese.case2.tau4, psiLj.ase.case2.tau4, psiLj.cp.case2.tau4,
#                  psiLj.effect.case3.tau4, psiLj.ese.case3.tau4, psiLj.ase.case3.tau4, psiLj.cp.case3.tau4,
#                  psiLj.effect.case4.tau4, psiLj.ese.case4.tau4, psiLj.ase.case4.tau4, psiLj.cp.case4.tau4,
#                  labels=c("(a)","","","","(b)","","","","(c)","","","","(d)","","",""), vjust = 1.5, ncol=4, nrow=4))

# ggsave("sim_tau4c.pdf", width = 12, height = 12,
#        plot_grid(psiLj.effect.case1.tau4, psiLj.effect.case2.tau4, psiLj.effect.case3.tau4, psiLj.effect.case4.tau4,
#                  psiLj.ese.case1.tau4, psiLj.ese.case2.tau4, psiLj.ese.case3.tau4, psiLj.ese.case4.tau4,
#                  psiLj.ase.case1.tau4, psiLj.ase.case2.tau4, psiLj.ase.case3.tau4, psiLj.ase.case4.tau4,
#                  psiLj.cp.case1.tau4, psiLj.cp.case2.tau4, psiLj.cp.case3.tau4, psiLj.cp.case4.tau4,
#                  labels=c("(a)","(b)","(c)","(d)","","","","","","","","","","","",""), vjust = 1.5, ncol=4, nrow=4))

## second pass
# Bai
load("survival/scenario 1/true_survival_1.RData")
load("survival/scenario 1/sim1.RData")
f.sim.case1 <- f.sim.result
rm(f.sim.result)

true.value.case1.tau1 <- true.value.tau1
true.value.case1.tau2 <- true.value.tau2
true.value.case1.tau3 <- true.value.tau3
true.value.case1.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("survival/scenario 1/sim1-Bai.RData")
f.sim.result <- lapply(f.sim.result, `rownames<-`, "Bai")
f.sim.case1 <- mapply(rbind, f.sim.result, f.sim.case1, SIMPLIFY=FALSE)
rm(f.sim.result)

# case 2
load("survival/scenario 2/true_survival_2.RData")
load("survival/scenario 2/sim2.RData")
f.sim.case2 <- f.sim.result
rm(f.sim.result)

true.value.case2.tau1 <- true.value.tau1
true.value.case2.tau2 <- true.value.tau2
true.value.case2.tau3 <- true.value.tau3
true.value.case2.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("survival/scenario 2/sim2-Bai.RData")
f.sim.result <- lapply(f.sim.result, `rownames<-`, "Bai")
f.sim.case2 <- mapply(rbind, f.sim.result, f.sim.case2, SIMPLIFY=FALSE)
rm(f.sim.result)

# case 3
load("competing/scenario 1/true_competing_7.RData")
load("competing/scenario 1/sim7.RData")
f.sim.case3 <- f.sim.result
rm(f.sim.result)

true.value.case3.tau1 <- true.value.tau1
true.value.case3.tau2 <- true.value.tau2
true.value.case3.tau3 <- true.value.tau3
true.value.case3.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("competing/scenario 1/sim7-Ozenne.RData")
f.sim.result <- lapply(f.sim.result, `rownames<-`, "Ozenne")
f.sim.case3 <- mapply(rbind, f.sim.result, f.sim.case3, SIMPLIFY=FALSE)
rm(f.sim.result)

# case 4
load("competing/scenario 2/true_competing_8.RData")
load("competing/scenario 2/sim8.RData")
f.sim.case4 <- f.sim.result
rm(f.sim.result)

true.value.case4.tau1 <- true.value.tau1
true.value.case4.tau2 <- true.value.tau2
true.value.case4.tau3 <- true.value.tau3
true.value.case4.tau4 <- true.value.tau4
rm(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4)

load("competing/scenario 2/sim8-Ozenne.RData")
f.sim.result <- lapply(f.sim.result, `rownames<-`, "Ozenne")
f.sim.case4 <- mapply(rbind, f.sim.result, f.sim.case4, SIMPLIFY=FALSE)
rm(f.sim.result)

## updated tau 4
nsim <- 1000
psiLj.case1.tau4 <- psiLj.case2.tau4 <- psiLj.case3.tau4 <- psiLj.case4.tau4 <- data.frame(matrix(NA, nrow=41, ncol=5))
colnames(psiLj.case1.tau4) <- colnames(psiLj.case2.tau4) <- colnames(psiLj.case3.tau4) <- colnames(psiLj.case4.tau4) <- c("true","effect","ese","ase","cp")
rownames(psiLj.case1.tau4) <- rownames(psiLj.case2.tau4) <-
  c("Bai","S.mc.ht.driptw", "S.mc.ht.driptw.t", "S.mc.ht.drow", "S.mc.ht.drmw", "S.mc.ht.drenw",
    "S.mc.ha.driptw", "S.mc.ha.driptw.t", "S.mc.ha.drow", "S.mc.ha.drmw", "S.mc.ha.drenw",
    "S.mt.ht.driptw", "S.mt.ht.driptw.t", "S.mt.ht.drow", "S.mt.ht.drmw", "S.mt.ht.drenw",
    "S.mt.ha.driptw", "S.mt.ha.driptw.t", "S.mt.ha.drow", "S.mt.ha.drmw", "S.mt.ha.drenw",
    "RMST.mc.ht.driptw", "RMST.mc.ht.driptw.t", "RMST.mc.ht.drow", "RMST.mc.ht.drmw", "RMST.mc.ht.drenw",
    "RMST.mc.ha.driptw", "RMST.mc.ha.driptw.t", "RMST.mc.ha.drow", "RMST.mc.ha.drmw", "RMST.mc.ha.drenw",
    "RMST.mt.ht.driptw", "RMST.mt.ht.driptw.t", "RMST.mt.ht.drow", "RMST.mt.ht.drmw", "RMST.mt.ht.drenw",
    "RMST.mt.ha.driptw", "RMST.mt.ha.driptw.t", "RMST.mt.ha.drow", "RMST.mt.ha.drmw", "RMST.mt.ha.drenw")
rownames(psiLj.case3.tau4) <- rownames(psiLj.case4.tau4) <-
  c("Ozenne","Fj.mc.ht.driptw", "Fj.mc.ht.driptw.t", "Fj.mc.ht.drow", "Fj.mc.ht.drmw", "Fj.mc.ht.drenw",
    "Fj.mc.ha.driptw", "Fj.mc.ha.driptw.t", "Fj.mc.ha.drow", "Fj.mc.ha.drmw", "Fj.mc.ha.drenw",
    "Fj.mt.ht.driptw", "Fj.mt.ht.driptw.t", "Fj.mt.ht.drow", "Fj.mt.ht.drmw", "Fj.mt.ht.drenw",
    "Fj.mt.ha.driptw", "Fj.mt.ha.driptw.t", "Fj.mt.ha.drow", "Fj.mt.ha.drmw", "Fj.mt.ha.drenw",
    "RMTLj.mc.ht.driptw", "RMTLj.mc.ht.driptw.t", "RMTLj.mc.ht.drow", "RMTLj.mc.ht.drmw", "RMTLj.mc.ht.drenw",
    "RMTLj.mc.ha.driptw", "RMTLj.mc.ha.driptw.t", "RMTLj.mc.ha.drow", "RMTLj.mc.ha.drmw", "RMTLj.mc.ha.drenw",
    "RMTLj.mt.ht.driptw", "RMTLj.mt.ht.driptw.t", "RMTLj.mt.ht.drow", "RMTLj.mt.ht.drmw", "RMTLj.mt.ht.drenw",
    "RMTLj.mt.ha.driptw", "RMTLj.mt.ha.driptw.t", "RMTLj.mt.ha.drow", "RMTLj.mt.ha.drmw", "RMTLj.mt.ha.drenw")
group.survival <- c("Bai","S.mc.ht.driptw", "S.mc.ha.driptw", "S.mt.ht.driptw", "S.mt.ha.driptw", "RMST.mc.ht.driptw", "RMST.mc.ha.driptw", "RMST.mt.ht.driptw", "RMST.mt.ha.driptw",
                    "S.mc.ht.driptw.t", "S.mc.ha.driptw.t", "S.mt.ht.driptw.t", "S.mt.ha.driptw.t", "RMST.mc.ht.driptw.t", "RMST.mc.ha.driptw.t", "RMST.mt.ht.driptw.t", "RMST.mt.ha.driptw.t",
                    "S.mc.ht.drow", "S.mc.ha.drow", "S.mt.ht.drow", "S.mt.ha.drow", "RMST.mc.ht.drow", "RMST.mc.ha.drow", "RMST.mt.ht.drow", "RMST.mt.ha.drow",
                    "S.mc.ht.drmw", "S.mc.ha.drmw", "S.mt.ht.drmw", "S.mt.ha.drmw", "RMST.mc.ht.drmw", "RMST.mc.ha.drmw", "RMST.mt.ht.drmw", "RMST.mt.ha.drmw",
                    "S.mc.ht.drenw", "S.mc.ha.drenw", "S.mt.ht.drenw", "S.mt.ha.drenw", "RMST.mc.ht.drenw", "RMST.mc.ha.drenw", "RMST.mt.ht.drenw", "RMST.mt.ha.drenw")
group.competing <- c("Ozenne","Fj.mc.ht.driptw", "Fj.mc.ha.driptw", "Fj.mt.ht.driptw", "Fj.mt.ha.driptw", "RMTLj.mc.ht.driptw", "RMTLj.mc.ha.driptw", "RMTLj.mt.ht.driptw", "RMTLj.mt.ha.driptw",
                     "Fj.mc.ht.driptw.t", "Fj.mc.ha.driptw.t", "Fj.mt.ht.driptw.t", "Fj.mt.ha.driptw.t", "RMTLj.mc.ht.driptw.t", "RMTLj.mc.ha.driptw.t", "RMTLj.mt.ht.driptw.t", "RMTLj.mt.ha.driptw.t",
                     "Fj.mc.ht.drow", "Fj.mc.ha.drow", "Fj.mt.ht.drow", "Fj.mt.ha.drow", "RMTLj.mc.ht.drow", "RMTLj.mc.ha.drow", "RMTLj.mt.ht.drow", "RMTLj.mt.ha.drow",
                     "Fj.mc.ht.drmw", "Fj.mc.ha.drmw", "Fj.mt.ht.drmw", "Fj.mt.ha.drmw", "RMTLj.mc.ht.drmw", "RMTLj.mc.ha.drmw", "RMTLj.mt.ht.drmw", "RMTLj.mt.ha.drmw",
                     "Fj.mc.ht.drenw", "Fj.mc.ha.drenw", "Fj.mt.ht.drenw", "Fj.mt.ha.drenw", "RMTLj.mc.ht.drenw", "RMTLj.mc.ha.drenw", "RMTLj.mt.ht.drenw", "RMTLj.mt.ha.drenw")

# ate
psiLj.case1.tau4$true <- as.numeric(true.value.case1.tau4[, c("rmst.diff" ,rep(c("rmst.diff","rmst.diff","rmst.ow.diff","rmst.mw.diff","rmst.enw.diff"), 8))])
psiLj.case1.tau4$effect <- apply(f.sim.case1$tau4[,grep("estimate.diff.v",colnames(f.sim.case1$tau4))], 1, mean) # - psiLj.case1.tau4$true
psiLj.case1.tau4$ese <- apply(f.sim.case1$tau4[,grep("estimate.diff.v",colnames(f.sim.case1$tau4))], 1, sd)
psiLj.case1.tau4$ase <- apply(f.sim.case1$tau4[,grep("estimate.diff.se",colnames(f.sim.case1$tau4))], 1, mean)
psiLj.case1.tau4$cp <- apply(f.sim.case1$tau4[,grep("estimate.diff.v",colnames(f.sim.case1$tau4))]+qnorm(0.05/2)*f.sim.case1$tau4[,grep("estimate.diff.se",colnames(f.sim.case1$tau4))] < psiLj.case1.tau4$true
                             & f.sim.case1$tau4[,grep("estimate.diff.v",colnames(f.sim.case1$tau4))]+qnorm(1-0.05/2)*f.sim.case1$tau4[,grep("estimate.diff.se",colnames(f.sim.case1$tau4))] > psiLj.case1.tau4$true, 1,sum)/nsim

psiLj.case2.tau4$true <- as.numeric(true.value.case2.tau4[, c("rmst.diff" ,rep(c("rmst.diff","rmst.diff","rmst.ow.diff","rmst.mw.diff","rmst.enw.diff"), 8))])
psiLj.case2.tau4$effect <- apply(f.sim.case2$tau4[,grep("estimate.diff.v",colnames(f.sim.case2$tau4))], 1, mean) # - psiLj.case2.tau4$true
psiLj.case2.tau4$ese <- apply(f.sim.case2$tau4[,grep("estimate.diff.v",colnames(f.sim.case2$tau4))], 1, sd)
psiLj.case2.tau4$ase <- apply(f.sim.case2$tau4[,grep("estimate.diff.se",colnames(f.sim.case2$tau4))], 1, mean)
psiLj.case2.tau4$cp <- apply(f.sim.case2$tau4[,grep("estimate.diff.v",colnames(f.sim.case2$tau4))]+qnorm(0.05/2)*f.sim.case2$tau4[,grep("estimate.diff.se",colnames(f.sim.case2$tau4))] < psiLj.case2.tau4$true
                             & f.sim.case2$tau4[,grep("estimate.diff.v",colnames(f.sim.case2$tau4))]+qnorm(1-0.05/2)*f.sim.case2$tau4[,grep("estimate.diff.se",colnames(f.sim.case2$tau4))] > psiLj.case2.tau4$true, 1,sum)/nsim

psiLj.case3.tau4$true <- as.numeric(true.value.case3.tau4[, c("rmtlj.diff" ,rep(c("rmtlj.diff","rmtlj.diff","rmtlj.ow.diff","rmtlj.mw.diff","rmtlj.enw.diff"), 8))])
psiLj.case3.tau4$effect <- apply(f.sim.case3$tau4[,grep("estimate.diff.v",colnames(f.sim.case3$tau4))], 1, mean) # - psiLj.case3.tau4$true
psiLj.case3.tau4$ese <- apply(f.sim.case3$tau4[,grep("estimate.diff.v",colnames(f.sim.case3$tau4))], 1, sd)
psiLj.case3.tau4$ase <- apply(f.sim.case3$tau4[,grep("estimate.diff.se",colnames(f.sim.case3$tau4))], 1, mean)
psiLj.case3.tau4$cp <- apply(f.sim.case3$tau4[,grep("estimate.diff.v",colnames(f.sim.case3$tau4))]+qnorm(0.05/2)*f.sim.case3$tau4[,grep("estimate.diff.se",colnames(f.sim.case3$tau4))] < psiLj.case3.tau4$true
                             & f.sim.case3$tau4[,grep("estimate.diff.v",colnames(f.sim.case3$tau4))]+qnorm(1-0.05/2)*f.sim.case3$tau4[,grep("estimate.diff.se",colnames(f.sim.case3$tau4))] > psiLj.case3.tau4$true, 1,sum)/nsim

psiLj.case4.tau4$true <- as.numeric(true.value.case4.tau4[, c("rmtlj.diff" ,rep(c("rmtlj.diff","rmtlj.diff","rmtlj.ow.diff","rmtlj.mw.diff","rmtlj.enw.diff"), 8))])
psiLj.case4.tau4$effect <- apply(f.sim.case4$tau4[,grep("estimate.diff.v",colnames(f.sim.case4$tau4))], 1, mean) # - psiLj.case4.tau4$true
psiLj.case4.tau4$ese <- apply(f.sim.case4$tau4[,grep("estimate.diff.v",colnames(f.sim.case4$tau4))], 1, sd)
psiLj.case4.tau4$ase <- apply(f.sim.case4$tau4[,grep("estimate.diff.se",colnames(f.sim.case4$tau4))], 1, mean)
psiLj.case4.tau4$cp <- apply(f.sim.case4$tau4[,grep("estimate.diff.v",colnames(f.sim.case4$tau4))]+qnorm(0.05/2)*f.sim.case4$tau4[,grep("estimate.diff.se",colnames(f.sim.case4$tau4))] < psiLj.case4.tau4$true
                             & f.sim.case4$tau4[,grep("estimate.diff.v",colnames(f.sim.case4$tau4))]+qnorm(1-0.05/2)*f.sim.case4$tau4[,grep("estimate.diff.se",colnames(f.sim.case4$tau4))] > psiLj.case4.tau4$true, 1,sum)/nsim

psiLj.case1.tau4 <- psiLj.case1.tau4[group.survival,]
psiLj.case2.tau4 <- psiLj.case2.tau4[group.survival,]
psiLj.case3.tau4 <- psiLj.case3.tau4[group.competing,]
psiLj.case4.tau4 <- psiLj.case4.tau4[group.competing,]

## plot
# effect
psiLj.effect.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=effect))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  scale_y_continuous(limits=c(-0.8, 0.8), breaks=c(-0.8, -0.4, 0, 0.4, 0.8))

psiLj.effect.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=effect))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  scale_y_continuous(limits=c(-0.8, 0.8), breaks=c(-0.8, -0.4, 0, 0.4, 0.8))

psiLj.effect.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=effect))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  scale_y_continuous(limits=c(-0.8, 0.8), breaks=c(-0.8, -0.4, 0, 0.4, 0.8))

psiLj.effect.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=effect))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_point(aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=true), shape=1, size=3)+ # yintercept=0, linetype="dotted", color = "black"
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated effect")+
  scale_y_continuous(limits=c(-0.8, 0.8), breaks=c(-0.8, -0.4, 0, 0.4, 0.8))

# ese
psiLj.ese.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=ese))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=ese))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=ese))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=ese))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# ase
psiLj.ase.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=ase))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=ase))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=ase))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=ase))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# cp
psiLj.cp.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), y=cp))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), y=cp))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), y=cp))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), y=cp))+
  geom_point(size=2, shape=c(16,16,rep(15,39)))+
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

ggsave("sim_tau4a.pdf", width = 20, height = 20,
       plot_grid(psiLj.effect.case1.tau4, psiLj.ese.case1.tau4, psiLj.ase.case1.tau4, psiLj.cp.case1.tau4,
                 psiLj.effect.case2.tau4, psiLj.ese.case2.tau4, psiLj.ase.case2.tau4, psiLj.cp.case2.tau4,
                 psiLj.effect.case3.tau4, psiLj.ese.case3.tau4, psiLj.ase.case3.tau4, psiLj.cp.case3.tau4,
                 psiLj.effect.case4.tau4, psiLj.ese.case4.tau4, psiLj.ase.case4.tau4, psiLj.cp.case4.tau4,
                 labels=c("(a)","","","","(b)","","","","(c)","","","","(d)","","",""), vjust = 1.5, ncol=4, nrow=4))

ggsave("sim_tau4b.pdf", width = 20, height = 20,
       plot_grid(psiLj.effect.case1.tau4, psiLj.effect.case2.tau4, psiLj.effect.case3.tau4, psiLj.effect.case4.tau4,
                 psiLj.ese.case1.tau4, psiLj.ese.case2.tau4, psiLj.ese.case3.tau4, psiLj.ese.case4.tau4,
                 psiLj.ase.case1.tau4, psiLj.ase.case2.tau4, psiLj.ase.case3.tau4, psiLj.ase.case4.tau4,
                 psiLj.cp.case1.tau4, psiLj.cp.case2.tau4, psiLj.cp.case3.tau4, psiLj.cp.case4.tau4,
                 labels=c("(a)","(b)","(c)","(d)","","","","","","","","","","","",""), vjust = 1.5, ncol=4, nrow=4))


ggsave("sim_tau4_effect.pdf", width = 8, height = 8,
       plot_grid(psiLj.effect.case1.tau4, psiLj.effect.case2.tau4, psiLj.effect.case3.tau4, psiLj.effect.case4.tau4,
                 labels=c("(a)","(b)","(c)","(d)"), vjust = 1.5, ncol=2, nrow=2))

ggsave("sim_tau4_ese.pdf", width = 8, height = 8,
       plot_grid(psiLj.ese.case1.tau4, psiLj.ese.case2.tau4, psiLj.ese.case3.tau4, psiLj.ese.case4.tau4,
                 labels=c("(a)","(b)","(c)","(d)"), vjust = 1.5, ncol=2, nrow=2))

ggsave("sim_tau4_ase.pdf", width = 8, height = 8,
       plot_grid(psiLj.ase.case1.tau4, psiLj.ase.case2.tau4, psiLj.ase.case3.tau4, psiLj.ase.case4.tau4,
                 labels=c("(a)","(b)","(c)","(d)"), vjust = 1.5, ncol=2, nrow=2))

ggsave("sim_tau4_cp.pdf", width = 8, height = 8,
       plot_grid(psiLj.cp.case1.tau4, psiLj.cp.case2.tau4, psiLj.cp.case3.tau4, psiLj.cp.case4.tau4,
                 labels=c("(a)","(b)","(c)","(d)"), vjust = 1.5, ncol=2, nrow=2))

# first pass
nsim <- 1000
psiLj.case1.tau1 <- psiLj.case2.tau1 <- psiLj.case3.tau1 <- psiLj.case4.tau1 <-
  psiLj.case1.tau2 <- psiLj.case2.tau2 <- psiLj.case3.tau2 <- psiLj.case4.tau2 <-
  psiLj.case1.tau3 <- psiLj.case2.tau3 <- psiLj.case3.tau3 <- psiLj.case4.tau3 <-
  psiLj.case1.tau4 <- psiLj.case2.tau4 <- psiLj.case3.tau4 <- psiLj.case4.tau4 <- data.frame(matrix(NA, nrow=41, ncol=5))
rownames(psiLj.case1.tau1) <- rownames(psiLj.case2.tau1) <- rownames(psiLj.case3.tau1) <- rownames(psiLj.case4.tau1) <-
  rownames(psiLj.case1.tau2) <- rownames(psiLj.case2.tau2) <- rownames(psiLj.case3.tau2) <- rownames(psiLj.case4.tau2) <-
  c("Bai","S.mc.ht.driptw", "S.mc.ht.driptw.t", "S.mc.ht.drow", "S.mc.ht.drmw", "S.mc.ht.drenw",
    "S.mc.ha.driptw", "S.mc.ha.driptw.t", "S.mc.ha.drow", "S.mc.ha.drmw", "S.mc.ha.drenw",
    "S.mt.ht.driptw", "S.mt.ht.driptw.t", "S.mt.ht.drow", "S.mt.ht.drmw", "S.mt.ht.drenw",
    "S.mt.ha.driptw", "S.mt.ha.driptw.t", "S.mt.ha.drow", "S.mt.ha.drmw", "S.mt.ha.drenw",
    "RMST.mc.ht.driptw", "RMST.mc.ht.driptw.t", "RMST.mc.ht.drow", "RMST.mc.ht.drmw", "RMST.mc.ht.drenw",
    "RMST.mc.ha.driptw", "RMST.mc.ha.driptw.t", "RMST.mc.ha.drow", "RMST.mc.ha.drmw", "RMST.mc.ha.drenw",
    "RMST.mt.ht.driptw", "RMST.mt.ht.driptw.t", "RMST.mt.ht.drow", "RMST.mt.ht.drmw", "RMST.mt.ht.drenw",
    "RMST.mt.ha.driptw", "RMST.mt.ha.driptw.t", "RMST.mt.ha.drow", "RMST.mt.ha.drmw", "RMST.mt.ha.drenw")
rownames(psiLj.case1.tau3) <- rownames(psiLj.case2.tau3) <- rownames(psiLj.case3.tau3) <- rownames(psiLj.case4.tau3) <-
  rownames(psiLj.case1.tau4) <- rownames(psiLj.case2.tau4) <- rownames(psiLj.case3.tau4) <- rownames(psiLj.case4.tau4) <-
  c("Bai","Fj.mc.ht.driptw", "Fj.mc.ht.driptw.t", "Fj.mc.ht.drow", "Fj.mc.ht.drmw", "Fj.mc.ht.drenw",
    "Fj.mc.ha.driptw", "Fj.mc.ha.driptw.t", "Fj.mc.ha.drow", "Fj.mc.ha.drmw", "Fj.mc.ha.drenw",
    "Fj.mt.ht.driptw", "Fj.mt.ht.driptw.t", "Fj.mt.ht.drow", "Fj.mt.ht.drmw", "Fj.mt.ht.drenw",
    "Fj.mt.ha.driptw", "Fj.mt.ha.driptw.t", "Fj.mt.ha.drow", "Fj.mt.ha.drmw", "Fj.mt.ha.drenw",
    "RMTLj.mc.ht.driptw", "RMTLj.mc.ht.driptw.t", "RMTLj.mc.ht.drow", "RMTLj.mc.ht.drmw", "RMTLj.mc.ht.drenw",
    "RMTLj.mc.ha.driptw", "RMTLj.mc.ha.driptw.t", "RMTLj.mc.ha.drow", "RMTLj.mc.ha.drmw", "RMTLj.mc.ha.drenw",
    "RMTLj.mt.ht.driptw", "RMTLj.mt.ht.driptw.t", "RMTLj.mt.ht.drow", "RMTLj.mt.ht.drmw", "RMTLj.mt.ht.drenw",
    "RMTLj.mt.ha.driptw", "RMTLj.mt.ha.driptw.t", "RMTLj.mt.ha.drow", "RMTLj.mt.ha.drmw", "RMTLj.mt.ha.drenw")

colnames(psiLj.case1.tau1) <- colnames(psiLj.case2.tau1) <- colnames(psiLj.case3.tau1) <- colnames(psiLj.case4.tau1) <-
  colnames(psiLj.case1.tau2) <- colnames(psiLj.case2.tau2) <- colnames(psiLj.case3.tau2) <- colnames(psiLj.case4.tau2) <-
  colnames(psiLj.case1.tau3) <- colnames(psiLj.case2.tau3) <- colnames(psiLj.case3.tau3) <- colnames(psiLj.case4.tau3) <-
  colnames(psiLj.case1.tau4) <- colnames(psiLj.case2.tau4) <- colnames(psiLj.case3.tau4) <- colnames(psiLj.case4.tau4) <- c("true","bias","ese","ase","cp")

## ate
# case 1
psiLj.case1.tau1$true <- as.numeric(true.value.case1.tau1[, c("rmst.diff" ,rep(c("rmst.diff","rmst.diff","rmst.ow.diff","rmst.mw.diff","rmst.enw.diff"), 8))])
psiLj.case1.tau1$bias <- apply(f.sim.case1$tau1[,grep("estimate.diff.v",colnames(f.sim.case1$tau1))], 1, mean) - psiLj.case1.tau1$true
psiLj.case1.tau1$ese <- apply(f.sim.case1$tau1[,grep("estimate.diff.v",colnames(f.sim.case1$tau1))], 1, sd)
psiLj.case1.tau1$ase <- apply(f.sim.case1$tau1[,grep("estimate.diff.se",colnames(f.sim.case1$tau1))], 1, mean)
psiLj.case1.tau1$cp <- apply(f.sim.case1$tau1[,grep("estimate.diff.v",colnames(f.sim.case1$tau1))]+qnorm(0.05/2)*f.sim.case1$tau1[,grep("estimate.diff.se",colnames(f.sim.case1$tau1))] < psiLj.case1.tau1$true
                             & f.sim.case1$tau1[,grep("estimate.diff.v",colnames(f.sim.case1$tau1))]+qnorm(1-0.05/2)*f.sim.case1$tau1[,grep("estimate.diff.se",colnames(f.sim.case1$tau1))] > psiLj.case1.tau1$true, 1,sum)/nsim

psiLj.case1.tau2$true <- as.numeric(true.value.case1.tau2[, c("rmst.diff" ,rep(c("rmst.diff","rmst.diff","rmst.ow.diff","rmst.mw.diff","rmst.enw.diff"), 8))])
psiLj.case1.tau2$bias <- apply(f.sim.case1$tau2[,grep("estimate.diff.v",colnames(f.sim.case1$tau2))], 1, mean) - psiLj.case1.tau2$true
psiLj.case1.tau2$ese <- apply(f.sim.case1$tau2[,grep("estimate.diff.v",colnames(f.sim.case1$tau2))], 1, sd)
psiLj.case1.tau2$ase <- apply(f.sim.case1$tau2[,grep("estimate.diff.se",colnames(f.sim.case1$tau2))], 1, mean)
psiLj.case1.tau2$cp <- apply(f.sim.case1$tau2[,grep("estimate.diff.v",colnames(f.sim.case1$tau2))]+qnorm(0.05/2)*f.sim.case1$tau2[,grep("estimate.diff.se",colnames(f.sim.case1$tau2))] < psiLj.case1.tau2$true
                             & f.sim.case1$tau2[,grep("estimate.diff.v",colnames(f.sim.case1$tau2))]+qnorm(1-0.05/2)*f.sim.case1$tau2[,grep("estimate.diff.se",colnames(f.sim.case1$tau2))] > psiLj.case1.tau2$true, 1,sum)/nsim

psiLj.case1.tau3$true <- as.numeric(true.value.case1.tau3[, c("rmst.diff" ,rep(c("rmst.diff","rmst.diff","rmst.ow.diff","rmst.mw.diff","rmst.enw.diff"), 8))])
psiLj.case1.tau3$bias <- apply(f.sim.case1$tau3[,grep("estimate.diff.v",colnames(f.sim.case1$tau3))], 1, mean) - psiLj.case1.tau3$true
psiLj.case1.tau3$ese <- apply(f.sim.case1$tau3[,grep("estimate.diff.v",colnames(f.sim.case1$tau3))], 1, sd)
psiLj.case1.tau3$ase <- apply(f.sim.case1$tau3[,grep("estimate.diff.se",colnames(f.sim.case1$tau3))], 1, mean)
psiLj.case1.tau3$cp <- apply(f.sim.case1$tau3[,grep("estimate.diff.v",colnames(f.sim.case1$tau3))]+qnorm(0.05/2)*f.sim.case1$tau3[,grep("estimate.diff.se",colnames(f.sim.case1$tau3))] < psiLj.case1.tau3$true
                             & f.sim.case1$tau3[,grep("estimate.diff.v",colnames(f.sim.case1$tau3))]+qnorm(1-0.05/2)*f.sim.case1$tau3[,grep("estimate.diff.se",colnames(f.sim.case1$tau3))] > psiLj.case1.tau3$true, 1,sum)/nsim

psiLj.case1.tau4$true <- as.numeric(true.value.case1.tau4[, c("rmst.diff" ,rep(c("rmst.diff","rmst.diff","rmst.ow.diff","rmst.mw.diff","rmst.enw.diff"), 8))])
psiLj.case1.tau4$bias <- apply(f.sim.case1$tau4[,grep("estimate.diff.v",colnames(f.sim.case1$tau4))], 1, mean) - psiLj.case1.tau4$true
psiLj.case1.tau4$ese <- apply(f.sim.case1$tau4[,grep("estimate.diff.v",colnames(f.sim.case1$tau4))], 1, sd)
psiLj.case1.tau4$ase <- apply(f.sim.case1$tau4[,grep("estimate.diff.se",colnames(f.sim.case1$tau4))], 1, mean)
psiLj.case1.tau4$cp <- apply(f.sim.case1$tau4[,grep("estimate.diff.v",colnames(f.sim.case1$tau4))]+qnorm(0.05/2)*f.sim.case1$tau4[,grep("estimate.diff.se",colnames(f.sim.case1$tau4))] < psiLj.case1.tau4$true
                             & f.sim.case1$tau4[,grep("estimate.diff.v",colnames(f.sim.case1$tau4))]+qnorm(1-0.05/2)*f.sim.case1$tau4[,grep("estimate.diff.se",colnames(f.sim.case1$tau4))] > psiLj.case1.tau4$true, 1,sum)/nsim

# case 2
psiLj.case2.tau1$true <- as.numeric(true.value.case2.tau1[, c("rmst.diff" ,rep(c("rmst.diff","rmst.diff","rmst.ow.diff","rmst.mw.diff","rmst.enw.diff"), 8))])
psiLj.case2.tau1$bias <- apply(f.sim.case2$tau1[,grep("estimate.diff.v",colnames(f.sim.case2$tau1))], 1, mean) - psiLj.case2.tau1$true
psiLj.case2.tau1$ese <- apply(f.sim.case2$tau1[,grep("estimate.diff.v",colnames(f.sim.case2$tau1))], 1, sd)
psiLj.case2.tau1$ase <- apply(f.sim.case2$tau1[,grep("estimate.diff.se",colnames(f.sim.case2$tau1))], 1, mean)
psiLj.case2.tau1$cp <- apply(f.sim.case2$tau1[,grep("estimate.diff.v",colnames(f.sim.case2$tau1))]+qnorm(0.05/2)*f.sim.case2$tau1[,grep("estimate.diff.se",colnames(f.sim.case2$tau1))] < psiLj.case2.tau1$true
                             & f.sim.case2$tau1[,grep("estimate.diff.v",colnames(f.sim.case2$tau1))]+qnorm(1-0.05/2)*f.sim.case2$tau1[,grep("estimate.diff.se",colnames(f.sim.case2$tau1))] > psiLj.case2.tau1$true, 1,sum)/nsim

psiLj.case2.tau2$true <- as.numeric(true.value.case2.tau2[, c("rmst.diff" ,rep(c("rmst.diff","rmst.diff","rmst.ow.diff","rmst.mw.diff","rmst.enw.diff"), 8))])
psiLj.case2.tau2$bias <- apply(f.sim.case2$tau2[,grep("estimate.diff.v",colnames(f.sim.case2$tau2))], 1, mean) - psiLj.case2.tau2$true
psiLj.case2.tau2$ese <- apply(f.sim.case2$tau2[,grep("estimate.diff.v",colnames(f.sim.case2$tau2))], 1, sd)
psiLj.case2.tau2$ase <- apply(f.sim.case2$tau2[,grep("estimate.diff.se",colnames(f.sim.case2$tau2))], 1, mean)
psiLj.case2.tau2$cp <- apply(f.sim.case2$tau2[,grep("estimate.diff.v",colnames(f.sim.case2$tau2))]+qnorm(0.05/2)*f.sim.case2$tau2[,grep("estimate.diff.se",colnames(f.sim.case2$tau2))] < psiLj.case2.tau2$true
                             & f.sim.case2$tau2[,grep("estimate.diff.v",colnames(f.sim.case2$tau2))]+qnorm(1-0.05/2)*f.sim.case2$tau2[,grep("estimate.diff.se",colnames(f.sim.case2$tau2))] > psiLj.case2.tau2$true, 1,sum)/nsim

psiLj.case2.tau3$true <- as.numeric(true.value.case2.tau3[, c("rmst.diff" ,rep(c("rmst.diff","rmst.diff","rmst.ow.diff","rmst.mw.diff","rmst.enw.diff"), 8))])
psiLj.case2.tau3$bias <- apply(f.sim.case2$tau3[,grep("estimate.diff.v",colnames(f.sim.case2$tau3))], 1, mean) - psiLj.case2.tau3$true
psiLj.case2.tau3$ese <- apply(f.sim.case2$tau3[,grep("estimate.diff.v",colnames(f.sim.case2$tau3))], 1, sd)
psiLj.case2.tau3$ase <- apply(f.sim.case2$tau3[,grep("estimate.diff.se",colnames(f.sim.case2$tau3))], 1, mean)
psiLj.case2.tau3$cp <- apply(f.sim.case2$tau3[,grep("estimate.diff.v",colnames(f.sim.case2$tau3))]+qnorm(0.05/2)*f.sim.case2$tau3[,grep("estimate.diff.se",colnames(f.sim.case2$tau3))] < psiLj.case2.tau3$true
                             & f.sim.case2$tau3[,grep("estimate.diff.v",colnames(f.sim.case2$tau3))]+qnorm(1-0.05/2)*f.sim.case2$tau3[,grep("estimate.diff.se",colnames(f.sim.case2$tau3))] > psiLj.case2.tau3$true, 1,sum)/nsim

psiLj.case2.tau4$true <- as.numeric(true.value.case2.tau4[, c("rmst.diff" ,rep(c("rmst.diff","rmst.diff","rmst.ow.diff","rmst.mw.diff","rmst.enw.diff"), 8))])
psiLj.case2.tau4$bias <- apply(f.sim.case2$tau4[,grep("estimate.diff.v",colnames(f.sim.case2$tau4))], 1, mean) - psiLj.case2.tau4$true
psiLj.case2.tau4$ese <- apply(f.sim.case2$tau4[,grep("estimate.diff.v",colnames(f.sim.case2$tau4))], 1, sd)
psiLj.case2.tau4$ase <- apply(f.sim.case2$tau4[,grep("estimate.diff.se",colnames(f.sim.case2$tau4))], 1, mean)
psiLj.case2.tau4$cp <- apply(f.sim.case2$tau4[,grep("estimate.diff.v",colnames(f.sim.case2$tau4))]+qnorm(0.05/2)*f.sim.case2$tau4[,grep("estimate.diff.se",colnames(f.sim.case2$tau4))] < psiLj.case2.tau4$true
                             & f.sim.case2$tau4[,grep("estimate.diff.v",colnames(f.sim.case2$tau4))]+qnorm(1-0.05/2)*f.sim.case2$tau4[,grep("estimate.diff.se",colnames(f.sim.case2$tau4))] > psiLj.case2.tau4$true, 1,sum)/nsim

# case 3
psiLj.case3.tau1$true <- as.numeric(true.value.case3.tau1[, c("rmtlj.diff" ,rep(c("rmtlj.diff","rmtlj.diff","rmtlj.ow.diff","rmtlj.mw.diff","rmtlj.enw.diff"), 8))])
psiLj.case3.tau1$bias <- apply(f.sim.case3$tau1[,grep("estimate.diff.v",colnames(f.sim.case3$tau1))], 1, mean) - psiLj.case3.tau1$true
psiLj.case3.tau1$ese <- apply(f.sim.case3$tau1[,grep("estimate.diff.v",colnames(f.sim.case3$tau1))], 1, sd)
psiLj.case3.tau1$ase <- apply(f.sim.case3$tau1[,grep("estimate.diff.se",colnames(f.sim.case3$tau1))], 1, mean)
psiLj.case3.tau1$cp <- apply(f.sim.case3$tau1[,grep("estimate.diff.v",colnames(f.sim.case3$tau1))]+qnorm(0.05/2)*f.sim.case3$tau1[,grep("estimate.diff.se",colnames(f.sim.case3$tau1))] < psiLj.case3.tau1$true
                             & f.sim.case3$tau1[,grep("estimate.diff.v",colnames(f.sim.case3$tau1))]+qnorm(1-0.05/2)*f.sim.case3$tau1[,grep("estimate.diff.se",colnames(f.sim.case3$tau1))] > psiLj.case3.tau1$true, 1,sum)/nsim

psiLj.case3.tau2$true <- as.numeric(true.value.case3.tau2[, c("rmtlj.diff" ,rep(c("rmtlj.diff","rmtlj.diff","rmtlj.ow.diff","rmtlj.mw.diff","rmtlj.enw.diff"), 8))])
psiLj.case3.tau2$bias <- apply(f.sim.case3$tau2[,grep("estimate.diff.v",colnames(f.sim.case3$tau2))], 1, mean) - psiLj.case3.tau2$true
psiLj.case3.tau2$ese <- apply(f.sim.case3$tau2[,grep("estimate.diff.v",colnames(f.sim.case3$tau2))], 1, sd)
psiLj.case3.tau2$ase <- apply(f.sim.case3$tau2[,grep("estimate.diff.se",colnames(f.sim.case3$tau2))], 1, mean)
psiLj.case3.tau2$cp <- apply(f.sim.case3$tau2[,grep("estimate.diff.v",colnames(f.sim.case3$tau2))]+qnorm(0.05/2)*f.sim.case3$tau2[,grep("estimate.diff.se",colnames(f.sim.case3$tau2))] < psiLj.case3.tau2$true
                             & f.sim.case3$tau2[,grep("estimate.diff.v",colnames(f.sim.case3$tau2))]+qnorm(1-0.05/2)*f.sim.case3$tau2[,grep("estimate.diff.se",colnames(f.sim.case3$tau2))] > psiLj.case3.tau2$true, 1,sum)/nsim

psiLj.case3.tau3$true <- as.numeric(true.value.case3.tau3[, c("rmtlj.diff" ,rep(c("rmtlj.diff","rmtlj.diff","rmtlj.ow.diff","rmtlj.mw.diff","rmtlj.enw.diff"), 8))])
psiLj.case3.tau3$bias <- apply(f.sim.case3$tau3[,grep("estimate.diff.v",colnames(f.sim.case3$tau3))], 1, mean) - psiLj.case3.tau3$true
psiLj.case3.tau3$ese <- apply(f.sim.case3$tau3[,grep("estimate.diff.v",colnames(f.sim.case3$tau3))], 1, sd)
psiLj.case3.tau3$ase <- apply(f.sim.case3$tau3[,grep("estimate.diff.se",colnames(f.sim.case3$tau3))], 1, mean)
psiLj.case3.tau3$cp <- apply(f.sim.case3$tau3[,grep("estimate.diff.v",colnames(f.sim.case3$tau3))]+qnorm(0.05/2)*f.sim.case3$tau3[,grep("estimate.diff.se",colnames(f.sim.case3$tau3))] < psiLj.case3.tau3$true
                             & f.sim.case3$tau3[,grep("estimate.diff.v",colnames(f.sim.case3$tau3))]+qnorm(1-0.05/2)*f.sim.case3$tau3[,grep("estimate.diff.se",colnames(f.sim.case3$tau3))] > psiLj.case3.tau3$true, 1,sum)/nsim

psiLj.case3.tau4$true <- as.numeric(true.value.case3.tau4[, c("rmtlj.diff" ,rep(c("rmtlj.diff","rmtlj.diff","rmtlj.ow.diff","rmtlj.mw.diff","rmtlj.enw.diff"), 8))])
psiLj.case3.tau4$bias <- apply(f.sim.case3$tau4[,grep("estimate.diff.v",colnames(f.sim.case3$tau4))], 1, mean) - psiLj.case3.tau4$true
psiLj.case3.tau4$ese <- apply(f.sim.case3$tau4[,grep("estimate.diff.v",colnames(f.sim.case3$tau4))], 1, sd)
psiLj.case3.tau4$ase <- apply(f.sim.case3$tau4[,grep("estimate.diff.se",colnames(f.sim.case3$tau4))], 1, mean)
psiLj.case3.tau4$cp <- apply(f.sim.case3$tau4[,grep("estimate.diff.v",colnames(f.sim.case3$tau4))]+qnorm(0.05/2)*f.sim.case3$tau4[,grep("estimate.diff.se",colnames(f.sim.case3$tau4))] < psiLj.case3.tau4$true
                             & f.sim.case3$tau4[,grep("estimate.diff.v",colnames(f.sim.case3$tau4))]+qnorm(1-0.05/2)*f.sim.case3$tau4[,grep("estimate.diff.se",colnames(f.sim.case3$tau4))] > psiLj.case3.tau4$true, 1,sum)/nsim

# case 4
psiLj.case4.tau1$true <- as.numeric(true.value.case4.tau1[, c("rmtlj.diff" ,rep(c("rmtlj.diff","rmtlj.diff","rmtlj.ow.diff","rmtlj.mw.diff","rmtlj.enw.diff"), 8))])
psiLj.case4.tau1$bias <- apply(f.sim.case4$tau1[,grep("estimate.diff.v",colnames(f.sim.case4$tau1))], 1, mean) - psiLj.case4.tau1$true
psiLj.case4.tau1$ese <- apply(f.sim.case4$tau1[,grep("estimate.diff.v",colnames(f.sim.case4$tau1))], 1, sd)
psiLj.case4.tau1$ase <- apply(f.sim.case4$tau1[,grep("estimate.diff.se",colnames(f.sim.case4$tau1))], 1, mean)
psiLj.case4.tau1$cp <- apply(f.sim.case4$tau1[,grep("estimate.diff.v",colnames(f.sim.case4$tau1))]+qnorm(0.05/2)*f.sim.case4$tau1[,grep("estimate.diff.se",colnames(f.sim.case4$tau1))] < psiLj.case4.tau1$true
                             & f.sim.case4$tau1[,grep("estimate.diff.v",colnames(f.sim.case4$tau1))]+qnorm(1-0.05/2)*f.sim.case4$tau1[,grep("estimate.diff.se",colnames(f.sim.case4$tau1))] > psiLj.case4.tau1$true, 1,sum)/nsim

psiLj.case4.tau2$true <- as.numeric(true.value.case4.tau2[, c("rmtlj.diff" ,rep(c("rmtlj.diff","rmtlj.diff","rmtlj.ow.diff","rmtlj.mw.diff","rmtlj.enw.diff"), 8))])
psiLj.case4.tau2$bias <- apply(f.sim.case4$tau2[,grep("estimate.diff.v",colnames(f.sim.case4$tau2))], 1, mean) - psiLj.case4.tau2$true
psiLj.case4.tau2$ese <- apply(f.sim.case4$tau2[,grep("estimate.diff.v",colnames(f.sim.case4$tau2))], 1, sd)
psiLj.case4.tau2$ase <- apply(f.sim.case4$tau2[,grep("estimate.diff.se",colnames(f.sim.case4$tau2))], 1, mean)
psiLj.case4.tau2$cp <- apply(f.sim.case4$tau2[,grep("estimate.diff.v",colnames(f.sim.case4$tau2))]+qnorm(0.05/2)*f.sim.case4$tau2[,grep("estimate.diff.se",colnames(f.sim.case4$tau2))] < psiLj.case4.tau2$true
                             & f.sim.case4$tau2[,grep("estimate.diff.v",colnames(f.sim.case4$tau2))]+qnorm(1-0.05/2)*f.sim.case4$tau2[,grep("estimate.diff.se",colnames(f.sim.case4$tau2))] > psiLj.case4.tau2$true, 1,sum)/nsim

psiLj.case4.tau3$true <- as.numeric(true.value.case4.tau3[, c("rmtlj.diff" ,rep(c("rmtlj.diff","rmtlj.diff","rmtlj.ow.diff","rmtlj.mw.diff","rmtlj.enw.diff"), 8))])
psiLj.case4.tau3$bias <- apply(f.sim.case4$tau3[,grep("estimate.diff.v",colnames(f.sim.case4$tau3))], 1, mean) - psiLj.case4.tau3$true
psiLj.case4.tau3$ese <- apply(f.sim.case4$tau3[,grep("estimate.diff.v",colnames(f.sim.case4$tau3))], 1, sd)
psiLj.case4.tau3$ase <- apply(f.sim.case4$tau3[,grep("estimate.diff.se",colnames(f.sim.case4$tau3))], 1, mean)
psiLj.case4.tau3$cp <- apply(f.sim.case4$tau3[,grep("estimate.diff.v",colnames(f.sim.case4$tau3))]+qnorm(0.05/2)*f.sim.case4$tau3[,grep("estimate.diff.se",colnames(f.sim.case4$tau3))] < psiLj.case4.tau3$true
                             & f.sim.case4$tau3[,grep("estimate.diff.v",colnames(f.sim.case4$tau3))]+qnorm(1-0.05/2)*f.sim.case4$tau3[,grep("estimate.diff.se",colnames(f.sim.case4$tau3))] > psiLj.case4.tau3$true, 1,sum)/nsim

psiLj.case4.tau4$true <- as.numeric(true.value.case4.tau4[, c("rmtlj.diff" ,rep(c("rmtlj.diff","rmtlj.diff","rmtlj.ow.diff","rmtlj.mw.diff","rmtlj.enw.diff"), 8))])
psiLj.case4.tau4$bias <- apply(f.sim.case4$tau4[,grep("estimate.diff.v",colnames(f.sim.case4$tau4))], 1, mean) - psiLj.case4.tau4$true
psiLj.case4.tau4$ese <- apply(f.sim.case4$tau4[,grep("estimate.diff.v",colnames(f.sim.case4$tau4))], 1, sd)
psiLj.case4.tau4$ase <- apply(f.sim.case4$tau4[,grep("estimate.diff.se",colnames(f.sim.case4$tau4))], 1, mean)
psiLj.case4.tau4$cp <- apply(f.sim.case4$tau4[,grep("estimate.diff.v",colnames(f.sim.case4$tau4))]+qnorm(0.05/2)*f.sim.case4$tau4[,grep("estimate.diff.se",colnames(f.sim.case4$tau4))] < psiLj.case4.tau4$true
                             & f.sim.case4$tau4[,grep("estimate.diff.v",colnames(f.sim.case4$tau4))]+qnorm(1-0.05/2)*f.sim.case4$tau4[,grep("estimate.diff.se",colnames(f.sim.case4$tau4))] > psiLj.case4.tau4$true, 1,sum)/nsim


##  bias
# case 1
psiLj.bias.case1.tau1 <- ggplot(psiLj.case1.tau1, aes(x=fct_inorder(rownames(psiLj.case1.tau1)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case1.tau2 <- ggplot(psiLj.case1.tau2, aes(x=fct_inorder(rownames(psiLj.case1.tau2)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case1.tau3 <- ggplot(psiLj.case1.tau3, aes(x=fct_inorder(rownames(psiLj.case1.tau3)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

# case 2
psiLj.bias.case2.tau1 <- ggplot(psiLj.case2.tau1, aes(x=fct_inorder(rownames(psiLj.case2.tau1)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case2.tau2 <- ggplot(psiLj.case2.tau2, aes(x=fct_inorder(rownames(psiLj.case2.tau2)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case2.tau3 <- ggplot(psiLj.case2.tau3, aes(x=fct_inorder(rownames(psiLj.case2.tau3)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

# case 3
psiLj.bias.case3.tau1 <- ggplot(psiLj.case3.tau1, aes(x=fct_inorder(rownames(psiLj.case3.tau1)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case3.tau2 <- ggplot(psiLj.case3.tau2, aes(x=fct_inorder(rownames(psiLj.case3.tau2)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case3.tau3 <- ggplot(psiLj.case3.tau3, aes(x=fct_inorder(rownames(psiLj.case3.tau3)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

# case 4
psiLj.bias.case4.tau1 <- ggplot(psiLj.case4.tau1, aes(x=fct_inorder(rownames(psiLj.case4.tau1)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case4.tau2 <- ggplot(psiLj.case4.tau2, aes(x=fct_inorder(rownames(psiLj.case4.tau2)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case4.tau3 <- ggplot(psiLj.case4.tau3, aes(x=fct_inorder(rownames(psiLj.case4.tau3)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

psiLj.bias.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), bias))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Bias")+
  scale_y_continuous(limits=c(-0.1, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1))

ggsave("sim_bias.pdf", width = 20, height = 20,
       plot_grid(psiLj.bias.case1.tau1, psiLj.bias.case1.tau2, psiLj.bias.case1.tau3, psiLj.bias.case1.tau4,
                 psiLj.bias.case2.tau1, psiLj.bias.case2.tau2, psiLj.bias.case2.tau3, psiLj.bias.case2.tau4,
                 psiLj.bias.case3.tau1, psiLj.bias.case3.tau2, psiLj.bias.case3.tau3, psiLj.bias.case3.tau4,
                 psiLj.bias.case4.tau1, psiLj.bias.case4.tau2, psiLj.bias.case4.tau3, psiLj.bias.case4.tau4,
                 labels=c("(a)","","","","(b)","","","","(c)","","","","(d)","","",""), vjust = 1.5, ncol=4, nrow=4))

# ese
# case 1
psiLj.ese.case1.tau1 <- ggplot(psiLj.case1.tau1, aes(x=fct_inorder(rownames(psiLj.case1.tau1)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case1.tau2 <- ggplot(psiLj.case1.tau2, aes(x=fct_inorder(rownames(psiLj.case1.tau2)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case1.tau3 <- ggplot(psiLj.case1.tau3, aes(x=fct_inorder(rownames(psiLj.case1.tau3)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# case 2
psiLj.ese.case2.tau1 <- ggplot(psiLj.case2.tau1, aes(x=fct_inorder(rownames(psiLj.case2.tau1)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case2.tau2 <- ggplot(psiLj.case2.tau2, aes(x=fct_inorder(rownames(psiLj.case2.tau2)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case2.tau3 <- ggplot(psiLj.case2.tau3, aes(x=fct_inorder(rownames(psiLj.case2.tau3)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# case 3
psiLj.ese.case3.tau1 <- ggplot(psiLj.case3.tau1, aes(x=fct_inorder(rownames(psiLj.case3.tau1)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case3.tau2 <- ggplot(psiLj.case3.tau2, aes(x=fct_inorder(rownames(psiLj.case3.tau2)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case3.tau3 <- ggplot(psiLj.case3.tau3, aes(x=fct_inorder(rownames(psiLj.case3.tau3)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# case 4
psiLj.ese.case4.tau1 <- ggplot(psiLj.case4.tau1, aes(x=fct_inorder(rownames(psiLj.case4.tau1)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case4.tau2 <- ggplot(psiLj.case4.tau2, aes(x=fct_inorder(rownames(psiLj.case4.tau2)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case4.tau3 <- ggplot(psiLj.case4.tau3, aes(x=fct_inorder(rownames(psiLj.case4.tau3)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ese.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), ese))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Monte Carlo SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

ggsave("sim_ese.pdf", width = 20, height = 20,
       plot_grid(psiLj.ese.case1.tau1, psiLj.ese.case1.tau2, psiLj.ese.case1.tau3, psiLj.ese.case1.tau4,
                 psiLj.ese.case2.tau1, psiLj.ese.case2.tau2, psiLj.ese.case2.tau3, psiLj.ese.case2.tau4,
                 psiLj.ese.case3.tau1, psiLj.ese.case3.tau2, psiLj.ese.case3.tau3, psiLj.ese.case3.tau4,
                 psiLj.ese.case4.tau1, psiLj.ese.case4.tau2, psiLj.ese.case4.tau3, psiLj.ese.case4.tau4,
                 labels=c("(a)","","","","(b)","","","","(c)","","","","(d)","","",""), vjust = 1.5, ncol=4, nrow=4))

# ase
# case 1
psiLj.ase.case1.tau1 <- ggplot(psiLj.case1.tau1, aes(x=fct_inorder(rownames(psiLj.case1.tau1)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case1.tau2 <- ggplot(psiLj.case1.tau2, aes(x=fct_inorder(rownames(psiLj.case1.tau2)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case1.tau3 <- ggplot(psiLj.case1.tau3, aes(x=fct_inorder(rownames(psiLj.case1.tau3)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# case 2
psiLj.ase.case2.tau1 <- ggplot(psiLj.case2.tau1, aes(x=fct_inorder(rownames(psiLj.case2.tau1)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case2.tau2 <- ggplot(psiLj.case2.tau2, aes(x=fct_inorder(rownames(psiLj.case2.tau2)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case2.tau3 <- ggplot(psiLj.case2.tau3, aes(x=fct_inorder(rownames(psiLj.case2.tau3)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# case 3
psiLj.ase.case3.tau1 <- ggplot(psiLj.case3.tau1, aes(x=fct_inorder(rownames(psiLj.case3.tau1)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case3.tau2 <- ggplot(psiLj.case3.tau2, aes(x=fct_inorder(rownames(psiLj.case3.tau2)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case3.tau3 <- ggplot(psiLj.case3.tau3, aes(x=fct_inorder(rownames(psiLj.case3.tau3)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

# case 4
psiLj.ase.case4.tau1 <- ggplot(psiLj.case4.tau1, aes(x=fct_inorder(rownames(psiLj.case4.tau1)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case4.tau2 <- ggplot(psiLj.case4.tau2, aes(x=fct_inorder(rownames(psiLj.case4.tau2)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case4.tau3 <- ggplot(psiLj.case4.tau3, aes(x=fct_inorder(rownames(psiLj.case4.tau3)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

psiLj.ase.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), ase))+
  geom_point() +
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Estimated SE")+
  scale_y_continuous(limits=c(0, 0.4), breaks=c(0, 0.1, 0.2, 0.3, 0.4))

ggsave("sim_ase.pdf", width = 20, height = 20,
       plot_grid(psiLj.ase.case1.tau1, psiLj.ase.case1.tau2, psiLj.ase.case1.tau3, psiLj.ase.case1.tau4,
                 psiLj.ase.case2.tau1, psiLj.ase.case2.tau2, psiLj.ase.case2.tau3, psiLj.ase.case2.tau4,
                 psiLj.ase.case3.tau1, psiLj.ase.case3.tau2, psiLj.ase.case3.tau3, psiLj.ase.case3.tau4,
                 psiLj.ase.case4.tau1, psiLj.ase.case4.tau2, psiLj.ase.case4.tau3, psiLj.ase.case4.tau4,
                 labels=c("(a)","","","","(b)","","","","(c)","","","","(d)","","",""), vjust = 1.5, ncol=4, nrow=4))

# cp
# case 1
psiLj.cp.case1.tau1 <- ggplot(psiLj.case1.tau1, aes(x=fct_inorder(rownames(psiLj.case1.tau1)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case1.tau2 <- ggplot(psiLj.case1.tau2, aes(x=fct_inorder(rownames(psiLj.case1.tau2)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case1.tau3 <- ggplot(psiLj.case1.tau3, aes(x=fct_inorder(rownames(psiLj.case1.tau3)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case1.tau4 <- ggplot(psiLj.case1.tau4, aes(x=fct_inorder(rownames(psiLj.case1.tau4)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

# case 2
psiLj.cp.case2.tau1 <- ggplot(psiLj.case2.tau1, aes(x=fct_inorder(rownames(psiLj.case2.tau1)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case2.tau2 <- ggplot(psiLj.case2.tau2, aes(x=fct_inorder(rownames(psiLj.case2.tau2)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case2.tau3 <- ggplot(psiLj.case2.tau3, aes(x=fct_inorder(rownames(psiLj.case2.tau3)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case2.tau4 <- ggplot(psiLj.case2.tau4, aes(x=fct_inorder(rownames(psiLj.case2.tau4)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

# case 3
psiLj.cp.case3.tau1 <- ggplot(psiLj.case3.tau1, aes(x=fct_inorder(rownames(psiLj.case3.tau1)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case3.tau2 <- ggplot(psiLj.case3.tau2, aes(x=fct_inorder(rownames(psiLj.case3.tau2)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case3.tau3 <- ggplot(psiLj.case3.tau3, aes(x=fct_inorder(rownames(psiLj.case3.tau3)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case3.tau4 <- ggplot(psiLj.case3.tau4, aes(x=fct_inorder(rownames(psiLj.case3.tau4)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

# case 4
psiLj.cp.case4.tau1 <- ggplot(psiLj.case4.tau1, aes(x=fct_inorder(rownames(psiLj.case4.tau1)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case4.tau2 <- ggplot(psiLj.case4.tau2, aes(x=fct_inorder(rownames(psiLj.case4.tau2)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case4.tau3 <- ggplot(psiLj.case4.tau3, aes(x=fct_inorder(rownames(psiLj.case4.tau3)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

psiLj.cp.case4.tau4 <- ggplot(psiLj.case4.tau4, aes(x=fct_inorder(rownames(psiLj.case4.tau4)), cp))+
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Estimator")+
  ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))

ggsave("sim_cp.pdf", width = 20, height = 20,
       plot_grid(psiLj.cp.case1.tau1, psiLj.cp.case1.tau2, psiLj.cp.case1.tau3, psiLj.cp.case1.tau4,
                 psiLj.cp.case2.tau1, psiLj.cp.case2.tau2, psiLj.cp.case2.tau3, psiLj.cp.case2.tau4,
                 psiLj.cp.case3.tau1, psiLj.cp.case3.tau2, psiLj.cp.case3.tau3, psiLj.cp.case3.tau4,
                 psiLj.cp.case4.tau1, psiLj.cp.case4.tau2, psiLj.cp.case4.tau3, psiLj.cp.case4.tau4,
                 labels=c("(a)","","","","(b)","","","","(c)","","","","(d)","","",""), vjust = 1.5, ncol=4, nrow=4))

# overlap plot
admin.cens <- 10
time.grid <- 0.01
taus <- c(1,2,3,4)

# case 1
dgp4_survival_ph <- function(n, admin.cens, time.grid){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=n, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(n,1,0.5)
  x5 <- rbinom(n,1,0.5)
  x6 <- rbinom(n,1,0.5)
  ps <- 1/(1+exp(-(0.3+0.2*x1+0.3*x2+0.3*x3-0.2*x4-0.3*x5-0.2*x6))) # 1/(1+exp(-(-1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6)))
  a <- rbinom(n, 1, ps) # -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6, ps)
  
  data.df$Ta0 <- simsurv(dist = "exponential", lambdas = 0.12, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))$eventtime
  data.df$Ta1 <- simsurv(dist = "exponential", lambdas = 0.15, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))$eventtime
  data.df$Ta <- (1-a)*data.df$Ta0 + a*data.df$Ta1
  
  data.df$C[data.df$a==0] <- simsurv(dist = "exponential", lambdas = 0.06, x = cbind(x0=1, data.df[data.df$a==0, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0.1, x1=0.4, x2=-0.7, x3=-0.4, x4=-0.5, x5=0.8, x6=-0.6))$eventtime
  data.df$C[data.df$a==1] <- simsurv(dist = "exponential", lambdas = 0.08, x = cbind(x0=1, data.df[data.df$a==1, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0, x1=0.5, x2=-0.6, x3=0.2, x4=0.6, x5=0.9, x6=-0.5))$eventtime
  data.df$C <- pmin(data.df$C, admin.cens, na.rm=TRUE)
  data.df$end.time <- ceiling(pmin(data.df$C, data.df$Ta)/time.grid)*time.grid # , data.df$Tj2, ceil(pmin(data.df$C,data.df$Tj1,data.df$Tj2)*365.25) # round(pmin(data.df$C,data.df$Tj1,data.df$Tj2),digits=2)
  data.df$event <- with(data.df,ifelse(Ta<=C, 1, 0))
  data.df <- data.df[order(data.df$end.time,-data.df$event),]
  data.df$id <- 1:n
  return(data.df)
}

n <- 2000
sim1.df <- dgp4_survival_ph(2*n, admin.cens, time.grid)
sum(sim1.df$a)/(2*n)
sum(sim1.df$event==0)/(2*n)
sum(sim1.df$event==1 & sim1.df$end.time < taus[1])/(2*n)
sum(sim1.df$event==1 & sim1.df$end.time < taus[2])/(2*n)
sum(sim1.df$event==1 & sim1.df$end.time < taus[3])/(2*n)
sum(sim1.df$event==1 & sim1.df$end.time < taus[4])/(2*n)

# case 2
dgp4_survival_ph <- function(n, admin.cens, time.grid){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=n, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(n,1,0.5)
  x5 <- rbinom(n,1,0.5)
  x6 <- rbinom(n,1,0.5)
  ps <- 1/(1+exp(-(-1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6)))
  a <- rbinom(n, 1, ps) # -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6, ps)
  
  data.df$Ta0 <- simsurv(dist = "exponential", lambdas = 0.12, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))$eventtime
  data.df$Ta1 <- simsurv(dist = "exponential", lambdas = 0.15, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))$eventtime
  data.df$Ta <- (1-a)*data.df$Ta0 + a*data.df$Ta1
  
  data.df$C[data.df$a==0] <- simsurv(dist = "exponential", lambdas = 0.06, x = cbind(x0=1, data.df[data.df$a==0, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0.1, x1=0.4, x2=-0.7, x3=-0.4, x4=-0.5, x5=0.8, x6=-0.6))$eventtime
  data.df$C[data.df$a==1] <- simsurv(dist = "exponential", lambdas = 0.08, x = cbind(x0=1, data.df[data.df$a==1, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0, x1=0.5, x2=-0.6, x3=0.2, x4=0.6, x5=0.9, x6=-0.5))$eventtime
  data.df$C <- pmin(data.df$C, admin.cens, na.rm=TRUE)
  data.df$end.time <- ceiling(pmin(data.df$C, data.df$Ta)/time.grid)*time.grid # , data.df$Tj2, ceil(pmin(data.df$C,data.df$Tj1,data.df$Tj2)*365.25) # round(pmin(data.df$C,data.df$Tj1,data.df$Tj2),digits=2)
  data.df$event <- with(data.df,ifelse(Ta<=C, 1, 0))
  data.df <- data.df[order(data.df$end.time,-data.df$event),]
  data.df$id <- 1:n
  return(data.df)
}

n <- 2000
sim2.df <- dgp4_survival_ph(2*n, admin.cens, time.grid)
sum(sim2.df$a)/(2*n)
sum(sim2.df$event==0)/(2*n)
sum(sim2.df$event==1 & sim2.df$end.time < taus[1])/(2*n)
sum(sim2.df$event==1 & sim2.df$end.time < taus[2])/(2*n)
sum(sim2.df$event==1 & sim2.df$end.time < taus[3])/(2*n)
sum(sim2.df$event==1 & sim2.df$end.time < taus[4])/(2*n)

# case 3
dgp4_competing <- function(n, admin.cens, time.grid){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=n, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(n,1,0.5)
  x5 <- rbinom(n,1,0.5)
  x6 <- rbinom(n,1,0.5)
  ps <- 1/(1+exp(-(0.3+0.2*x1+0.3*x2+0.3*x3-0.2*x4-0.3*x5-0.2*x6))) # 1/(1+exp(-(-1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6)))
  a <- rbinom(n, 1, ps) # -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6, ps)
  
  data.df$Tj1a0 <- simsurv(dist = "exponential", lambdas = 0.12, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))$eventtime
  data.df$Tj1a1 <- simsurv(dist = "exponential", lambdas = 0.15, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))$eventtime
  data.df$Tj1 <- (1-a)*data.df$Tj1a0 + a*data.df$Tj1a1
  
  data.df$Tj2a0 <- simsurv(dist = "exponential", lambdas = 0.1, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.12, x1=-0.1, x2=0.3, x3=0.1, x4=0.2, x5=-0.4, x6=0.5))$eventtime
  data.df$Tj2a1 <- simsurv(dist = "exponential", lambdas = 0.08, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=-0.2, x2=-0.1, x3=0.2, x4=0.3, x5=0.3, x6=-0.3))$eventtime
  data.df$Tj2 <- (1-a)*data.df$Tj2a0 + a*data.df$Tj2a1
  
  data.df$C[data.df$a==0] <- simsurv(dist = "exponential", lambdas = 0.12, x = cbind(x0=1, data.df[data.df$a==0, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0.1, x1=0.4, x2=-0.7, x3=-0.4, x4=-0.5, x5=0.8, x6=-0.6))$eventtime
  data.df$C[data.df$a==1] <- simsurv(dist = "exponential", lambdas = 0.14, x = cbind(x0=1, data.df[data.df$a==1, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0, x1=0.5, x2=-0.6, x3=0.2, x4=0.6, x5=0.9, x6=-0.5))$eventtime
  data.df$C <- pmin(data.df$C, admin.cens, na.rm=TRUE)
  data.df$end.time <- ceiling(pmin(data.df$C, data.df$Tj1, data.df$Tj2)/time.grid)*time.grid # , data.df$Tj2, ceil(pmin(data.df$C,data.df$Tj1,data.df$Tj2)*365.25) # round(pmin(data.df$C,data.df$Tj1,data.df$Tj2),digits=2)
  data.df$event <- with(data.df, ifelse(pmin(Tj1,Tj2)>C,0, ifelse(Tj1<Tj2,1,2)))
  data.df <- data.df[order(data.df$end.time,-data.df$event),]
  data.df$id <- 1:n
  return(data.df)
}

n <- 2000
sim3.df <- dgp4_competing(2*n, admin.cens, time.grid)
sum(sim3.df$a)/(2*n)
sum(sim3.df$event==0)/(2*n)
sum(sim3.df$event==1 & sim3.df$end.time < taus[1])/(2*n)
sum(sim3.df$event==1 & sim3.df$end.time < taus[2])/(2*n)
sum(sim3.df$event==1 & sim3.df$end.time < taus[3])/(2*n)
sum(sim3.df$event==1 & sim3.df$end.time < taus[4])/(2*n)

# case 4
dgp4_competing <- function(n, admin.cens, time.grid){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=n, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(n,1,0.5)
  x5 <- rbinom(n,1,0.5)
  x6 <- rbinom(n,1,0.5)
  ps <- 1/(1+exp(-(-1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6))) # 1/(1+exp(-(0.3+0.2*x1+0.3*x2+0.3*x3-0.2*x4-0.3*x5-0.2*x6)))
  a <- rbinom(n, 1, ps) # -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6, ps)
  
  data.df$Tj1a0 <- simsurv(dist = "exponential", lambdas = 0.12, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))$eventtime
  data.df$Tj1a1 <- simsurv(dist = "exponential", lambdas = 0.15, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))$eventtime
  data.df$Tj1 <- (1-a)*data.df$Tj1a0 + a*data.df$Tj1a1
  
  data.df$Tj2a0 <- simsurv(dist = "exponential", lambdas = 0.1, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.12, x1=-0.1, x2=0.3, x3=0.1, x4=0.2, x5=-0.4, x6=0.5))$eventtime
  data.df$Tj2a1 <- simsurv(dist = "exponential", lambdas = 0.08, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=-0.2, x2=-0.1, x3=0.2, x4=0.3, x5=0.3, x6=-0.3))$eventtime
  data.df$Tj2 <- (1-a)*data.df$Tj2a0 + a*data.df$Tj2a1
  
  data.df$C[data.df$a==0] <- simsurv(dist = "exponential", lambdas = 0.12, x = cbind(x0=1, data.df[data.df$a==0, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0.1, x1=0.4, x2=-0.7, x3=-0.4, x4=-0.5, x5=0.8, x6=-0.6))$eventtime
  data.df$C[data.df$a==1] <- simsurv(dist = "exponential", lambdas = 0.14, x = cbind(x0=1, data.df[data.df$a==1, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0, x1=0.5, x2=-0.6, x3=0.2, x4=0.6, x5=0.9, x6=-0.5))$eventtime
  data.df$C <- pmin(data.df$C, admin.cens, na.rm=TRUE)
  data.df$end.time <- ceiling(pmin(data.df$C, data.df$Tj1, data.df$Tj2)/time.grid)*time.grid # , data.df$Tj2, ceil(pmin(data.df$C,data.df$Tj1,data.df$Tj2)*365.25) # round(pmin(data.df$C,data.df$Tj1,data.df$Tj2),digits=2)
  data.df$event <- with(data.df, ifelse(pmin(Tj1,Tj2)>C,0, ifelse(Tj1<Tj2,1,2)))
  data.df <- data.df[order(data.df$end.time,-data.df$event),]
  data.df$id <- 1:n
  return(data.df)
}

n <- 2000
sim4.df <- dgp4_competing(2*n, admin.cens, time.grid)
sum(sim4.df$a)/(2*n)
sum(sim4.df$event==0)/(2*n)
sum(sim4.df$event==1 & sim4.df$end.time < taus[1])/(2*n)
sum(sim4.df$event==1 & sim4.df$end.time < taus[2])/(2*n)
sum(sim4.df$event==1 & sim4.df$end.time < taus[3])/(2*n)
sum(sim4.df$event==1 & sim4.df$end.time < taus[4])/(2*n)

## plot
# case 1
sim1.ps <- ggplot(sim1.df, aes(x=ps, group=factor(a)))+
  geom_density(aes(linetype=ifelse(a, "dotted", "solid")))+
  xlab("True propensity score")+
  ylab("Density")+
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1))+
  labs(fill = "A:")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = guide_legend(title.position = "left"))+
  scale_y_continuous(limits=c(0, 15), breaks=c(0, 5, 10,15))+
  ggtitle("(a) Setting 1 and Setting 3")

# case 2
sim2.ps <- ggplot(sim2.df, aes(x=ps, group=factor(a)))+
  geom_density(aes(linetype=ifelse(a, "dotted", "solid")))+
  xlab("True propensity score")+
  ylab("Density")+
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1))+
  labs(fill = "A:")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = guide_legend(title.position = "left"))+
  scale_y_continuous(limits=c(0, 15), breaks=c(0, 5, 10,15))+
  ggtitle("(b) Setting 2 and Setting 4")

# case 3
sim3.ps <- ggplot(sim3.df, aes(x=ps, group=factor(a)))+
  geom_density(aes(linetype=ifelse(a, "dotted", "solid")))+
  xlab("True propensity score")+
  ylab("Density")+
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1))+
  labs(fill = "A:")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = guide_legend(title.position = "left"))+
  scale_y_continuous(limits=c(0, 15), breaks=c(0, 5, 10,15))

# case 4
sim4.ps <- ggplot(sim4.df, aes(x=ps, group=factor(a)))+
  geom_density(aes(linetype=ifelse(a, "dotted", "solid")))+
  xlab("True propensity score")+
  ylab("Density")+
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1))+
  labs(fill = "A:")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = guide_legend(title.position = "left"))+
  scale_y_continuous(limits=c(0, 15), breaks=c(0, 5, 10,15))

# ggsave("sim_ps.pdf", width = 12, height = 3,
#        plot_grid(sim1.ps, sim2.ps, sim3.ps, sim4.ps, labels=c("(a)","(b)","(c)","(d)"),ncol=4,nrow=1))

ggarrange(sim1.ps, sim2.ps, nrow = 1, ncol = 2)
ggsave("sim_ps.pdf", width=8, height=4)





psiLj.cp.case2.tau3 <- ggplot(psiLj.case2.tau3, aes(x=fct_inorder(rownames(psiLj.case2.tau3)), cp))+
  # geom_point() +
  # geom_hline(yintercept=0.95, linetype="dotted", color = "black")+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
  #       plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
  #       panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  # xlab("Estimator")+
  # ylab("Coverage probability")+
  scale_y_continuous(labels = scales::percent, limits=c(0.7,1),
                     breaks=c(0.70, 0.80, 0.90, 0.95, 1))



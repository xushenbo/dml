library(ggplot2)
library(ggpubr)

followup <- 15 # years
tau <- 10

### curves plot
## cancer incidence
load("cancer_incidence_rv.RData")
f.sim$dr.direct.mc.ha.driptw <- f.sim$dr.direct.mc.ha.driptw/365.25
f.sim$dr.direct.mc.ha.drow <- f.sim$dr.direct.mc.ha.drow/365.25
# ATE
cancer.incidence.ate.diff <- ggplot(f.sim$dr.direct.mc.ha.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(-0.6, 0.2), breaks = c(-0.6, -0.4, -0.2, 0, 0.2)) +  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a) ATE, competing") +
  theme_classic()

# ATO
cancer.incidence.ato.diff <- ggplot(f.sim$dr.direct.mc.ha.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(-0.6, 0.2), breaks = c(-0.6, -0.4, -0.2, 0, 0.2)) +  xlab("Follow-up time (years)") +
  xlab("Follow-up time (years)") +
  ylab("Difference in RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b) ATO, competing") +
  theme_classic()

## treatment-specific
# IPTW
cancer.incidence.ate.a0 <- ggplot(f.sim$dr.direct.mc.ha.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a0), color="black") +
  geom_line(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = c(0, 0.5, 1, 1.5, 2, 2.5)) +
  xlab("Follow-up time (years)") +
  ylab("RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(a) Sulphonylurea (A=0)") +
  theme_classic()

cancer.incidence.ate.a1 <- ggplot(f.sim$dr.direct.mc.ha.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a1), color="black") +
  geom_line(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = c(0, 0.5, 1, 1.5, 2, 2.5)) +
  xlab("Follow-up time (years)") +
  ylab("RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(b) Metformin (A=1)") +
  theme_classic()

# OW
cancer.incidence.ato.a0 <- ggplot(f.sim$dr.direct.mc.ha.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a0), color="black") +
  geom_line(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = c(0, 0.5, 1, 1.5, 2, 2.5)) +
  xlab("Follow-up time (years)") +
  ylab("RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(c)") +
  theme_classic()

cancer.incidence.ato.a1 <- ggplot(f.sim$dr.direct.mc.ha.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a1), color="black") +
  geom_line(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = c(0, 0.5, 1, 1.5, 2, 2.5)) +
  xlab("Follow-up time (years)") +
  ylab("RMTL (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(d)") +
  theme_classic()

## cancer mortality
load("cancer_mortality_rv.RData")
f.sim$dr.direct.mc.ha.driptw <- f.sim$dr.direct.mc.ha.driptw/365.25
f.sim$dr.direct.mc.ha.drow <- f.sim$dr.direct.mc.ha.drow/365.25
# ATE
cancer.mortality.ate.diff <- ggplot(f.sim$dr.direct.mc.ha.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4)) +  xlab("Follow-up time (years)") +
  ylab("Difference in RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(c) ATE, survival") +
  theme_classic()

# ATO
cancer.mortality.ato.diff <- ggplot(f.sim$dr.direct.mc.ha.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.diff.lower.cb, ymax = estimate.diff.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.diff), color="black") +
  geom_line(aes(y=estimate.diff.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.diff.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4)) +  xlab("Follow-up time (years)") +
  ylab("Difference in RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(d) ATO, survival") +
  theme_classic()

## treatment-specific
# IPTW
cancer.mortality.ate.a0 <- ggplot(f.sim$dr.direct.mc.ha.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a0), color="black") +
  geom_line(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
  xlab("Follow-up time (years)") +
  ylab("RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(e)") +
  theme_classic()

cancer.mortality.ate.a1 <- ggplot(f.sim$dr.direct.mc.ha.driptw, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a1), color="black") +
  geom_line(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
  xlab("Follow-up time (years)") +
  ylab("RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(f)") +
  theme_classic()

# OW
cancer.mortality.ato.a0 <- ggplot(f.sim$dr.direct.mc.ha.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a0.lower.cb, ymax = estimate.a0.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a0), color="black") +
  geom_line(aes(y=estimate.a0.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a0.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
  xlab("Follow-up time (years)") +
  ylab("RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(g)") +
  theme_classic()

cancer.mortality.ato.a1 <- ggplot(f.sim$dr.direct.mc.ha.drow, aes(x=t)) +
  geom_ribbon(aes(x = t, ymin = estimate.a1.lower.cb, ymax = estimate.a1.upper.cb), fill = "gray") +
  geom_line(aes(y=estimate.a1), color="black") +
  geom_line(aes(y=estimate.a1.lower.ci), color = "black", linetype="dashed") +
  geom_line(aes(y=estimate.a1.upper.ci), colour = "black", linetype="dashed") +
  scale_x_continuous(limits = c(0, followup), breaks = c(0, followup/5, 2*followup/5, 3*followup/5, 4*followup/5, followup)) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
  xlab("Follow-up time (years)") +
  ylab("RMST (years)") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(title="(h)") +
  theme_classic()

## arrange
ggarrange(cancer.incidence.ate.diff, cancer.incidence.ato.diff, cancer.mortality.ate.diff, cancer.mortality.ato.diff,
          nrow = 2, ncol = 2)
ggsave("app_effect.pdf", width=8, height=8)

ggarrange(cancer.incidence.ate.a0, cancer.incidence.ate.a1,
          cancer.incidence.ato.a0, cancer.incidence.ato.a1,
          cancer.mortality.ate.a0, cancer.mortality.ate.a1,
          cancer.mortality.ato.a0, cancer.mortality.ato.a1, nrow = 4, ncol = 2)
ggsave("app_treatment_specific.pdf", width=8, height=16)


# defense slides
ggarrange(cancer.incidence.ate.diff, cancer.incidence.ato.diff, nrow = 1, ncol = 2)
ggsave("app_effect_incidence.pdf", width=8, height=4)


ggarrange(cancer.mortality.ate.diff, cancer.mortality.ato.diff, nrow = 1, ncol = 2)
ggsave("app_effect_mortality.pdf", width=8, height=4)




ggarrange(cancer.incidence.ate.a0, cancer.incidence.ate.a1, cancer.incidence.ate.diff,
          cancer.incidence.ato.a0, cancer.incidence.ato.a1, cancer.incidence.ato.diff,
          cancer.mortality.ate.a0, cancer.mortality.ate.a1, cancer.mortality.ate.diff,
          cancer.mortality.ato.a0, cancer.mortality.ato.a1, cancer.mortality.ato.diff,
          nrow = 4, ncol = 3)
ggsave("application.pdf", width=12, height=16)

### table
## cancer incidence
load("cancer_incidence_rv.RData")
f.sim$dr.direct.mc.ha.driptw <- f.sim$dr.direct.mc.ha.driptw/365.25
f.sim$dr.direct.mc.ha.drow <- f.sim$dr.direct.mc.ha.drow/365.25
# ATE
paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")

# ATO
paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")

## cancer mortality
load("cancer_mortality_rv.RData")
f.sim$dr.direct.mc.ha.driptw <- f.sim$dr.direct.mc.ha.driptw/365.25
f.sim$dr.direct.mc.ha.drow <- f.sim$dr.direct.mc.ha.drow/365.25

# ATE
paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.driptw$t, f.sim$dr.direct.mc.ha.driptw$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")

# ATO
paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a0.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.a1.upper.cb, xout=tau)$y, digits=3), ")")

paste0(round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff, xout=tau)$y, digits=3), " (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.se, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.lower.ci, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.upper.ci, xout=tau)$y, digits=3), ")",
       " & (", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.lower.cb, xout=tau)$y, digits=3), ", ", round(approx(f.sim$dr.direct.mc.ha.drow$t, f.sim$dr.direct.mc.ha.drow$estimate.diff.upper.cb, xout=tau)$y, digits=3), ")")

### ps density plot
# cancer incidence
a.sl.lib <- c("SL.mean","SL.glm","SL.nnet","SL.rpartPrune","SL.xgboost","SL.ranger","SL.step","SL.glmnet","SL.earth") # a.sl.lib <- c("SL.mean","SL.glm") # ,"SL.kernelKnn","SL.gam"
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>% # too slow
  filter(is.na(cancer_eventdate)|(prescription_date<cancer_eventdate)) %>%
  filter(!(cancer==1 & is.na(cancer_eventdate))) %>%
  mutate(event1 = cancer) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(cancer_eventdate, deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  mutate(event2 = ifelse(cancer, 0, death)) %>%
  mutate(end.time2 = end.time1) %>%
  mutate(event3 = ifelse((!is.na(deathdate_n) & cancer==0), 2, cancer)) %>%
  mutate(end.time3 = end.time1) %>%
  mutate(event5 = cancer) %>%
  mutate(end.time5 = as.numeric(difftime(if_else(as.logical(cancer), cancer_eventdate, if_else(as.logical(death), tod_n, pmin(tod_n, lcd_n, na.rm=TRUE))), prescription_date, units="days"))) %>%
  mutate(event6 = if_else(as.logical(cancer), 0, if_else(as.logical(death),1,0))) %>%
  mutate(end.time6 = as.numeric(difftime(if_else(as.logical(cancer), tod_n, if_else(as.logical(death), deathdate_n, pmin(tod_n, lcd_n, na.rm=TRUE))), prescription_date, units="days"))) %>%
  mutate(event7 = 1*(death|cancer)) %>%
  mutate(end.time7 = end.time1) %>%
  filter(end.time1>0) %>%
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]
n <- nrow(data.df)
data.df$a <- data.df$ref_sulf

train.ind <- sample(1:n, ceiling(n/2), replace=FALSE) # sample training set
train.df <- data.df[train.ind,]
train.df <- train.df[order(train.df$patid),]
estimation.df <- data.df[setdiff(1:n, train.ind),]
estimation.df <- estimation.df[order(estimation.df$patid),] # estimation.df <- estimation.df[order(estimation.df$end.time,-estimation.df$event),]

train.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = estimation.df$a, X = estimation.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
estimation.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = train.df$a, X = train.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred) # estimation.df$ps <- predict(glm(a~x1+x2+x3+x4,data=estimation.df,family=binomial()),type="response")
data.df <- rbind(train.df, estimation.df)
cancer.incidence.ps <- ggplot(data.df, aes(x=ps, group=factor(a)))+
  geom_density(aes(linetype=ifelse(a, "dotted", "solid")))+
  xlab("Propensity score")+
  ylab("Density")+
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1))+
  labs(fill = "A:")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = guide_legend(title.position = "left"))+
  scale_y_continuous(limits=c(0, 22), breaks=c(0, 5, 10, 15, 20))+
  ggtitle("(a) Cancer incidence")

# cancer mortality
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>%
  filter(prescription_date > cancer_eventdate) %>%
  mutate(event1 = death) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]
n <- nrow(data.df)
data.df$a <- data.df$ref_sulf

train.ind <- sample(1:n, ceiling(n/2), replace=FALSE) # sample training set
train.df <- data.df[train.ind,]
train.df <- train.df[order(train.df$patid),]
estimation.df <- data.df[setdiff(1:n, train.ind),]
estimation.df <- estimation.df[order(estimation.df$patid),] # estimation.df <- estimation.df[order(estimation.df$end.time,-estimation.df$event),]

train.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = estimation.df$a, X = estimation.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                           cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
estimation.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = train.df$a, X = train.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                                cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred) # estimation.df$ps <- predict(glm(a~x1+x2+x3+x4,data=estimation.df,family=binomial()),type="response")

data.df <- rbind(train.df, estimation.df)
cancer.mortality.ps <- ggplot(data.df, aes(x=ps, group=factor(a)))+
  geom_density(aes(linetype=ifelse(a, "dotted", "solid")))+
  xlab("Propensity score")+
  ylab("Density")+
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1))+
  labs(fill = "A:")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = guide_legend(title.position = "left"))+
  scale_y_continuous(limits=c(0, 22), breaks=c(0, 5, 10, 15, 20))+
  ggtitle("(b) Mortality after cancer diagnosis")

ggarrange(cancer.incidence.ps, cancer.mortality.ps, nrow = 1, ncol = 2)
ggsave("application_ps.pdf", width=8, height=4)


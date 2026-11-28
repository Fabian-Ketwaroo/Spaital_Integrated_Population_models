# Compare the two case-study runs produced by GRCA_SV_BayesNSGP.R (run from the repository root):
#   Case study/sipm_RW.rds       (NIMBLE scalar RW sampler on the spatial effects)
#   Case study/sipm_RW_NN_GP.rds (BayesNSGP sampler RW_NN_GP)
#
# What is compared, in order of importance:
#   1. run time and effective sample size per second
#   2. the identifiable quantities: survival (sa), recruitment (f), hyperparameters, fixed effects
#   3. the survival spatial field w.as. This field is only weakly identified (it trades off with
#      the intercept and the longitude/latitude terms), so two chains of the SAME sampler already
#      differ; the between-run agreement is therefore judged against the within-run agreement.
# A summary table is written to results/comparison_summary.csv.

library(coda)
library(ggplot2)

res_dir <- "Case Study"
rw   <- readRDS(file.path(res_dir, "sipm_RW.rds"))
nngp <- readRDS(file.path(res_dir, "sipm_RW_NN_GP.rds"))
stopifnot("the two runs must use identical MCMC settings" = identical(rw$settings, nngp$settings))
out <- list()   # collects the numbers for the summary table

# helpers ------------------------------------------------------------------------------------
site_means <- function(s, idx) sapply(s, function(ch) colMeans(ch[, idx, drop = FALSE]))  # nodes x chains
cols    <- function(s, pat) grep(pat, colnames(s[[1]]))
rhat    <- function(s) gelman.diag(s, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1]

# 1. Run time and effective sample size per second -------------------------------------------
hours <- c(RW = rw$time[["elapsed"]], RW_NN_GP = nngp$time[["elapsed"]]) / 3600
cat("1. Run time (hours)\n"); print(round(hours, 2))
cat("   Speed-up (RW / RW_NN_GP):", round(hours[["RW"]] / hours[["RW_NN_GP"]], 1), "\n")

ess_rate <- function(res) effectiveSize(res$samples) / res$time[["elapsed"]]
r_rw <- ess_rate(rw); r_nngp <- ess_rate(nngp); ratio <- r_nngp / r_rw
cat("\n   ESS/second ratio (RW_NN_GP / RW), all monitored parameters\n"); print(round(summary(ratio), 1))
hyper <- c("rho.as", "alpha2.as", "rho.af", "alpha2.af")
cat("   Hyperparameters:\n")
print(round(data.frame(RW = r_rw[hyper], RW_NN_GP = r_nngp[hyper], ratio = ratio[hyper]), 4))
out$hours_RW <- hours[["RW"]]; out$hours_RW_NN_GP <- hours[["RW_NN_GP"]]
out$speedup <- hours[["RW"]] / hours[["RW_NN_GP"]]
out$ess_per_sec_ratio_median <- median(ratio); out$ess_per_sec_ratio_min <- min(ratio); out$ess_per_sec_ratio_max <- max(ratio)

# 2. Identifiable quantities ------------------------------------------------------------------
cat("\n2. Agreement of identifiable quantities\n")
# 2a. survival (sa) and recruitment (f): posterior means over all nodes (pooled chains)
for (q in c("sa", "f")) {
  i  <- cols(rw$samples2, paste0("^", q, "\\["))
  i2 <- match(colnames(rw$samples2[[1]])[i], colnames(nngp$samples2[[1]]))
  a <- rowMeans(site_means(rw$samples2, i)); b <- rowMeans(site_means(nngp$samples2, i2))
  ok <- is.finite(a) & is.finite(b)
  cat(sprintf("   %-2s correlation of posterior means: %.3f  (%d nodes); median |difference| = %.4f\n",
              q, cor(a[ok], b[ok]), sum(ok), median(abs(a[ok] - b[ok]))))
  out[[paste0("cor_", q)]] <- cor(a[ok], b[ok])
}

# 2b. hyperparameters and fixed effects (first monitor set)
summ <- function(res) {
  s <- res$samples
  m <- do.call(rbind, s)
  data.frame(mean = colMeans(m), sd = apply(m, 2, sd), n.eff = effectiveSize(s), Rhat = rhat(s))
}
sum_rw <- summ(rw); sum_nngp <- summ(nngp)
cat("\n   Max R-hat, all monitored parameters: RW =", round(max(sum_rw$Rhat, na.rm = TRUE), 2),
    " RW_NN_GP =", round(max(sum_nngp$Rhat, na.rm = TRUE), 2), "\n")
d_sd <- abs(sum_nngp$mean - sum_rw$mean) / pmax(sum_rw$sd, 1e-12)
cat("   Difference in posterior means, in posterior sd units: median =", round(median(d_sd, na.rm = TRUE), 3),
    " max =", round(max(d_sd, na.rm = TRUE), 2), "\n")
mcse <- function(s) s$sd / sqrt(pmax(s$n.eff, 1))
z <- (sum_nngp$mean - sum_rw$mean) / sqrt(mcse(sum_nngp)^2 + mcse(sum_rw)^2)
cat("   Share of parameters differing by < 2 Monte Carlo SE:", round(mean(abs(z) < 2, na.rm = TRUE), 3),
    "\n   (the Monte Carlo SEs are too small where R-hat is large, because ESS ignores between-chain\n",
    "    differences; read this number together with the R-hat above)\n")
out$maxRhat_RW <- max(sum_rw$Rhat, na.rm = TRUE); out$maxRhat_RW_NN_GP <- max(sum_nngp$Rhat, na.rm = TRUE)
out$median_diff_in_sd <- median(d_sd, na.rm = TRUE); out$share_within_2mcse <- mean(abs(z) < 2, na.rm = TRUE)

cat("\n   Posterior means of the hyperparameters, per chain:\n")
hp <- c(hyper, "beta.s.lat", "beta.s.long")
for (nm in c("RW", "RW_NN_GP")) {
  m <- sapply(get(if (nm == "RW") "rw" else "nngp")$samples, function(ch) colMeans(ch[, hp]))
  colnames(m) <- paste0(nm, "_ch", seq_len(ncol(m))); print(round(m, 3))
}

p <- ggplot(data.frame(RW = sum_rw$mean, RW_NN_GP = sum_nngp$mean), aes(RW, RW_NN_GP)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50") + geom_point(alpha = 0.7) +
  theme_bw(base_size = 14) + labs(x = "Posterior mean (RW)", y = "Posterior mean (RW_NN_GP)")
ggsave(file.path(res_dir, "posterior_means_RW_vs_RWNNGP.png"), p, width = 5, height = 5, dpi = 300)

# 3. Survival spatial field w.as: limit of the model, so compare with the within-run baseline --
cat("\n3. Survival spatial field w.as (weakly identified)\n")
idx <- cols(rw$samples2, "^w\\.as\\[")
cm  <- list(RW = site_means(rw$samples2, idx), RW_NN_GP = site_means(nngp$samples2, idx))
w_in_rw <- cor(cm$RW[, 1], cm$RW[, 2]); w_in_nn <- cor(cm$RW_NN_GP[, 1], cm$RW_NN_GP[, 2])
w_btw   <- cor(rowMeans(cm$RW), rowMeans(cm$RW_NN_GP))
cat(sprintf("   correlation of posterior means: within RW (chain 1 vs 2) = %.3f | within RW_NN_GP = %.3f | between runs = %.3f\n",
            w_in_rw, w_in_nn, w_btw))
cat("   A between-run value inside the range of the within-run values means the two samplers differ\n",
    "  no more than two chains of the same sampler.\n")
for (nm in c("RW", "RW_NN_GP")) {
  s <- get(if (nm == "RW") "rw" else "nngp")$samples2[, idx]
  ess <- effectiveSize(s); rh <- rhat(s)
  cat(sprintf("   %-9s median ESS = %4.0f | share R-hat > 1.1 = %.2f | max R-hat = %.2f\n",
              nm, median(ess), mean(rh > 1.1, na.rm = TRUE), max(rh, na.rm = TRUE)))
  out[[paste0("w_ess_median_", nm)]] <- median(ess)
  out[[paste0("w_share_rhat_gt_1.1_", nm)]] <- mean(rh > 1.1, na.rm = TRUE)
  out[[paste0("w_max_rhat_", nm)]] <- max(rh, na.rm = TRUE)
}
out$w_cor_within_RW <- w_in_rw; out$w_cor_within_RW_NN_GP <- w_in_nn; out$w_cor_between <- w_btw

png(file.path(res_dir, "w_as_agreement.png"), 1500, 500, res = 130)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot(cm$RW[, 1], cm$RW[, 2], pch = 16, cex = .5, main = "RW: chain 1 vs 2", xlab = "chain 1", ylab = "chain 2"); abline(0, 1)
plot(cm$RW_NN_GP[, 1], cm$RW_NN_GP[, 2], pch = 16, cex = .5, main = "RW_NN_GP: chain 1 vs 2", xlab = "chain 1", ylab = "chain 2"); abline(0, 1)
plot(rowMeans(cm$RW), rowMeans(cm$RW_NN_GP), pch = 16, cex = .5, main = "RW vs RW_NN_GP", xlab = "RW", ylab = "RW_NN_GP"); abline(0, 1)
dev.off()

# 4. Summary table -------------------------------------------------------------------------------
tab <- data.frame(quantity = names(out), value = signif(unlist(out), 4))
write.csv(tab, file.path(res_dir, "comparison_summary.csv"), row.names = FALSE)
cat("\nSummary written to", file.path(res_dir, "comparison_summary.csv"), "\n")

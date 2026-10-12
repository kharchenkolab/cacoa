S <- Sys.getenv("S"); res <- readRDS(file.path(S, "bench/tail/scc_r1_B19999.rds"))
countGE <- function(v, x) { sv <- sort(v); length(sv) - findInterval(x - 1e-12, sv) }
calib <- function(p) { al <- c(1e-2, 1e-3, 1e-4, 1e-5); setNames(vapply(al, function(a) mean(p < a) / a, numeric(1)), paste0("r", al)) }
perm <- res$heat$perm; perm <- perm[, colSums(is.finite(perm)) > 0.9 * nrow(perm)]; B <- nrow(perm); m <- ncol(perm)
set.seed(7); held <- sort(sample(seq(4000, B), 4000)); H <- perm[held, ]
full <- perm[-held, ]; mu.f <- colMeans(full); sd.f <- apply(full, 2, sd)
zH <- abs(sweep(sweep(H, 2, mu.f), 2, sd.f, "/"))
pp <- function(pool, z) (countGE(pool, as.numeric(z)) + 1) / (length(pool) + 1)
poolOf <- function(rows) sort(abs(as.numeric(sweep(sweep(perm[rows, ], 2, mu.f), 2, sd.f, "/"))))
cat("true moments throughout; pools from different relabeling sets\n")
out <- rbind("rows 1..1999" = calib(pp(poolOf(1:1999), zH)), "rows 2000..3998" = calib(pp(poolOf(2000:3998), zH)),
             "random 1999 non-held rows" = calib(pp(poolOf(sample(setdiff(seq_len(B), held), 1999)), zH)),
             "random 3998" = calib(pp(poolOf(sample(setdiff(seq_len(B), held), 3998)), zH)),
             "random 7996" = calib(pp(poolOf(sample(setdiff(seq_len(B), held), 7996)), zH)),
             "all non-held (15999)" = calib(pp(poolOf(setdiff(seq_len(B), held)), zH)))
print(round(out, 3))
# how many distinct relabelings and cells form the pooled tail of a 1999-row pool
Z <- abs(sweep(sweep(perm[1:1999, ], 2, mu.f), 2, sd.f, "/")); q <- quantile(Z, 1 - 1e-4); idx <- which(Z > q, arr.ind = TRUE)
cat(sprintf("pool of 1999 rows: %d values above the 1e-4 quantile from %d distinct relabelings and %d distinct cells; top 10 relabelings carry %d\n",
            nrow(idx), length(unique(idx[, 1])), length(unique(idx[, 2])), sum(sort(table(idx[, 1]), decreasing = TRUE)[1:10])))
q3 <- quantile(Z, 1 - 1e-3); idx3 <- which(Z > q3, arr.ind = TRUE)
cat(sprintf("                   %d values above the 1e-3 quantile from %d relabelings and %d cells\n", nrow(idx3), length(unique(idx3[, 1])), length(unique(idx3[, 2]))))
# same for the held-out draws: which relabelings carry the exceedances of the pooled 1e-4 threshold
ex <- which(zH > q, arr.ind = TRUE); cat(sprintf("held-out: %d exceedances (expected %d) from %d relabelings; top 10 relabelings carry %d\n", nrow(ex), round(length(zH) * 1e-4), length(unique(ex[, 1])), sum(sort(table(ex[, 1]), decreasing = TRUE)[1:10])))

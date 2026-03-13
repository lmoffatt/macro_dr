## Moment Poisson Bootstrap

The idea is to replace `Moment_statistics` with a `Moment_poisson_bootstrap`.

`Moment_poisson_bootstrap` would maintain (and update) a collection of bootstrap replicates, i.e. a vector of `Moment_statistics`.

The output of the bootstrap must be the distribution of the moments, i.e. the empirical cumulative distribution of the statistic across replicates, together with requested confidence intervals (e.g. 5%–95%).

## Design Notes

- If you only have the final aggregated `(count, mean, variance)` summary, you cannot reconstruct an exact nonparametric bootstrap of the original data; you need to update the bootstrap online while streaming samples (or store the samples).
- For a dataset of size `N`, the classical bootstrap corresponds to multinomial replicate weights `(k_1..k_N)` with `sum k_i = N`. Marginally, each `k_i` is `Binomial(N, 1/N)`, but the weights are not independent.
- A standard streaming alternative is the **Poisson bootstrap**: draw per-sample weights `k_i ~ Poisson(1)` (independent). Then the total weight `K = sum k_i` varies between replicates (`K ~ Poisson(N)`), which is usually acceptable for large `N`.
- Implementation-wise, `Moment_bootstrap` can update each replicate’s `Moment_statistics` using the sampled weight for the incoming observation (i.e., add the same observation `k` times, or implement a weighted-update variant).
- For confidence intervals, compute the statistic of interest from each replicate and then use empirical quantiles (e.g. 5% and 95%) of the replicate distribution.

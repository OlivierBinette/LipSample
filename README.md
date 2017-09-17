# Sampling-Lipschitz-continuous-densities
MATLAB code to sample arbitrary Lipschitz continuous densities on the interval.

## Example

Sample from the normal distribution restricted to [-5, 5].

```matlab
sample = lipsample(@normpdf, 0.25, [-5 5], 10000000);
```

Plot the result.

```matlab
  pretty_hist(sample, [-5 5]);
```

![](norm_sample.png)

# Sampling Lipschitz continuous densities
MATLAB code to sample arbitrary Lipschitz continuous densities on the interval.

## Example

Sample from the normal distribution restricted to [-5, 5].

```matlab
sample = lipsample(@normpdf, 0.25, [-5 5], 10000000);
```
Here 0.25 is the Lipschitz continuity order of the normal distribution, [-5, 5] is the interval to which it is restricted, and 10000000 is the number of samples.

Plot the result.

```matlab
pretty_hist(sample, [-5 5]);
```

![](norm_sample.png)

## Usage and functionalities

`sample = lipsample(@f, L, [a b], m);` Draws _m_ samples from the probability density function _f_ which is Lipschitz continuous of order _L_ on _[a,b]_. For _m_ larger than 1000 * (b-a) * L, samples are drawn from a piecewise linear approximation of _f_ at distance at most 0.001 in the supremum norm. Exact samples or different error tolerances can be obtained by using the 'Tolerance' key.

`sample = lipsample(@f, L, [a b], m, 'Tolerance', epsilon);` ... If _epsilon_ = 0, exact samples are drawn for all sample sizes. Otherwise, for large _m_, samples are drawn from a piecewise linear approximation of _f_ at distance less than _epsilon_ in the supremum norm.

## Notes
The density function _f_ should be roughly normalized. It would be straightforward to adapt the method to unnormalized densities, but that would require the user to also specify the number of components to be used for the mixture envelope of _f_. 

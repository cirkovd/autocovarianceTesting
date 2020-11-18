autocovarianceTesting: Test for Equality of Autocovariance Functions in
Time Series
================
Daniel Cirkovic

The package ‘autocovarianceTesting’ implements multiple hypothesis tests
for testing equality of autocovariance functions for two mean zero,
stationary time series. The time series are allowed to be multivariate
and linearly dependent. We implement methods assuming independence of
the series from [“Testing equality of stationary autocovariances.” Lund
et. al (2009)](https://doi.org/10.1111/j.1467-9892.2009.00616.x) and new
extensions to linearly dependent times series. We also implement a
weighting mechanism to the former tests by borrowing ideas from [“New
Weighted Portmanteau Statistics for Time Series Goodness of Fit Testing”
Fisher and Gallagher
(2012)](https://doi.org/10.1080/01621459.2012.688465). The package also
implements bootstrap methods from [“A computational bootstrap procedure
to compare two dependent time series” Jin et. al
(2019)](https://doi.org/10.1080/00949655.2019.1639704), offering
extensions to multivariate series.

‘autocovarianceTesting’ includes one function: `autocovarianceTest`.
`autocovarianceTest` allows the the user to supply two time series and
select from any combination of the previously mentioned methods to test
for equality of autocovariance functions. Choices of methods are made
through the `test` argument. See function documentation for more
details.

## Installation

*To be filled in once package completed*

## Example

*To be filled in once package completed*

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Sampling.R
\name{sample_one_point}
\alias{sample_one_point}
\title{Sampling one point}
\usage{
sample_one_point(s)
}
\arguments{
\item{s}{A matrix with each row (h'(x), h(x_i),x_i,xmin,xmax,cumulative probability, normalized cumulative probability) at each node. This can be calculated using envelop_density(h,x,domain)}
}
\description{
Can sample one point from the given piece-wise exponential density using the upper hull. The Inverse-CDF method is used here. A random number is generated, then it's used with the help of the calculated inverse CDF function of the piece-wise exponential density to generate an x.
}


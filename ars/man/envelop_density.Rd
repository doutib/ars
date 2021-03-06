% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Envelop.R
\name{envelop_density}
\alias{envelop_density}
\title{Envelope density}
\usage{
envelop_density(h,x,domain)
}
\arguments{
\item{h}{The log of the density function}

\item{x}{The nodes, found by abscissae()}

\item{domain}{A vector of length 2 giving the domain (left bound, right bound) of the density function}
}
\description{
Returns a normalized density function using the upper hull. A matrix with each row (h'(x_i), h(x_i),x_i,xmin,xmax,cumulative probability, normalized cumulative probability) at each node is returned. With this infomation we can sample with sk.
}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/carSignal-v15.9.10.R
\name{get.normalised.ABP.stimulus}
\alias{get.normalised.ABP.stimulus}
\title{This function creates an normalised ABP step stimulus.}
\usage{
get.normalised.ABP.stimulus(sampling.time = 0.1, time.until.release = 10,
  time.after.release = 20, smooth.step.stimulus = FALSE, filter.order = 2,
  cutoff.frequency = 0.2,
  left.stabilisation.time = ifelse(smooth.step.stimulus, 30, 0),
  time.rounding.digits = format.info(sampling.time)[2],
  time.tol = sampling.time/100)
}
\arguments{
\item{sampling.time:}{sampling time for the ABP stimulus (in seconds).
Default value: 0.1 s}

\item{time.until.release:}{time included in the resulting signal until the
time of cuff release (in seconds).
Default value: 10 s}

\item{time.after.release:}{time included in the resulting signals after the
time of cuff release.
Default value: 20 s}

\item{smooth.step.stimulus:}{logical, whether the ABP step stimulus should
be smoothed (filtered).
Default value: FALSE}

\item{filter.order:}{the order of the low-pass Butterworth filter used to
smooth the ABP step stimulus.
Default value: 2}

\item{cutoff.frequency:}{the cutoff frequency for the low-pass Butterworth
filter used to smooth the ABP step stimulus (in Hz).
Default value: 0.20 Hz}

\item{left.stabilisation.time:}{the time to be used for stabilisation of the
step. This time is included in the step when
a filter is applied, but removed after that.
Default value: 30 s when a filter will be
applied, 0 s otherwise.}

\item{time.rounding.digits:}{the number of decimal when time instants are
used.
Default value: the number of decimals in the
sampling time}
}
\value{
The function's answer is a list with the following:
.$time.instants: the time instants in which the signal is sampled.
.$ABP.normalised: the normalised ABP stimulus generated
.$sampling.time: the sampling time used
.$time.release: the time instant in which the thigh cuffs are supposedly released. Currently at the second 0.0.
}
\description{
This function creates an normalised ABP step stimulus.
}


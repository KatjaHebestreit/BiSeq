\name{logisticRegression}
\alias{logisticRegression}
\alias{logisticRegression,formula,character,BSrel,missing-method}
\alias{logisticRegression,formula,character,BSrel,numeric-method}


\title{
Estimates and tests a group factor
}
\description{
It is used to fit a linear model on the log odds of each (smoothed) methylation level. The first independent variable in \code{formula} is tested to be unequal to zero.
}
\usage{
logisticRegression(formula, link, object, mc.cores)
}

\arguments{
  \item{formula}{An object of class \code{formula}. For the first independent variable the p-value and the effect on methylation is returned.
  }
  \item{link}{A character specifying the link function. Currently,
  \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"log"}, \code{"loglog"} are supported.
}
  \item{object}{A \code{BSrel} object.
  }
    \item{mc.cores}{
    Passed to \code{mclapply}.
}
}

\value{
A \code{data.frame} containing the position, chromosome, P value, estimated
methylation level in group 1 and group 2 and methylation difference of
group 1 and group 2.
}

\author{
Katja Hebestreit
}

\seealso{
  \code{mclapply}, \code{glm}
  }

\examples{
data(predictedMeth)

logisticResults <- logisticRegression(formula = ~group, link = "logit",
                                      object = predictedMeth) 
}

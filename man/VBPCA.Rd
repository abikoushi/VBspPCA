\name{VBPCA}
\alias{VBPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variational Bayes method for probablistic pricipal component analysis with sparase matrix
}
\description{
Variational Bayes method for probablistic pricipal component analysis with sparase matrix
}
\usage{
VBPCA(Y, rank, iter = 10, prior_prec = 1, a = 1, b = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
 input data matrix to decompose
}
  \item{rank}{
  rank of the output matices
}
  \item{iter}{
  the number of the iterations
}
  \item{prior_prec}{
precision paremeter of Gaussian prior distribution
}
  \item{a}{
shape paremeter of the gamma prior distribution
}
  \item{b}{
  rate paremeter of the gamma prior distribution
}
}
\value{
  \item{mean_row }{shape parameters of the posterir gamma distribution for row-axis}
  \item{prec_row }{rate parameters of the posterir gamma distribution for row-axis}
  \item{mean_col }{shape parameters of the posterir gamma distribution for col-axis}
  \item{prec_col }{rate parameters of the posterir gamma distribution for col-axis}
  \item{logprob }{lower bound of the log-joint probability}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (Y, rank, iter = 10, prior_prec = 1, a = 1, b = 1) 
{
    if (!any(class(Y) == "dgTMatrix")) {
        Y <- as(Y, "TsparseMatrix")
    }
    doVB_norm(y = Y@x, rowi = Y@i, coli = Y@j, Nr = Y@Dim[1], 
        Nc = Y@Dim[2], L = rank, iter = iter, prior_prec = prior_prec, 
        a = a, b = b)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

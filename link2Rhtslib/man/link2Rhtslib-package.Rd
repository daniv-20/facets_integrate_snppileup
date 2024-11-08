\name{link2Rhtslib-package}
\alias{link2Rhtslib-package}
\alias{link2Rhtslib}
\docType{package}
\title{
  
  Package with native code that links to Rhtslib, for confirming proper
  linkage and include directives.

}
\description{
  
  Package imports a header in HTSlib in Rhtslib pkg and has a routine
  that calls a function in the libhts library.

}
\details{

  A successful call \code{htsVersion} means the correct header was
  found, and linking was successful.

  \code{callLocal} calls a native routine in the package without external
  dependencies.
  
}
\author{
Maintainer: Nathaniel Hayden <nhayden@fredhutch.org>
}

\examples{

  htsVersion()
  callLocal()

}

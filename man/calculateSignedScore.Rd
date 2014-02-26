\name{calculateSignedScore}
\alias{calculateSignedScore}
\title{Calculates pairwise score between every gene in data1 and every gene in data2}
\usage{
calculateSignedScore(label, data1, data2)
}
\arguments{
  \item{label}{group labels}

  \item{data1}{first gene matrix}

  \item{data2}{second gene matrix}
}
\value{
signedscore Signed socres
}
\description{
Calculates pairwise score between every gene in data1 and
every gene in data2. It is signed because
score(i,j)=P(X_i<X_j|1)-P(X_i<X_j|0) - P(X_i>X_j|1)
+P(X_i>X_j|0) + C The third and the forth terms are for
avoiding the equlities. The C term is proportion to the
secondary score to avoid the ties.  C =
(E(X_j-X_i|1)-E(X_j-X_i|0))/K where K is big enough to make
sure that the secondary score does not intervene with the
primary score. a,M,k,N are usable for fisher exact test. M
is the size of the class situation == 0 and N is the total
number of samples. a_ij = #(X_i<X_j|0) k_ij = #(X_i<X_j)
P_ij = P(X_i<X_j|0) and Q(X_i<X_j|1)
}
\keyword{Calculate}
\keyword{score}
\keyword{signed}


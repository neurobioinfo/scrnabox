\name{split_up_down}
\alias{split_up_down}
\title{split contrast
}
\description{
This function splits the contrast to up and down according the cut_pvalue, cut_log_low, and cut_log_up
}
\usage{
split_up_down(file_path,dest_path,cut_pvalue,cut_log_low,cut_log_up,n_row)
}
\arguments{
  \item{file_path}{path of directory where contrasts are}
  \item{dest_path}{path of directory where you want to save new files}
  \item{cut_pvalue}{threshold of p-value}
  \item{cut_log_low}{threshold of low up}
  \item{cut_log_up}{threshold of up low}
  \item{n_row}{no of genes}

}
\details{
This function splits the contrast to up and down
}

\author{ Saeid Amiri}

\examples{
#file_path='~/contrasts'
#dest_path='~/contrasts_up_down'
#split_up_down(file_path,dest_path,cut_pvalue=0.5,cut_log_low=-0.2,cut_log_up=0.2,n_row=200)
}

\keyword{utils}

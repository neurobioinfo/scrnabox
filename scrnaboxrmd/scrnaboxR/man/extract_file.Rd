\name{extract_file}
\alias{extract_file}
\title{extract file
}
\description{
This function imports the files from a given path and save the contents as list.
}
\usage{
extract_file(file_path)
}
\arguments{
  \item{file_path}{path of folder}
}
\details{
The output of this function is a vector and list; a vector of names of files, and the list of imported valued from files in the defined folder.
}

\author{Saeid Amiri}

\examples{
#file_path<-"/Users/sam/Desktop/final_scrna_dge/test"
#file0<-extract_file(file_path)
#file0[[2]]
#file0[[1]]
}

\keyword{utils}

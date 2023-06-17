\name{csv_to_xlsx}
\alias{csv_to_xlsx}
\title{csv_to_xlsx
}
\description{
This function imports the csv files from a given path and save them as sheet in xlsx file.
}
\usage{
csv_to_xlsx(files_path,dest_path,name='sheet')
}
\arguments{
  \item{files_path}{folders where csv files are}
  \item{dest_path}{folder where xlsx will saved}
  \item{sheet}{name of sheet with strt with what you define here}
}
\details{
The output of this function is an xlsx where there sheets are  contents of csv files.
}

\author{Saeid Amiri}

\examples{
#files_path<-"/Users/sam/Desktop/final_scrna_dge/csvs"
#dest_path<-"/Users/sam/Desktop/final_scrna_dge/dest"
#csv_to_xlsx(files_path,dest_path,name='sheet')
}

\keyword{utils}

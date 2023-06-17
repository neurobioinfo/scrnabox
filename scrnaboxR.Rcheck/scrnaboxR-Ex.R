pkgname <- "scrnaboxR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "scrnaboxR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('scrnaboxR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("annotate")
### * annotate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: annotation
### Title: generate annotation
### Aliases: annotation
### Keywords: Annotation

### ** Examples

#level_cluster='integrated_snn_res.0.7'
#PWD='~/Desktop/annot/annotation_d/'
#PSUE='~/Desktop/annot/annotation_d/step6/objs/seu_int_clu.rds'
#top_sel=5
# db <- c('Descartes_Cell_Types_and_Tissue_2021',
#'CellMarker_Augmented_2021','Azimuth_Cell_Types_2021')
#ClusterMarkers<-read.csv(file = "~/Desktop/step7/objs7/ClusterMarkers.csv",row.names = 1)
# annotation(level_cluster,ClusterMarkers,PWD,PSUE,top_sel,db)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("annotate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_index")
### * compute_index

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Compute_index
### Title: Different indices clustering
### Aliases: Compute_index
### Keywords: Annotation

### ** Examples

# X<-sample(c('a1', 'a2'), 1000, replace=TRUE)
# Y<-sample(c('b1', 'b2', 'b3'), 1000, replace=TRUE)
# aa<-compute_index(X,Y)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_index", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("csv_to_xlsx")
### * csv_to_xlsx

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: csv_to_xlsx
### Title: csv_to_xlsx
### Aliases: csv_to_xlsx
### Keywords: utils

### ** Examples

#files_path<-"/Users/sam/Desktop/final_scrna_dge/csvs"
#dest_path<-"/Users/sam/Desktop/final_scrna_dge/dest"
#csv_to_xlsx(files_path,dest_path,name='sheet')



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("csv_to_xlsx", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dist_contrast")
### * dist_contrast

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dist_contrast
### Title: distance matrix of contrasts
### Aliases: dist_contrast
### Keywords: distance

### ** Examples

#file_path<-"/Users/sam/Desktop/final_scrna_dge/test"
#file0<-extract_file(file_path)
#con<-file0[[2]]
#dst<-dist_contrast(con)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dist_contrast", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dittplot")
### * dittplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dittplot
### Title: dittplot
### Aliases: dittplot
### Keywords: plot

### ** Examples

con<-list()
#con[[1]]<-read.csv('/Users/sam/Desktop/sa_int/intract_contrast/all/des1.csv')
#con[[2]]<-read.csv('/Users/sam/Desktop/sa_int/intract_contrast/all/des2.csv')
#dist_contrast<-function(con)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dittplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("extract_file")
### * extract_file

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: extract_file
### Title: extract file
### Aliases: extract_file
### Keywords: utils

### ** Examples

#file_path<-"/Users/sam/Desktop/final_scrna_dge/test"
#file0<-extract_file(file_path)
#file0[[2]]
#file0[[1]]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("extract_file", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("split_up_down")
### * split_up_down

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: split_up_down
### Title: split contrast
### Aliases: split_up_down
### Keywords: utils

### ** Examples

#file_path='~/contrasts'
#dest_path='~/contrasts_up_down'
#split_up_down(file_path,dest_path,cut_pvalue=0.5,cut_log_low=-0.2,cut_log_up=0.2,n_row=200)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("split_up_down", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

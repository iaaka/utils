#!/software/R-4.1.0/bin/Rscript

# Pavel Mazin, 2022
# pm19@sanger.ac.uk

library(argparser,quietly = TRUE)

p = arg_parser("Background ATAC quantification.\nProgramm quantifies ATAC in all avaliable barcodes and saves results as three gziped text files (colnames, rownames and sparce matrix in MM format)",hide.opts=TRUE)
p = add_argument(p, "fragments", help="input fragmnent file (atac_fragments.tsv.gz). File should be indexed.")
p = add_argument(p, "peaks", help="input peak file in bed format")
p = add_argument(p, "out", help="folder to save results. Script will create the folder, it will throw an exception if folder already exists.")

pars = parse_args(p)

suppressMessages({
  library(GenomicRanges,quietly = TRUE)
  library(Signac,quietly = TRUE)
  library(Matrix,quietly = TRUE)
})
#' Save sparce matrix as three gzipped files
#'
#'
#' @param mtx sparce matrix to be saved
#' @param path output folder. Function will create the folder, it will throw an exception if folder already exists.
#' @param row.fname,col.fname,mtx.fname file name to save rownames,colnames, and matrix (without 'gz' at the end)
#'
#' @details all output files will be gziped.
#'
#' @return
#' @export
#'
#' @examples
saveMM = function(mtx,path,row.fname='features.tsv',col.fname='barcodes.tsv',mtx.fname='matrix.mtx'){
  dir.create(path)
  r = gzfile(paste0(path,'/',row.fname,'.gz'), "w")
  c = gzfile(paste0(path,'/',col.fname,'.gz'), "w")
  
  writeLines(rownames(mtx), r)
  writeLines(colnames(mtx), c)
  Matrix::writeMM(mtx,paste0(path,'/',mtx.fname))
  system(paste0("gzip ",path,'/',mtx.fname))
  close(r)
  close(c)
}

#' Quantifies atac in all cells
#'
#' @param peaks dataframe: chr,start,end
#' @param frags.fname path to fragnment file
#' @param cells optional character vector of barcodes to consider. If NUll, all barcodes will be used.
#'
#' @return sparce matrix
#' @export
#'
#' @examples
atac.bkg.quant = function(peaks,frags.fname,cells=NULL){
  peaks = makeGRangesFromDataFrame(peaks[,1:3])
  frags = CreateFragmentObject(path = frags.fname,cells = cells)
  
  gratac = FeatureMatrix(fragments = frags,
                         features = peaks,
                         cells = cells)
}

peaks = read.table(pars$peaks)
colnames(peaks)[1:3]= c('chr','start','end')
r = atac.bkg.quant(peaks,pars$fragments)
saveMM(r,path=pars$out)

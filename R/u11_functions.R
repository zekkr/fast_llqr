# Formal research entry: same R adapters and core as fastllqr 0.2.0.
# Usage: u11 <- load_u11(); u11$llqr_seq_ppro(x,y,...)
load_u11 <- function(root=getwd(), build_dir=file.path(root,"tmp/u11-build"),
                     rebuild=FALSE, profile="portable") {
  root <- normalizePath(root)
  builder <- new.env(parent=baseenv())
  sys.source(file.path(root,"reproduction/build.R"),envir=builder)
  if(rebuild || !file.exists(file.path(build_dir,"ssqr_optimized.so")))
    builder$build_u11(build_dir,profile,root)
  env <- new.env(parent=globalenv())
  for(file in list.files(file.path(root,"R/u11"),pattern="[.]R$",full.names=TRUE))
    sys.source(file,envir=env)
  dll <- dyn.load(file.path(normalizePath(build_dir),"ssqr_optimized.so"))
  env$C_ssqr_kernel_path <- getNativeSymbolInfo("ssqr_kernel_path",dll)
  list(llqr_seq_ppro=env$llqr_seq_ppro,tvcqr_seq_ppro=env$tvcqr_seq_ppro,
       dll=dll,backend="unified_u11")
}

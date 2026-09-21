# Build a frozen experiment in an isolated directory, without editing its sources.
build_u11 <- function(destination, profile="portable", root=getwd()) {
  dir.create(destination,recursive=TRUE,showWarnings=FALSE)
  destination <- normalizePath(destination)
  source_dir <- file.path(root,"src/fortran/u11")
  frozen <- file.path(root,"reproduction/frozen/fortran")
  fc <- Sys.getenv("FC",Sys.which("gfortran"))
  if (!nzchar(fc)) stop("A Fortran compiler is required; set FC.")
  flags <- switch(profile, portable=c("-O2"),
    checked=c("-O0","-fcheck=all","-fbacktrace","-ffpe-trap=invalid,zero,overflow"),
    `paper-hpc`=c("-O3","-march=native","-funroll-loops","-ffast-math"),
    stop("Unknown build profile"))
  rlib <- file.path(R.home(),"lib")
  base <- c("-shared","-fPIC","-ffree-line-length-none",flags)
  if (Sys.info()[["sysname"]]=="Darwin")
    base <- c(base,paste0("-Wl,-syslibroot,",system2("xcrun","--show-sdk-path",stdout=TRUE)))
  commands <- character()
  compile <- function(output,files) {
    args <- c(base,shQuote(files),paste0("-L",shQuote(rlib)),"-lRlapack","-lRblas","-o",shQuote(output))
    commands <<- c(commands,paste(shQuote(fc),paste(args,collapse=" ")))
    if(system2(fc,args)!=0L) stop("Fortran build failed")
  }
  old <- setwd(destination); on.exit(setwd(old),add=TRUE)
  compile("ssqr_optimized.so",file.path(source_dir,c("weighted_qr_core.f90","kernel_entry.f90")))
  for(name in c("llqr_seq_lean_sortskip","tvcqr_seq_lean_nohistory"))
    compile(paste0(name,".so"),file.path(frozen,paste0(name,".f90")))
  writeLines(commands,"compile_commands.txt")
  writeLines(system2(fc,"--version",stdout=TRUE),"compiler_version.txt")
  invisible(destination)
}

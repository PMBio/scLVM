.onAttach <- function(libname, pkgname){
  scLVM_path = system.file(package="scLVM")
  python.assign('sclvm_path', scLVM_path)
  
  #tell the user what to do if necessary python dependencies are missing.
  limix_error = "limix is not in the python path. Please specify the limix path by calling the configLimix function after loading the package."
  package_error_base = c("It seems python package","is not installed. Try to install it using 'pip install")
  package_error_end  = "If the package is installed already, make sure rPyhton is using the correct python version. If you have several versions of python installed, have a look at our installation guide on github."
  
  python_error_fun <- function(e){
    if(pmatch("No module named ",e[[1]])==1){
      miss_package = strsplit(e[[1]]," ")[[1]][4]
      if(miss_package=="limix"){
        print(limix_error)}else{
        print(paste(package_error_base[1],miss_package,package_error_base[2],
                    miss_package,"'",package_error_end,sep = " "))
      }
    }else{
      print(e)
    }
  }
  
  
  tryCatch(python.load(system.file("pysrc","init_data.py",package="scLVM")),
           error = function(e) python_error_fun(e))
}

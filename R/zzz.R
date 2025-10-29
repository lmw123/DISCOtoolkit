.onAttach <- function(libname, pkgname) {
  tryCatch({
    if(packageVersion("DISCOtoolkit") == fromJSON("https://www.immunesinglecell.com/disco_v3_api/getToolkitVersion")$version) {
      packageStartupMessage(paste0("Welcome to DISCOtoolkit. version ", packageVersion("DISCOtoolkit")))
    } else {
      packageStartupMessage("A new version of DISCOtoolkit is available. Please update to ensure all functions work properly.")
    }
  }, error = function(e){
    packageStartupMessage("Warning: Unable to check the latest version of DISCOtoolkit. Please check https://github.com/lmw123/DISCOtoolkit")
  })
}

.onLoad <- function(libname, pkgname) {
  op <- options()

  op.disco <-  list(
    disco_url = "https://www.immunesinglecell.com/disco_v3_api/",
    timeout = 6000
  )

  options(op.disco)
  invisible()
}


GetJson <- function(url, info.msg, error.msg){
  tryCatch({
    message(info.msg)
    return(fromJSON(url))
  }, error = function(e){
    stop(error.msg)
  })
}


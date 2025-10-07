.onAttach <- function(libname, pkgname) {
<<<<<<< HEAD
  tryCatch({
    if(packageVersion("DISCOtoolkit") == fromJSON("https://www.immunesinglecell.com/disco_v3_api/getToolkitVersion")$version) {
      packageStartupMessage("Welcome to DISCOtoolkit. version 1.2.0 (Oct 7, 2025)")
=======

  tryCatch({
    if(packageVersion("DISCOtoolkit") == fromJSON("https://disco.bii.a-star.edu.sg/disco_v3_api/getToolkitVersion")$version) {
      packageStartupMessage("Welcome to DISCOtoolkit.")
>>>>>>> 055018bc3c98e747bebd54708bc0da167a57f0e2
    } else {
      packageStartupMessage("A new version of DISCOtoolkit is available. Please update to ensure all functions work properly.")
    }
  }, error = function(e){
<<<<<<< HEAD
    packageStartupMessage("Warning: Unable to check the latest version of DISCOtoolkit. Please check https://github.com/lmw123/DISCOtoolkit")
=======
    packageStartupMessage("Welcome to DISCOtoolkit.")
>>>>>>> 055018bc3c98e747bebd54708bc0da167a57f0e2
  })
}

.onLoad <- function(libname, pkgname) {
  op <- options()

  op.disco <-  list(
<<<<<<< HEAD
    disco_url = "https://www.immunesinglecell.com/disco_v3_api/",
    timeout = 6000
  )

=======
    disco_url = "https://disco.bii.a-star.edu.sg/disco_v3_api/",
    timeout = 6000
  )

  # tryCatch({
  #   op.disco <- list(
  #     disco_url = fromJSON("http://www.immunesinglecell.org/api/vishuo/getToolkitUrl")$url,
  #     timeout = 6000
  #   )
  # }, error = function(e){
  #   message("Fail to get url prefix, use default value")
  # })

  # toset <- !(names(op.disco) %in% names(op))
  # if (any(toset)) options(op.disco[toset])
>>>>>>> 055018bc3c98e747bebd54708bc0da167a57f0e2
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


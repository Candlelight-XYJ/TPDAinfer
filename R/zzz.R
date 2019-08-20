##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0(pkgname, " v", pkgVersion, "  ",
                "For help: https://github.com/Candlelight-XYJ/", pkgname, "\n\n")

  packageStartupMessage(paste0(msg))
}

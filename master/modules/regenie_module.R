# =============================================================================
# modules/regenie_module.R
#
# Loads a standalone Shiny app (app.R) from a directory without launching it.
# Strategy: read the file, strip the shinyApp() call, source into an isolated
# environment, and return that environment so the master app can access
# $ui and $server directly.
# =============================================================================

.load_shiny_app_from_dir <- function(app_dir) {
  app_file <- file.path(app_dir, "app.R")

  if (!file.exists(app_file)) {
    message("[load] app.R not found at: ", app_file)
    return(NULL)
  }

  # Read file and remove the shinyApp() call so sourcing does not launch the app
  lines <- readLines(app_file, warn = FALSE)
  lines <- lines[!grepl("shinyApp", lines)]

  tmp <- tempfile(fileext = ".R")
  writeLines(lines, tmp)
  on.exit(unlink(tmp), add = TRUE)

  env <- new.env(parent = globalenv())
  env$REGENIE_APP_DIR <- app_dir
  tryCatch({
    old_wd <- setwd(app_dir)
    on.exit(setwd(old_wd), add = TRUE)
    source(tmp, local = env)
  }, error = function(e) {
    message("[load] Error sourcing app: ", e$message)
    return(NULL)
  })

  if (!exists("ui", envir = env) || !exists("server", envir = env)) {
    message("[load] Could not find 'ui' or 'server' after sourcing ", app_file)
    return(NULL)
  }

  env
}

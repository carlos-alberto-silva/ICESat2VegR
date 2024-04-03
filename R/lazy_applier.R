lazy_apply_dt_call <- function(dt, call, group.by = "") {
  inner.env <- new.env()
  func.names.args <- all.names(call)
  fn.name <- func.names.args[1]
  dt.name <- paste0("var", as.character(runif(1)), collapse = "")

  fn.is.function <- tryCatch(is.function(get(fn.name)), error = function(e) FALSE)
  if (!fn.is.function) {
    inner.env[[fn.name]] <- parent.frame(2)[[fn.name]]
    stopifnot("Could not find function" = is.function(inner.env[[fn.name]]))
  }

  call_string <-  as.character(list(call))
  print(call_string)
  entire_call <-  paste0(dt.name, "[, ", call_string, ", ", group.by, "]")
  print(entire_call)

  inner.env[[dt.name]] <- as.data.table(dt)
  result <- eval(parse(text = entire_call), inner.env)
  return(result)
}

seg_position <- function(lon_ph, lat_ph) {
  metrics <- list(
    lat_seg = min(lat_ph) + diff(range(lat_ph)) / 2,
    lon_seg = min(lon_ph) + diff(range(lon_ph)) / 2
  )
  return(metrics)
}

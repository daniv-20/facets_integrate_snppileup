callLocal <- function() {
    .Call("localfunc")
}

htsVersion <- function() {
    .Call("extvers")
}

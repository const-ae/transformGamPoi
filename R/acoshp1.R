


acoshp1 <- function(x){
  acosh(x + 1)
}

setGeneric("acoshp1", function(x){
  standardGeneric("acoshp1")
})

setMethod("acoshp1", signature = "sparseMatrix", function(x){
  acoshp1(as(x, "CsparseMatrix"))
})


setMethod("acoshp1", signature = "CsparseMatrix", function(x){
  # acosh(0 + 1) == 0
  x@x <- acosh(x@x + 1)
  x
})

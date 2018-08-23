# Auxiliary functions
left.inv <- function(ys, Fs) {
  ys <- sort(ys)
  Fs <- sort(Fs)
  stats::stepfun(Fs, c(ys,max(ys)), right=TRUE)
}

wecdf <- function(y,outcome,w) weighted.mean((outcome <= y),w)

#' Provides a function to plot the directional means of a mixture of von Mises-Fisher distribution.
#' @import ggplot2
#' @import tibble
#' @import Matrix
#' @param movMFresul  results from the function movMF_EM.
#' @param with_color add colors to the plot.
#' @param map_value a function to transform the value of representation
#' @md
#' @return returns a plot of the directional means.
#'
#' @export




plot_proto <- function(movMFresult, with_color=FALSE, map_values=NULL) {
  class_order <- order(movMFresult$alpha)
  proto_pre_order <- order_proto(movMFresult$mu[class_order,])
  proto_order <- proto_pre_order$order_bin
  xs <- seq(0,1, length.out=length(proto_order)+1)
  cs <- cumsum(movMFresult$alpha[class_order])
  ys <- (c(0,cs[-length(cs)])+cs)/2
  xmins <- xs[-length(xs)]
  xmaxs <- xs[-1]
  ymins <- c(0,cs[-length(cs)])
  ymaxs <- cs
  positions <- cbind(expand.grid(xmin=xmins, ymin=ymins),
                     expand.grid(xmax=xmaxs, ymax=ymaxs))
  positions$value <-
    as.vector(t(as.matrix(movMFresult$mu[class_order,proto_order])))
  if(!is.null(map_values)) {
    positions$value <- map_values(positions$value)
  }
  if(with_color) {
    plot_mat(positions, cs, proto_pre_order$change_index/length(xmins))
  } else {
    plot_mat(positions, cs)
  }
}

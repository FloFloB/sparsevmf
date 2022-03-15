#' Provides a function to plot the data set according to the result of a mixture of von Mises-Fisher distribution.
#' @import ggplot2
#' @import tibble
#' @import Matrix
#' @param movMFresul  results from the function movMF_EM.
#' @param X dataset used by the function movMF_EM.
#' @param with_color add colors to the plot.
#' @md
#' @return returns a plot of the dataset.
#'
#' @export




plot_data <- function(movMFresult, X, with_color=FALSE) {
  class_order <- order(movMFresult$alpha)
  proto_pre_order <- order_proto(movMFresult$mu[class_order,])
  proto_order <- proto_pre_order$order_bin
  data_cl_mapper <- rep(NA,length(class_order))
  data_cl_mapper[class_order] <- 1:length(class_order)
  data_order <- order(data_cl_mapper[movMFresult$cluster])
  xs <- seq(0, 1, length.out=length(proto_order)+1)
  ys <- seq(0, 1, length.out=nrow(X)+1)
  xmins <- xs[-length(xs)]
  xmaxs <- xs[-1]
  ymins <- ys[-length(ys)]
  ymaxs <- ys[-1]
  positions <- cbind(expand.grid(xmin=xmins, ymin=ymins),
                     expand.grid(xmax=xmaxs, ymax=ymaxs))
  positions$value <- as.vector(t(as.matrix(X[data_order, proto_order])))
  positions$value <- as.integer(positions$value>0)
  cs <- c(0,ys[1+which(diff(sort(data_cl_mapper[movMFresult$cluster]))!=0)],1)
  if(with_color) {
    plot_mat(positions, cs, proto_pre_order$change_index/length(xmins))
  } else {
    plot_mat(positions, cs)
  }
}

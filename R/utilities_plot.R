firstIndex <- function(vec) {
  unq <- unique(vec)
  sapply(unq, function(x) {min(which(vec == x))})
}

order_proto <- function(proto, prec=sqrt(.Machine$double.eps)) {
  proto_bin <- ifelse(as.matrix(abs(proto))<prec,0,1)
  order_bin <- order(decreasing=T, colSums(proto_bin),
                     apply(proto_bin,2,paste,collapse=''),colSums(abs(proto)))
  change_index <- firstIndex(colSums(proto_bin)[order_bin])-1
  change_index[1] <-0
  change_index[length(change_index)+1] <-length(order_bin)
  list(order_bin=order_bin,change_index=change_index)
}

plot_mat <- function(positions, cs=NULL, change_index=NULL) {
  plt <- ggplot(positions, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=value)) +
    geom_rect()+ scale_fill_gradient(low="white",high="black",limits=c(0,1)) +
    theme(legend.position = "none",
          axis.text  = element_blank(),
          panel.grid = element_blank(),
          axis.line  = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank()) +
    scale_x_continuous(expand = expansion(0)) +
    scale_y_continuous(expand = expansion(0))
  if(!is.null(change_index)) {
    colors <- RColorBrewer::brewer.pal(12, "Set3")
    colors <- rep(colors,length.out=length(change_index)-1)
    plt <- plt + annotate("rect",
                          xmin=change_index[-length(change_index)],
                          xmax=change_index[-1],
                          ymin=rep(0, length(change_index)-1),
                          ymax=rep(1, length(change_index)-1),
                          fill=colors,
                          alpha=0.2)
  }
  if(!is.null(cs)) {
    plt <- plt + geom_hline(yintercept=c(0,cs),col="grey")
  }
  plt
}

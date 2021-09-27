#' @title Plot a cancer progression network
#' 
#' @description This function plots a cancer progression network inferred 
#' in plotCPN. 
#' 
#' @param fit An object obtained from plotCPN. 
#' 
#' @author Phillip B. Nicol <philnicol740@gmail.com>
#' 
#' @examples 
#' ## Load the example data
#' data("example")
#' 
#' ## Fit a conjunctive Bayesian network (CBN)
#' out <- fitCPN(example, epsilon=0.01)
#' 
#' ## Plot the graph
#' plotCPN(out)
#' 
plotCPN <- function(fit) {
  G <- fit$graph
  xy <- igraph::layout.auto(G)
  V(G)$x <- xy[,1]
  V(G)$y <- xy[,2]
  
  p <- ggraph::ggraph(G, "manual", x=V(G)$x,y=V(G)$y)
  p <- p + ggraph::geom_edge_link(arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(3, 'mm'))
  p <- p + ggraph::geom_node_point(size=5)
  p <- p+ggraph::geom_node_label(aes(label=name), repel=TRUE)
  plot(p)
}


plotBootstrapCPN <- function(bs_out, threshold = 0.5) {
  bestEdges <- bs_out$Edges[bs_out$counts/bs_out$N > threshold,]

  G <- igraph::make_directed_graph(c(t(bestEdges)))
  xy <- igraph::layout.auto(G)
  V(G)$x <- xy[,1]
  V(G)$y <- xy[,2]
  E(G)$weight <- bs_out$counts[bs_out$counts/bs_out$N > threshold]/bs_out$N
  E(G)$color <- rgb(1-E(G)$weight,1-E(G)$weight,1-E(G)$weight)
  
  p <- ggraph::ggraph(G, "manual", x=V(G)$x,y=V(G)$y)
  p <- p + ggraph::geom_edge_link(aes(colour=E(G)$color),
                                  arrow = arrow(
    length = unit(4, 'mm')), end_cap = circle(3, 'mm'))
  p <- p + ggraph::scale_edge_color_manual(values=E(G)$color)
  p <- p + ggraph::geom_node_point(size=5)
  p <- p+ggraph::geom_node_label(aes(label=name), repel=TRUE)
  p <- p + theme(legend.position="none")
  plot(p)
}
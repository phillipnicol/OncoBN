
plotCPN <- function(fit) {
  g <- make_directed_graph(fit$edgelist)
  V(g)$color <- "white"
  V(g)$shape <- "none"
  plot(g, label.color = "blue", label.cex = 1.6, edge.arrow.size = 0.5)
}

plotBootstrapCPN <- function(bs_out, threshold = 0.5) {
  bestEdges <- bs_out$Edges[bs_out$counts/bs_out$N > threshold,]
  g <- make_directed_graph(c(t(bestEdges)))
  V(g)$color <- "white"
  V(g)$shape <- "none"
  E(g)$weight <- bs_out$counts[bs_out$counts/bs_out$N > threshold]/bs_out$N
  E(g)$color <- rgb(1-E(g)$weight,1-E(g)$weight,1-E(g)$weight)
  plot(g, label.color = "blue", label.cex = 1.6, edge.color = E(g)$color, edge.arrow.size = 0.5)
}
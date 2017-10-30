#' Threshold plots for posterior latent networks
#'
#' @param w_adj Weighted adjacency matrix for latent network (symmetric).
#' @param thresh Threshold (betweeen 0 and 1) for an edge to be plotted. Defaults to 0.1 (i.e., edges appearing in at least 10\% of posterior latent networks will be shown).
#' @param obs The adjacency matric of the observed network on the same vertices as \code{w_adj}. Default is NULL. If provided, plots the latent and observed networks side by side.
#' @return A threshold plot
#' @examples
#'   data(smYr1)
#'   data(smYr1_hire)
#'   year1 <- CAgibbs(smYr1, smYr1_hire, c(2,8), c(2,8), burn=100, thin=5, numpostdraws=100)
#'     # in practice run for longer
#'
#'   threshplot(year1$Theta_prop) #latent network only
#'   threshplot(year1$Theta_prop, obs=smYr1) #latent network side by side with observed network

threshplot <- function(w_adj,
                       thresh=0.1,
                       obs=NULL) {

  stopifnot(isSymmetric(w_adj))

  #setup for latent plot
  adj.thresh <- ifelse(w_adj >= thresh, 1, 0)
  g <- network(adj.thresh)
  el <- data.frame(as.edgelist(g))
  el$weight <- rep(0,nrow(el))
  for (i in 1:nrow(el)) {
    el$weight[i] <- w_adj[el$X1[i], el$X2[i]]
  }
  el$wcolor <- paste("gray", 100-100*el$weight, sep="")
  g %e% "weight" <- el$weight
  g %e% "wcolor" <- el$wcolor

  # fix coordinates (will be used for both plots if obs != NULL)
  coord.net <- network.layout.fruchtermanreingold(g, layout.par=NULL)
  g %v% "x" <- coord.net[, 1]
  g %v% "y" <- coord.net[, 2]
  b <- theme(panel.background = element_rect(color = "grey50"))
  t1 <- paste("Latent Network, Tie Threshold = ", thresh)

  # latent plot
  lat.plot <- ggnet2(g,
                     mode=c("x","y"),
                     node.size=6,
                     node.color="cornflowerblue",
                     edge.size=1.5,
                     edge.color="wcolor")

  if (is.null(obs) == FALSE) {

    # setup for observed plot
    g.obs <- network(obs)
    g.obs %v% "x" <- coord.net[, 1]
    g.obs %v% "y" <- coord.net[, 2]

    # observed plot
    obs.plot <- ggnet2(g.obs,
                       mode=c("x","y"),
                       node.size=6,
                       node.color="cornflowerblue",
                       edge.size=1.5,
                       edge.color="black",
                       arrow.size = 6, arrow.gap = 0.02)

    # put together plots
    grid.arrange(lat.plot + ggtitle(t1) + b,
                 obs.plot + ggtitle("Observed Network") + b,
                 ncol=2)

  } else {

    lat.plot + ggtitle(t1) + b

  }

}

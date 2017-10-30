## Use roxygen to make documentation for package

#' Gibbs Sampler
#'
#' @param obs_network The adjacency matrix A of the observed directed network. A is an n by n matrix, where n is the number of nodes. A_ij equals 1 if a tie exists between i and j. A_ij equals 0 if there is no tie between i and j. A_ii equals 0 by default. If fitting collaboration models on multiple time points for the same set of people, n should be the full number of people. Values of \code{obs_network} for people who have not been hired by that time will be 0.
#' @param hired An indicator vector of length n indicating whether each person in \code{obs_network} had been hired by this point in time.
#' @param p_prior Vector of hyperparameters for the prior Beta distribution of p. Specify both alpha and beta for the Beta(alpha,beta) distribution, e.g. \code{c(2,8)}.
#' @param q_prior Vector of hyperparameters for the prior Beta distribution of q. Specify both alpha and beta for the Beta(alpha,beta) distribution, e.g. \code{c(2,8)}.
#' @param d The density used for Theta_ij if i and j both hired by the observed time. If i and/or j not hired, Theta_ij is zero. Defaults to 0.15, a reasonable network density.
#' @param burn Burn-in period for Gibbs sampler. Defaults to 1000 iterations.
#' @param thin Thinning interval for Gibbs sampler. Defaults to 25 iterations.
#' @param numpostdraws Number of postrior draws to keep. Sampler will run \code{numpostdraws*thin + burn} iterations in total, and keep just nonpostdraws. Defaults to 1000 iterations.
#' @return A list consisting of
#'   \item{p}{Matrix (dimension numpostdraws x n) where each row is vector of p_i's for a posterior draw.}
#'   \item{q}{Matrix (dimension numpostdraws x n) where each row is vector of q_i's for a posterior draw.}
#'   \item{Theta_prop}{Proportion of times each tie appears in latent adjacency matrices over all draws. The sum of \code{Theta_all} divided by \code{numpostdraws}.}
#'   \item{Theta_all}{List of all latent adjacency matrices.}
#'   \item{dens_vec}{Vector consisting of the network density for each posterior draw.}
#' @examples
#'   # This example uses data from Year 1
#'   data(smYr1)
#'   data(smYr1_hire)
#'   year1 <- CAgibbs(smYr1, smYr1_hire, c(2,8), c(2,8), burn=100, thin=5, numpostdraws=100)
#'     # in practice run for longer
#'
#'   colMeans(year1$p_mat) #posterior mean for p (false negative) for each person at time 1
#'   colMeans(year1$q_mat) #posterior mean for q (false positive) for each person at time 1
#'
#'   threshplot(year1$Theta_prop) #latent network only
#'   threshplot(year1$Theta_prop, obs=smYr1) #latent network side by side with observed network

CAgibbs <- function(obs_network,
                    hired,
                    p_prior,
                    q_prior,
                    d = 0.15,
                    burn = 1000,
                    thin = 25,
                    numpostdraws = 1000) {

  # Total number of Gibbs iterations
  steps <- numpostdraws * thin + burn

  # Set up
  n <- nrow(obs_network)

  delta_mat <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (hired[i] == 1 & hired[j] == 1) {
    		delta_mat[i,j] <- d
    	} else {
    		delta_mat[i,j] <- 0
    	}
    }
  }
  diag(delta_mat) <- 0

  p <- rep(0, n)
  q <- rep(0, n)

  # Create storage for posterior draws
  p_mat <- matrix(0, numpostdraws, n)
  q_mat <- matrix(0, numpostdraws, n)
  Theta_sum <- matrix(0, n, n)
  Theta_all <- list()
  dens_vec <- rep(0, numpostdraws)

  # Initial Theta
  # NOTES we start Theta at all 0, which means the initial draws of p_i and q_i are from their prior dists
  Theta <- matrix(0, n, n)

  # Set up counter of kept iterations
  count <- 0

  # Gibbs sampler
  for (l in 1:steps) {

    # Draw new p_i's and q_i's
    for (i in 1:n) {
      p[i] <- rbeta(1, (as.numeric((1-obs_network[i,]) %*% Theta[i,]) + p_prior[1]),
                     (as.numeric(obs_network[i,] %*% Theta[i,]) + p_prior[2]))
      q[i] <- rbeta(1, (as.numeric(obs_network[i,] %*% (1-Theta[i,])) + q_prior[1]),
                     (as.numeric((1-obs_network[i,]) %*% (1-Theta[i,]) - 1) + q_prior[2]))
    }

    # Draw new Theta
    for (i in 2:n) {
  		for (j in 1:(i-1)) {

  			# Calculate probabilities needed to get new r_ij
  			prob1 <- delta_mat[i,j] * p[i]^(1-obs_network[i,j]) * (1-p[i])^(obs_network[i,j]) * (p[j])^(1-obs_network[j,i]) * (1-p[j])^(obs_network[j,i])
  			prob0 <- (1-delta_mat[i,j]) * q[i]^(obs_network[i,j]) * (1-q[i])^(1-obs_network[i,j]) * (q[j])^(obs_network[j,i]) * (1-q[j])^(1-obs_network[j,i])

  			## NOTES
  			# so r_ij is prob1 / (prob1 + prob 0)

  			# normalize
  			temp <- prob1 + prob0
  			prob1 <- prob1 / temp
  			prob0 <- prob0 / temp

  			# sample Theta[i,j] from Bernoulli(r_ij) and enforce mutuality
  			Theta[i,j] <- sample((0:1), 1, prob=c(prob0, prob1))
  			Theta[j,i] <- Theta[i,j]
  		}
    }

    # Keep only thinned iterations
    if (l > burn & (l-burn) %% thin == 0) {
      # print counter
      count <- count + 1
      if (count %% 10 == 0) {
        print(count)
      }

      # save and update
      p_mat[count, ] <- p
      q_mat[count, ] <- q
      Theta_sum <- Theta_sum + Theta
      Theta_all[[count]] <- Theta
      dens_vec[count] <- network.density(as.network(Theta))
    }

  }

  Theta_prop <- Theta_sum / numpostdraws

  return(list(p_mat=p_mat,
              q_mat=q_mat,
              Theta_prop=Theta_prop,
              Theta_all=Theta_all,
              dens_vec=dens_vec))

}

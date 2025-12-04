#' Generate random samples from a Dirichlet distribution
#'
#' @param n Integer, number of samples to generate.
#' @param alpha Numeric vector of concentration parameters (must be > 0).
#'
#' @return A matrix with n rows and length(alpha) columns,
#'         where each row sums to 1.
#'
#' @examples
#' set.seed(123)
#' alpha <- c(2, 3, 5)
#' samples <- rdirichlet(5, alpha)
#' rowSums(samples)  # should all be 1
#'
#' # Example in your code
#' prop <- c(0.4, 0.3, 0.2, 0.1)
#' s <- 50
#' nsamples <- 10
#' ct <- c("A", "B", "C", "D")
#' pdir <- rdirichlet(nsamples, prop * s)
#' pure_scenario_dataframe <- data.frame(
#'   matrix(pdir, ncol = length(ct),
#'          dimnames = list(paste0("dirichlet", 1:nsamples), ct))
#' )
#' @export
rdirichlet <- function(n, alpha) {
  if (any(alpha <= 0)) stop("All elements of alpha must be > 0.")
  k <- length(alpha)
  x <- matrix(
    stats::rgamma(n * k, shape = rep(alpha, each = n), rate = 1),
    ncol = k,
    byrow = TRUE
  )
  x / rowSums(x)
}



#' Select top pseudo-bulk marker genes from MiCBuS DE results
#'
#' Ranks and returns the top markers.
#'
#' @param DEpseudo.df Data frame returned by \code{MiCBuS()}.
#' @param ntop Integer. Number of top marker genes to return.
#'
#' @return Character vector of top-ranked marker genes.
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom stats median
#' @importFrom dplyr mutate arrange slice_head pull desc
psMarker_Top <- function(DEpseudo.df, ntop = 200) {

  # prevent R CMD check "no visible binding for global variable"
  baseMean <- FC_exp_score <- gene <- NULL

  if (!all(c("log2FoldChange", "baseMean", "gene") %in% colnames(DEpseudo.df)))
    stop("DEpseudo.df must contain log2FoldChange, baseMean, and gene columns.")

  median50 <- stats::median(DEpseudo.df$baseMean) * 50

  DEpseudo.df <- DEpseudo.df %>%
    dplyr::mutate(baseMean_adj = pmin(baseMean, median50))

  DEpseudo.df$FC_exp_score <-
    (DEpseudo.df$log2FoldChange)^2 *
    log10(DEpseudo.df$baseMean_adj + 1)

  psMarkerTop <- DEpseudo.df %>%
    dplyr::arrange(dplyr::desc(FC_exp_score)) %>%
    dplyr::slice_head(n = ntop) %>%
    dplyr::pull(gene)

  return(psMarkerTop)
}


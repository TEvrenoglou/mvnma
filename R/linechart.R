#' Line chart showing the VIKOR results
#'
#' @description
#' A line chart showing the results of VišeKriterijumska Optimizacija I
#' Kompromisno Rešenje (VIKOR) across the three metrics Q, S and R for each
#' treatment. Within each metric, lower values represent better treatment
#' performance.
#' 
#' @param x An object of class \code{\link{vikor}}.
#' @param sort A character specifying the order of treatments on the x-axis.
#'   By default, the order is according to the Q-metric (\code{"Q"}). Ordering
#'   treatments according to the S (\code{"S"}) and R (\code{"R"}) metrics is
#'   also possible.
#' @param exclude A character specifying a metric that will not be displayed in
#'   the graph. By default, all metrics are displayed (\code{"none"}).
#'   Alternative options are to exclude either the Q-metric (\code{"Q"}),
#'   the S-metric (\code{"S"}), or the R-metric (\code{"R"}).
#' @param n A numeric value indicating the number of treatments to be plotted
#'   in the graph. By default, all treatments are displayed. If specified,
#'   only the first \code{"n"} treatments according to the hierarchy specified
#'   from \code{"sort"} argument are plotted.
#' @param linewidth A numeric value specifying the width of the lines
#'   (default: 1.1).
#' @param size A numeric value specifying the size of the points
#'   (default: 2).
#' @param \dots Additional arguments passed to \code{\link[ggplot2]{ggplot}}
#'   function.
#'
#' @return
#' A \code{ggplot} object.
#' 
#' @examples
#' # Locate file "mvnma_examples.rda" with mvnma() results
#' .fname <- system.file("extdata/mvnma_examples.rda", package = "mvnma")
#' load(.fname)
#' 
#' # Rank treatments using SUCRAs (default)
#' ranks12 <- mvrank(mvnma12, small.values = c("undes", "undes"))
#' ranks12
#' 
#' # Get the best compromise solution across the efficacy outcomes
#' vk12 <- vikor(ranks12)
#' 
#' # Visualize the results with default settings
#' linechart(vk12)
#' 
#' # Sort by the "R" metric
#' linechart(vk12, sort = "R")
#' 
#' # Sort by the "R" metric and include only the first 3 treatments
#' linechart(vk12, sort = "R", n = 3)
#' 
#' # Exclude the "R" metric
#' #
#' linechart(vk12, exclude = "R")
#' 
#' @export linechart

linechart <- function(x, 
                      sort = "Q",
                      exclude = "none",
                      n = nrow(x),
                      linewidth = 1.1,
                      size = 2,
                      ...) {
  
  chkclass(x, "vikor")
  #
  sort <- setchar(sort, c("Q" , "S", "R"))
  exclude <- setchar(exclude, val = c("none", "Q", "S", "R"))
  #
  chknumeric(n, min = 1, max = nrow(x), length = 1)
  chknumeric(linewidth, min = 0, zero = TRUE, length = 1)
  chknumeric(size, min = 0, zero = TRUE, length = 1)
  #
  if (sort == exclude)
    stop("The excluded list cannot coincide with the list used for sorting.",
         call. = FALSE)
  #
  chknumeric(n, min = 0)
  
  treat <- rep(row.names(x), 3)
  values <- c(x$Q, x$S, x$R)
  
  type <- rep(c("Q", "S", "R"), each = length(treat) / 3)
  
  dat <- data.frame(treat, type, values)
  #
  dat.sort <- dat %>% 
    filter(type == sort) %>% 
    arrange((values))
  
  dat.sort <- dat.sort[seq_len(n), ]
  sort.treat <- dat.sort$treat
  levels.list <- c("Q", "S", "R")
  
  if (exclude == "none") {
    dat <- dat %>% 
      filter(treat %in% sort.treat) %>% 
      mutate(type = factor(type, levels = levels.list)) %>% 
      mutate(treat = factor(treat, levels = sort.treat))
  }
  else {
    l <- which(levels.list == exclude)
    #
    levels.list <- levels.list[-l]
    #
    dat <- dat %>%
      filter(type != exclude) %>% 
      filter(treat %in% sort.treat) %>% 
      mutate(type = factor(type, levels = levels.list)) %>% 
      mutate(treat = factor(treat, levels = sort.treat))
  }
  
  p <- ggplot(dat, aes(x = treat, y = values, color = type, group = type)) +
    geom_line(linewidth = linewidth) +
    geom_point(size = size) +
    theme_minimal() +
    xlab("") +
    ylab("Value") +
    ylim(c(0, 1)) +
    guides(color = guide_legend(title = "Metric"))
  #
  attr(p, "data") <- dat
  #
  p
}

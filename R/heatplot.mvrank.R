#' A heatplot to visualize the ranking in a multivariate network meta-analysis
#' 
#' @description
#' This function produces a heatplot displaying the results of function
#' \code{\link{mvrank}}. The graph can be used to visualize the ranking output
#' when the method used to rank the treatments is either the surface under the
#' cumulative ranking curve (SUCRA) or the probability of best value (pBV).
#' 
#' @param x An object of class \code{\link{mvrank}}.
#' @param sort An optional argument to define an outcome to be used as a
#'   reference when sorting the order of treatments on the x-axis.
#'   If \code{NULL}, then the first outcome, as it appears in
#'   \code{\link{mvrank}}, is used as the reference. The user can either
#'   provide a single number corresponding to the position of the outcome in
#'   the \code{\link{mvrank}} object or a character string with the name of
#'   the outcome; can be abbreviated.
#' @param col.num The color of the numbers in the squares (default: "white").
#' @param num.size The size of the numbers in the squares (default: 4.5).
#' @param col.low The color for the low end of the gradient
#'   (default: "lightblue").
#' @param col.high The color for the high end of the gradient
#'   (default: "darkblue").
#' @param width.bar The width of the bars in the graph (default: 1.2).
#' @param angle.x The angle (in [0, 360]) to rotate the text appearing on
#'   the x-axis (default: 0).
#' @param angle.y The angle (in [0, 360]) to rotate the text appearing on the
#'   y-axis (default: 0).
#' @param hjust.x Horizontal justification for the text on the x-axis
#'   (default: 0.5).
#' @param hjust.y Horizontal justification for the text on the y-axis
#'   (default: 1).
#' @param hjust.legend Horizontal justification for the text in the legend
#'   title (default: 0.5).
#' @param size.x Text size (in pts) on the x-axis (default: 8).
#' @param size.y Text size (in pts) on the y-axis (default: 10).
#' @param legend.position Position of the legend. Options are "top", "bottom",
#'   "right" (default), and "left"; can be abbreviated.
#' @param legend.direction Direction of the legend. Options are "vertical"
#'   (default) and "horizontal"; can be abbreviated.
#' @param \dots Additional arguments (ignored).
#'
#' @return
#' A \code{ggplot} object.
#' 
#' @examples
#' # Locate file "mvnma_example.rda" with mvnma() results
#' .fname <- system.file("extdata/mvnma_examples.rda", package = "mvnma")
#' load(.fname)
#' 
#' # Rank the treatments
#' ranks_sucra <- mvrank(mvnma_all, 
#'   small.values = c("undes", "undes", "des", "des", "des"),
#'   method = "SUCRA")
#' ranks_sucra
#' 
#' # Create a heatplot sorting the results according to the first outcome
#' # (default), i.e., Early_Response
#' heatplot(ranks_sucra)
#' 
#' # Create a heatplot sorting the results according to the second outcome
#' # appearing in ranks_sucra (i.e Early_Remission)
#' heatplot(ranks_sucra, sort = 2)
#' 
#' # Create a heatplot sorting the results by explicitly mentioning the
#' # name of the outcome Early_Remission
#' heatplot(ranks_sucra, sort = "Early_Remission")
#' 
#' @method heatplot mvrank 
#' @export

heatplot.mvrank <- function(x,
                            sort = NULL,
                            col.num = "white",
                            num.size = 4.5,
                            col.low = "lightblue",
                            col.high = "darkblue",
                            width.bar = 1.2,
                            angle.x = 0,
                            angle.y = 0, 
                            hjust.x = 0.5,
                            hjust.y = 1,
                            hjust.legend = 0.5,
                            size.x = 8,
                            size.y = 10,
                            legend.position = "right",
                            legend.direction = "vertical",
                            ...) {
  
  chkclass(x, "mvrank")
  #
  method <- attr(x, "method")
  #
  if (!(method %in% c("SUCRA", "pBV")))
    stop("Heatplot can only be produced for mvnma object created with ",
         "'method = \"SUCRA\"' or \"pBV\".",
         call. = FALSE)
  
  outcome <- treatment <- val <- NULL
  #
  legend.position <-
    setchar(legend.position, c("top", "bottom", "right", "left"))
  #
  legend.direction <- setchar(legend.direction, c("horizontal", "vertical"))
  
  # Create a data frame with all ranks
  #
  dat.ranks <- vector("list")
  
  for (i in 1:length(x)) {
    dat.ranks[[i]] <- x[[i]] 
    dat.ranks[[i]]$outcome <- names(x)[i]
  }
  #
  dat.ranks <- bind_rows(dat.ranks)
  
  # Make sure that argument sort is either a character (can be abbreviated) or
  # a numeric value
  #
  if (is.null(sort)) {
    sort <-  names(x)[1]
    order <- dat.ranks %>% filter(outcome==sort)
  }
  else {
    if (length(sort) > 1) {
      stop("Argument 'sort' must be of length 1.", call. = FALSE)
    }
    #
    if ((!is.numeric(sort)) & (!is.character(sort))) {
      stop("Argument 'sort' must be either a character string or a ",
           "numeric value of length 1.", call. = FALSE)
    }
    else if (is.character(sort)) {
      sort <- setchar(sort, names(x))
      order <- dat.ranks %>% filter(outcome == sort)
    }
    else if (is.numeric(sort)) {
      if (sort > length(x)) {
        stop("Argument 'sort' must be a number between ", 1," and ", length(x),
             ".", call. = FALSE)
      }
      #
      sort <- names(x)[sort] 
      order <- dat.ranks %>% filter(outcome==sort)
    }
  }
  
  # Create order for treatments and outcomes
  #
  order <- order$treatment
  #
  dat.ranks$treatment <- factor(dat.ranks$treatment, levels = order)
  #
  dat.ranks$outcome <- factor(dat.ranks$outcome, levels = names(x))
  #
  dat.ranks$outcome <- relevel(dat.ranks$outcome, ref = sort)
  #
  names(dat.ranks) <- c("treatment","val","outcome")
  
  # Create the heatplot
  #
  p <- ggplot(dat.ranks, aes(treatment, fct_rev(outcome), fill = val)) + 
    geom_tile(color = "black") +
    geom_text(aes(label = paste(format(round(val,digits = 2), nsmall = 2))),
              size = num.size, color = col.num) +
    scale_fill_gradient(low = col.low, high = col.high) + 
    guides(fill = guide_colourbar(label = FALSE, ticks = FALSE,
                                  barwidth = width.bar))+
    xlab("") +
    ylab("") +
    labs(fill = method) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_void() +
    theme(
      # bold font for legend text
      legend.text = element_text(face = "bold"),
      # set thickness of axis ticks
      axis.ticks.y = element_blank(),
      axis.text = element_text(face = "bold"),
      #remove plot background
      plot.background = element_blank(),
      #remove plot border
      panel.border = element_blank(),
      axis.text.x = element_text(angle = angle.x, hjust = hjust.x,
                                 size = size.x),
      axis.text.y = element_text(angle = angle.y, hjust = hjust.y,
                                 size = size.y),
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = element_text(hjust = hjust.legend,)
    )
  #
  attr(p, "dat.ranks") <- dat.ranks
  attr(p, "method") <- method
  #
  p
}

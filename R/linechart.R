#' Linechart showing the results of \code{\link{vikor}} across the three metrics Q, S and R.
#'
#' @description
#' A linechart showing the results of \code{\link{vikor}} across the three metrics Q, S and R for each treatment. Within each metric, lower values represent better treatment performance.
#' 
#' @param x An object of class \code{\link{vikor}}.
#' @param sort.by A character specifying the order of treatments on the x-axis. By default, the order
#' is according to the Q-metric (\code{"Q"}). Ordering treatments according to the S (\code{"S"}) and R (\code{"R"}) metrics is also possible.
#' @param exclude A character specifying a metric that will not be displayed in the graph. By default, all metrics are displayed (\code{"none"}). 
#' Alternative options are to exlude either the Q-metric (\code{"Q"}), or the S-metric (\code{"S"}) or the R-metric (\code{"R"}).
#' @param k A numeric value indicating the number of treatments to be plotted in the graph. By default, all treatments are displayed. If specified,
#' only the first \code{"k"} treatments according to the hierarchy specified from \code{"sort.by"} argument are plotted.
#' @param linewidth A numeric value specifying the width of the lines  (default: 1.1).
#' @param point.size A numeric value specifying the size of the points (default: 2).
#' @param \dots Additional arguments passed to \code{\link[ggplot2]{ggplot}} function.
#'
#' @return
#' A \code{ggplot} object.
#' 
#' @examples
#' \donttest{
#' library("netmeta")
#' 
#' # Use 'pairwise' to obtain contrast based data for the first two outcomes
#' data("Linde2015")
#' # Early response
#' p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#' # Early remissions
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Define outcome labels
#' outcomes <- c("Early_Response", "Early_Remission")
#'  
#' # Fit the model combining only the two efficacy outcomes
#' set.seed(1909)
#' mvnma12 <- mvnma(p1, p2,
#'   reference.group = "Placebo", outclab = outcomes,
#'   n.iter = 1000, n.burnin = 100)
#' mvnma12
#' 
#' # Rank treatments using SUCRAs
#' ranks12 <- mvrank(mvnma12, small.values = c("und", "und"), method = "sucra")
#' ranks12
#' 
#' # Get the best compromise solution across the efficacy outcomes
#' vk <- vikor(ranks12)
#' 
#' # Visualize the results with default settings
#' linechart(vk)
#' # sort by the "R" metric
#' linechart(vk,sort.by = "R")
#' # sort by the "R" metric and include only the first 3 treatments
#' linechart(vk,sort.by = "R",k=3)
#' # Exclude the "R" metric
#' linechart(vk,exclude = "R")
#' }
#' @export linechart

linechart <- function(x, 
                       sort.by = "Q",
                       exclude = "none",
                       k=nrow(x),
                       linewidth = 1.1,
                       point.size = 2,
                       ...){
  
chkclass(x, "vikor")
  
sort.by <- setchar(sort.by,val = c("Q","S","R"))

exclude <- setchar(exclude,val = c("none","Q","S","R"))

if(sort.by == exclude){
  
  stop("The excluded list cannot coincide with the list sorting by.")

}

chknumeric(k, min = 0)
  
treat <- rep(row.names(x),3)

values <- c(x$Q,x$S,x$R)

type <- rep(c("Q","S","R"),each = length(treat)/3)

dat <- data.frame("treat"= treat,
                  "type" = type,
                  "values" = values)

dat.sort <- dat %>% 
filter(type==sort.by) %>% 
arrange((values))

dat.sort <- dat.sort[1:k,]

sort.treat <- dat.sort$treat

levels.list <- c("Q","S","R")

if(exclude == "none"){

dat <- dat %>% 
  filter(treat %in% sort.treat) %>% 
  mutate(type = factor(type, levels = levels.list)) %>% 
  mutate(treat = factor(treat, levels = sort.treat))

}else{
  
  l <- which(levels.list == exclude)
  
  levels.list <- levels.list[-l]
  
  dat <- dat %>%
    filter(type != exclude) %>% 
    filter(treat %in% sort.treat) %>% 
    mutate(type = factor(type, levels = levels.list)) %>% 
    mutate(treat = factor(treat, levels = sort.treat))
    
}


graph <- ggplot(dat, aes(x = treat, y = values, color = type,group = type)) +
  geom_line(linewidth = linewidth) +
  geom_point(size = point.size) +
  theme_minimal() +
  xlab("") +
  ylab("Value") +
  ylim(c(0, 1)) +
  guides(color = guide_legend(title = "Metric"))


attr(graph, "data") <- dat

graph

}



#' A heatplot to visualize the output of \code{\link{mvrank}} across all outcomes
#' 
#' @description
#' This function produces a heatplot displaying the results of function \code{\link{mvrank}}. 
#' The graph can be used to visualize the ranking output when the method used to rank the treatments
#' is either the SUCRA or the pBV method.
#' 
#' @param x An object of class \code{\link{mvrank}}.
#' @param sortby An optional argument to define an outcome to be used as a reference when sorting
#' the order of treatments on the x-axis. If \code{NULL}, then the first outcome, as it appears in \code{\link{mvrank}},
#' is used as the reference. If not \code{NULL}, the user can specify this argument using a single numeric value indicating
#' the position of the outcome in the order of \code{\link{mvrank}} (i.e., setting it equal to 1 will use the first outcome as the reference, as it appears in \code{\link{mvrank}}). Alternatively, the argument can be specified
#' using a character string explicitly stating the name of the outcome. In this case, the name can be abbreviated.
#' @param col.num The color of the numbers in the squares (default: "white").
#' @param num.size The size of the numbers in the squares (default: 4.5).
#' @param col.low The color for the low end of the gradient (default: "lightblue").
#' @param col.high The color for the high end of the gradient (default: "darkblue").
#' @param width.bar The width of the bars in the graph (default: 1.2).
#' @param angle.x The angle (in [0, 360]) to rotate the text appearing on the x-axis (default: 0).
#' @param angle.y The angle (in [0, 360]) to rotate the text appearing on the y-axis (default: 0).
#' @param hjust.x Horizontal justification for the text on the x-axis (default: 0.5).
#' @param hjust.y Horizontal justification for the text on the y-axis (default: 1).
#' @param hjust.legend Horizontal justification for the text in the legend title (default: 0.5).
#' @param size.x Text size (in pts) on the x-axis (default: 8).
#' @param size.y Text size (in pts) on the y-axis (default: 10).
#' @param legend.position Position of the legend. Options are "top", "bottom", "right" (default), and "left". Can be abbreviated.
#' @param legend.direction Direction of the legend. Options are "vertical" (default) and "horizontal". Can be abbreviated.
#' @param \dots Additional arguments passed to \code{\link[ggplot2]{ggplot}}.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' library("netmeta")
#' 
#' data("Linde2015")
#' 
#' # Use 'pairwise' to obtain contrast based data for each one of the five
#' # available outcomes 
#'
#' # Early response
#' p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Early remissions
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Adverse events
#' p3 <- pairwise(treat = list(treatment1, treatment2,treatment3),
#'   event = list(ae1, ae2, ae3),  n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Loss to follow-up
#' p4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss1, loss2, loss3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Loss_to_follow_up_(AE)
#' p5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss.ae1, loss.ae2, loss.ae3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Define outcome labels
#' outcomes <- c("Early_Response", "Early_Remission",
#'   "Adverse_events", "Loss_to_follow_up", "Loss_to_follow_up_AE")
#' 
#' # Fit the model combining all five outcomes
#' set.seed(1909)
#' mvnma_all <- mvnma(p1, p2, p3, p4, p5,
#'   reference.group = "Placebo", outclab = outcomes,
#'   n.iter = 1000, n.burnin = 100)
#' # Rank the treatments
#'   ranks_sucra <- mvrank(mvnma_all, 
#'   small.values = c("undes", "undes", "des", "des", "des"),
#'   method = "SUCRA")
#' ranks_sucra
#' # Create a heatplot sorting the results according to the first outcome appearing in ranks_sucra (i.e Early_Response)
#'   heatplot(ranks_sucra,sortby = 1)
#' # Create a heatplot sorting the results by explicitly mentioning the name of the outcome Early_Response
#'   heatplot(ranks_sucra,sortby = "Early_Response") 
#'   }
#' @method heatplot mvrank 
#' @export

heatplot.mvrank <- function(x,
                    sortby=NULL,
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
                    ...
                    ){
  
  chkclass(x, "mvrank")
  
  if (!(attr(x, "method") %in% c("SUCRA", "pBV")))
    stop("Heatplot can only be produced for ",
         "'method = \"SUCRA\"' and 'method = \"pBV\"'.",
         call. = FALSE)
  
  legend.position <- setchar(legend.position,c("top","bottom","right","left"))
  
  legend.direction <- setchar(legend.direction,c("horizontal","vertical"))
  
  method <- attributes(x)$method
  
  ## create a dataframe with all ranks
  dat.ranks <- vector("list")
  
  for(i in 1:length(x)){
  
  dat.ranks[[i]] <- x[[i]] 
  
  dat.ranks[[i]]$outcome <- names(x)[i]
  
  }
  
  dat.ranks <- bind_rows(dat.ranks)
  
  # make sure that argument sortby is either a character (can be abbreviated) or a numeric value
  
  if(is.null(sortby)){
    
    sortby <-  names(x)[1]
    
    order <- dat.ranks %>% filter(outcome==sortby)
    
  }else{
  
  if(length(sortby)>1){
    
    stop("Argument 'sortby' should be of length 1.")
  }
  
  if((!is.numeric(sortby)) & (!is.character(sortby))){
    
   stop("Argument 'sortby' should be either a character string or a numeric value of length 1.") 
  
    }else if(is.character(sortby)){  
  
    sortby <- setchar(sortby,names(x))
    
    order <- dat.ranks %>% filter(outcome==sortby)
    
    }else if(is.numeric(sortby)){
    
      sortby <- sortby

      if(sortby>length(x)){
        stop(paste("Argument 'sortby' should be a number between ",1," and ",length(x),sep = ""))
      }
      
      sortby <- names(x)[sortby] 
      order <- dat.ranks %>% filter(outcome==sortby)
    }
  }
  
  ## create order for treatments and outcomes
  order <- order$treatment
  
  dat.ranks$treatment <- factor(dat.ranks$treatment,levels = order)
  
  dat.ranks$outcome <- factor(dat.ranks$outcome,levels = names(x))
  
  dat.ranks$outcome <- relevel(dat.ranks$outcome,ref = sortby)
  
  names(dat.ranks) <- c("treatment","val","outcome")
  
  # create the heatplot
  g <- ggplot(dat.ranks, aes(treatment, fct_rev(outcome), fill= val)) + 
    geom_tile(color = "black") +
    geom_text(aes(label = paste(format(round(val,digits = 2),nsmall=2) )),
              size=num.size,color=col.num) +
    
    scale_fill_gradient(low=col.low, high=col.high) + 
  
    guides(fill = guide_colourbar(label = FALSE,
                                  ticks = FALSE,
                                  barwidth = width.bar))+
    xlab("")+
    ylab("")+
    labs(fill= method)+
    scale_y_discrete(expand=c(0,0))+
    theme_void()+
    theme(
      #bold font for legend text
      legend.text = element_text(face="bold"),
      #set thickness of axis ticks
      axis.ticks.y = element_blank(),
      axis.text = element_text(face = "bold"),
      #remove plot background
      plot.background=element_blank(),
      #remove plot border
      panel.border=element_blank(),
      axis.text.x=element_text(angle=angle.x,hjust=hjust.x,size=size.x),
      
      axis.text.y=element_text(angle = angle.y,hjust=hjust.y,size=size.y),
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = element_text(hjust = hjust.legend,)
      
    )
  
  attr(g,"dat.ranks") <- dat.ranks
  attr(g,"method") <- method
  
 g
  
}


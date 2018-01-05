#
# Author: Gael Yvert, CNRS
# Gael.Yvert@ens-lyon.fr
#
#---------------------------------------------
#' Indexes bins.
#' @param x A value between 0 and 1
#' @param n An integer indicating the number of bins between 0 and 1
#' @return The index of the bin containing x, or NA if x > 1 or x > 0
#' @export

x2i <- function(x,n){
  seg <- seq(0,1,length.out=n+1);
  if (x == 0)
    out = 1
  else if (x > 1)
    out = NA
  else if (x < 0)
    out = NA
  else
    out = max(which(seg<x))
  return(out);
}
#---------------------------------------------
#' Estimates ordinate value of the density at position x.
#' @param d a density function, as returned by density()
#' @param x a numerical value
#' @return the estimated ordinate value of the density at position x
#' @export

density_value <- function(x, d){
  n = length(d$x);
  if (x > max(d$x) | x < min(d$x) ) #out of bounds
    out = 0
  else if (x == min(d$x))
    out = d$y[which.min(d$x)]
  else{
    i = max(which(d$x < x))
    out = 0.5*(d$y[i+1] + d$y[i])
  }
  return(out);
}

#---------------------------------------------
#' Translate density into colour, for graphix representation.
#' @param d A density function, as returned by density()
#' @param ylim The limite of the drawned box on the Y axis, default is range(d$x)
#' @param xlim The limite of the drawned box on the X axis, default is c(0,1)
#' @param nstripes (integer) The number of horizontal stripes to plot per experiment, default is 1000
#' @param colors The colors to use to construct the ramp, default is c("white", "darkblue")
#' @param ncolors (integer) The number of gradual colors of the ramp, default is 100
#' @return The colors and coordinates to use to plot a given density
#' @examples
#' d = density(x=1:1000)
#' ylim = range(d$x)
#' xlim = c(0,1)
#' nstripes = 1000
#' colors = c("white", "darkblue"),
#' ncolors = 100
#' plot(1:2, xlim = xlim, ylim = ylim, type = "n")
#' draw.one.density (d, ylim, xlim, nstripes, colors, ncolors)
#' @export

draw.one.density <- function(d, ylim = range(d$x), xlim = c(0,1), nstripes = 1000, colors = c("white", "darkblue"), ncolors = 100){

  # coordinates of rectangles to draw
  xleft  = rep(xlim[1], times = nstripes)
  xright = rep(xlim[2], times = nstripes)
  y.all = seq(from = ylim[1], to = ylim[2], length.out = nstripes + 1)
  ybottom = y.all[1:nstripes]
  ytop = y.all[2:(nstripes + 1)]

  # prepare the palette
  colR <- colorRampPalette(colors)
  col <- colR(ncolors)

  # translate density into color
  ymid = 0.5*(ybottom+ytop) # the GFP intensity value to consider for every stripe
  yval <- sapply(ymid, FUN = density_value, d) # the density fo cells at these GFP values
  yval <- yval/max(yval) # scaling to have densities from 0 to 1
  coloridx <- sapply(yval, FUN = x2i, n = ncolors) # index of color corresponding to the density value
  rect.col = col[coloridx]

  # draw
  rect(xleft, ybottom, xright, ytop, border = NA, col = rect.col);
}

#---------------------------------------------
#' Draws a list of densities according to an index of color.
#' @param ld A list of density functions
#' @param nstripes (integer) The number of horizontal stripes to plot per experiment, default is 1000
#' @param colors The colors to use to construct the ramp, default is c("white", "darkblue")
#' @param ncolors (integer) The number of gradual colors of the ramp, default is 100
#' @param xpos A vector of positions on the x-axis to plot each experiment.
#       If NULL (default), then experiments are plotted side by side in their index order
#' @param ... Additional arguments to pass to the plot() function
#' @return A plot of a serie of density profiles according to an index of color
#' @examples
#' ld = list(density(x=1:1000), density(x=1:100), density(x=1:10))
#' ylim = range(d$x)
#' xlim = c(0,1)
#' nstripes = 1000
#' colors = c("white", "darkblue"),
#' ncolors = 100
#' dyn.profile (ld)
#' @export

dyn.profile <- function(ld, colors = c("white", "darkblue"), ncolors = 100, nstripes = 1000, xpos = NULL, ylim = NULL, ...){
  n = length(ld);
  # get the y-axis limits
  if (is.null(ylim)){
    ranges <- NULL
    for (i in 1:n)
      ranges <- c(ranges, range(ld[[i]]$x))
    ylim = range(ranges)
  }
  # get the x-axis limits
  if (is.null(xpos))
    xlim = c(0,n)
  else{
    if ( length(unique(xpos)) != length(xpos))
      stop("the following vector of positions has tied values:\n", paste(xpos, collapse = ";"), "\n dyn.profile() aborted", sep = " ");
    #determine width of vertical bands
    sx <- sort(xpos)
    x1 = c(NA, sx)
    x2 = c(sx, NA)
    bw = 0.5*min(x2-x1, na.rm = TRUE);
    #set limits
    xlim = c(min(xpos) - bw, max(xpos) + bw);
  }

  # draw
  if (is.null(xpos)){
    # initialize plot
    plot(1:2, xlim = xlim, ylim = ylim, type = "n", ...)
    for (i in 1:n)
      draw.one.density(ld[[i]], ylim = ylim, xlim = c(i-1,i), nstripes = nstripes, colors = colors, ncolors = ncolors)
  }
  else{
    # initialize plot
    plot(1:2, xlim = xlim, ylim = ylim, type = "n", ...)
    for (i in 1:n){
      draw.one.density(ld[[i]], ylim = ylim, xlim = c(xpos[i] - bw, xpos[i] + bw), nstripes = nstripes, colors = colors, ncolors = ncolors)
    }
  }
}

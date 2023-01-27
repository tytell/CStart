require(dplyr)
require(gsignal)

#' Estimates the location of the center of mass in a single frame
#' 
#' `get_com_in_frame` returns a single value, which is the estimated location
#' of the center of mass of a fish, given the location of the midline and the
#' width, height, or mass of each segment along the body.
#' 
#' You can pass in only the body width, in which case it will estimate COM
#' based on the horizontal projected area (assuming constant body height), or the
#' width and height, in which case it assumes that the body is constant density
#' and made up of segments of cones with oval cross-sections.
#' 
#' @param mx x coordinate of the midline
#' @param my y coordinate of the midline
#' @param width Width of the body (in fractions of the body length) (optional).
#' @param height Height of the body (in fractions of body length) (optional).
#' @param segmentmass Mass of each body segment. Should have one fewer element
#' than the number of coordinates along the midline.
#' 
#' @returns A list (`comx` and `comy`) with the x and y coordinates of the 
#' center of mass.
get_com_in_frame <- function(mx, my, width = NULL, height = NULL, segmentmass = NULL)
{
  if (any(is.na(mx)) || any(is.na(my))) {
    return(list(comx = NA,
                comy = NA))
  }

  if (is.null(segmentmass) && is.null(width)) {
    stop("get_com_in_frame requires at least the width of the body or the segmentmass")
  }
  
  if (any(width > 0.5) || any(height > 1.0)) {
    warning("width and height should be in fractions of body length.")
  }
  
  npt <- length(mx)
  
  if (is.null(segmentmass)) {
    width <- width / 2
    
    if (is.null(height)) {
      height <- rep_along(width, mean(width))
    }
    else if (length(height) == 1) {
      height <- rep_along(width, height)
    }
    else {
      height <- height / 2
    }
    
    dw <- diff(width)
    dh <- diff(height)
    dl <- 1 / (npt-1)
    
    width <- width[1:(npt-1)]
    height <- height[1:(npt-1)]
    
    segmentmass <- pi * (width * height * dl +
                           0.5 * (dw / dl * height + dh / dl * width) * dl^2 +
                           1/3 * dw * dh * dl)
  }
  else if (length(segmentmass) != npt-1) {
    stop("get_com_in_frame: segmentmass should have one less than the number of body points")
  }
  
  s0 <- c(0, cumsum(diff(mx)^2 + diff(my)^2))
  
  # get the midpoint of each segment
  s1 <- (s0[1:(npt-1)] + s0[2:npt]) / 2
  
  # interpolate the location of the midpoint
  q <- approx(s0, mx, s1)
  mx1 <- q$y
  
  q <- approx(s0, my, s1)
  my1 <- q$y
  
  totalmass <- sum(segmentmass)
  
  return(list(comx = sum(mx1 * segmentmass) / totalmass,
              comy = sum(my1 * segmentmass) / totalmass))
}

#' Estimates the location of the center of mass in a grouped tibble
#' 
#' Estimates the location of the center of mass for many midlines in a tibble.
#' 
#' `get_com` is meant to be called using `group_modify`
#' 
#' @param df The tibble or data frame
#' @param xvar The name of the x coordinate in the data frame
#' @param yvar The name of the y coordinate in the data frame
#' @param width Width of the body (in fractions of the body length) (optional).
#' @param height Height of the body (in fractions of body length) (optional).
#' @param segmentmass Mass of each body segment. Should have one fewer element
#' than the number of coordinates along the midline.
#'
#' @return A data frame with additional columns `comx` and `comy`.
#' 
#' @examples 
#' data |> 
#'   group_by(t) |> 
#'   group_modify(~ get_com(.x, mx,my, width=width))
#'
#' @seealso get_com_in_frame
get_com <- function(df, xvar, yvar, width=NULL, height=NULL, segmentmass=NULL)
{
  a <- pull(df, {{xvar}})
  b <- pull(df, {{yvar}})
  
  comxy <- get_com_in_frame(a, b, width = width, height = height,
                            segmentmass = segmentmass)
  as_data_frame(comxy)
}

#' Estimates the location of the stretched-straight center of mass
#' 
#' `get_stretched_straight_com` estimates the location of the center of mass
#' of a fish when the body is straight, assuming constant body density.
#' 
#' Uses the width and height of the body (or just the width) to estimate the
#' distance of the COM from the snout when the body is straight.
#' 
#' @param width Width of the body (in fractions of the body length).
#' @param height Height of the body (in fractions of body length).
#' 
#' @return The distance of the COM from the snout.
get_stretched_straight_com <- function(width, height)
{
  npt <- length(width)
  mx <- (seq_along(width)-1) / (npt-1)
  my <- rep_along(width, 0)
  
  comxy <- get_com_in_frame(mx,my, width=width, height=height)
  return(comxy$comx)
}

#' Interpolates the position of the stretched-straight center of mass
#' 
#' `interp_stretched_straight_com_in_frame` interpolates the x,y location
#' of the stretched-straight center of mass, given a single curved midline.
#' 
#' @param s Arc length along the body
#' @param mx x coordinate of the midline
#' @param my y coordinate of the midline
#' @param ssdist Distance of the stretched-straight center of mass from the
#' snout, in the same units as `s`
#' 
#' @return A list containing `comx` and `comy`, the x and y coordinates of COM.
interp_stretched_straight_com_in_frame <- function(s, mx,my, ssdist)
{
  x1 <- approx(s,mx, ssdist)
  y1 <- approx(s,my, ssdist)
  
  return(list(comx = x1$y, comy = y1$y))
}

#' Interpolates the stretched-straight location of the center of mass in a grouped tibble
#' 
#' Interpolates the location of the stretched-straight center of mass for 
#' many midlines in a tibble.
#' 
#' `interp_stretched_straight_com` is meant to be called using `group_modify`
#' 
#' @param df The tibble or data frame
#' @param xvar The name of the x coordinate in the data frame
#' @param yvar The name of the y coordinate in the data frame
#' @param ssdist Distance of the stretched-straight center of mass from the
#' snout, in the same units as `s`
#'
#' @return A data frame with additional columns `comx` and `comy`.
#' 
#' @examples 
#' data |> 
#'   group_by(t) |> 
#'   group_modify(~ interp_stretched_straight_com(.x, mx,my, ssdist=0.36))
#'
#' @seealso interp_stretched_straight_com_in_frame
interp_stretched_straight_com <- function(df, xvar, yvar, ssdist)
{
  a <- pull(df, {{xvar}})
  b <- pull(df, {{yvar}})

  if (any(is.na(a)) || any(is.na(b))) {
    return(data.frame(comx = NA,
                      comy = NA))
  }
  
  s <- c(0, cumsum(diff(a)^2 + diff(b)^2))
  
  comxy <- interp_stretched_straight_com_in_frame(s, a,b, ssdist*s[length(s)])
  
  return(as_data_frame(comxy))
}

#' Smooths one position variable
#' 
#' Smooths one position variable using a low-pass Butterworth filter.
#' 
#' The cutoff frequency of the filter is defined by the filter duration, where
#' fcut = 1/filterdur. If `stage1dur` is given, the filter duration is 25% of
#' the stage 1 duration.
#' 
#' Requires one of `smoothdur`, `stage1dur`, or `filt`
#' 
#' @param x Variable to smooth
#' @param samplerate Sample rate for the data set (Hz)
#' @param smoothdur Duration of the filter (sec) (optional)
#' @param stage1dur Duration of stage 1 for a C-start. Filter duration will be 
#'   1/4 of the stage 1 duration. (sec) (optional)
#' @param filt Filter structure from filter design functions, such as 
#'   `butter` (optional)
#' @param endmode Method for smoothing data at the ends of the time series.
#'   Could be "none" or "rev", which means that a reversed version of the data set
#'   will be appended at the beginning and end before filtering.
#'   
#' @return The smoothed data vector
smooth_1_position <- function(x, samplerate=NULL, 
                            smoothdur=NULL, stage1dur=NULL,
                            filt=NULL, endmode = "none")
{
  if (is.null(smoothdur) && is.null(stage1dur) && is.null(filt)) {
    stop("smooth_position requires either smoothdur or stage1dur or a filter structure")
  }
  else if (!is.null(stage1dur)) {
    smoothdur <- stage1dur / 4
  }

  if (is.null(filt)) {
    if (is.null(samplerate)) {
      stop("smooth_1_position requires samplerate to design a filter")
    }
    filt <- gsignal::butter(9, (1/smoothdur)/(samplerate/2), type="low",
                            output="Sos")
  }
  
  good <- !is.na(x)
  n <- sum(good)
  
  xg <- x[good]
  
  if (endmode == "none") {
    xlong <- c(xg)
    xslong <- gsignal::filtfilt(filt, xlong)
    
    smoothrange <- seq(1, n)
  } else if (endmode == "rev") {
    xlong <- c(rev(xg), xg, rev(xg))
    xslong <- gsignal::filtfilt(filt, xlong)

    smoothrange <- seq((n+1), (2*n))
  }
  else {
    stop("endmode must be one of \"none\" or \"rev\"")
  }
  
  xs <- rep_along(x, NA_real_)
  xs[good] <- xslong[smoothrange]
  
  return(xs)
}

#' Smooths x and y positions together
#' 
#' Smooths x and y positions in a data frame using a Butterworth filter in a
#' grouped tibble. Meant to be called using `group_modify`
#' 
#' `smooth_positions` smooths x and y positions using a Butterworth filter.
#' 
#' The cutoff frequency of the filter is defined by the filter duration, where
#' fcut = 1/filterdur. If `stage1dur` is given, the filter duration is 25% of
#' the stage 1 duration.
#' 
#' Requires one of `smoothdur`, `stage1dur`, or `filt`
#' 
#' @param xvar Name of the x variable to smooth
#' @param yvar Name of the y variable to smooth
#' @param samplerate Sample rate for the data set (Hz)
#' @param smoothdur Duration of the filter (sec) (optional)
#' @param stage1dur Duration of stage 1 for a C-start. Filter duration will be 
#'   1/4 of the stage 1 duration. (sec) (optional)
#' @param filt Filter structure from filter design functions, such as 
#'   `butter` (optional)
#'   
#' @return The tibble with new variables called "{xvar}.s" and "{yvar}.s", where
#'   "{xvar}" and "{yvar}" are whatever the original names for the x and y
#'   positions.
smooth_positions <- function(df, xvar,yvar, samplerate, 
                             smoothdur=NULL, stage1dur=NULL,
                             filt=NULL)
{
  if (is.null(smoothdur) && is.null(stage1dur) && is.null(filt)) {
    stop("smooth_position requires either smoothdur or stage1dur or a filter structure")
  }
  else if (!is.null(stage1dur)) {
    smoothdur <- stage1dur / 4
  }
  
  if (is.null(filt)) {
    if (is.null(samplerate)) {
      stop("smooth_1_position requires samplerate to design a filter")
    }
    filt <- gsignal::butter(9, (1/smoothdur)/(samplerate/2), type="low",
                            output="Sos")
  }

  df |> 
    mutate("{{xvar}}.s" := smooth_1_position({{xvar}}, filt=filt),
           "{{yvar}}.s" := smooth_1_position({{yvar}}, filt=filt))
}

#' Numerical first derivative with a central difference algorithm
#' 
#' Takes a first derivative of the `y` variable using a simple central difference
#' algorithm. Does not assume equal spacing in `t`
#' 
#' @param t Time, or the x variable
#' @param y Y variable to differentiate
#' @returns The first derivative dy/dt as a vector with the same size as y.
nderivative <- function(t, y)
{
  (lead(y) - lag(y)) / (lead(t) - lag(t))
}

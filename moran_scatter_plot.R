#'
#' @author Lukas Dargel, modified by Thibault Laurent
#'
#' @param X a numeric vector
#' @param geo a list of polygon geographies
#' @param my_listw a listw object
#' @param X_name a character indicating the name of X
#' @param only_scatter a logical if `TRUE` only the scatter plot is shown
#' @param only_mapa a logical if `TRUE` only the map is shown
#'
#' @return a ggplot graphic
moran_scatter_plot <- function(
    X,
    geo,
    my_listw,
    X_name = "X",
    only_scatter,
    only_map,
    cols = 2,
    ...) {
  
  stopifnot(require("colorspace", quietly = TRUE),
            require("ggplot2", quietly = TRUE),
            require("patchwork", quietly = TRUE),
            require("spdep", quietly = TRUE),
            require("sf", quietly = TRUE))
  
  
  # W <- nb2listw(poly2nb(geo))
  X <- scale(X)
  WX <- lag.listw(my_listw, X)
  lonlat <- st_coordinates(st_centroid(geo))
  
  mp <- data.frame(X, WX, lon = lonlat[,1], lat = lonlat[,2])
  #colnames(mp) <- c(paste0(X_name, c(""," (lag)")),"Latitude", "Longitude")
  
  mp_scatter <- ggplot(mp) +
    geom_point(aes(
      x = X, y = WX, fill = X, col = WX,
      shape = if_else(0 > sign(X * WX), "Pos", "Neg"),
      size = sqrt(X ^ 2 + WX ^ 2),
      stroke = sqrt(X ^ 2 + WX ^ 2))) +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    geom_abline(slope = cor(X,WX), col = "grey65", lty = "dashed") +
    scale_shape_manual(name = X_name, values = c(21,22)) +
    scale_color_continuous_diverging("Blue-Red", name = X_name, p1 = .25, p2 = 1) +
    scale_fill_continuous_diverging("Blue-Red", name = X_name, p1 = .25, p2 = 1) +
    coord_fixed() +
    ylim(range(X, WX)) +
    xlim(range(X, WX)) +   
    themeLCD() +
    labs(x = X_name, y = paste0(X_name, " (lag)")) +
    theme(legend.position = "none")
  if (!missing(only_scatter) && isTRUE(only_scatter))
    return(mp_scatter)
  
  mp_map <- ggplot(mp) +
    geom_sf(data = geo, fill = NA) +
    geom_point(aes(
      x = lon, y = lat, fill = X, col = WX,
      shape = if_else(0 > sign(X * WX), "Pos", "Neg"),
      size = sqrt(X^2 + WX^2),
      stroke = sqrt(X^2 + WX^2))) +
    scale_shape_manual(name = X_name, values = c(21, 22)) +
    scale_color_continuous_diverging("Blue-Red", name = X_name, 
                                     p1 = .25, p2 = 1) +
    scale_fill_continuous_diverging("Blue-Red", name = X_name, 
                                    p1 = .25, p2 = 1) +
    coord_sf() +
    theme_map() +
    labs(x = "Longitude", y = "Latitude") +
    theme(legend.position = "none")
  if (!missing(only_map) && isTRUE(only_map))
    return(only_map)
  
  mp_scatter + mp_map
}


#' @author Lukas Dargel
themeLCD <- function(b_size = 15) {
  theme_bw(base_size = b_size) %+replace%
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "gray80", colour = "black"),
      strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
      axis.text.y = element_text(hjust = 0, color = "black", face = "bold"),
      axis.text.x = element_text(hjust = 0.5, color = "black", face = "bold")
    )
}

#' @author Lukas Dargel
theme_map <- function(b_size = 15, ...) {
  theme_bw(base_size = b_size) %+replace%
    theme(
      text = element_text(family = "Ubuntu Regular", color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_line(color = "#ebebe5", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "#FFFFFF", color = NA),
      panel.background = element_rect(fill = "#FFFFFF", color = NA),
      legend.background = element_rect(fill = "#FFFFFF", color = NA),
      panel.border = element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm"),
      ...
    )
}
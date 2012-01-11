# find the best position for a legend given a set of points
bestLegendPos <- function(x, y) {
  x1 <- min(x) + (max(x) - min(x))/3
  x2 <- min(x) + (max(x) - min(x))*2/3
  y1 <- min(y) + (max(y) - min(y))/3
  y2 <- min(y) + (max(y) - min(y))*2/3
  n.bottomleft <- sum(x < x1 & y < y1)
  n.bottom <- sum(x1 < x & x < x2 & y < y1)
  n.bottomright <- sum(x2 < x & y < y1)
  n.left <- sum(x < x1 & y1 < y & y < y2)
  n.center <- sum(x1 < x & x < x2 & y1 < y & y < y2)
  n.right <- sum(x2 < x & y1 < y & y < y2)
  n.topleft <- sum(x < x1 & y2 < y)
  n.top <- sum(x1 < x & x < x2 & y2 < y)
  n.topright <- sum(x2 < x & y2 < y)

  minpts <- min(n.bottomleft, n.bottom, n.bottomright,
                n.left, n.center, n.right,
                n.topleft, n.top, n.topright)

  if (minpts == n.topleft) return("topleft")
  if (minpts == n.topright) return("topright")
  if (minpts == n.bottomleft) return("bottomleft")
  if (minpts == n.bottomright) return("bottomright")
  if (minpts == n.top) return("top")
  if (minpts == n.bottom) return("bottom")
  if (minpts == n.left) return("left")
  if (minpts == n.right) return("right")
  if (minpts == n.center) return("center")
}

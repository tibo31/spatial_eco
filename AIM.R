# Author of these codes :
# thibault.laurent@tse-fr.eu
# 1. DAW function 
daw <- function(sources, targets, y, nature = "extensive", 
                st.inter = NULL, scaling = TRUE, onebyone = FALSE) {
  # computational conditions on the sources and the targets
  # - sources must be a SpatialPolygonsDataFrame including the Y variable
  # - targets may be either a SpatialPolygons or a SpatialPolygonsDataFrame
  stopifnot(class(sources) == "SpatialPolygonsDataFrame", 
            class(targets)%in%c("SpatialPolygonsDataFrame", "SpatialPolygons"), 
            all(y%in%names(sources)), nature%in%c("extensive","intensive"))
  # CRS should be the same ....
  if (!all(is.na(c(proj4string(sources), proj4string(targets))))) {
    if (any(is.na(c(proj4string(sources), proj4string(targets))))) {
      warning("Sources and Targets do not have the same CRS")
    } else {
      if (!all.equal(proj4string(sources),proj4string(targets)))
        warning("Sources and Targets do not have the same CRS")
    }
  }
  # if Targets are not a spatialpolygonsdataframe, we prepare it
  if (class(targets)=="SpatialPolygons") {
    targets <- SpatialPolygonsDataFrame(targets, data.frame(ID = sapply(slot(targets, "polygons"), slot, "ID"),
                                                            row.names = "ID"))
  }
  
  # Spatial units were y is not observed are removed
  if (length(y) == 1) {
    res.NA <- is.na(sources@data[, y])
  } else {
    res.NA <- apply(sources@data[, y], 1, function(x) any(is.na(x)))
    if(any(res.NA)) 
      print(paste("Warnings: In the sources, we removed spatial unit number",
                  which(res.NA), "because it contains at least one missing value among Y variables"))
  }
  
  sources <- sources[!res.NA,]
  
  # initialization
  targets@data[, paste(y, "daw", sep = "")] <- NA
  ns <- nrow(sources)
  nt <- nrow(targets)
  p <- length(y)
  
  areaInter <- matrix(0, nrow = ns, ncol = nt) # matrix containing the area the intersections
  
  # vector containing the area of the sources or the targets
  if (nature == "extensive")
    areaSource <- gArea(sources, byid = TRUE)
  
  # areaTarget has to be computed for nature=extensive and if scale=TRUE
  areaTarget <- gArea(targets, byid = TRUE)
  
  dimnames(areaInter)[[1]] <- row.names(sources)
  dimnames(areaInter)[[2]] <- row.names(targets)
  
  if (is.null(st.inter)) {# Polygons intersected between sources and targets
    i.s.t <- intersect.spdf(sources, targets, out = "spdf", onebyone = onebyone)
    # we check in the intersection where are the sources and the targets
    fips.i.s.t  <- strsplit(sapply(slot(i.s.t, "polygons"), slot, "ID"), "_", perl = TRUE) 
    fipsTarget2 <- unlist(lapply(fips.i.s.t, function(x) x[2])) 
    fipsSource2 <- unlist(lapply(fips.i.s.t, function(x) x[1])) 
    # we fill the matrix of area 
    areaInter[cbind(fipsSource2, fipsTarget2)] <- sapply(slot(i.s.t, "polygons"), slot, "area")
  } else {
    areaInter[cbind(st.inter$fipsSource, st.inter$fipsTarget)] <- st.inter$area
  }
  
  
  if(nature == "extensive") {
    # Estimation of Y on the intersections :
    # Compute Weight matrix
    W <- apply(areaInter, MARGIN = 2, function(x) x/areaSource)
    # Estimated values on the intersections for the variable hh.owner
    estimtargets.daw <- t(W) %*% as.matrix(sources@data[, y], ncol = p)
  } else {
    W <- t(apply(areaInter, MARGIN = 1, function(x) x/areaTarget))
    estimtargets.daw <- t(W) %*% as.matrix(sources@data[, y], ncol = p)
  }
  
  # if scale=TRUE, Yt is multiplied by: |Tt|/sum|A_st|
  if(scaling)
  {estimtargets.daw<-estimtargets.daw*areaTarget/apply(areaInter,2,sum)
  }
  
  # the results
  # If targets do not contains any sources, value will be NA
  targets@data[apply(areaInter,2,sum)>0,paste(y,"daw",sep="")]<-estimtargets.daw[apply(areaInter,2,sum)>0,]
  
  return(targets)
}

# 2. DAX function 
dax <- function(sources, targets, y, st.df, x, Z.y = NULL, Z.x = NULL, scaling = TRUE) {
  # y are the names of the variables to be treated (i.e. transformed from sources to targets).
  # They must be necessarly all extensive (Z.y=NULL) or intensive (Z.y non NULL).
  # In the case of intensivity, the correction on y will be done with the same variable Z.y 
  
  # st.df: data.frame which contains a colum for the sources id (fipsSource), the targets (fipsTarget) id 
  # and the a column called "x" (auxiliary information at the intersection level)
  
  # computational conditions on the sources and the targets
  # - sources must be a SpatialPolygonsDataFrame including the Y variable
  # - targets may be either a SpatialPolygons or a SpatialPolygonsDataFrame
  stopifnot(class(sources) == "SpatialPolygonsDataFrame",
            class(targets)%in%c("SpatialPolygonsDataFrame", "SpatialPolygons"),
            all(y%in%names(sources)))
  
  # verifications on the data.frame of the intersections
  stopifnot(c("fipsTarget","fipsSource")%in%colnames(st.df),
            all(st.df[, "fipsTarget"]%in%row.names(targets)),
            all(st.df[, "fipsSource"]%in%row.names(sources)),
            x%in%names(st.df))
  
  # verification on the projection
  if (!all(is.na(c(proj4string(sources), proj4string(targets))))) {
    if (any(is.na(c(proj4string(sources), proj4string(targets))))) { 
      warning("Sources and Targets do not have the same CRS")
    } else {
      if (!all.equal(proj4string(sources),proj4string(targets)))
        warning("Sources and Targets do not have the same CRS")
    }
  }
  
  # if targets is not a spatialpolygonsdataframe, we prepare it
  if (class(targets) == "SpatialPolygons") {
    targets <- SpatialPolygonsDataFrame(targets,
                                        data.frame(ID = sapply(slot(targets, "polygons"), slot, "ID"),
                                                   row.names = "ID"))
  }
  # we prepare the variables to be created
  y.label <- paste(y, "dax", sep = "")
  targets@data[, y.label] <- NA
  
  # If x is not available at the source level, we can build it with the
  # auxiliary information of the intersection
  if (!x%in%names(sources) & is.null(Z.x)) {
    sources@data[, x] <- NA
    sources@data[, x] <- sapply(row.names(sources),
                                function(z) sum(st.df[st.df$fipsSource == z, x]) )
  }
  
  # if Z.y=NULL, y is supposed to be extensive. Otherwise, user give the name of the variable Z 
  # which permits to transform y as an extensive variable. 
  # condition on Z.y : if non NULL, it should exist at least on intersection. 
  # If Z.y does exist on target and sources, it will be used, otherwise, it will be built with intersection
  
  if (!is.null(Z.y)) {
    stopifnot(length(Z.y) == 1, (Z.y%in%names(targets) & Z.y%in%names(sources))|Z.y%in%names(st.df))
    if (!Z.y%in%names(sources)) {
      Z.y.source <- by(st.df[,Z.y], st.df$fipsSource, sum)
      sources@data[, Z.y] <- Z.y.source[row.names(sources)]
    }
    sources@data[, y] <- sources@data[, y]*sources@data[, Z.y]
    if (!Z.y%in%names(targets)) {
      Z.y.target <- by(st.df[,Z.y], st.df$fipsTarget, sum)
      Z.y.target <- Z.y.target[row.names(targets)]
    } else {
      Z.y.target<-targets@data[,Z.y]
    } 
  } 
  # if Z.x=NULL, x is supposed to be extensive. Otherwise, user give the name of the variable Z 
  # which permits to transform x as an extensive variable. 
  # condition on Z.x : if non NULL, it should exist both on source and intersection
  if (!is.null(Z.x)) {
    stopifnot(length(Z.x) == 1, Z.x%in%names(st.df))
    if (!Z.x%in%names(sources)) {
      Z.x.source <- by(st.df[,Z.x], st.df$fipsSource, sum)
      sources@data[, Z.x] <- Z.x.source[row.names(sources)]
    }
    
    sources@data[, x] <- sources@data[, x]*sources@data[, Z.x]
    st.df[, x] <- st.df[, x]*st.df[, Z.x] 
  }   
  
  # initialization
  fipsTarget <- row.names(targets)
  fipsSource <- row.names(sources)
  ns <- nrow(sources)
  nt <- nrow(targets)
  p <- length(y)
  
  # matrix of the intersections
  interX <- matrix(0, nrow = ns, ncol = nt) # matrix containing the X on the intersections
  dimnames(interX)[[1]] <- row.names(sources)
  dimnames(interX)[[2]] <- row.names(targets)
  interX[cbind(as.character(st.df[, "fipsSource"]), as.character(st.df[, "fipsTarget"]))] <- st.df[, x]
  interX <- interX/matrix(sources@data[, x], ns, nt)
  interX[is.na(interX)] <- 0
  
  # we compute the Y on the targets  
  estimtargets.dax <- t(interX)%*%as.matrix(sources@data[, y], ncol = p)
  
  if(scaling) {
    if(!x%in%names(targets)) {
      x.targets <- sapply(row.names(targets), function(z) sum(st.df[st.df$fipsTarget == z, x]))
    } else {
      x.targets <- targets@data[, x]
    }
    
    sum.X <- tapply(st.df[, x], st.df[, "fipsTarget"], sum, na.rm = TRUE)
    estimtargets.dax <- as.matrix(estimtargets.dax)*as.vector(as.vector(x.targets)/sum.X[row.names(targets)]) 
  }
  
  if (!is.null(Z.y)) {
    estimtargets.dax <- estimtargets.dax/Z.y.target
  }
  
  # res.inter<-aggregate(st.df[,"area"],list(st.df[,"fipsTarget"]),sum)
  # targets@data[,y.label]<-ifelse(row.names(targets)%in%as.character(res.inter[,"Group.1"]),estimtargets.dax,NA)
  targets@data[, y.label] <- estimtargets.dax
  
  return(targets)
}

# 3. intersection method 
intersect.spdf<-function(sources, targets, out = "df", onebyone = FALSE) {
  # this function creates either a data.frame with the variables 
  # fipsSource and fipsTarget and the variable area
  # either a SpatialPolygonsDataFrame if out="spdf"
  
  stopifnot(out%in%c("df","spdf"))
  
  fipsTarget <- row.names(targets)
  fipsSource <- row.names(sources)
  
  # we compute the area of the intersection  
  if (onebyone) {
    ns <- length(sources)
    nt <- length(targets)
    res <- sources[1, ]
    for (i in 1:ns) {
      for (j in 1:nt) {
        res.b <- gIntersection(sources[i, ], targets[j, ], byid = TRUE)
        if (class(res.b) == "SpatialCollections") {
          res <- spRbind(res, res.b@polyobj)
        } else {
          if (class(res.b) == "SpatialPolygons") {
            res <- spRbind(res, res.b)
          }
        }
      }
    }
    res <- res[-1, ]
  } else {
    res <- gIntersection(sources, targets, byid = TRUE, drop_lower_td = TRUE)
  }
  
  fipsSource <- strsplit(sapply(slot(res, "polygons"), slot, "ID"),
                         paste(paste("", row.names(targets)), collapse = "|"), perl = TRUE) 
  fipsSource <- unlist(lapply(fipsSource, function(x) x[1])) 
  fipsTarget <- strsplit(sapply(slot(res, "polygons"), slot, "ID"),
                         paste(paste(row.names(sources), ""), collapse = "|"), perl = TRUE) 
  fipsTarget <- unlist(lapply(fipsTarget, function(x) x[length(x)]))  
  
  # fipsSource=gsub(paste(paste("",fipsTarget),collapse="|"),"",sapply(slot(res, "polygons"), slot, "ID"))
  # fipsTarget=gsub(paste(paste(fipsSource,""),collapse="|"),"",sapply(slot(res, "polygons"), slot, "ID"))
  
  res.df <- data.frame(fipsSource = fipsSource, fipsTarget = fipsTarget,
                       area = gArea(res, byid = TRUE),
                       row.names = paste(fipsSource, fipsTarget, sep = "_"))
  
  if (out == "df") {
    return(res.df)
  } else {
    row.names(res) <- row.names(res.df)
    return(SpatialPolygonsDataFrame(res, res.df))
  }
}
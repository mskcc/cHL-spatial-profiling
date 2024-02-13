#' Get a bounding box for a single FOV
#'
#' Given a table of X/Y coordinates of all cells, get the list of 
#' min/max coordinates that form a box enclosing all cells in a FOV. Optionally
#' crob the FOV by shrinking the bounding box by a fixed amount on all
#' sides.
#' 
#' @param coords  tibble containing columns X and Y representing cell coordinates
#' @param crop    default = 0; a fixed value to cut from all sides of FOV. WARNING:
#'                this value MUST be in the same unit as X and Y values. 
#'
#' @return list containing values X0, X1, Y0, Y1
getBoundingBox <- function(coords, crop = 0){
    list(X0 = min(coords$X) + crop, 
         X1 = max(coords$X) - crop,
         Y0 = min(coords$Y) + crop, 
         Y1 = max(coords$Y) - crop)
}


#' Calculate total area of bounding box
#'
#' Simple area calculation on FOV bounding box
#' 
#' @param bb  list representing FOV bounding box, containing
#'            values for keys X0, X1, Y0, Y1
#'
#' @return area in same unit as bounding box values
areaBB<-function(bb){
    (bb$X1-bb$X0)*(bb$Y1-bb$Y0)
}


#' Determine if a single point falls inside a polygon
#' 
#' Determine if a single point falls inside a polygon
#'
#' @param pts    table containing columns X and Y representing point coordinates
#' @param poly   table containing columns X and Y representing a single polygon
pointsInPolygon <- function(pts,poly){
    point.in.polygon(pts$X,pts$Y,poly$X,poly$Y)==1
}


#' Get the list of points that fall inside a polygon
#' 
#' Given a table of X/Y coordinates and a list of polygon boundaries,
#' return the coordinates table filtered for points that are inside the polygons
#'
#' @param pts    table containing columns X and Y representing point coordinates
#' @param bList  a list containing one or more tables of X/Y coordinates
#'               each representing a boundary of a single polygon
#' 
#' @return filtered table of points inside any of the given polygons
pointsInsidePolygonList <- function(pts,bList) {

    inside=bList %>%
        map(function(bn){point.in.polygon(pts$X,pts$Y,bn$X,bn$Y)>0}) %>%
        as.data.frame %>%
        apply(.,1,any)

    inside

}


#' Generate data frame of random points (is it random???)
#' 
#' Generate a data frame of random points that fall within the plot area
#' 
#' @param   nGrid   number of points to generate
#' @param   bbG     a list containing X0,X1,Y0,Y1 representing the boundary
#'                  of the plot area
#' @return data frame of X and Y values of the random points
generateSobolGrid <- function(nGrid,bbG) {
    gg=sobol(nGrid,2,scrambling=2)
    data.frame(
        X=gg[,1]*(bbG$X1-bbG$X0)+bbG$X0,
        Y=gg[,2]*(bbG$Y1-bbG$Y0)+bbG$Y0
    )
}


#' Calculate total tissue area in a single FOV
#'
#' Calculate total tissue area in a single FOV
#' 
#' @param fovDat           tibble of all halo data for a single FOV
#' @param boundaries       list of one or more tables, each containing X/Y coordinates
#'                         representing a single exclusion boundary
#' @param maxG             ???? default=5
#' @param unit             ['px'|'um'|'mm'] value will be returned in these square units; default = 'um' 
#'
#' @return total FOV area in unit of choice 
#' @export
calculateAreaTotalFOV <- function(fovDat, boundaries, maxG=5, unit = "um"){

    bb   <- getBoundingBox(fovDat)
    gg   <- generateSobolGrid(10^maxG, bb)

    if(!is.null(boundaries) && length(boundaries) > 0){
        gg <- gg[!pointsInsidePolygonList(gg, boundaries),]
    }

    ## area is in square pixels
    px2 <- areaBB(bb) * nrow(gg)/10^maxG

    areas <- list(px = px2,
                  um = px2 * unique(fovDat$PixelScaling)^2,
                  mm = px2 * unique(fovDat$PixelScaling)^2 * (0.001^2))

    areas[[unit]]
}

#' Calculate area for all FOVs in a single sample
#' 
#' Get a table of CellDive_ID, FOV_number, Area and Area_unit
#' including total area for each FOV in one sample
#' 
#' @param sHaloFile    path to Halo object file for a single sample
#' @param unit         ['px'|'um'|'mm'] areas will be returned in these square units; default = 'mm'
#' @param maxG        
#' @param trimFOVs     default = FALSE; when true, data will be filtered for cells with value in column
#'                     'Exclude.Boundary' set to FALSE; this effectively calculates the area
#'                     of the FOV AFTER cropping to eliminate edge artifacts
#' @param outFile      default = NULL; RDA file that already contains or will contain recalculated areas
#' @param forceRecalc  default = FALSE; recalculate areas even if outFile exists and contains data 
#'
#' @return table with columns CellDive_ID, FOV_number, Area and Area_unit
getAllSampleFOVAreas <- function(sHaloFile, unit = "mm", maxG = 5,
                                  trimFOVs = FALSE, outFile = NULL, forceRecalc = FALSE){

    if(fileDone(outFile) && !forceRecalc){
        log_info(paste0("Loading pre-computed Sample FOV areas from file ",outFile))
        readRDS(outFile)
    }

    log_info("No pre-computed Sample FOV areas file found and/or forced rerun turned ON. Calculating areas now...")

    dat <- readRDS(sHaloFile)
    if(trimFOVs){
        log_warn("TRIMMING FOVs using column 'Exclude.Boundary' before calculating areas.")
        dat$geom.data <- dat$geom.data %>% filter(!Exclude.Boundary)
    }

    res <- tibble()
    for(fov in unique(dat$geom.dat$SPOT)){
        log_debug("Excluding ", length(dat$exclusionRegions[[paste0("spot_", fov)]]), " regions from FOV number ", fov, ".")
        area <- calculateAreaTotalFOV(dat$geom.dat %>% filter(SPOT == fov),
                                      dat$exclusionRegions[[paste0("spot_", fov)]],
                                      maxG = maxG,
                                      unit = unit)
        res <- res %>%
               bind_rows(tibble(CellDive_ID = unique(dat$geom.dat$Sample),
                                FOV_number = fov,
                                Area = area,
                                Area_unit = paste0(unit, "^2")))
    }

    if(!is.null(outFile)){
        saveRDS(res, outFile)
    }

    res
}



#' Get densities of cell populations
#' 
#' Given a tibble of cell population counts and one of areas per calculation unit, calculate
#' density for each population 
#' 
#' @param counts      tibble of counts for all cell populations for which densities are 
#'                    to be calculated
#' @param areas       tibble of areas including a single value for each calculation unit
#' @param calcUnit    data column name containing values for which a single density value
#'                    should be calculated (default: FOV_ID)
#' @param outfile     optional; path to RDA file where results should be saved
#'
#' @return tibble of cell population densities
#' @export
getPopulationDensities <- function(counts, areas, calcUnit = "FOV_ID", outFile = NULL, forceRecalc = FALSE){

    if(fileDone(outFile)){
        log_info(paste0("Loading pre-computed densities from file ",outFile))
        return(readRDS(outFile))
    }

    dens <- tibble()

    log_info("No pre-computed densities file found. Generating now...")
    areasPer <- areas %>%
                group_by_at(calcUnit) %>%
                summarize(TotalArea = sum(Area, na.rm = T))

    # remove NA only if summing multiple FOV values (i.e., calcUnit == "CellDive_ID")
    unitCount <- counts %>% 
                 filter(is.na(Count)) %>% 
                 group_by(Population) %>%
                 summarize(unitCount = n()) %>%
                 pull(unitCount) %>%
                 unique 
    rmNA <- ifelse(unitCount == 1 || calcUnit == "FOV_ID", FALSE, TRUE)

    dens <- counts %>%
            group_by_at(c(calcUnit,"Population")) %>%
            mutate(TotalCounts = sum(Count, na.rm = rmNA)) %>%
            left_join(areasPer, by=c(calcUnit)) %>%
            mutate(Density = TotalCounts/TotalArea)

    if(!is.null(outFile)){
        saveRDS(dens, file = outFile)
    }

    dens
}


#' Remove all data points that fall outside of the bounding box
#' 
#' Filter out rows in a data frame that fall outside of a preset 
#' FOV bounding box
#' 
#' @param df  data frame to be filtered, including columns X and Y
#' @param bb  list of four points named X0, X1, Y0, Y1 with X and Y coordinates
#'            of a preset FOV bounding box
#' 
#' @return filtered data frame
trimDFToBB<-function(df,bb) {

    df %>%
    filter(X>=bb$X0, X<=bb$X1, Y>=bb$Y0, Y<=bb$Y1)

}


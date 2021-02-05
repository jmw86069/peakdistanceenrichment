
#' Count peak-feature overlap with increasing distance
#'
#' Count peak-feature overlap with increasing distance
#'
#' This function counts features overlaps across peaks,
#' while progressively expanding the peak width.
#' It also tabulates counts by feature category, using
#' the feature annotations in `category_colname`.
#'
#' @param peaksL list named by experiment, containing a
#'    list of GRanges, where each set of GRanges represents
#'    a set of peaks.
#' @param featuresL list named by experiment, containing
#'    GRanges representing features, with `category_colname`
#'    containing category values.
#' @param randomsL list named by experiment, containing
#'    GRanges from which random peaks should be drawn, using
#'    the same number of peaks as each corresponding entry
#'    in `peaksL`.
#' @param do_random logical indicating whether to choose
#'    random peaks from `randomsL` using the same number of
#'    peaks as the corresponding entry in `peaksL`.
#' @param category_colname character string indicating the
#'    colname present in the GRanges for each `featuresL`
#'    element. If the values are factors, the factor
#'    levels should be preserved in order. When the
#'    column `category_colname` does not exist, or for
#'    any NA entries, `default_category` is used as a
#'    replacement.
#' @param default_category vector length=1 used when
#'    `category_colname` is not present in
#'    `colnames(values(featuresL[[i]]))`.
#' @param expand_kb numeric vector of distances in kilobases (kb)
#'    to use when expanding the peaks in `peaksL`. When
#'    `expand_kb` is `NULL` then the default values are used,
#'    within the range defined by `expand_range_kb`.
#' @param expand_range_kb numeric range of distances to allow
#'    in kilobases (kb) for `expand_kb`.
#' @param expand_factor numeric value multiplied by the expand width
#'    during `GenomicRanges::resize()`, since the resize happens
#'    at peak center and defines a fixed width. For example,
#'    use `expand_factor=2` when this distance should represent
#'    distance from the peak in either direction, for which the
#'    peaks are expanded from center by this distance.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' hg19.sizes.txt
#'
#' @export
peak_feature_ranged_counter <- function
(peaksL=NULL,
 featuresL=NULL,
 randomsL=NULL,
 category_colname="hit",
 default_category="hit",
 expand_kb=NULL,
 expand_range_kb=c(5, 5000),
 expand_factor=2,
 do_random=FALSE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take a set of peaks (ATACseq peaks intended)
   ## a set of features (gene TSSes)
   ## a vector of expand distances, and
   ## expand the peaks by this distance on each side,
   ## each time count the feature overlaps, subdivide
   ## overlap features by category, using category_colname
   ##
   ## Define Experiments
   if ("GRanges" %in% class(peaksL)) {
      peaksL <- list(Experiment=peaksL);
      if ("GRanges" %in% class(featuresL)) {
         featuresL <- list(Experiment=featuresL);
      }
      if (length(randomsL) > 0 && "GRanges" %in% class(randomsL)) {
         randomsL <- list(Experiment=randomsL);
      }
   }
   Experiments <- jamba::nameVectorN(peaksL);
   if (any(!Experiments %in% names(featuresL))) {
      stop(paste0("All names(peaksL) must be present in names(featuresL), ",
         "where names(peaksL)=",
         paste0("'", names(peaksL), "'", collapse=", ")));
   }
   if (length(randomsL) == 0) {
      do_random <- FALSE;
   }
   if (do_random &&
         length(randomsL) > 0 &&
         any(!Experiments %in% names(randomsL))) {
      stop(paste0("All names(peaksL) must be present in names(randomsL), ",
         "where names(peaksL)=",
         paste0("'", names(peaksL), "'", collapse=", ")));
   }

   ## Default expand_kb
   if (length(expand_kb) == 0) {
      expand_kb <- get_expanded_ranges(x=5,
         expand_range_kb);
   }
   expand_range_kb <- range(expand_range_kb, na.rm=TRUE);
   expand_kb <- expand_kb[expand_kb >= min(expand_range_kb) &
         expand_kb <= max(expand_range_kb)]
   if (length(expand_kb) == 0) {
      stop("expand_kb length is zero, after applying expand_range_kb.");
   }
   ## Check experiment names
   if (!all(names(peaksL) %in% names(featuresL))) {
      stop("All names(peaksL) must be present in names(featuresL).")
   }

   ## Detect category levels across all input data
   cat_levelsL <- lapply(jamba::nameVectorN(featuresL), function(Experiment){
      iGR <- featuresL[[Experiment]];
      if ("list" %in% class(iGR)) {
         iGR <- GenomicRanges::GRangesList(iGR)@unlistData;
      }
      if (!category_colname %in% colnames(values(iGR))) {
         if (verbose) {
            jamba::printDebug("peak_feature_ranged_counter(): ",
               "assigning default_category:",
               default_category,
               " to missing category_colname:",
               category_colname);
         }
         GenomicRanges::values(iGR)[[category_colname]] <- default_category;
      }
      ## Get unique values, replacing NA values with default_category
      cat_levels <- jamba::rmNA(
         naValue=default_category,
         unique(GenomicRanges::values(iGR)[[category_colname]]));
   });
   cat_levels <- unique(unlist(cat_levelsL));
   if ("factor" %in% class(cat_levels)) {
      cat_levels <- levels(cat_levels);
   } else {
      cat_levels <- as.character(jamba::mixedSort(cat_levels, keepNegative=TRUE));
   }
   if (verbose) {
      jamba::printDebug("peak_feature_ranged_counter(): ",
         "cat_levels:", cat_levels);
   }

   ##################################
   ## Iterate Experiments
   ## - each Experiment contains one or more sets of peaks
   ## - iterate each set of peaks
   ##   -
   atacRangesL <- lapply(Experiments, function(Experiment){
      if (verbose) {
         printDebug("peak_feature_enrich(): ",
            "Processing Experiment:", Experiment);
      }
      bedL <- peaksL[[Experiment]];
      if ("list" %in% class(featuresL[[Experiment]])) {
         if (!all(names(bedL) %in% names(featuresL[[Experiment]]))) {
            stop(paste0(
               "When featuresL is a list of lists, all child names must match peaksL.",
               "This failed for Experiment '", Experiment, "'."));
         }
      }

      ##################################
      ## Iterate each PeakSet
      atacRangesDF <- jamba::rbindList(lapply(jamba::nameVectorN(bedL), function(PeakSet){
         iGR <- bedL[[PeakSet]];
         ## Optionally create random peaks
         if (do_random) {
            if ("list" %in% class(featuresL[[Experiment]])) {
               randomGR <- randomsL[[Experiment]][[PeakSet]];
            } else {
               randomGR <- randomsL[[Experiment]];
            }
            iGR <- randomGR[sample(seq_along(randomGR),
               size=length(iGR))];
            if (verbose) {
               jamba::printDebug("peak_feature_enrich(): ",
                  "Selected ",
                  length(iGR),
                  " random peaks from ",
                  length(randomGR),
                  " provided peaks.");
            }
         }
         if ("list" %in% class(featuresL[[Experiment]])) {
            tssGR <- featuresL[[Experiment]][[PeakSet]];
         } else {
            tssGR <- featuresL[[Experiment]];
         }
         if (verbose) {
            jamba::printDebug("peak_feature_enrich(): ",
               "      Processing PeakSet:", PeakSet,
               ", length(PeakSet):", length(iGR),
               ", length(featuresGR):", length(tssGR));
         }

         ##################################
         ## Iterate each expand_kb
         atacExDF <- jamba::rbindList(lapply(expand_kb, function(iExpand){
            expand_width_use <- iExpand*1000*expand_factor;
            ## expand the regions
            ## reduce them so overlapping regions are combined
            iGRex <- GenomicRanges::reduce(
               GenomicRanges::resize(iGR,
                  width=expand_width_use,
                  fix="center"));
            ol <- GenomicRanges::findOverlaps(iGRex, tssGR);
            ## get the tssGR sites that overlap the expanded regions
            sol <- unique(S4Vectors::subjectHits(ol));
            unsol <- setdiff(seq_along(tssGR), sol);

            ## Tabulate overlaps by cat_levels
            if (length(sol) == 0) {
               hitV <- rep(0, length(cat_levels));
               nonHitV <- table(GenomicRanges::values(tssGR)[[category_colname]])[cat_levels];
            } else if (length(sol) == length(tssGR)) {
               hitV <- table(GenomicRanges::values(tssGR)[[category_colname]])[cat_levels];
               nonHitV <- rep(0, length(cat_levels));
            } else {
               hitV <- table(GenomicRanges::values(tssGR[sol])[[category_colname]])[cat_levels];
               nonHitV <- table(GenomicRanges::values(tssGR[unsol])[[category_colname]])[cat_levels];
            }
            hitV <- jamba::rmNA(naValue=0, hitV);
            nonHitV <- jamba::rmNA(naValue=0, nonHitV);
            names(hitV) <- cat_levels;
            names(nonHitV) <- paste0(cat_levels, "_non");
            iDF <- data.frame(check.names=FALSE,
               Experiment=Experiment,
               PeakSet=PeakSet,
               Width=sum(width(iGRex)),
               Range=iExpand,
               as.list(hitV),
               as.list(nonHitV));
            iDF;
         }));
         atacExDF;
      }));
   });
   atacRangesLDF <- jamba::rbindList(atacRangesL);
   atacRangesDF <- tidyr::gather(atacRangesLDF,
      key="Category",
      value="Count",
      `-1`:`1_non`);
   atacRangesDF$Subclass <- ifelse(grepl("_non", atacRangesDF[["Category"]]),
      "Non-Overlap",
      "Overlap");
   return(atacRangesDF);
}

#' Test hypergeometric enrichment of range count categories
#'
#' Test hypergeometric enrichment of range count categories
#'
#' This function takes output from `peak_feature_ranged_counter()`
#' and runs hypergeometric tests for each corresponding category
#' value.
#'
#' Specifically, the argument `test_cats` defines the category
#' values to test in each enrichment test. For example, the
#' category values may include `c("-1", "0", "1")`, and one may
#' want to test `c("1")`. This test will count the number of
#' peaks with category `"1"` that overlap features,
#' compared to the number of peaks with other category
#' values, in this case `c("-1", "0")`, that overlap features.
#'
#' It is possible to combine category values, for example one
#' may test `c("-1", "1")` and in this scenario, any peak
#' with the category `"1"` or `"-1"` will be tested against
#' all other peaks, in this case category `"0"`.
#'
#'
#' @param range_counts `data.frame` output from
#'    `peak_feature_ranged_counter()`.
#' @param test_cats `list` of `character` vectors, where each
#'    vector contains one or more category values to test for
#'    enrichment against all other category values. When
#'    `test_cats` is `NULL` then each category value is tested
#'    individually against the other category values.
#' @param return_type `character` string indicating the type
#'    of output to return: `"phyper"` will return the results
#'    of the hypergeometric test using `phyper()`; and `"wide"`
#'    will return a wide `data.frame` that contains the exact
#'    numbers used for this test.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
peak_feature_ranged_enrichment <- function
(range_counts,
 test_cats=NULL,
 return_type=c("phyper", "wide"),
 verbose=FALSE,
 ...)
{
   ## Purpose is to take counts from peak_feature_ranged_counter()
   ## and run hypergeometric enrichment
   return_type <- match.arg(return_type);
   cat_levels <- unique(range_counts$Category);
   all_cats <- jamba::unvigrep("_non$", cat_levels);
   if (length(test_cats) == 0) {
      test_cats <- as.list(jamba::unvigrep("_non$", cat_levels));
   } else if (!"list" %in% class(test_cats)) {
      test_cats <- as.list(test_cats);
   }
   if (length(names(test_cats)) == 0) {
      names(test_cats) <- jamba::cPaste(test_cats);
   }
   if (verbose) {
      printDebug("peak_feature_ranged_enrichment(): ",
         "test_cats:");
      print(test_cats);
   }
   ## make wider form
   range_counts$Category <- factor(range_counts$Category,
      levels=cat_levels);
   #head(dplyr::select(range_counts, -Subclass))
   range_counts_wide <- tidyr::spread(
      dplyr::select(range_counts, -Subclass),
      key="Category",
      value="Count");
   if ("wide" %in% return_type) {
      return(range_counts_wide);
   }

   paramColnames <- intersect(
      c("Experiment", "PeakSet", "Range",
         "Iteration", "Test"),
      colnames(range_counts_wide));
   paramsDF <- jamba::rbindList(lapply(jamba::nameVectorN(test_cats), function(test_cat_name){
      test_cat <- test_cats[[test_cat_name]];
      drawn1 <- test_cats[[test_cat_name]];
      undrawn1 <- paste0(drawn1, "_non");
      drawn0 <- setdiff(all_cats, drawn1);
      undrawn0 <- paste0(drawn0, "_non");
      if (verbose) {
         printDebug("peak_feature_ranged_enrichment(): ",
            "drawn1:",drawn1,
            ", drawn0:", drawn0,
            ", undrawn1:", undrawn1,
            ", undrawn0:", undrawn0);
      }
      q <- rowSums(range_counts_wide[,c(drawn1),drop=FALSE]);
      m <- rowSums(range_counts_wide[,c(drawn1, undrawn1),drop=FALSE]);
      n <- rowSums(range_counts_wide[,c(drawn0, undrawn0),drop=FALSE]);
      k <- rowSums(range_counts_wide[,all_cats]);
      data.frame(q_drawn1=q,
         m_drawn1_undrawn1=m,
         n_drawn0_undrawn0=n,
         k_drawn0_drawn1=k,
         test_cat=rep(test_cat_name, length.out=length(q)),
         range_counts_wide[,paramColnames,drop=FALSE]);
   }));
   paramsDF$test_cat <- factor(paramsDF$test_cat,
      levels=names(test_cats));
   paramsDF$phyper.p <- phyper(q=paramsDF$q_drawn1,
      m=paramsDF$m_drawn1_undrawn1,
      n=paramsDF$n_drawn0_undrawn0,
      k=paramsDF$k_drawn0_drawn1,
      lower.tail=FALSE);
   paramsDFcolnames <- unique(c(
      intersect(colnames(range_counts),colnames(paramsDF)),
      colnames(paramsDF)));
   return(paramsDF[,paramsDFcolnames,drop=FALSE]);
   #return(paramsDF);
}

#' Get vector of expanded ranges for peak widths
#'
#' Get vector of expanded ranges for peak widths
#'
#' This function takes a range of distances in kilobases,
#' and a default initial step size,
#' and returns a reasonable set of human-readable distances
#' that spans this range in using roughly log-scale steps.
#'
#' @param x `numeric` value representing the step size.
#' @param expand_range_kb `numeric` range of distances to allow
#'    in kilobases (kb) for `expand_kb`.
#' @param ... additional arguments are ignored.
#'
#' @return `numeric` vector of distances in kilobases.
#'
#' @examples
#' opar <- par("mfrow"=c(2,2))
#' for (x in c(2, 5, 10, 20)) {
#'    y <- get_expanded_ranges(x);
#'    plot(y,
#'       ylim=c(0, 5300),
#'       pch=20,
#'       xlab="step",
#'       main=paste("step size", x))
#' }
#' par(opar)
#'
#' @export
get_expanded_ranges <- function
(x=5,
 expand_range_kb=c(5, 5000),
 ...)
{
   #
   step_from <- c(5,10,20,50,100,500,1000,2000,5000,10000,20000);
   step_to <- c(10,20,50,100,500,1000,2000,5000,10000,20000,50000);
   step_by <- x * c(1,1,1,2,5,10,20,50,100,200,500);
   data.frame(step_from, step_to, step_by)
   expand_kbL <- lapply(seq_along(step_from), function(i){
      i_from <- step_from[i];
      i_to <- step_to[i];
      i_by <- step_by[i];
      seq(from=i_from, to=i_to, by=i_by);
   });
   expand_df <- data.frame(size=unlist(expand_kbL),
      step=rep(step_by,
         lengths(expand_kbL)));
   expand_df <- expand_df[match(unique(expand_df$size), expand_df$size),,drop=FALSE];
   expand_df$diff <- c(Inf, diff(expand_df[,1]));
   ## Remove steps closer together than x
   expand_df <- subset(expand_df, diff >= step);
   expand_kb <- expand_df$size;

   # restrict range of values
   expand_range_kb <- range(expand_range_kb, na.rm=TRUE);
   expand_kb <- expand_kb[expand_kb >= min(expand_range_kb) &
         expand_kb <= max(expand_range_kb)];
   return(expand_kb);
}

#' Plot ranged enrichment hypergeometric results
#'
#' Plot ranged enrichment hypergeometric results
#'
#' This function takes `"phyper"` output from
#' `peak_feature_ranged_enrichment()` and produces a plot
#' using `ggplot()`.
#'
#' @param paramsDF `data.frame` produced by
#'    `peak_feature_ranged_enrichment()` with
#'    `return_type="phyper"`.
#' @param p_cutoff `numeric` threshold used to draw a horizontal
#'    line indicating a threshold for statistical significance.
#' @param p_color `character` with R color, used to color the
#'    horizontal line defined by `p_cutoff`.
#' @param draw_points `logical` indicating whether each point
#'    should be drawn, or when `draw_points=FALSE` only the
#'    line is drawn.
#' @param lwd `numeric` indicating the line width.
#' @param pt_cex `numeric` indicating point size, used when
#'    `draw_points=TRUE`.
#' @param color_sub `character` vector of R colors, whose names
#'    correspond to category values tested for enrichment.
#' @param fill_sub `character` vector, same as `color_sub` but
#'    used when a fill color is required.
#' @param facet_formula `formula` sufficient to inform
#'    `facet_grid()` used to create separate ggplot panels.
#' @param shape_sub `integer` vector of R point shapes, whose
#'    names are peak names tested for enrichment.
#' @param alpha_sub `numeric` vector of alpha transparency values
#'    ranging between 0 and 1, whose names are `Test` and
#'    `Random`. These values are only used when random permutation
#'    tests have been performed, and the entries that resulted from
#'    random permutations are stored in column `"Test"` and have
#'    value `"Random"`.
#' @param y_include `numeric` vector used to extend the y-axis
#'    range to include at least these values. This effect is useful
#'    in requiring the y-axis includes y=0 for example, and y=10
#'    so that the y-axis range spans at least a reasonable viewing
#'    range for interpreting the data. The y-axis units are
#'    `-log10(pvalue)` values.
#' @param x_include `numeric` vector used to extend the x-axis to
#'    ensure the range includes this value.
#' @param test_cat_title `character` string used as the title for the
#'    category legend.
#' @param PeakSet_title `character` string used as a legend for point
#'    shape, only relevant when `draw_points=TRUE`.
#' @param main `character` string used as the title, called
#'    with `ggtitle()`.
#' @param xlab,ylab `character` string used for the x-axis and y-axis
#'    label, respectively.
#' @param scales `character` string passed to `facet_grid()` which
#'    defines the y-axis and x-axis range across multiple facet
#'    panels: `"free_y"` will allow each panel to have independent
#'    y-axis range; `"fixed"` will use one y-axis range across
#'    all facet panels.
#' @param ... additional arguments are ignored.
#'
#' @return `ggplot` object, which can be viewed by calling
#'    `print()`, or further modified.
#'
#' @export
plot_ranged_enrichment <- function
(paramsDF,
 p_cutoff=0.01,
 p_color="grey40",
 draw_points=TRUE,
 lwd=1,
 pt_cex=3,
 color_sub=NULL,
 fill_sub=color_sub,
 facet_formula=Experiment~PeakSet,
 shape_sub=c(`higher`=24, `lower`=25),
 alpha_sub=c(`Test`=1, `Random`=0.075),
 y_include=c(0,10),
 x_include=c(5),
 test_cat_title="Category",
 PeakSet_title="PeakSet",
 main="Enrichment by distance from peaks",
 xlab="Distance from peak center (kb)",
 ylab="-log10(hypergeometric P-value)",
 scales="free_y",
 ...)
{
   ## Purpose is to create a ggplot plot object
   shape_sub <- shape_sub[as.character(unique(paramsDF$PeakSet))];
   if (length(shape_sub) < length(unique(paramsDF$PeakSet))) {
      shape_sub <- NULL;
   }

   ## Handle random data
   if (!c("Iteration") %in% colnames(paramsDF)) {
      paramsDF[,"Iteration"] <- 0;
   }

   ## Define colors as needed
   if (length(color_sub) == 0) {
      color_sub <- colorjam::group2colors(unique(paramsDF$test_cat));
   }
   if (length(fill_sub) == 0) {
      fill_sub <- color_sub;
   }

   if ("Test" %in% colnames(paramsDF)) {
      if (!"factor" %in% class(paramsDF[,"Test"])) {
         paramsDF[,"Test"] <- factor(paramsDF[,"Test"],
            levels=jamba::provigrep(c("test", "random"),
               unique(paramsDF[,"Test"])));
      }
      gg_enrich <- ggplot2::ggplot(paramsDF,
         ggplot2::aes(x=Range,
            y=-log10(phyper.p),
            group=paste0(PeakSet, "_", test_cat, "_", Iteration),
            shape=PeakSet,
            alpha=Test,
            color=test_cat,
            fill=test_cat));
   } else {
      gg_enrich <- ggplot2::ggplot(paramsDF,
         ggplot2::aes(x=Range,
            y=-log10(phyper.p),
            group=paste0(PeakSet, "_", test_cat, "_", Iteration),
            shape=PeakSet,
            color=test_cat,
            fill=test_cat));
   }

   if (length(p_cutoff) > 0) {
      gg_enrich <- gg_enrich +
         ggplot2::geom_hline(yintercept=-log10(p_cutoff),
            color=p_color,
            lty="dashed",
            size=1);
   }
   if (any(!is.na(c(y_include, x_include)))) {
      gg_enrich <- gg_enrich +
         ggplot2::expand_limits(y=y_include,
            x=x_include);
   }
   #   geom_ribbon(aes(ymin=0, ymax=-log10(phyper.p)), alpha=0.4) +
   if (draw_points) {
      gg_enrich <- gg_enrich +
         ggplot2::geom_point(size=pt_cex) +
         ggplot2::scale_shape_manual(values=c(`higher`=24, `lower`=25),
            name=PeakSet_title);
   }
   gg_enrich <- gg_enrich +
      ggplot2::geom_line(size=lwd) +
      colorjam::theme_jam() +
      ggplot2::scale_color_manual(name=test_cat_title,
         values=color_sub) +
      ggplot2::scale_alpha_manual(values=alpha_sub);
   if (length(fill_sub) == 0) {
      fill_sub <- sapply(as.character(unique(paramsDF$test_cat)), function(i)
         "transparent");
   }
   gg_enrich <- gg_enrich +
      ggplot2::scale_fill_manual(name=test_cat_title,
         values=fill_sub);
   gg_enrich <- gg_enrich +
      ggplot2::ylab(ylab) +
      ggplot2::xlab(xlab) +
      ggplot2::facet_grid(facet_formula,
         scales=scales);
   if (length(main) > 0 && sum(nchar(main)) > 0) {
      gg_enrich <- gg_enrich +
         ggplot2::ggtitle(main);
   }
   gg_enrich;
}

#' Plot feature-peak overlap counts by distance
#'
#' Plot feature-peak overlap counts by distance
#'
#' This function is a small utility function intended
#' to display the individual counts produced by
#' `peak_feature_ranged_counter()`.
#'
#' @param counts_df `data.frame` produced by
#'    `peak_feature_ranged_counter()`.
#' @param color_sub `character` vector of R colors, named
#'    by values in the `Category` column of `counts_df`.
#' @param draw_points `logical` indicating if individual points
#'    should be drawn, or when `draw_points=FALSE` only the
#'    line is drawn.
#' @param lwd `numeric` indicating the line width.
#' @param ... additional arguments are ignored.
#'
#' @return `ggplot` object
#'
#' @export
plot_ranged_counts <- function
(range_counts,
 color_sub=NULL,
 draw_points=TRUE,
 lwd=1,
 ...)
{
   ## If random permutations were performed, use alpha transparency
   if (all(c("Iteration", "Test") %in% colnames(range_counts))) {
      gg_counts <- ggplot2::ggplot(range_counts,
         ggplot2::aes(x=Range,
            y=Count,
            alpha=Test,
            group=paste0(Experiment, "_", PeakSet, "_", Category, "_",
               Test, "_", Iteration),
            color=Category)) +
         ggplot2::scale_alpha_manual(values=c(`Test`=1, `Random`=0.1));
   } else {
      gg_counts <- ggplot2::ggplot(range_counts,
         ggplot2::aes(x=Range,
            y=Count,
            group=paste0(Experiment, "_", PeakSet, "_", Category),
            color=Category));
   }
   if (draw_points) {
      gg_counts <- gg_counts +
         ggplot2::geom_point();
   }
   gg_counts <- gg_counts +
      ggplot2::geom_line(size=lwd) +
      colorjam::theme_jam();
   if (length(color_sub) > 0) {
      gg_counts <- gg_counts +
         ggplot2::scale_color_manual(values=color_sub);
   } else {
      gg_counts <- gg_counts +
         colorjam::scale_color_jam();
   }
   gg_counts <- gg_counts +
      ggplot2::facet_grid(PeakSet~Experiment+Subclass) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10();
   gg_counts;
}

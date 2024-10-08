#' Plots expression profile of a gene set
#' @param pathway Gene set to plot.
#' @param E matrix with gene expression values
#' @param center a logical value indicating whether the gene expression should be centered to have zero mean before the analysis takes place.
#' The default is TRUE. The value is passed to \link[base]{scale}.
#' @param scale a logical value indicating whether the gene expression should be scaled to have unit variance before the analysis takes place.
#' The default is FALSE. The value is passed to \link[base]{scale}.
#' @param titles sample titles to use for labels
#' @param conditions sample grouping to use for coloring
#' @return ggplot object with the coregulation profile plot
#' @import data.table ggplot2
#' @export
plotCoregulationProfile <- function(pathway, E,
                                    center=TRUE,
                                    scale=FALSE,
                                    titles=colnames(E),
                                    conditions=NULL) {
    E <- t(base::scale(t(E), center=center, scale = scale))

    genes <- pathway

    dt <- as.data.table(E[rownames(E) %in% genes, , drop=FALSE], keep.rownames = TRUE)

    colnames(dt) <- c("gene", titles)

    dt[, id := seq_len(.N)]

    mdt <- melt(dt, measure.vars = colnames(dt)[2:(ncol(dt) - 1)], value.name = "expressionValue",
                variable.name = "sample",id.vars = c("id", "gene"))
    mdt[, gene := as.factor(gene)]
    mdt[, sample := factor(sample, levels=titles)]

    if (!is.null(conditions)) {
        if (is.character(conditions)) {
            conditions <- factor(conditions, levels=unique(conditions))
        }
    }

    pointDt <- data.table(x = seq_len(ncol(E)),
                          y = colSums(E[rownames(E) %in% genes, , drop=FALSE]) / sum(rownames(E) %in% genes),
                          condition=if (!is.null(conditions)) { conditions  } else "x")


    profilePlot <- ggplot(mdt, aes(x=sample, y=expressionValue, group=gene, color=gene),
                          show.legend=FALSE) +
        scale_color_discrete(guide="none") +
        geom_point(alpha = 0.1) +
        geom_path(alpha = 0.2) +
        geom_line(data = pointDt, aes(x = x, y = y),
                  group = "mean", color = "#13242a", linewidth = 1.5) +
        geom_hline(yintercept = min(pointDt$y), color = "#495057", linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = max(pointDt$y), color = "#495057", linetype = "dashed", linewidth = 1) +
        (if (!is.null(conditions)) {
            geom_point(shape=21, size=4,
                       data = pointDt,
                       aes(x = x, y = y, fill=condition),
                       group="mean", color="black")
        } else {
            geom_point(shape=21, size=4,
                       data = pointDt,
                       aes(x = x, y = y),
                       fill="black",
                       group="mean", color="black")
        }) +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ylab("expression") +
        NULL

    profilePlot
}

#' Plots table of gene set profiles.
#' @param gesecaRes Table with geseca results.
#' @param pathways Pathways to plot table, as in `geseca` function.
#' @param E gene expression matrix, as in `geseca` function.
#' @param center a logical value indicating whether the gene expression should be centered to have zero mean before the analysis takes place.
#' The default is TRUE. The value is passed to \link[base]{scale}.
#' @param scale a logical value indicating whether the gene expression should be scaled to have unit variance before the analysis takes place.
#' The default is FALSE. The value is passed to \link[base]{scale}.
#' @param colwidths Vector of five elements corresponding to column width for
#'      grid.arrange. Can be both units and simple numeric vector, in latter case
#'      it defines proportions, not actual sizes. If column width is set to zero, the column is not drawn.
#' @param titles sample titles to use an axis labels. Default to `colnames(E)`
#' @param colors vector of three colors to use in the color scheme
#' @param pathwayLabelStyle list with style parameter adjustments for pathway labels.
#'      For example, `list(size=10, color="red")` set the font size to 10 and color to red.
#'      See `cowplot::draw_text` for possible options.
#' @param headerLabelStyle similar to `pathwayLabelStyle` but for the table header.
#' @param valueStyle similar to `pathwayLabelStyle` but for pctVar and p-value columns.
#' @param axisLabelStyle list with style parameter adjustments for sample labels.
#'      See `ggplot2::element_text` for possible options.
#' @param axisLabelHeightScale height of the row with axis labels compared to other rows.
#'      When set to `NULL` the value is determined automatically.
#' @return ggplot object with gene set profile plots
#' @import ggplot2
#' @import cowplot
#' @export
plotGesecaTable <- function(gesecaRes,
                            pathways,
                            E,
                            center=TRUE,
                            scale=FALSE,
                            colwidths=c(5, 3, 0.8, 1.2, 1.2),
                            titles=colnames(E),
                            colors=c("blue", "white", "red"),
                            pathwayLabelStyle=NULL,
                            headerLabelStyle=NULL,
                            valueStyle=NULL,
                            axisLabelStyle=NULL,
                            axisLabelHeightScale=NULL){

    pathwayLabelStyleDefault <- list(size=12, hjust=1, x=0.95, vjust=0)
    pathwayLabelStyle <- modifyList(pathwayLabelStyleDefault, as.list(pathwayLabelStyle))

    headerLabelStyleDefault <- list(size=12)
    headerLabelStyle <- modifyList(headerLabelStyleDefault, as.list(headerLabelStyle))

    valueStyleDefault <- list(size=12, vjust=0)
    valueStyle <- modifyList(valueStyleDefault, as.list(valueStyle))

    axisLabelStyleDefault <- list(angle = 90, hjust = 1, size=10)
    axisLabelStyle <- modifyList(axisLabelStyleDefault, as.list(axisLabelStyle))

    if (is.null(axisLabelHeightScale)) {
        axisLabelHeightScale <- max(sapply(titles, nchar))/4*
                                axisLabelStyle$size/pathwayLabelStyle$size
    }

    gesecaRes <- gesecaRes[pathway %in% names(pathways)]
    pathways <- pathways[gesecaRes$pathway]
    # ^^ works with #40, as there can't be no empty pathways in the results


    E <- t(base::scale(t(E), center=center, scale = scale))
    colnames(E) <- titles

    pathways <- lapply(pathways, function(p) {
        unname(as.vector(na.omit(fmatch(p, rownames(E)))))
    })


    prjs <- t(do.call(cbind, lapply(pathways, function(p){
        scale(colSums(E[p, , drop=FALSE]))
    })))


    rownames(prjs) <- names(pathways)
    prjspd <- as.data.table(prjs, keep.rownames = "pathway")

    prjspd <- copy(melt(prjspd, id.vars = "pathway",
                        measure.vars = colnames(prjspd)[2:ncol(prjspd)],
                        variable.name = "sample", variable.factor=FALSE))
    prjspd[, pathway := factor(pathway, levels = rev(rownames(prjs)))]
    prjspd[, sample := factor(sample, levels=titles)]

    maxValue <- max(prjspd$value)
    minValue <- min(prjspd$value)

    color_legend <- cowplot::get_legend(
        ggplot(prjspd[pathway %fin% names(pathways)[[1]]],
               aes(x=sample, y=pathway, fill=value)) +
            geom_tile() +
            scale_fill_gradientn(limits=c(-3, 3), breaks=c(-3, 0, 3),
                                 oob=scales::squish,
                                 colors=colors,
                                 # guide = guide,
                                 name = "z-score"
            ) +
            theme(legend.position = "bottom"))

    ps <- lapply(names(pathways), function(pn) {
        p <- pathways[[pn]]
        annotation <- gesecaRes[match(pn, gesecaRes$pathway), ]
        list(
            cowplotText(pn, pathwayLabelStyle),
            ggplot(prjspd[pathway %fin% pn],
                   aes(x=sample, y=pathway, fill=value)) +
                geom_tile(color = "black", size = min(10/ncol(E), 0.5)) +
                scale_fill_gradientn(limits=c(-3, 3), breaks=c(-3, 0, 3),
                                     oob=scales::squish,
                                     colors=colors,
                                     # guide = guide,
                                     name = "z-score"
                ) +
                # scale_fill_gradient2(low = "blue",
                #                      high = "red",
                #                      mid = "white",
                #                      limit = c(minValue, maxValue),
                #                      space = "Lab") +
                theme(panel.background = element_blank(),
                      axis.line = element_blank(),
                      axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid = element_blank(),
                      axis.title = element_blank(),
                      plot.margin = unit(c(0.05,0,0.05,0), "npc"),
                      panel.spacing = unit(c(0.05,0,0.05,0), "npc"),
                      legend.position = "none") +
                # coord_equal() +
                NULL,
            cowplotText(sprintf("%.3f", annotation$pctVar), valueStyle),
            cowplotLabel(valueToExpExpression(annotation$pval), valueStyle),
            cowplotLabel(valueToExpExpression(annotation$padj), valueStyle)
        )
    })

    sampleTitle <- ggplot(data = data.table(sample=unique(prjspd$sample)), aes(x=sample)) +
        theme(panel.background = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              axis.title = element_blank(),
              plot.margin = unit(c(0,0,0,0), "npc"),
              panel.spacing = unit(c(0,0,0,0), "npc"),
              axis.title.x = element_blank(),
              axis.text.x = do.call(element_text, as.list(axisLabelStyle)))

    grobs <- c(
        list(nullGrob(),
             color_legend,
             nullGrob(),
             nullGrob(),
             nullGrob()),
        list(cowplotText("Pathway",
                         modifyList(headerLabelStyle, pathwayLabelStyle[c("hjust", "x")])
        )),
        lapply(c("Profile", "pctVar", "pval", "padj"), cowplotText, style=headerLabelStyle),
        unlist(ps, recursive = FALSE),
        list(nullGrob(),
             sampleTitle,
             nullGrob(),
             nullGrob(),
             nullGrob())
        )

    # not drawing column if corresponding colwidth is set to zero
    grobsToDraw <- rep(as.numeric(colwidths) != 0, length(grobs)/length(colwidths))

    p <- cowplot::plot_grid(plotlist=grobs[grobsToDraw],
                     ncol=sum(as.numeric(colwidths) != 0),
                     rel_widths=colwidths[as.numeric(colwidths) != 0],
                     rel_heights=c(1, 1, rep(1, length(pathways)), axisLabelHeightScale))

    p
}

#' Plot a spatial expression profile of a gene set
#' @param pathway Gene set to plot or a list of gene sets (see details below)
#' @param object Seurat object
#' @param title plot title
#' @param assay assay to use for obtaining scaled data, preferably with
#' the same universe of genes in the scaled data
#' @param colors vector of three colors to use in the color scheme
#' @param guide option for `ggplot2::scale_color_gradientn` to control for presence of the color legend
#' the same universe of genes in the scaled data
#' @param image.alpha adjust the opacity of the background images
#' @param ... optional arguments for \link[Seurat]{SpatialFeaturePlot}
#' @return ggplot object (or a list of objects) with the coregulation profile plot
#'
#' When the input is a list of pathways, pathway names are used for titles.
#' A list of ggplot objects a returned in that case.
#
#' @import ggplot2
#' @export
plotCoregulationProfileSpatial <- function(pathway,
                                           object,
                                           title=NULL,
                                           assay=DefaultAssay(object),
                                           colors=c("darkblue", "lightgrey", "darkred"),
                                           guide="colourbar",
                                           image.alpha = 0,
                                           ...) {
    stopifnot(requireNamespace("Seurat"))
    # TODO duplicated code with plotCoregulationProfileReduction
    if (is.list(pathway)) {
        if (is.null(title)) {
            titles <- names(pathway)
        }
        else {
            if (length(title) != length(pathway)) {
                stop("Length of the specified titles does not match count of pathways")
            }
            titles <- title
        }
        ps <- lapply(seq_along(pathway), function(i)
            plotCoregulationProfileSpatial(pathway[[i]],
                                           object = object,
                                           title = titles[i],
                                           assay = assay,
                                           colors = colors,
                                           image.alpha = image.alpha,
                                           ...))
        names(ps) <- names(pathway)
        ps <- unlist(ps, recursive = FALSE)
        return(ps)
    }


    obj2 <- addGesecaScores(list(pathway = pathway), object,
                                    assay = assay, scale = TRUE)
    ps <- Seurat::SpatialFeaturePlot(obj2, features = "pathway",
                                     combine = FALSE, image.alpha = image.alpha, ...)
    # suppress message of replacing existing color palette
    suppressMessages(ps <- lapply(ps, function(p){
        res <- p + scale_fill_gradientn(limits = c(-3, 3),
                                        breaks = c(-3, 0, 3),
                                        oob = scales::squish,
                                        colors = colors,
                                        guide = guide,
                                        name = "z-score")
        res <- res + theme(legend.position = theme_get()$legend.position)
        return(res)
    }))

    if (!is.null(title)){
        ps <- lapply(ps, function(p) p + ggtitle(title))
    }
    return(ps)
}

addGesecaScores <- function(pathways,
                            object,
                            assay=DefaultAssay(object),
                            prefix="",
                            scale=FALSE) {
    x <- GetAssay(object, assay)
    E <- x@scale.data

    res <- object


    for (i in seq_along(pathways)) {
        pathway <- pathways[[i]]
        pathway <- intersect(unique(pathway), rownames(E))

        score <- colSums(E[pathway, , drop=FALSE])/sqrt(length(pathway))
        score <- scale(score, center=TRUE, scale=scale)
        res@meta.data[[paste0(prefix, names(pathways)[i])]] <- score
    }


    return(res)
}

#' Plot a spatial expression profile of a gene set
#' @param pathway Gene set to plot or a list of gene sets (see details below)
#' @param object Seurat object
#' @param title plot title
#' @param assay assay to use for obtaining scaled data, preferably with
#' @param reduction reduction to use for plotting (one of the `Seurat::Reductions(object)`)
#' @param colors vector of three colors to use in the color scheme
#' @param guide option for `ggplot2::scale_color_gradientn` to control for presence of the color legend
#' the same universe of genes in the scaled data
#' @param ... additional arguments for Seurat::FeaturePlot
#' @return ggplot object (or a list of objects) with the coregulation profile plot
#'
#' When the input is a list of pathways, pathway names are used for titles.
#' A list of ggplot objects a returned in that case.
#
#' @import ggplot2
#' @export
plotCoregulationProfileReduction <- function(pathway, object, title=NULL,
                                             assay=DefaultAssay(object),
                                             reduction=NULL,
                                             colors=c("darkblue", "lightgrey", "darkred"),
                                             guide="colourbar",
                                             ...) {
    stopifnot(requireNamespace("Seurat"))

    if (is.list(pathway)) {
        if (is.null(title)) {
            titles <- names(pathway)
        } else {
            if (length(title) != length(pathway)) {
                stop("Length of the specified titles does not match count of pathways")
            }
            titles <- title
        }
        ps <- lapply(seq_along(pathway), function(i)
            plotCoregulationProfileReduction(pathway[[i]],
                                           object=object,
                                           title=titles[i],
                                           assay=assay,
                                           reduction=reduction,
                                           colors=colors,
                                           guide=guide,
                                           ...))
        names(ps) <- names(pathway)
        return(ps)
    }

    obj2 <- addGesecaScores(list(pathway=pathway), object, assay=assay,
                            scale=TRUE)

    p <- Seurat::FeaturePlot(obj2, features = "pathway",
                                    combine = FALSE, reduction=reduction, ...)[[1]]
    p <- p + coord_fixed()
    p$scales$scales[p$scales$find("color")] <- NULL

    # suppress message of replacing existing color palette
    suppressMessages(p2 <- p +
        scale_color_gradientn(limits=c(-3, 3), breaks=c(-3, 0, 3),
                             colors=colors,
                             oob=scales::squish,
                             guide=guide,
                             name = "z-score"
        ))

    if (!is.null(title)) {
        p2 <- p2 + ggtitle(title)
    }
    p2
}

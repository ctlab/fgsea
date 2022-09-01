#' @import data.table ggplot2
plotCoregulationProfile <- function(pathway, E,
                                    titles=colnames(E),
                                    conditions=NULL) {
    E <- t(scale(t(E), scale = FALSE))

    genes <- pathway

    dt <- as.data.table(E[rownames(E) %in% genes, ], keep.rownames = TRUE)

    colnames(dt) <- c("gene", titles)

    dt[, id := seq_len(.N)]

    mdt <- melt(dt, measure.vars = colnames(dt)[2:(ncol(dt) - 1)], value.name = "expressionValue",
                variable.name = "sample",id.vars = c("id", "gene"))
    mdt[, gene := as.factor(gene)]
    mdt[, sample := factor(sample, levels=titles)]

    pointDt <- data.table(x = seq_len(ncol(E)),
                          y = colSums(E[rownames(E) %in% genes, ]) / sum(rownames(E) %in% genes),
                          condition=if (!is.null(conditions)) { conditions  } else "x")


    profilePlot <- ggplot(mdt, aes(x=sample, y=expressionValue, group=gene, color=gene),
                          show.legend=FALSE) +
        scale_color_discrete(guide="none") +
        geom_point(alpha = 0.1) +
        geom_path(alpha = 0.2) +
        geom_line(data = pointDt, aes(x = x, y = y),
                  group = "mean", color = "#13242a", size = 1.5) +
        geom_hline(yintercept = min(pointDt$y), color = "#495057", linetype = "dashed", size = 1) +
        geom_hline(yintercept = max(pointDt$y), color = "#495057", linetype = "dashed", size = 1) +
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


plotGesecaTable <- function(pathways,
                            E,
                            gesecaRes,
                            colwidths=c(5, 3, 0.8, 1.2, 1.2),
                            titles=colnames(E),
                            render = TRUE){

    E <- t(scale(t(E), scale = FALSE))
    colnames(E) <- titles

    pathways <- lapply(pathways, function(p) {
        unname(as.vector(na.omit(fmatch(p, rownames(E)))))
    })

    # fixes #40
    pathways <- pathways[sapply(pathways, length) > 0]

    prjs <- t(do.call(cbind, lapply(pathways, function(p){
        colSums(E[p, ])/(length(p))
    })))


    rownames(prjs) <- names(pathways)
    prjspd <- as.data.table(prjs, keep.rownames = "pathway")

    prjspd <- copy(melt(prjspd, id.vars = "pathway",
                        measure.vars = colnames(prjspd)[2:ncol(prjspd)],
                        variable.name = "sample"))
    prjspd[, pathway := factor(pathway, levels = rev(rownames(prjs)))]

    maxValue <- max(prjspd$value)
    minValue <- min(prjspd$value)


    ps <- lapply(names(pathways), function(pn) {
        p <- pathways[[pn]]
        annotation <- gesecaRes[match(pn, gesecaRes$pathway), ]
        list(
            textGrob(pn, just=c("right", "centre"), x=unit(0.95, "npc")),
            ggplot(prjspd[pathway %fin% pn],
                   aes(x=sample, y=pathway, fill=value)) +
                geom_tile(color = "black", size = 0.75) +
                scale_fill_gradient2(low = "blue",
                                     high = "red",
                                     mid = "white",
                                     limit = c(minValue, maxValue),
                                     space = "Lab") +
                theme(panel.background = element_blank(),
                      axis.line = element_blank(),
                      axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid = element_blank(),
                      axis.title = element_blank(),
                      plot.margin = rep(unit(0, "null"), 4),
                      panel.spacing = rep(unit(0, "null"), 4),
                      legend.position = "none") +
                # coord_equal() +
                NULL,
            textGrob(sprintf("%.3f", annotation$pctVar)),
            textGrob(sprintf("%.1e", annotation$pval)),
            textGrob(sprintf("%.1e", annotation$padj))
        )
    })

    sampleTitle <- ggplot(data = prjspd, aes(x = sample)) +
        theme(panel.background = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              axis.title = element_blank(),
              plot.margin = rep(unit(0, "null"), 4),
              panel.spacing = rep(unit(0, "null"), 4),
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90))

    grobs <- c(
        list(textGrob("Pathway", just="right", x=unit(0.95, "npc"))),
        lapply(c("Projection", "pctVar", "pval", "padj"), textGrob),
        unlist(ps, recursive = FALSE),
        list(nullGrob(),
             sampleTitle,
             nullGrob(),
             nullGrob(),
             nullGrob()),
        list(nullGrob(),
             nullGrob(),
             nullGrob(),
             nullGrob(),
             nullGrob()))


    grobsToDraw <- rep(as.numeric(colwidths) != 0, length(grobs)/length(colwidths))

    p <- arrangeGrob(grobs=grobs[grobsToDraw],
                     ncol=sum(as.numeric(colwidths) != 0),
                     widths=colwidths[as.numeric(colwidths) != 0])
    if (render) {
        grid.draw(p)
    } else {
        as_ggplot(p)
    }
}

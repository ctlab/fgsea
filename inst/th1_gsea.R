library(reactome.db)
library(org.Mm.eg.db)
library(data.table)
library(fgsea)

rnk.file <- "./inst/extdata/naive.vs.th1.rnk"
gmt.file <- "./inst/extdata/mouse.reactome.gmt"

##### Loading pathways from Reactome

mouse.universe <- keys(org.Mm.eg.db, "ENTREZID")

# Selecting reactome gene sets
pathways <- na.omit(select(reactome.db, keys=mouse.universe, c("PATHID"),
                           keytype = 'ENTREZID'))
pathways <- split(pathways$ENTREZID, pathways$PATHID)

pathway2name <- as.data.table(na.omit(select(reactome.db, names(pathways),
                                             c("PATHNAME"), 'PATHID')))
# Remove organism prefix
pathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)]
pathway2name <- setNames(pathway2name$PATHNAME, names=pathway2name$PATHID)
pathway2name <- iconv(pathway2name, "latin1", "ASCII", sub="")

pathway.lines <- sapply(names(pathways), function(p) {
    link <- p
    name <- paste0(p, "_", pathway2name[p])
    name <- gsub("[ ()/]+", "_", name)
    sprintf("%s\t%s\t%s", name, link, paste0(pathways[[p]], collapse="\t"))
})
write(pathway.lines, file=gmt.file)

##### Running Broad GSEA

tag <- "naive_vs_th1"
res.dir <- "./inst/broad_gsea"
system.time(
system2("java",
        c(
            "-cp",
            "~/Dropbox/gsea2-2.1.0.jar",
            "-Xmx512m",
            "xtools.gsea.GseaPreranked",
            "-gmx", gmt.file,
            "-collapse", "false",
            "-mode", "Max_probe",
            "-norm", "meandiv",
            "-nperm", "1000",
            "-rnk", rnk.file,
            "-scoring_scheme", "weighted",
            "-rpt_label", tag,
            "-include_only_symbols", "true",
            "-make_sets", "true",
            "-plot_top_x", "20",
            "-rnd_seed", "timestamp",
            "-set_max", "500",
            "-set_min", "15",
            "-rnd_seed", "42",
            "-out", res.dir,
            "-gui", "false")))


##### Running fast GSEA

ranks <- read.table(rnk.file,
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- structure(ranks$t, names=ranks$ID)

set.seed(42)
system.time(
fgseaRes <- fgsea(pathways = pathways, stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=1000,
                  nproc=1)
)
print(sum(fgseaRes$padj < 1e-2)) # 0

set.seed(42)
system.time(
fgseaRes <- fgsea(pathways = pathways, stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=1e4,
                  nproc=1)
)
print(sum(fgseaRes$padj < 1e-2)) # 78

set.seed(42)
system.time(
fgseaRes <- fgsea(pathways = pathways, stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=100000,
                  nproc=1)
)
print(sum(fgseaRes$padj < 1e-2)) # 77

fgseaRes <- fgseaRes[order(pval)]
write.table(fgseaRes, file="./inst/extdata/fgsea_res.tsv", sep="\t", quote = FALSE, row.names = FALSE)

system.time(x <- fgsea(pathways=pathways["5990988"], stats=ranks, nperm=1e6, nproc=4))

plotGsea <- function(pathway, stats) {
    stats <- sort(stats, decreasing = TRUE)
    S <- match(pathway, names(stats))
    r <- stats
    p <- 1

    S <- sort(S)

    m <- length(S)
    N <- length(r)
    NR <- (sum(abs(r[S])^p))
    rAdj <- abs(r[S])^p
    rCumSum <- cumsum(rAdj) / NR

    tops <- rCumSum - (S - seq_along(S)) / (N - m)
    bottoms <- tops - rAdj / NR
    maxP <- max(tops)

    xs1 <- S - seq_along(S)
    xs2 <- xs1
    x <- c(as.vector(rbind(xs1, xs2)), N-m)

    ys1 <- cumsum(rAdj)
    ys2 <- ys1 - rAdj
    y <- c(as.vector(rbind(ys2, ys1)), NR)

    qplot(x, y, geom="line") + theme_bw() +
        xlab("rank") + ylab("stat") +
        geom_abline(intercept=0, slope=NR/(N-m), linetype="dashed") +
        geom_abline(intercept=NR*maxP, slope=NR/(N-m),
                    linetype="dotted") +
        geom_vline(x=xs1[which.max(tops)],
                   linetype="dotted")
}

p <- plotGsea(pathways[["5991022"]], stats=ranks)

plotGseaUpdate <- function(rsample, k, stats) {
    stats <- sort(stats, decreasing = TRUE)
    r <- stats
    p <- 1

    prevMax <- -1

    while (TRUE) {
        S <- rsample[1:k]
        S <- sort(S)

        m <- length(S)
        N <- length(r)
        NR <- (sum(abs(r[S])^p))
        rAdj <- abs(r[S])^p
        rCumSum <- cumsum(rAdj) / NR

        tops <- rCumSum - (S - seq_along(S)) / (N - m)
        bottoms <- tops - rAdj / NR
        maxP <- max(tops)

        xs1 <- S - seq_along(S)
        xs2 <- xs1
        x <- c(as.vector(rbind(xs1, xs2)), N-m)

        ys1 <- cumsum(rAdj)
        ys2 <- ys1 - rAdj
        y <- c(as.vector(rbind(ys2, ys1)), NR)

        message(which.max(tops))
        if (prevMax - which.max(tops) > 3) {
            p2 <- p1 +
                geom_line(data=data.frame(x=x, y=y), aes(x=x, y=y), color="black") +
                geom_abline(intercept=0, slope=NR/(N-m),
                            linetype="dashed", color="black") +
                geom_abline(intercept=NR*maxP, slope=NR/(N-m),
                            linetype="dotted", color="black") +
                geom_vline(x=xs1[which.max(tops)],
                           linetype="dotted", color="black")
            break
        }
        prevMax <- which.max(tops)
        p1 <- ggplot() +
            geom_line(data=data.frame(x=x, y=y), aes(x=x, y=y), color="#777777") +
            geom_abline(intercept=0, slope=NR/(N-m),
                        linetype="dashed", color="#777777") +
            geom_abline(intercept=NR*maxP, slope=NR/(N-m),
                        linetype="dotted", color="#777777") +
            geom_vline(x=xs1[which.max(tops)],
                       linetype="dotted", color="#777777")
        k <- k-1
    }

    message(sprintf("k - 1  = %s", k))
    p2 + theme_bw() + xlab("rank") + ylab("stat")
}

set.seed(42)
p <- plotGseaUpdate(rsample=sample(seq_along(ranks), size=400), k=200, stats=ranks)
pp <- p + xlim(0, 2500) + ylim(0, 450)


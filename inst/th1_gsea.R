library(reactome.db)
library(org.Mm.eg.db)
library(data.table)
library(fgsea)

rnk.file <- "./inst/naive.vs.th1.rnk"
gmt.file <- "./inst/mouse.reactome.gmt"

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
pathway2name <- structure(pathway2name$PATHNAME, names=pathway2name$PATHID)

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

ranks <- read.table("./inst/naive.vs.th1.rnk",
                    header=T, colClasses = c("character", "numeric"))
ranks <- structure(ranks$t, names=ranks$ID)

set.seed(42)
system.time(
fgseaRes <- fgsea(pathways = pathways, stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=2,
                  nproc=1)
)
print(sum(fgseaRes$padj < 1e-2)) # 0

set.seed(42)
system.time(
fgseaRes <- fgsea(pathways = pathways, stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=20,
                  nproc=1)
)
print(sum(fgseaRes$padj < 1e-2)) # 76

set.seed(42)
system.time(
fgseaRes <- fgsea(pathways = pathways, stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=200,
                  nproc=1)
)
print(sum(fgseaRes$padj < 1e-2)) # 78

fgseaRes <- fgseaRes[order(pval)]
write.table(fgseaRes, file="./inst/fgsea_res.tsv", sep="\t", quote = F, row.names = F)

system.time(x <- fgsea(pathways=pathways["5992313"], stats=ranks, nperm=1e6, nproc=4))

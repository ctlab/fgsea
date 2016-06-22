library(reactome.db)
library(org.Mm.eg.db)
library(data.table)

mouse.universe <- keys(org.Mm.eg.db, "ENTREZID")

# Selecting reactome gene sets
pathways <- na.omit(select(reactome.db, keys=mouse.universe, c("PATHID"),
                           keytype = 'ENTREZID'))
pathways <- split(pathways$ENTREZID, pathways$PATHID)

pathway2name <- as.data.table(na.omit(select(reactome.db, names(pathways),
                                             c("PATHNAME"), 'PATHID')))
# Remove organism prefix
pathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)]
pathway2name <- setNames(pathway2name$PATHNAME, pathway2name$PATHID)
pathway2name <- iconv(pathway2name, "latin1", "ASCII", sub="")

pathway.lines <- sapply(names(pathways), function(p) {
    link <- p
    name <- paste0(p, "_", pathway2name[p])
    name <- gsub("[ ()/]+", "_", name)
    sprintf("%s\t%s\t%s", name, link, paste0(pathways[[p]], collapse="\t"))
})
write(pathway.lines, file="./inst/extdata/mouse.reactome.gmt")

#!/usr/bin/env Rscript
rnk.file <- "./inst/extdata/naive.vs.th1.rnk"
gmt.file <- "./inst/extdata/mouse.reactome.gmt"

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
                "-nperm", "100",
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

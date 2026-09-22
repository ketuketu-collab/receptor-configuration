suppressPackageStartupMessages({library(Seurat); library(Matrix)})
d <- Sys.getenv("SCP3755_DIR", file.path(Sys.getenv("DATA_ROOT", "/Volumes/M4_SSD/projects"), "scp3755_uc_atlas"))
o <- readRDS(file.path(d,"scp3755.RDS"))
cat("class:",class(o),"\n"); md <- o@meta.data
cat("cells:",nrow(md),"assays:",names(o@assays),"default:",DefaultAssay(o),"\n")
a <- o@assays[[DefaultAssay(o)]]
cat("assay slots:",slotNames(a),"\n")
cnt <- NULL
for (s in c("counts","data")) if (s %in% slotNames(a)) { m <- slot(a,s); if (!is.null(m) && nrow(m)>0 && is.null(cnt)) { cnt <- m; cat("using slot:",s," dim:",dim(m)," colsum1:",sum(m[,1])," max:",max(m[,1]),"\n") } }
stopifnot(!is.null(cnt))
grp <- paste(md$SAMID, md$Annotations_Level_2, sep="|")
idx <- split(seq_along(grp), grp); ug <- names(idx)
cat("pseudobulk groups:",length(ug),"\n")
PB <- sapply(ug, function(g) Matrix::rowSums(cnt[, idx[[g]], drop=FALSE]))
rownames(PB) <- rownames(cnt)
info <- do.call(rbind, strsplit(ug,"\\|", perl=TRUE))
meta <- data.frame(group=ug, SAMID=info[,1], celltype=info[,2],
                   n_cells=sapply(ug,function(g) length(idx[[g]])), stringsAsFactors=FALSE)
write.csv(as.matrix(PB), file.path(d,"pb_counts.csv"), quote=FALSE)
write.csv(meta, file.path(d,"pb_meta.csv"), row.names=FALSE, quote=FALSE)
cat("done. genes:",nrow(PB),"\n")

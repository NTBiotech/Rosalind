dna = as.character(read.table("rosalind_revc.txt"))

c = c()
bases = c("G", "A", "T", "C")
comp = data.frame(bases, rev(bases))
names(comp) = c("bases", "cobases")
for (b in as.list(strsplit(dna, split = ""))[[1]]) {
    c = append(c, comp[["cobases"]][which(comp[["bases"]] == b)])
    }

revc = paste(as.character(rev(c)), collapse = "")

cat(revc, file = "reverse_complement_dna.txt")

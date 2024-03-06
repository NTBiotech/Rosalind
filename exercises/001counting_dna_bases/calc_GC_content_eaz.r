getwd()
seq = as.character(read.table("rosalind_dna.txt", ))

bcont = function(base, sequ) {
    v = 0
    for (x in strsplit(sequ, "")) {
        for (b in x) {
            if (b == base) {
                v = v+1
    }}
}
    return(v)
}

bases <- c("A", "C", "G", "T")

df = data.frame(row.names = bases)

basecontent = c(bcont("A", seq), bcont("C", seq), 
    bcont("G", seq), bcont("T", seq))
cbind(df, basecontent)

barplot(basecontent, names.arg = bases)

dna = as.character(read.table("rosalind_rna.txt"))

transcribe = function(seq) {
    rseq <- c()
    for (x in as.list(strsplit(seq, "")[[1]])) {
        print(x)
        if (x == "T") {rseq = append(rseq, "U")
        }   else {rseq = append(rseq, x)}}
    return(paste(rseq, collapse = ""))
}
rna = transcribe(dna)
cat(rna, file = "trancribed_rna.txt")

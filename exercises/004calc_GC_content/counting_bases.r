#Problem
#
#A string is simply an ordered collection of symbols selected from some alphabet and formed into a word; the length of a string is the number of symbols that it contains.
#
#An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."
#
#Given: A DNA string s of length at most 1000 nt.
#
#Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s

#setwd("calc_GC_content")
getwd()

seq = read.table("rosalind_gc.txt", )

cols = apply(seq, 1, function(seq) {return(grepl(">Rosalind", seq))})
numcols <- which(cols)
headers = seq[cols, ]
cseq = list()
G_content = c()
for (n in numcols){
    c = match(n, numcols)
    print(c)
    if(c==length(numcols)){
        x = seq[n:nrow(seq), 1]
    }else{
        x = seq[n:numcols[c+1]-1,1]
        if(c!=1){x = x[2:length(x)]}
        }
        x = x[2:length(x)]
        x1 = ""
        for (y in x){
            x1 = paste(x1, y, sep = "")
        }
        cseq = append(cseq, x1)
    print(x1)
    
}
seqlist = setNames(cseq, headers)

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

df = data.frame(row.names = c("G", "A", "T", "C"))
for (s in seqlist) {
    cont = c(bcont("G", s),bcont("A", s), bcont("T", s),bcont("C", s))
    df = cbind(df, cont)
}

colnames(df) = headers
df
rbind(df, colSums(df))

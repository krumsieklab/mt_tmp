preD <- D

View(metadata(D))
View(metadata(preD))

all.equal(assay(D), assay(preD))
all.equal(as.data.frame(colData(D)), as.data.frame(colData(preD)))
all.equal(as.data.frame(rowData(D)), as.data.frame(rowData(preD)))

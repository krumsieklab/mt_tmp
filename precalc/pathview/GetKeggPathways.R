require(AnnotationDbi)
require(org.Hs.eg.db)
require(graphite)

# download kegg pathways
pwdb <- pathways(species = "hsapiens", database = "kegg")
pwdb <- mcmapply(function(pwname) {
  pw <- pwdb[[pwname]]
  pw %>% 
    graphite::edges(which = "mixed") %>% # "metabolites", "proteins", or "mixed"
    mutate(name = pwname,
           ID = pathwayId(pw))
},
names(pwdb),
SIMPLIFY = F,
mc.cleanup = T,
mc.cores = 3)

save(pwdb,file = codes.makepath("snippets/packages/metabotools_external/pathview/KeggPathways.Rds"))


# download kegg pathways
pwdb <- pathways(species = "mmusculus", database = "kegg")
pwdb <- mcmapply(function(pwname) {
  pw <- pwdb[[pwname]]
  pw %>% 
    graphite::edges(which = "mixed") %>% # "metabolites", "proteins", or "mixed"
    mutate(name = pwname,
           ID = pathwayId(pw))
},
names(pwdb),
SIMPLIFY = F,
mc.cleanup = T,
mc.cores = 3)

save(pwdb,file = codes.makepath("snippets/packages/metabotools_external/pathview/KeggPathways_mouse.Rds"))

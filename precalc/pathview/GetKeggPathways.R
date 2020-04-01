library(AnnotationDbi)
library(org.Hs.eg.db)
library(graphite)

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

save(pwdb,file = data.makepath("MT_precalc/pathview/KeggPathways.Rds"))


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

save(pwdb,file = data.makepath("MT_precalc/pathview/KeggPathways_mouse.Rds"))

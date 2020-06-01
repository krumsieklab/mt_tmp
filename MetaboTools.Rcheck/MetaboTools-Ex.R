pkgname <- "MetaboTools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "MetaboTools-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('MetaboTools')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("mt_anno_metabolites_HMDB2KEGG")
### * mt_anno_metabolites_HMDB2KEGG

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_anno_metabolites_HMDB2KEGG
### Title: Create KEGG identifiers from HMDB identifiers.
### Aliases: mt_anno_metabolites_HMDB2KEGG

### ** Examples

# Load data, sheet with sample annotations, sheet with metabolite annotations and add KEGG identifiers
## Not run: 
##D D <-
##D   # load raw data
##D   mt_files_data_xls(file=file, sheet="data", samplesInRows=T, ID="SAMPLE_NAME") %>%
##D   # sample annotations from metabolomics run
##D   mt_files_anno_xls(file=file, sheet="sampleinfo", annosfor="samples", IDanno = "SAMPLE_NAME") %>%
##D   # metabolite annotations
##D   mt_files_anno_xls(file=file, sheet="metinfo", annosfor="metabolites", IDanno="BIOCHEMICAL", IDdata = "name") %>%
##D   # add KEGG identifiers
##D   mt_anno_metabolites_HMDB2KEGG(in_col="HMDb_id", out_col="KEGG") %>%
##D   ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_anno_metabolites_HMDB2KEGG", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_anno_pathways_HMDB")
### * mt_anno_pathways_HMDB

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_anno_pathways_HMDB
### Title: Add pathway information.
### Aliases: mt_anno_pathways_HMDB

### ** Examples

# annotate metabolites using smp_db
## Not run: 
##D ... %>%
##D   mt_anno_pathways_HMDB(in_col = "HMDb_ID", out_col = "smp_db",
##D                         pwdb_name = "SMP",
##D                         db_dir = codes.makepath("snippets/packages/metabotools_external/hmdb")) %>%
##D ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_anno_pathways_HMDB", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_anno_pathways_from_file")
### * mt_anno_pathways_from_file

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_anno_pathways_from_file
### Title: Add pathway information.
### Aliases: mt_anno_pathways_from_file

### ** Examples

## Not run: 
##D # annotate metabolites using the SMP column of the pathway database flat file
##D ... %>%
##D       mt_anno_pathways_from_file(in_col = "HMDb_ID",
##D                                  out_col = "janpw",
##D                                  file = codes.makepath("packages/metabotools_external/hmdb/hmdb_preprocessed_4.0.csv"),
##D                                  met_ID_col = "HMDB_id",
##D                                  pw_ID_col = "SMP",
##D                                  pw_name_col = "pathway_name") %>%
##D ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_anno_pathways_from_file", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_anno_pathways_graphite")
### * mt_anno_pathways_graphite

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_anno_pathways_graphite
### Title: Add pathway information.
### Aliases: mt_anno_pathways_graphite

### ** Examples

# annotate metabolites using the humancyc database from the graphite package
## Not run: 
##D ... %>%
##D     mt_add_pathways(in_col = "KEGG",
##D                     out_col = "humancyc_db",
##D                     pwdb_species = "hsapiens",
##D                     pwdb_name = "humancyc",
##D                     n_cpus = 5) %>%
##D ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_anno_pathways_graphite", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_anno_pathways_remove_redundant")
### * mt_anno_pathways_remove_redundant

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_anno_pathways_remove_redundant
### Title: Remove redundant pathway annotations.
### Aliases: mt_anno_pathways_remove_redundant

### ** Examples

# first annotate metabolites using smp_db and then remove redundant pathways
## Not run: 
##D ... %>%
##D   mt_anno_pathways_HMDB(in_col = "HMDb_ID", out_col = "smp_db",
##D   pwdb_name = "SMP", db_dir = codes.makepath("snippets/packages/metabotools_external/hmdb")) %>%
##D   mt_anno_pathways_remove_redundant(met_ID_col = "HMDb_ID", pw_ID_col = "smp_db") %>%
##D ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_anno_pathways_remove_redundant", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_files_anno_xls")
### * mt_files_anno_xls

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_files_anno_xls
### Title: Load annotations from Excel file.
### Aliases: mt_files_anno_xls

### ** Examples

# Load data, two sheets with sample annotations, and one sheet with metabolite annotations from the same file
## Not run: 
##D D <- 
##D   # load raw data
##D   mt_files_data_xls(file=file, sheet="data", samplesInRows=T, ID="SAMPLE_NAME") %>% 
##D   # sample annotations from metabolomics run
##D   mt_files_anno_xls(file=file, sheet="sampleinfo", anno_type="samples", anno_ID = "SAMPLE_NAME") %>% 
##D   # sample annotations from clinical table
##D   mt_files_anno_xls(file=file, sheet="clinicals", anno_type="samples", anno_ID="SAMPLE_NAME") %>% 
##D   # metabolite annotations`
##D   mt_files_anno_xls(file=file, sheet="metinfo", anno_type="metabolites", anno_ID="BIOCHEMICAL", data_ID = "name") %>% 
##D   ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_files_anno_xls", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_files_checksum")
### * mt_files_checksum

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_files_checksum
### Title: Validate MD5 checksum of file.
### Aliases: mt_files_checksum

### ** Examples

# first call, to get the checksum (will crash, deliberately)
## Not run: ... %>% mt_files_checksum(file="input.xlsx", checksum="") %>% ...   

# copy-paste the correct ('actual') checksum from the error message into the call:
## Not run: ... %>% mt_files_checksum(file="input.xlsx", checksum="688048bd1eb9c771be0eb49548a6f947") %>% ...   




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_files_checksum", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_files_data_xls")
### * mt_files_data_xls

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_files_data_xls
### Title: Load data matrix from Excel file.
### Aliases: mt_files_data_xls

### ** Examples

# Load data, two sheets with sample annotations, and one sheet with metabolite annotations from the same file
## Not run: 
##D D <- 
##D   # load raw data
##D   mt_files_data_xls(file=file, sheet="data", samples_in_rows=T, ID_col="SAMPLE_NAME") %>% 
##D   # sample annotations from metabolomics run
##D   mt_files_anno_xls(file=file, sheet="sampleinfo", annosfor="samples", IDanno = "SAMPLE_NAME") %>% 
##D   # sample annotations from clinical table
##D   mt_files_anno_xls(file=file, sheet="clinicals", annosfor="samples", IDanno="SAMPLE_NAME") %>% 
##D   # metabolite annotations`
##D   mt_files_anno_xls(file=file, sheet="metinfo", annosfor="metabolites", IDanno="BIOCHEMICAL", IDdata = "name") %>% 
##D   ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_files_data_xls", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_files_load_UCD")
### * mt_files_load_UCD

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_files_load_UCD
### Title: Load UC Davis-format data.
### Aliases: mt_files_load_UCD

### ** Examples

## Not run: 
##D D <- 
##D   # load data
##D   mt_files_load_UCD(codes.makepath("Mt/sampledata.xlsx"), "OrigScale") %>%
##D   ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_files_load_UCD", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_files_load_WCM")
### * mt_files_load_WCM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_files_load_WCM
### Title: Load WCM core-format data.
### Aliases: mt_files_load_WCM

### ** Examples

## Not run: 
##D D <-
##D   # load data in WCM format
##D   mt_files_load_WCM(
##D     file='DM369_Metabolites separated by nucleotides 12_6_18_altered.xlsx',
##D     sheet='peakHeight_metabolite_intensiti') %>%
##D   ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_files_load_WCM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_files_load_metabolon")
### * mt_files_load_metabolon

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_files_load_metabolon
### Title: Load Metabolon-format data.
### Aliases: mt_files_load_metabolon

### ** Examples

## Not run: 
##D D <- 
##D   # load data
##D   mt_files_load_metabolon(codes.makepath("Mt/sampledata.xlsx"), "OrigScale") %>%
##D   ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_files_load_metabolon", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_files_write_pathwayannos")
### * mt_files_write_pathwayannos

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_files_write_pathwayannos
### Title: Write out all pathway annotations, possible redundant pathways,
###   and metabolite annotations.
### Aliases: mt_files_write_pathwayannos

### ** Examples

## Not run: 
##D ... %>% 
##D mt_add_pathways(in_col = "KEGG",
##D     out_col = "kegg_db",
##D     pw_species = "hsapiens",
##D     pw_name = "kegg") %>%
##D mt_anno_pathways_remove_redundant(met_ID_col = "KEGG", pw_col = "kegg_db") %>%
##D mt_files_write_pathwayannos(pwfield='kegg_db', file='pwannos.xlsx')
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_files_write_pathwayannos", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_files_write_stats")
### * mt_files_write_stats

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_files_write_stats
### Title: Export statistical results from pipeline into Excel file.
### Aliases: mt_files_write_stats

### ** Examples

## Not run: 
##D # Write out all results
##D ... %>% mt_files_write_stats(file="results.xlsx") %>%
##D # Write out specific result]
##D ... %>% mt_files_write_stats(file="results.xlsx", compnames="comp1") %>%
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_files_write_stats", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_files_write_xls")
### * mt_files_write_xls

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_files_write_xls
### Title: Outputs assay, colData and rowData into an Excel file.
### Aliases: mt_files_write_xls

### ** Examples

## Not run: %>% mt_files_write_xls(file = "out.xlsx") %>%



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_files_write_xls", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_modify_aggPW")
### * mt_modify_aggPW

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_modify_aggPW
### Title: Generate aggregated pathway values.
### Aliases: mt_modify_aggPW

### ** Examples

## Not run: 
##D %>% mt_modify_aggPW(pw_col="SUB_PATHWAY", method="aggmean") %>%  # subpathways from metabolon
##D 
##D # add KEGG pathways and use those
##D %>%
##D   mt_anno_pathways_HMDB(in_col = "HMDB", out_col = "kegg_db", pwdb_name = "KEGG", db_dir = codes.makepath("snippets/packages/metabotools_external/hmdb")) %>%
##D   mt_anno_pathways_remove_redundant(met_ID_col = "HMDB", pw_col = "kegg_db") %>%
##D   mt_modify_aggPW(pw_col="kegg_db", method="aggmean") %>%
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_modify_aggPW", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_modify_applytoanno")
### * mt_modify_applytoanno

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_modify_applytoanno
### Title: Apply function to metabolite or sample annotation column.
### Aliases: mt_modify_applytoanno

### ** Examples

 ## Not run: 
##D ... %>%
##D  # ensure factor for casecontrol variable
##D  mt_modify_applytoanno(
##D    anno_type='samples',
##D    col_name='casecontrol',
##D    fun=as.factor) %>%
##D  ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_modify_applytoanno", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_modify_averagesample")
### * mt_modify_averagesample

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_modify_averagesample
### Title: mt-modify_averagesample
### Aliases: mt_modify_averagesample

### ** Examples

## Not run: ... %>% mt_modify_averagesample(group_by = "RID") %>% ...




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_modify_averagesample", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_modify_cv")
### * mt_modify_cv

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_modify_cv
### Title: mt_modify cv
### Aliases: mt_modify_cv

### ** Examples

## Not run: ... %>% mt_modify_cv(qc_samples=="PQC", col_lab = "PQC_cv") %>% ...




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_modify_cv", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_modify_filter_metabolites")
### * mt_modify_filter_metabolites

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_modify_filter_metabolites
### Title: Filter metabolites.
### Aliases: mt_modify_filter_metabolites

### ** Examples

## Not run: ... %>% mt_modify_filter_metabolites(SUPER_PATHWAY=="Nucleotide") %>% ...




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_modify_filter_metabolites", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_modify_filter_samples")
### * mt_modify_filter_samples

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_modify_filter_samples
### Title: Filter samples.
### Aliases: mt_modify_filter_samples

### ** Examples

# filter to two specific groups of samples
## Not run: ... %>% mt_modify_filter_samples(sample_filter = GROUP %in% c("FL","ctrl")) %>% ...




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_modify_filter_samples", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_modify_mutate")
### * mt_modify_mutate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_modify_mutate
### Title: Create new variable by expression
### Aliases: mt_modify_mutate

### ** Examples

## Not run: 
##D # Convert numeric sample annotation field 'num' to factor
##D ...  %>%
##D  mt_modify_mutate(anno_type="samples", col_name="num", term = as.factor(num)) %>% ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_modify_mutate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_modify_ratios")
### * mt_modify_ratios

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_modify_ratios
### Title: Generate metabolite ratios.
### Aliases: mt_modify_ratios

### ** Examples

## Not run: 
##D # Transform dataset to ratios
##D ... %>%  mt_modify_ratios() %>% ... # proceed with statistical analysis
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_modify_ratios", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_modify_sampleanno_reorder")
### * mt_modify_sampleanno_reorder

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_modify_sampleanno_reorder
### Title: Reorder a sample annotation.
### Aliases: mt_modify_sampleanno_reorder

### ** Examples

## Not run: %>% mt_modify_sampleanno_reorder("Group", c('WT_Norm_F','WT_Hyp_F','KO_Norm_F','KO_Hyp_F','WT:KO2_F','WT:WT2_F')) %>%




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_modify_sampleanno_reorder", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_PCA")
### * mt_plots_PCA

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_PCA
### Title: 2D PCA of samples.
### Aliases: mt_plots_PCA

### ** Examples

## Not run: 
##D ## PCA on scaledata, color and shape by "Group" variable in colData
##D ... $>$ mt_plots_PCA(scaledata=T, color=Group, shape=Group, title="PCA - scaled data") %>% ...
##D ## PCA scores plot on non-scaled data, with ellipse and extra explained variance plot, and two ggadds (white background and centering of title)
##D mt_plots_PCA(title="PCA scores", show = 'scores', scaledata=F, PCa=1, PCb=2, ellipse=0.95, expvarplot=T, ggadd = theme_bw() + theme(plot.title=element_text(hjust=0.5)))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_PCA", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_UMAP")
### * mt_plots_UMAP

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_UMAP
### Title: 2D UMAP of samples.
### Aliases: mt_plots_UMAP

### ** Examples

## Not run: 
##D ## UMAP on scaledata, color and shape by "Group" variable in colData
##D ... $>$ mt_plots_UMAP(scaledata=T, color=Group, shape=Group, title="UMAP - scaled data") %>% ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_UMAP", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_boxplot")
### * mt_plots_boxplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_boxplot
### Title: Boxplots
### Aliases: mt_plots_boxplot

### ** Examples

## Not run: 
##D # boxplots as overview of results with a result already in 'comp'
##D # color by "Group" variable in colData
##D mt_plots_boxplot(x                  = Group,
##D                  statname           = "comp",
##D                  correct_confounder = ~BATCH_MOCK,
##D                  metab_filter       = p.value<0.01,
##D                  metab_sort         = p.value,
##D                  annotation         = "{sprintf('P-value: %.1e', p.value)}\nStatistic: {sprintf('%.2f', statistic)}",
##D                  rows               = 2,
##D                  cols               = 2) %>%
##D                  ...
## End(Not run)

 @author Jonas Zierer

 @import ggplot2
 @importFrom magrittr %>% %<>%
 @import SummarizedExperiment

 @export



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_boxplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_compare2stats")
### * mt_plots_compare2stats

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_compare2stats
### Title: Comparative plot between two comparisons.
### Aliases: mt_plots_compare2stats

### ** Examples

## Not run: 
##D ## compare two stats from inside the same pipeline
##D ... %>%
##D mt_plots_compare2stats(stat1='WT',
##D   filter1= p.adj<0.1,
##D   stat2='KO',
##D   filter2= p.adj<0.1,
##D   filterop = 'OR'
##D ) %>% ...
##D 
##D ## compare two stats from different pipelines, as part of the pipeline of the second
##D # 'comp' is a string that contains the name of a comparison (here both SEs have the same comparison on two datasets)
##D .. %>% mt_plots_compare2stats(
##D   stat1 = comp, filter1 = p.adj<0.1,
##D   D2 = firstPipeSE, stat2 = comp, filter2 = p.adj<0.1,
##D   filterop = "OR") %>% ...
##D 
##D ## compare two stats from different pipelines, output as plot object
##D ## not part of the actual MT pipelines, but separate call
##D # 'comp' is a string that contains the name of a comparison (here both SEs have the same comparison on two datasets)
##D gg <- mt_plots_compare2stats(
##D   D1 = D1, stat1 = comp, filter1 = p.adj<0.1,
##D   D2 = D2, stat2 = comp, filter2 = p.adj<0.1,
##D   filterop = "OR", return.plot.only=T)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_compare2stats", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_equalizer")
### * mt_plots_equalizer

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_equalizer
### Title: 'Equalizer' plots.
### Aliases: mt_plots_equalizer

### ** Examples

## Not run: 
##D # super-pathway / sub-pathway equalizer
##D # sub-pathway analysis must already be stored in D_sub, and this is part of the super-pathway pipeline, with a result already in 'comp'
##D  ... %>%
##D  mt_plots_equalizer(
##D   comp1='comp',
##D   D2=D_sub,
##D   comp2=='comp',
##D   legend.fine="sub pathway",
##D   legend.coarse='super pathway',
##D   vertline.fine = p.adj < 0.1,
##D   vertline.coarse = p.adj < 0.1) %>%
##D ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_equalizer", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_net")
### * mt_plots_net

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_net
### Title: mt_plots_net
### Aliases: mt_plots_net

### ** Examples

#' # in the context of a SE pipeline
... %>% mt_plots_net(statsname = "xxx") %>% ...    # standard call
... %>% mt_plots_net(statsname = "xxx", corr_filter = p.adj < 0.5, node_coloring="Li's", save.html="Network.html", height=800) %>% ...    # filters only significant correlations and colors the nodes according to the results in the indicated test, saves visnetwork to file




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_net", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_pathview")
### * mt_plots_pathview

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_pathview
### Title: mt_plots_pathview
### Aliases: mt_plots_pathview

### ** Examples

## Not run: 
##D # plot all pathways with at least one significant metabolite from the statistical comparison "comp" in them
##D mt_plots_pathview(D = D,
##D                   met.id="KEGG_mapped",
##D                   statname = "comp",
##D                   color.scale = -sign(fc)*log10(p.adj),
##D                   color.range = -log10(0.01),
##D                   metab.filter = p.adj < 0.05,
##D                   show.only.filtered = TRUE,
##D                   path.database = "./Pathview_database",
##D                   path.output = "./results/pathview",
##D                   same.layer = F,
##D                   add.pwname.suffix = T
##D                   ) %>%
##D                   ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_pathview", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_pheatmap")
### * mt_plots_pheatmap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_pheatmap
### Title: Heatmap plot by pheatmap::pheatmap
### Aliases: mt_plots_pheatmap

### ** Examples

## Not run: 
##D D %>%
##D mt_plots_pheatmap(annotation_row = c("SUPER_PATHWAY", "PLATFORM", "RI"),
##D                   annotation_col = c("GROUP_DESC","BATCH_MOCK","gender"),
##D                   fD = function(x) scale(exp(scale(x))),
##D                   clustering_distance_cols =  "correlation",
##D                   clustering_distance_rows = "minkowski")
##D  
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_pheatmap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_pvalhist")
### * mt_plots_pvalhist

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_pvalhist
### Title: Generate p-value histogram.
### Aliases: mt_plots_pvalhist

### ** Examples

... %>% mt_plots_pvalhist() %>% ...                  # for all
... %>% mt_plots_pvalhist(statnames='comp') %>% ...  # for one




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_pvalhist", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_pvalqq")
### * mt_plots_pvalqq

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_pvalqq
### Title: Generate p-value qq plot
### Aliases: mt_plots_pvalqq

### ** Examples

## Not run: ... %>% mt_plots_pvalhist(statnames='comp') %>% ...  # for one





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_pvalqq", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_qc_dilutionplot")
### * mt_plots_qc_dilutionplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_qc_dilutionplot
### Title: Correlate variable with dilution factors from quotient
###   normalization.
### Aliases: mt_plots_qc_dilutionplot

### ** Examples

%>% mt_plots_qc_dilutionplot(comp="group") %>% # compare with 'group' sample annotation




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_qc_dilutionplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_qc_missingness")
### * mt_plots_qc_missingness

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_qc_missingness
### Title: Visualize missing value structure of dataset.
### Aliases: mt_plots_qc_missingness

### ** Examples

%>%  mt_plots_qc_missingness() %>% # without horizontal line
%>%  mt_plots_qc_missingness(metMax=0.5) %>% # with horizontal line at 50%




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_qc_missingness", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_sampleboxplot")
### * mt_plots_sampleboxplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_sampleboxplot
### Title: Boxplot of samples.
### Aliases: mt_plots_sampleboxplot

### ** Examples

## Not run: 
##D ## sample boxplot, color by colData 'group' variable, with specific title, on log scale,
##D ... %>% mt_plots_sampleboxplot(color=group, plottitle='after quotient normalization', logged=T) %>% ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_sampleboxplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_scatter")
### * mt_plots_scatter

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_scatter
### Title: Scatter plots
### Aliases: mt_plots_scatter

### ** Examples

## Not run: 
##D # scatter plots as overview of results with a result already in 'comp'
##D # color by "age" variable in colData
##D mt_plots_boxplot(x                  = age,
##D                  statname           = "comp",
##D                  correct_confounder = ~BATCH_MOCK,
##D                  metab_filter       = p.value<0.01,
##D                  metab_sort         = p.value,
##D                  annotation         = "{sprintf('P-value: %.1e', p.value)}\nStatistic: {sprintf('%.2f', statistic)}",
##D                  rows               = 2,
##D                  cols               = 2) %>%
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_scatter", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_plots_volcano")
### * mt_plots_volcano

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_plots_volcano
### Title: mt_plots_volcano
### Aliases: mt_plots_volcano

### ** Examples

## Not run: 
##D # Volcano plot as overview of results with a result already in 'comp'
##D ... %>%
##D mt_plots_volcano(statname     = "comp",
##D  metab_filter = p.adj < 0.1,
##D  colour       = p.value < 0.05) %>%
##D  ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_plots_volcano", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_post_addFC")
### * mt_post_addFC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_post_addFC
### Title: Compute fold-change
### Aliases: mt_post_addFC

### ** Examples

## Not run: 
##D # add fold-changes to the result table of the statistical comparison called "comparison1", after correcting for variable "age"
##D ... %>%
##D  mt_post_addFC(stat_name="comparison1", correct_confounder=age) %>% ...
##D  
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_post_addFC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_post_multTest")
### * mt_post_multTest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_post_multTest
### Title: Multiple testing correction
### Aliases: mt_post_multTest

### ** Examples

## Not run: 
##D # correct the statistical comparison called "Li's" using Benjamini-Hochberg
##D ... %>%
##D  mt_post_multTest(stat_name="Li's", method="BH") %>% ...
##D  
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_post_multTest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_post_pgain")
### * mt_post_pgain

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_post_pgain
### Title: Compute p-gain from metabolite ratio test
### Aliases: mt_post_pgain

### ** Examples

## Not run: 
##D # add p-gains to the result table of the statistical comparison called "comparison1"
##D ... %>%
##D  mt_post_pgain(stat_name="comparison1") %>% ...
##D  
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_post_pgain", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_batch_COMBAT")
### * mt_pre_batch_COMBAT

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_batch_COMBAT
### Title: COMBAT batch correction.
### Aliases: mt_pre_batch_COMBAT

### ** Examples

## Not run: ... %>% mt_pre_batch_COMBAT(batches="BATCH") %>% ...





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_batch_COMBAT", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_batch_median")
### * mt_pre_batch_median

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_batch_median
### Title: Median batch correction.
### Aliases: mt_pre_batch_median

### ** Examples

## Not run: ... %>% mt_pre_batch_median(batches="BATCH") %>% ...





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_batch_median", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_confounding_correction")
### * mt_pre_confounding_correction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_confounding_correction
### Title: Confounding correction
### Aliases: mt_pre_confounding_correction

### ** Examples


 ## Not run: 
##D # not run
##D  #... %>% mt_pre_confounding_correction( formula = ~batch + age, strata = "RUN_DAY" )
##D  
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_confounding_correction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_confounding_correction_stepwise_aic")
### * mt_pre_confounding_correction_stepwise_aic

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_confounding_correction_stepwise_aic
### Title: Confounding correction using stepwise AIC
### Aliases: mt_pre_confounding_correction_stepwise_aic

### ** Examples

 ## Not run: 
##D #... %>% mt_pre_confounding_correction_stepwise_aic(id_col = "RID", med_file = meds,
##D                          to_remove = c( "Med.Anticholinesterases", "Med.NMDAAntag"), n_cores = 10)
##D  
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_confounding_correction_stepwise_aic", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_filtermiss")
### * mt_pre_filtermiss

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_filtermiss
### Title: Filter by missingness.
### Aliases: mt_pre_filtermiss

### ** Examples

## Not run: 
##D # first remove samples with >10% missingness, then metabolites with >20% missingness
##D ... %>%
##D   mt_pre_filtermiss(sample_max=0.1) %>%
##D   mt_pre_filtermiss(met_max=0.2) %>%
##D ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_filtermiss", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_impute_knn")
### * mt_pre_impute_knn

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_impute_knn
### Title: Impute using kNN method.
### Aliases: mt_pre_impute_knn

### ** Examples

## Not run: 
##D # in the context of a SE pipeline
##D ... %>% mt_pre_impute_knn() %>% ...    # standard call
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_impute_knn", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_impute_knn_multicore")
### * mt_pre_impute_knn_multicore

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_impute_knn_multicore
### Title: Impute using kNN method with multi-core functionality.
### Aliases: mt_pre_impute_knn_multicore

### ** Examples

## Not run: 
##D # in the context of a SE pipeline
##D ... %>% mt_pre_impute_knn_multicore() %>% ...    # standard call
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_impute_knn_multicore", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_norm_external")
### * mt_pre_norm_external

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_norm_external
### Title: Normalize by an external sample annotation
### Aliases: mt_pre_norm_external

### ** Examples

## Not run: 
##D #' # in the context of a SE pipeline
##D ... %>% mt_pre_norm_quot(field='DNA') %>% ...    # normalize by values in field DNA
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_norm_external", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_norm_quot")
### * mt_pre_norm_quot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_norm_quot
### Title: Quotient normalization
### Aliases: mt_pre_norm_quot

### ** Examples

## Not run: 
##D #' # in the context of a SE pipeline
##D ... %>% mt_pre_norm_quot() %>% ...    # standard call
##D ... %>% mt_pre_norm_quot(ref_samples = GROUP=="ctrl") %>% ...    # use reference samples where 'GROUP' field in colData is 'ctrl'
##D #' ... %>% mt_pre_norm_quot(met_max = 0.2) %>% ...    # use only metabolites with <= 20% missing values to compute the reference used for normalization
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_norm_quot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_outlier")
### * mt_pre_outlier

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_outlier
### Title: Identifies sample outliers.
### Aliases: mt_pre_outlier

### ** Examples

## Not run: 
##D # first identify samples that have more than 50% univariate outliers, then identify multivariate outliers with a leverage >4m/n
##D ... %>%
##D   mt_pre_outlier(method="univariate", thresh=4, perc=0.5) %>%
##D   mt_pre_outlier(method="leverage", thresh=4) %>%
##D   mt_pre_outlier(method="mahalanobis", pval=0.01) %>%
##D ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_outlier", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_outliercorrection")
### * mt_pre_outliercorrection

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_outliercorrection
### Title: Identifies single outliers in samples.
### Aliases: mt_pre_outliercorrection

### ** Examples

## Not run: 
##D ... %>%
##D   mt_pre_outliercorrection(threshold=3, sample_num_correction=F) %>%
##D ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_outliercorrection", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_trans_exp")
### * mt_pre_trans_exp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_trans_exp
### Title: Exponentiate, base 2 by default.
### Aliases: mt_pre_trans_exp

### ** Examples

## Not run: 
##D # in the context of a SE pipeline
##D ... %>% mt_pre_trans_exp() %>% ...    # standard call, base 2
##D ... %>% mt_pre_trans_exp(base=10) %>% ...    # base 10
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_trans_exp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_trans_log")
### * mt_pre_trans_log

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_trans_log
### Title: Log, base 2 by default.
### Aliases: mt_pre_trans_log

### ** Examples

## Not run: 
##D # in the context of a SE pipeline
##D ... %>% mt_pre_trans_log() %>% ...    # standard call, base 2
##D ... %>% mt_pre_trans_log(base=10) %>% ...    # base 10
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_trans_log", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_trans_relative")
### * mt_pre_trans_relative

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_trans_relative
### Title: Scale each sample relative to the mean of a given group of
###   samples.
### Aliases: mt_pre_trans_relative

### ** Examples

## Not run: 
##D # normalize to control group, data not logged
##D ... %>% mt_pre_trans_relative(ref_samples = GROUP=="ctrl", is_logged=F) %>% ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_trans_relative", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_trans_scale")
### * mt_pre_trans_scale

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_trans_scale
### Title: Scale data, mean 0 / sd 1 by default.
### Aliases: mt_pre_trans_scale

### ** Examples

## Not run: 
##D # in the context of a SE pipeline
##D ... %>% mt_pre_trans_scale() %>% ...    # standard call, center and scale
##D ... %>% mt_pre_trans_scale(scale=F) %>% ...    # only mean centering
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_trans_scale", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_pre_zeroToNA")
### * mt_pre_zeroToNA

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_pre_zeroToNA
### Title: Set zeros in dataset to NA.
### Aliases: mt_pre_zeroToNA

### ** Examples

## Not run: 
##D # in the context of a SE pipeline
##D ... %>% mt_pre_zeroToNA() %>% ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_pre_zeroToNA", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_stats_multiv_net_GeneNet")
### * mt_stats_multiv_net_GeneNet

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_stats_multiv_net_GeneNet
### Title: Computes partial correlation matrix using the GeneNet estimator.
### Aliases: mt_stats_multiv_net_GeneNet

### ** Examples

## Not run: 
##D ... %>%
##D   mt_stats_multiv_net_GeneNet(stat_name ="pcor") %>%
##D ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_stats_multiv_net_GeneNet", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_stats_pathway_enrichment")
### * mt_stats_pathway_enrichment

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_stats_pathway_enrichment
### Title: Pathway enrichment script using statistical analysis results
###   from metabotools pipeline.
### Aliases: mt_stats_pathway_enrichment

### ** Examples

## Not run: 
##D %>% mt_stats_pathway_enrichment("kegg_db",
##D                                 grouping_var = "Group",
##D                                 control_grp_name = "Vehicle",
##D                                 case_grp_name = c("treatment2", "treatment1") %>%
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_stats_pathway_enrichment", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_stats_pathway_enrichment_with_differential_analysis")
### * mt_stats_pathway_enrichment_with_differential_analysis

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_stats_pathway_enrichment_with_differential_analysis
### Title: Pathway enrichment using different methods.
### Aliases: mt_stats_pathway_enrichment_with_differential_analysis

### ** Examples

## Not run: 
##D %>% mt_stats_pathway_enrichment("kegg_db",
##D                                 grp_col = "Group",
##D                                 ctrl_grp = "Vehicle",
##D                                 case_grp = c("treatment2", "treatment1") %>%
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_stats_pathway_enrichment_with_differential_analysis", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_stats_univ_lm")
### * mt_stats_univ_lm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_stats_univ_lm
### Title: Univariate GLMs.
### Aliases: mt_stats_univ_lm

### ** Examples

## Not run: 
##D # run lm with no confounders, "Group" as outcome
##D # filter to groups "Li_2" and "Li_5"
##D # name the comparison "Li's"
##D ... %>%
##D  mt_stats_univ_lm(
##D    formula      = ~ Group,
##D    sample_filter = (Group %in% c("Li_2","Li_5")),
##D    stat_name         = "Li's"
##D  ) %>% ...
##D  
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_stats_univ_lm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_stats_univ_missingness")
### * mt_stats_univ_missingness

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_stats_univ_missingness
### Title: Perform missingness significance analysis.
### Aliases: mt_stats_univ_missingness

### ** Examples

# run on sample field 'Group', name output stats object 'miss'
... %>% mt_stats_univ_missingness(comp_col = 'Group', stat_name='miss') %>% ...




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_stats_univ_missingness", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_stats_univ_tau")
### * mt_stats_univ_tau

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_stats_univ_tau
### Title: Computes Kendall's rank correlation. If present, NAs will be
###   omitted.
### Aliases: mt_stats_univ_tau

### ** Examples

## Not run: 
##D ... %>%
##D   mt_stats_univ_tau(var = "Stage", sample_filter = (GROUP %in% "Tumor"), name = "tau") %>%
##D ...
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_stats_univ_tau", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

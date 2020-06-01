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

%>% mt_stats_pathway_enrichment("kegg_db",
                                grouping_var = "Group",
                                control_grp_name = "Vehicle",
                                case_grp_name = c("treatment2", "treatment1") %>%





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

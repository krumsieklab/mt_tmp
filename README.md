# **Maplet**: **M**etabolomics **A**nalysis **P**ipe**L**inE **T**oolbox

Maplet is an R package for statistical data analysis with a special focus on metabolomics datasets. It allows users to create self-contained analytical pipelines. The toolbox builds upon the bioconductor package [SummarizedExperiment (SE)](https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html), which serves as a central repository for each pipeline’s data, analysis steps, and results. Maplet provides a suite of functions for interacting with this container including but not limited to data loading, annotation, statistical analysis, visualization, and reporting. Maplet is designed to work with the pipe operator (%>%) from the [magrittr](https://magrittr.tidyverse.org/) package. This operator allows for smooth connections between pipeline steps, without the need for temporary variables or multiple assignments. The combination of these elements allows for the creation of pipelines which are simple to follow, highly modular, and easily reproducible.

## Installation
The latest stable version of Maplet can be easily installed using the following command:
```{r}
devtools::install_github(repo="krumsieklab/maplet@v1.0.0", subdir="Maplet")
```

To install from the latest commit:
```{r}
devtools::install_github(repo="krumsieklab/maplet", subdir="Maplet")
```

**Note:** Maplet is in active development. Any commit without a release tag is not guaranteed to be stable. 

## Getting Started
Users should review the available examples in the [mt/examples](https://gitlab.com/krumsieklab/mt/-/tree/master/examples) folder. A large example pipeline demonstrating the general setup of a Maplet pipeline and most of the functions provided by Maplet is provided [here](https://gitlab.com/krumsieklab/mt/-/blob/master/examples/MT_example_pipeline.R). The mt/example folder also contains several stand-alone examples for more specialized functions not included in the main example.  
  
Users are also encouraged to review the Maplet Reference Guide sections 1 and 2.

## Want to Get Involved?
Anyone is welcome to contribute to the Maplet R package. For users wishing to contribute a new function, we recommend first reviewing the code and documentation of a few existing functions. We also recommend reviewing section 3 of the Maplet Reference Guide.  

**Note**: Maplet functions follow strict naming conventions (refer to section 3.2.1 of the Maplet Reference Guide for details). A submitted function may initially be rejected if it does not follow the function naming rules. We recommend reaching out to the maintainer (kpc4002@med.cornell.edu) if you have questions about naming and function development.

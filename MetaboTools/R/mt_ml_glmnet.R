#' GLMNET with outer and inner cross-validation
#'
#' Runs cv.glment (inner cross-validation) on an outer cross-validation loop. Then runs glmnet on all data.
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_name name under which this comparison will be stored, must be unique to all other statistical results
#' @param label_col name of column in colData(D) to be used as the response variable (y)
#' @param outer_k number of folds in the outer cross-validation loop
#' @param inner_k number of folds in the inner cross-validation loop
#' @param family response type; default: binomial
#' @param alpha elasticnet mixing parameter
#' @param lambda regularization parameter
#' @param rand_seed integer to set RNG state so results (including CV partitioning) can be exactly reproduced
#'
#' @return results$ouput: a list containing the foloowing:
#' \itemize{
#'  \item{table}{the coefficient matrix returned for glmnet trained on all data}
#'  \item{lstobj}{the glmnet model trained using on all data}
#' }
#' @return results$output2: a list containing the following:
#' \itemize{
#'  \item{idx}{list of indices used to create the train:test partitions for the outer cross-validation loop}
#'  \item{yhat}{list of all predictions from the outer cross-validation loop}
#'  \item{labels}{list of test labels from the outer cross-validation loop}
#'  \item{models}{list of all models from the outer cross-validation loop}
#'  \item{params}{list containing the original parameters}
#' }
#'
#' @examples
#' # example of how to run function
#' \dontrun{}
#'
#' @author KC
#'
#' @export

mt_ml_glmnet <- function(D,
                            stat_name,
                            label_col,
                            outer_k,
                            inner_k,
                            family="binomial",
                            alpha,
                            rand_seed){

  # validate arguments
  if(any(missing(stat_name), missing(label_col), missing(outer_k), missing(inner_k), missing(alpha), missing(rand_seed))){
    stop(paste0("Values must be provided for ALL of the following variables:\n",
                paste0(c("stat_name", "label_col", "outer_k", "inner_k", "alpha", "rand_seed"),collapse = "\n")))
  }
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(length(label_col)==1)

  if(label_col %in% colnames(colData(D))==F){
    stop("label_col is not a colData column name")
  }

  # get x and y
  x <- t(assay(D))
  y <- colData(D)[,label_col]

  # get outer cv folds
  set.seed(rand_seed)
  idxList <- caret::createFolds(factor(y), k=outer_k)

  # create training and testing dataset partitions
  sets <- lapply(1:outer_k, FUN=function(i){
    trainData <- x[unlist(idxList[-i]), ]
    trainLabel <- y[unlist(idxList[-i])]
    testData <- x[unlist(idxList[i]), ]
    testLabel <- y[unlist(idxList[i])]
    sets <- list(trainData, trainLabel, testData, testLabel)
    names(sets) <- c("trainData", "trainLabel", "testData", "testLabel")
    return(sets)
  })

  modList <- predList <- trueList <- list()

  # train glmnets using outer and inner cross-validation
  for(i in 1:outer_k){

    ## !! NEED UPDATE to allow user to select s (leave as iMod$lambda.min for now)
    set.seed(rand_seed)
    iMod <- glmnet::cv.glmnet(x=sets[[i]]$trainData, y=sets[[i]]$trainLabel, alpha=alpha, nfolds=inner_k, family=family)
    iPred <- predict(iMod, sets[[i]]$testData, s=iMod$lambda.min, type="response")

    modList[[i]] <- iMod
    predList[[i]] <- iPred
    trueList[[i]] <- sets[[i]]$testLabel

  }

  # train glmnet using all data
  mod <- glmnet::glmnet(x=x, y=y, alpha=alpha, family = family)
  coef_tab <- data.matrix(coef(mod, s=mod$lambda))

  # add status information
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = paste0('ran glmnet; outer k: ' , outer_k, ', inner k: ', inner_k),
      output = list(
        table = coef_tab,
        name = stat_name,
        lstobj = mod,
        outcome = NA
      ),
      output2 = list(
        idx = idxList,
        yhat = predList,
        labels = trueList,
        models = modList,
        params = list(
          x = x,
          y = y,
          outer_k = outer_k,
          inner_k = inner_k,
          alpha = alpha,
          rand_seed = rand_seed
        )
      )
    )

  # return
  D


}

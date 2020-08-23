#' GLMNET with outer and inner cross-validation
#'
#' Runs cv.glment (inner cross-validation) on an outer cross-validation loop. Then runs glmnet on all data.
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_name name under which this comparison will be stored, must be unique to all other statistical results
#' @param x input matrix
#' @param y response variable
#' @param outer_k number of folds in the outer cross-validation loop
#' @param inner_k number of folds in the inner cross-validation loop
#' @param alpha
#' @param lambda
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
#'
#' @export

mt_stats_glmnet <- function(D,
                            stat_name,
                            x,
                            y,
                            outer_k,
                            inner_k,
                            alpha,
                            lambda = NULL,
                            rand_seed){

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

  modList <- predList <- list()

  # train glmnets using outer and inner cross-validation
  for(i in 1:outer_k){

    set.seed(rand_seed)
    iMod <- glmnet::cv.glmnet(x=sets[[i]]$trainData, y=sets[[i]]$trainLabel, alpha=alpha, lambda=lambda, nfolds=inner_k, family="binomial")
    iPred <- predict(iMod, sets[[i]]$testData, s=iMod$lambda.min, type="response")

    modList[[i]] <- iMod
    predList[[i]] <- iPred

  }

  # train glmnet using all data
  mod <- glmnet::glmnet(x=x, y=y, alpha=alpha, lambda = lambda, family = "binomial")
  coef_tab <- data.matrix(coef(mod, s=mod$lambda))

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = 'glmnet',
      output = list(
        table = coef_tab,
        name = stat_name,
        lstobj = mod,
        outcome = NA
      ),
      output2 = list(
        idx = idxList,
        yhat = predList,
        models = modList,
        params = list(
          x = x,
          y = y,
          outer_k = outer_k,
          inner_k = inner_k,
          alpha = alpha,
          lambda = lambda,
          rand_seed = rand_seed
        )
      )
    )

  # return
  D


}

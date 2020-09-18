#' Plot Evaluation Measures for glmnet (binary) Predictions
#'
#' This function creates four plots: (1) ROC curve, (2) evaluation measures across thresholds (folds ignored), (3) AUC scores per fold,
#' (4) evaluation measures per fold (if cutoff provided, single values; else across thresholds).
#'
#' @param D \code{SummarizedExperiment} input
#' @param result_name name under which the glmnet results are stored, must be unique to all other glmnet results
#' @param plot_measures a vector of measures to include; can be one or more of the following: spec, sens, f1, acc, and/or ppv; default c("spec", "sens", "ppv")
#' @param cutoff threshold cutoff for prediction probability class assignment; function will select threshold closest to user provided value
#' @param pos_class which label to set as positive class
#'
#' @return results$ouput: list of four or more plots
#'
#' @examples
#' # example of how to run function
#' \dontrun{D %<>% mt_plots_glmnet_roc_auc(result_name = "glmnet result", plot_measures = c("spec", "sens") cutoff = 0.75)}
#'
#' @author KC
#'
#' @import ggplot2
#'
#' @export

mt_plots_glmnet_roc_auc <- function(D,
                                    result_name,
                                    plot_measures =  c("spec", "sens", "ppv", "f1", "acc"),
                                    cutoff,
                                    pos_class){

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  if(missing(result_name)){
    stop("A value must be provided for result_name")
  }

  if(missing(pos_class)){
    stop("pos_class is missing! The positive class must be specified.")
  }

  res_list <- MetaboTools:::mti_get_ml_res_by_name(D, result_name)

  pred_list <- res_list$output2$yhat
  true_list <- res_list$output2$labels

  names(pred_list) <- seq(1:length(pred_list))
  pred_df <- purrr::map_df(pred_list, ~as.data.frame(.x), .id="fold")
  names(true_list) <- seq(1:length(true_list))
  true_df <- purrr::map_df(true_list, ~as.data.frame(.x), .id="fold")
  all_folds_df <- cbind(pred_df, label=true_df$.x)
  colnames(all_folds_df)[2] <- "prediction"

  # check pos_class is present in true class labels
  if(pos_class %in% all_folds_df$label==F){
    stop("Provided positive class not found in true class labels.")
  }

  # convert true labels as boolean
  all_folds_df$label <- all_folds_df$label == pos_class

  # order predictions
  all_folds_df <- all_folds_df[order(all_folds_df$prediction),]

  # get roc plot data
  measures_list = MetaboTools:::mti_get_measures_list(trueclass=all_folds_df$label, predprob=all_folds_df$prediction)
  roc_data=data.frame(FPR=rev(1-measures_list$specvals),
                      sensitivity=rev(measures_list$sensvals))

  # get measures plot (folds ignored) data
  measures_data <- data.frame(sens=measures_list$sensvals,
                              spec=measures_list$specvals,
                              ppv=measures_list$ppvvals,
                              acc=measures_list$accvals,
                              f1=measures_list$f1vals,
                              thresholds=measures_list$thresholds)
  measures_data <- measures_data[,c(plot_measures, "thresholds")]
  measures_plot_data <- reshape2::melt(measures_data, id=c("thresholds"))

  # get evaluation measures per fold
  folds <- length(unique(all_folds_df$fold))
  fold_aucs <- c()
  fold_measures <- list()
  for(i in 1:folds){
    f_df <- all_folds_df[all_folds_df$fold==i,]
    fold_result <- MetaboTools:::mti_get_measures_list(f_df$label, f_df$prediction)

    # results for measures plot
    if(!missing(cutoff)){
      cut_val <- abs(cutoff - fold_result$thresholds)
      ct <- which(cut_val == min(cut_val))
      fms <- c(sens=fold_result$sensvals[ct],
               spec=fold_result$specvals[ct],
               ppv=fold_result$ppvvals[ct],
               acc=fold_result$accvals[ct],
               f1=fold_result$f1vals[ct],
               fold=i,
               cutval=fold_result$thresholds[ct],
               cutidx = ct)

    }else{
      fms <- data.frame(sens=fold_result$sensvals,
                        spec=fold_result$specvals,
                        ppv=fold_result$ppvvals,
                        acc=fold_result$accvals,
                        f1=fold_result$f1vals,
                        thresholds=fold_result$thresholds,
                        fold = i)
    }
    fold_measures[[i]] <- fms

    # results for AUC plot
    fold_aucs <- c(fold_aucs, fold_result$AUC)

  }

  if(!missing(cutoff)){
    fold_measures_df <- do.call(rbind, fold_measures)
    fold_measures_df <- as.data.frame(fold_measures_df[,c(plot_measures, "fold")])
    fold_measures_plot_df <- reshape2::melt(fold_measures_df, id=c("fold"))
  }

  fold_auc_df <- data.frame(folds=1:folds, AUC=fold_aucs)


  plot_list <- list()
  plot_idx <- 1
  # (1) ROC Plot
  roc_plot <- ggplot(data=roc_data, aes(x=FPR, y=sensitivity)) +
    geom_line() + geom_abline(intercept=0, slope=1, colour="blue") +
    ggtitle(sprintf("AUC: %.3f (All Folds)", measures_list$AUC)) + xlab("FPR") + ylab("TPR") +
    xlim(0,1) + ylim(0,1)
  plot_list[[plot_idx]] <- roc_plot
  plot_idx <- plot_idx + 1

  # (2) Evaluation Measures Line Plot (fold ignored)
  # threshold values plot
  measures_plot <- ggplot(measures_plot_data, aes(x=thresholds, y=value, color=variable)) + geom_line() +
    xlab("Thresholds (values)") +
    ylab("Evaluation Measure Values") +
    labs(color='Measures')
  if(!missing(cutoff)){
    measures_plot <- measures_plot + geom_vline(xintercept = cutoff)
    mp_title <- paste0("Evaluaiton Measures (All Folds), cutoff: ", round(cutoff,2))
  }else{
    mp_title <- "Evaluation Measures (All Folds)"
  }
  measures_plot <- measures_plot + ggtitle(mp_title)
  plot_list[[plot_idx]] <- measures_plot
  plot_idx <- plot_idx + 1

  # (3) AUC Per Fold Plot
  auc_plot <- ggplot(fold_auc_df, aes(x = factor(folds), y = AUC)) + geom_point(size=3) + ylim(0,1) +
    xlab("Fold") +
    ggtitle("AUC (Per Fold)")
  plot_list[[plot_idx]] <- auc_plot
  plot_idx <- plot_idx + 1

  # (4) Evaluation Measures Per Fold Plot(s)
  if(!missing(cutoff)){
    # plot evaluation measures at specified cutoff
    measure_pf_plot <- ggplot(fold_measures_plot_df, aes(x= factor(fold), y=value, color=variable)) + geom_point(size=3) +
      xlab("Fold") +
      ylab("Evaluation Measure Value") +
      ggtitle(paste0("Evaluation Measures (Per Fold), cutoff: ", round(cutoff,2)))
    plot_list[[plot_idx]] <- measure_pf_plot
  }else{
    for(i in 1:length(fold_measures)){
      fms <- fold_measures[[i]]
      fold <- unique(fms$fold)
      fms <- fms[,c(plot_measures, "thresholds")]
      fms <- reshape2::melt(fms, id=c("thresholds"))

      # thresholds values plot
      measures_pf_plot <- ggplot(fms, aes(x=thresholds, y=value, color=variable)) + geom_line() +
        xlab("Thresholds") +
        ylab("Evaluation Measure Values") +
        labs(color='Measures') +
        ggtitle(paste0("Evaluation Measures: Fold ", fold))
      plot_list[[plot_idx]] <- measures_pf_plot
      plot_idx <- plot_idx + 1
    }
  }

  ## add status information & plot
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = "glmnet evaluation measure plots",
      output = plot_list
    )
  ## return
  D


}

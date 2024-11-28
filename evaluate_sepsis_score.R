library(data.table)

evaluate_sepsis_score <-function(data, want_auc=F){
  # This helps with performance! (about factor of 2 quicker)
  data <- as.data.table(data)
  
  # Set parameters.
  dt_early <- -8
  dt_optimal <- -4
  dt_late <- 3
  max_u_tp <- 1
  min_u_fn <- -2
  u_fp <- -0.05
  u_tn <- 0
  
  # % Load labels and predictions.
  num_files <- length(unique(data$patient))
  labels <- data$SepsisLabel
  predictions <- data$PredictedLabel
  probabilities <- data$PredictedProbability
  
  num_records <- nrow(data)
  if (!all(labels==0 | labels==1)){
    print('Labels must satisfy label == 0 or label == 1.')    
    return(NA)
  }
  
  if (!all(predictions==0 | predictions==1)){
    print('Predictions must satisfy prediction == 0 or prediction == 1.');
    return(NA)
  }

  if (!all(probabilities>=0 | probabilities<=1)){
    print('Probabilities do not satisfy 0 <= probability <= 1.');
    return(NA)
  }

  min_probability_positive <- min(probabilities[predictions == 1])
  max_probability_negative <- max(probabilities[predictions == 0])
  if (min_probability_positive <= max_probability_negative){
    print('Predictions are inconsistent with probabilities, i.e., a positive prediction has a lower (or equal) probability than a negative prediction.');
  }

  # Compute AUC, accuracy, and F-measure.
  if (want_auc==T){
    auc <- compute_auc(labels, probabilities)
  }else{
    auc <- c(0,0)
  }
  auroc <- auc[1]
  auprc <- auc[2]
  accuracy_f <- compute_accuracy_f_measure(labels, predictions)
  accuracy <- accuracy_f [1]
  f_measure <- accuracy_f [2]

  # Compute utility.
  observed_utilities <- numeric(num_files)
  best_utilities <- numeric(num_files)
  worst_utilities <- numeric(num_files)
  inaction_utilities <- numeric(num_files)

  patients<-unique(data$patient)
  for (k in 1:length(patients)){
    pdata <- data[data$patient==patients[k],]
    p_labels <- pdata$SepsisLabel
    num_records <- length(p_labels)
    observed_predictions <- pdata$PredictedLabel
    best_predictions <- numeric(num_records)
    worst_predictions <- numeric(num_records)
    inaction_predictions <- numeric(num_records);
    
    if (sum(p_labels)>0){
      t_sepsis <- which(p_labels == 1)[1] - dt_optimal
      best_predictions[max(1, t_sepsis + dt_early) : min(t_sepsis + dt_late, num_records)] <- 1
    }else{
      best_predictions[] <- 0
    }
    worst_predictions <- (1 - best_predictions)

    observed_utilities[k] <- compute_prediction_utility(p_labels, observed_predictions, dt_early, dt_optimal, dt_late, max_u_tp, min_u_fn, u_fp, u_tn)
    best_utilities[k] <- compute_prediction_utility(p_labels, best_predictions, dt_early, dt_optimal, dt_late, max_u_tp, min_u_fn, u_fp, u_tn)
    worst_utilities[k] <- compute_prediction_utility(p_labels, worst_predictions, dt_early, dt_optimal, dt_late, max_u_tp, min_u_fn, u_fp, u_tn)
    inaction_utilities[k] <- compute_prediction_utility(p_labels, inaction_predictions, dt_early, dt_optimal, dt_late, max_u_tp, min_u_fn, u_fp, u_tn)
  }
  
  unnormalized_observed_utility <- sum(observed_utilities);
  unnormalized_best_utility <- sum(best_utilities);
  unnormalized_worst_utility <- sum(worst_utilities);
  unnormalized_inaction_utility <- sum(inaction_utilities);

  normalized_observed_utility <- (unnormalized_observed_utility - unnormalized_inaction_utility) / 
                                 (unnormalized_best_utility - unnormalized_inaction_utility);
  res <- c(AUROC=auroc,AUPRC=auprc,Accuracy=accuracy,Fmeasure=f_measure,Utility=normalized_observed_utility)
  #print(res)
  return(res)
}

# The compute_auc function computes AUROC and AUPRC as well as other summary
# statistics (TP, FP, FN, TN, TPR, TNR, PPV, NPV, etc.) that can be exposed
# from this function.
compute_auc <- function(labels, predictions){
  # Check inputs for errors.
  if (length(predictions) != length(labels)){
    print('Numbers of predictions and labels must be the same.')
    return(NA)
  }

  n <- length(labels)
  if (!all(labels==0 | labels==1)){
    print('Labels must satisfy label == 0 or label == 1.')    
    return(NA)
  }

  if (!all(predictions>=0 && predictions<=1)){
    print('Predictions do not satisfy 0 <= prediction <= 1.')
    return(NA)
  }

  # Find prediction thresholds.
  thresholds <- rev(sort(unique(predictions)));

  if (thresholds[1] != 1){
    thresholds <- c(1, thresholds)
  }

  if (thresholds[length(thresholds)] != 0){
    thresholds <- c(thresholds, 0);
  }

  m <- length(thresholds);

  # Populate contingency table across prediction thresholds.
  tp <- numeric(m);
  fp <- numeric(m);
  fn <- numeric(m);
  tn <- numeric(m);
  
  # Find indices that sort predicted probabilities from largest to smallest.
  idx <- order(predictions, decreasing=T)
  
  i <- 1;
  for (j in 1:m){
    # Initialize contingency table for j-th prediction threshold.
    if (j == 1){
      tp[j] = 0;
      fp[j] = 0;
      fn[j] = sum(labels);
      tn[j] = n - fn[j];
    }else{
      tp[j] = tp[j-1];
      fp[j] = fp[j-1];
      fn[j] = fn[j-1];
      tn[j] = tn[j-1];
    }

    # Update contingency table for i-th largest prediction probability.
    while (i <= n && predictions[idx[i]] >= thresholds[j]){
      if (labels[idx[i]] == 1){
        tp[j] <- tp[j] + 1;
        fn[j] <- fn[j] - 1;
      }else{
        fp[j] <- fp[j] + 1;
        tn[j] <- tn[j] - 1;
      }
      i <- i + 1;
    }
  }
  
  # Summarize contingency table.
  tpr <- numeric(m)
  tnr <- numeric(m)
  ppv <- numeric(m)
  npv <- numeric(m)

  for (j in 1:m){
    if ((tp[j] + fn[j]) > 0){
      tpr[j] <- tp[j] / (tp[j] + fn[j]);
    }else{
      tpr[j] <- 1;
    }
    
    if ((fp[j] + tn[j]) > 0){
      tnr[j] <- tn[j] / (fp[j] + tn[j]);
    }else{
      tnr[j] <- 1;
    }
    
    if ((tp[j] + fp[j]) > 0){
      ppv[j] <- tp[j] / (tp[j] + fp[j]);
    }else{
      ppv[j] <- 1;
    }
    
    if ((fn[j] + tn[j]) > 0){
      npv[j] <- tn[j] / (fn[j] + tn[j]);
    }else{
      npv[j] <- 1;
    }
  }

  # Compute AUROC as the area under a piecewise linear function of TPR /
  # sensitivity (x-axis) and TNR / specificity (y-axis) and AUPRC as the area
  # under a piecewise constant of TPR / recall (x-axis) and PPV / precision
  # (y-axis).
  auroc <- 0;
  auprc <- 0;
  for (j in 1:(m-1)){
    auroc <- auroc + 0.5 * (tpr[j + 1] - tpr[j]) * (tnr[j + 1] + tnr[j]);
    auprc <- auprc + (tpr[j + 1] - tpr[j]) * ppv[j + 1];
  }
  return(c(auroc, auprc))
}
#labels <- c(0,0,0,0,1,1)
#predictions <- c(0.3,0.4,0.6,0.7,0.8,0.8)
#auroc <- 1
#auprc <- 1
#compute_auc(labels, predictions)

compute_accuracy_f_measure <- function (labels, predictions){
  # Check inputs for errors.
  if (length(predictions) != length(labels)){
    print('Numbers of predictions and labels must be the same.')
    return(NA)
  }
  
  n <- length(labels)
  if (!all(labels==0 | labels==1)){
    print('Labels must satisfy label == 0 or label == 1.')    
    return(NA)
  }
  
  if (!all(predictions==0 | predictions==1)){
    print('Predictions must satisfy prediction == 0 or prediction == 1.');
    return(NA)
  }

  #Populate contingency table.
  tp <- sum(labels == 1 & predictions == 1); 
  fp <- sum(labels == 0 & predictions == 1);
  fn <- sum(labels == 1 & predictions == 0); 
  tn <- sum(labels == 0 & predictions == 0);

  # Summarize contingency table.
  if ((tp + fp + fn + tn) > 0){
    accuracy <- (tp + tn) / (tp + fp + fn + tn);
  }else{
    accuracy <- 1;
  }

  if ((2 * tp + fp + fn) > 0){
    f_measure <- 2 * tp / (2 * tp + fp + fn);
  }else{
    f_measure <- 1;
  }
  return(c(accuracy, f_measure))
}
#labels <- c(0,0,0,0,1,1)
#predictions <- c(0,0,1,1,1,1)
#accuracy <- 0.66667
#f_measure <- 0.66667
#compute_accuracy_f_measure(labels, predictions)

compute_prediction_utility <- function(labels, predictions, dt_early=-8, dt_optimal=-4,
                                       dt_late=3, max_u_tp=1, min_u_fn=-2, u_fp=-0.05, u_tn=0){
  # Check inputs for errors.
  if (length(predictions) != length(labels)){
    print('Numbers of predictions and labels must be the same.')
    return(NA)
  }

  n <- length(labels)
  if (!all(labels==0 | labels==1)){
    print('Labels must satisfy label == 0 or label == 1.')    
    return(NA)
  }
  
  if (!all(predictions==0 | predictions==1)){
    print('Predictions must satisfy prediction == 0 or prediction == 1.');
    return(NA)
  }

  if (dt_early >= dt_optimal){
    print('The earliest beneficial time for predictions must be before the optimal time.')
    return(NA)
  }

  if (dt_optimal >= dt_late){
    error('The optimal time for predictions must be before the latest beneficial time.')
    return(NA)
  }
  # Does the patient eventually have sepsis?
  if (sum(labels)>0){
    is_septic <- T
    t_sepsis <- which(labels == 1)[1] - dt_optimal;
  }else{
    is_septic <- F;
    t_sepsis <- Inf;
  }

  # Define slopes and intercept points for affine utility functions of the
  # form u = m * t + b.
  m_1 <- max_u_tp / (dt_optimal - dt_early);
  b_1 <- -m_1 * dt_early;
  m_2 <- -max_u_tp / (dt_late - dt_optimal);
  b_2 <- -m_2 * dt_late;
  m_3 <- min_u_fn / (dt_late - dt_optimal);
  b_3 <- -m_3 * dt_optimal;

  # Compare predicted and true conditions.
  u <- integer(n)

  for (t in 1:n){
    if (t <= (t_sepsis + dt_late)){
    # TP
      if (is_septic && predictions[t]){
        if (t <= (t_sepsis + dt_optimal)){
          u[t] <- max(m_1 * (t - t_sepsis) + b_1, u_fp)
        }else if (t <= (t_sepsis + dt_late)){
          u[t] <- m_2 * (t - t_sepsis) + b_2
        }
      # FN
      }else if (is_septic && !predictions[t]){
        if (t <= (t_sepsis + dt_optimal)){
          u[t] <- 0
        }else if (t <= (t_sepsis + dt_late)){
          u[t] <- m_3 * (t - t_sepsis) + b_3
        }
      # FP
      }else if (!is_septic && predictions[t]){
        u[t] <- u_fp
      # TN
      }else if (!is_septic && !predictions[t]){
        u[t] <- u_tn
      }
    }
  }

  # Find total utility for patient.
  utility <- sum(u)
  return(utility)
}
#labels <- c(0,0,0,0,1,1)
#predictions <- c(0,0,1,1,1,1)
#utility <- 3.3889
#utility <- compute_prediction_utility(labels, predictions)
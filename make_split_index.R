set.seed(216)

cv_index <- sample(ggplot2::cut_number(1:253,5)|>as.numeric())
save(cv_index, file = 'cv_split.rda')
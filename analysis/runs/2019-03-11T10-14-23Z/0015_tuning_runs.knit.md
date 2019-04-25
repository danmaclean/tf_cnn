---
title: "tuned_model.Rmd"
author: "Dan MacLean"
date: "11/03/2019"
output: html_document
---





```r
arab_data <- readRDS("../data/scaled_arab_TF_one_channel.RDS")
new_order <- sample(1:8704)

shuffled_data <- list(
  x = arab_data$x[new_order,,,, drop = F],
  y = as.numeric(arab_data$y[new_order])
)


train_i <- 1:6900
val_i <- 6901:7802
test_i <- 7803:8704

x_train <- shuffled_data$x[train_i,,, ,drop = F]
x_val <- shuffled_data$x[val_i,,, ,drop = F]
x_test <- shuffled_data$x[test_i,,, ,drop = F]

y_train <- shuffled_data$y[train_i]
y_val <- shuffled_data$y[val_i]
y_test <- shuffled_data$y[test_i]
```



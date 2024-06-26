library(dplyr)
library(pbapply)
library(scales)
library(latex2exp)
library(randomForest)
library(glmnet)
library(e1071)
library(gbm3)

source("pair_plot.R")

islrblue <- "#3333b2"
islrred <- "#990000"
islrgreen <- "#009900"

dube <- read.csv("proteomics_dube_2023.csv", stringsAsFactors = TRUE) |>
  select(-SubjectID)
n <- nrow(dube)
p <- ncol(dube)

#########################
# Quantitative variable #
#########################

set.seed(42)

# Train / test split
prop_test <- 0.75
s <- sample(1:n, 0.75 * n)

train <- dube[s, ]
test <- dube[setdiff(1:n, s), ]

# Pick a quantitative variable for modelling

# LASSO regression

x_train <- model.matrix(EHBP1 ~ ., data = train)[, -1]
x_test <- model.matrix(EHBP1 ~ ., data = test)[, -1]
y_train <- train[, "EHBP1"]
y_test <- test[, "EHBP1"]

cv_lasso <- cv.glmnet(x_train, y_train,
  standardise = TRUE, nfolds = 30,
  thresh = 1e-12
)
lambda_opt <- cv_lasso$lambda.min

pdf("lasso_fitting.pdf", width = 12)
par(mfrow = c(1, 2))
plot(
  log(cv_lasso$lambda),
  cv_lasso$cvm,
  type = "l",
  xlab = expression(log(lambda)),
  ylab = "LOOCV error",
  ylim = c(0, 1)
)
lines(log(cv_lasso$lambda), cv_lasso$cvup, lty = 2)
lines(log(cv_lasso$lambda), cv_lasso$cvlo, lty = 2)
abline(v = log(lambda_opt), lty = 2)
plot(cv_lasso$glmnet.fit, "lambda",
  xlab = expression(lambda),
  ylab = "Standardised Coefficients"
)
abline(v = log(lambda_opt), lty = 2)
dev.off()

coefs_lasso <- predict(cv_lasso, s = lambda_opt, type = "coefficients")
coefs_lasso_nz <- coefs_lasso[coefs_lasso[, 1] != 0, ]

coefs_lasso_nz |>
  sort(decreasing = TRUE) |>
  barplot(names.arg = names(coefs_lasso_nz), cex.names = 0.75, col = NA)

lasso_pred <- predict(cv_lasso, s = lambda_opt, newx = x_test)

mean((lasso_pred - y_test)^2)

# Pair plot of relevant variables
names_lasso <- c("EHBP1", names(coefs_lasso_nz)[-1])

pdf("lasso_nz_cov.pdf")
pair_plot(dube[, names_lasso], cex = 0.2)
dev.off()

# Random forest

tuneRF(x_train, y_train, ntreeTry = 500, trace = FALSE, plot = FALSE)

mtry <- c(seq(10, p - 1, by = 100)) |> sort()

rf_fits <- pblapply(mtry, function(mt) {
  rf_fit <- randomForest(
    EHBP1 ~ .,
    data = train,
    xtest = x_test |> data.frame(),
    ytest = y_test,
    mtry = mt,
    ntree = 500,
    importance = TRUE
  )
})

rf_perf <- pblapply(rf_fits, function(fit) {
  list(mse = fit$test$mse, balanced = fit$test$mse[500])
})

best_model <- sapply(rf_perf, function(i) i$balanced) |> which.min()

plot(rf_perf[[best_model]]$mse,
  lwd = 3.5,
  col = islrred,
  type = "l",
  ylim = range(rf_perf),
  xlab = "No. trees",
  ylab = "Test set MSE"
)

for (i in 1:length(rf_perf)) {
  if (i != best_model) {
    lines(
      rf_perf[[i]]$mse,
      col = ifelse(best_model == i, "red", alpha("black", 0.25)),
      lwd = ifelse(best_model == i, 3.0, 1.0),
      type = "l"
    )
  }
}

varImpPlot(rf_fits[[best_model]], main = NA)

inc_pur <- importance(rf_fits[[best_model]]) |>
  data.frame() |>
  arrange(IncNodePurity) |>
  tail(30) |>
  rownames()

inc_mse <- importance(rf_fits[[best_model]]) |>
  data.frame() |>
  arrange(X.IncMSE) |>
  tail(30) |>
  rownames()

pair_plot(dube[, c("EHBP1", inc_mse)])

#########################
# Categorical variable  #
#########################

pr_out <- prcomp(dube[, 3:p], scale = TRUE)

pr_var <- pr_out$sdev^2
pve <- pr_var / sum(pr_var)

par(mfrow = c(1, 2))
plot(pve, xlab = "No. PC", ylab = "PVE", type = "b")
abline(v = 12, lty = 2, col = "red")
plot(cumsum(pve), xlab = "No. PC", ylab = "Cum. PCE", type = "b")
abline(v = 12, lty = 2, col = "red")

s <- sample(1:n, 0.5 * n, replace = FALSE)

dube_pr <- cbind(dube[, 1:2], pr_out$x[, 1:12])
train_pr <- dube_pr[s, ]
test_pr <- dube_pr[setdiff(1:n, s), ]

# SVM with radial kernel

set.seed(42)

svm_fit_pr_lin <- tune(
  svm,
  temp ~ .,
  data = train_pr,
  kernel = "linear",
  ranges = list(
    cost = 10^seq(-3, 3, length.out = 20)
  )
)

svm_lin_opt <- summary(svm_fit_pr_lin)$best.model
svm_lin_pred <- predict(svm_lin_opt, test_pr)
table(svm_lin_pred, test_pr$temp)

svm_fit_pr_poly <- tune(
  svm,
  temp ~ .,
  data = train_pr,
  kernel = "polynomial",
  ranges = list(
    cost = 10^seq(-3, 3, length.out = 20),
    degree = 1:8
  )
)

svm_poly_opt <- summary(svm_fit_pr_poly)$best.model
svm_poly_pred <- predict(svm_poly_opt, test_pr)
table(svm_poly_pred, test_pr$temp)

svm_fit_pr_radial <- tune(
  svm,
  temp ~ .,
  data = train_pr,
  kernel = "radial",
  ranges = list(
    gamma = 10^seq(-3, 2, length.out = 20),
    cost = 10^seq(-3, 3, length.out = 20)
  ),
  validation.x = test_pr
)

svm_radial_opt <- summary(svm_fit_pr_radial)$best.model
svm_radial_pred <- predict(svm_radial_opt, test_pr)
table(svm_radial_pred, test_pr$temp)

# KNN
library(class)

knn_tune <- tune.knn(
  train_pr[, -c(1, 2)],
  train_pr[, 1],
  k = 1:16,
  tunecontrol = tune.control(
    nrepeat = 10,
    cross = 5
  )
)

plot(knn_tune$performances$k, knn_tune$performances$error,
  type = "b",
  xlab = "K", ylab = "CV error"
)

knn_opt <- knn_tune$best.model
table(knn_opt$cl, test_pr$temp)

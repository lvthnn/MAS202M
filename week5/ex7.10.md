---
title: "Exercise 7.10"
author: "Kári Hlynsson"
date: "2024-02-13"
output: html_document
---

This question relates to the `College` data set.

## Part (a)
Split the data into a training and test set. Using out-of-state tuition as the response and the other variables
as the predictors, perform forward stepwise selection on the training set in order to identify a satisfactory model
that uses just a subset of the predictors.


```r
library(ISLR2)
library(leaps)
data(College)
attach(College)

set.seed(1)

n <- nrow(College)
train <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.75, 0.25))

reg.fit <- regsubsets(Outstate ~ ., data = College[train, ], nvmax = 17, method = "forward")
reg.summary <- summary(reg.fit)
```

Take a look at some statistics:


```r
par(mfrow = c(1, 3))

plot(reg.summary$bic, type = "b", xlab = "No. predictors", ylab = "BIC")
plot(reg.summary$cp, type = "b", xlab = "No. predictors", ylab = expression(C[p]))
plot(reg.summary$adjr2, type = "b", xlab = "No. predictors", ylab = expression("adj R"^2))
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

We will use the model with 6 predictors since it is within 0.2 standard deviations of each of the optimum statistics.
This is also done just to simplify the model (I don't want to deal with > 10 predictors).


```r
coef.opt <- coef(reg.fit, id = 6)
vars <- names(coef.opt)

coef.opt
```

```
##   (Intercept)    PrivateYes    Room.Board           PhD   perc.alumni 
## -3116.5965010  2680.7956450     0.9706120    31.7946951    55.5388924 
##        Expend     Grad.Rate 
##     0.2174337    25.8072183
```

## Part (b)
Fit a GAM on the training data, using out-of-state tuition as the response and the features selected in the previous step
as the predictors. Plot the results, and explain your findings.


```r
library(gam)
```

```
## Error in library(gam): there is no package called 'gam'
```

```r
gam.fit <- gam(Outstate ~ Private + s(Room.Board, 2) + s(PhD, 2) + s(perc.alumni, 2) + s(Grad.Rate, 4) + s(Expend, 5), data = College[train, ])
```

```
## Error in gam(Outstate ~ Private + s(Room.Board, 2) + s(PhD, 2) + s(perc.alumni, : could not find function "gam"
```

```r
par(mfrow = c(2, 3))
plot(gam.fit)
```

```
## Error in eval(expr, envir, enclos): object 'gam.fit' not found
```

We see that `Room.Board`, `PhD` and `perc.alumni` all seem to be linear at first glance.

## Part (c)
Evaluate the model obtained on the test set, and explain the results obtained.


```r
gam.pred <- predict(gam.fit, College[!train, ])
```

```
## Error in eval(expr, envir, enclos): object 'gam.fit' not found
```

```r
mse <- mean((gam.pred - College[!train, ]$Outstate)^2)
```

```
## Error in eval(expr, envir, enclos): object 'gam.pred' not found
```

```r
mse
```

```
## Error in eval(expr, envir, enclos): object 'mse' not found
```
The MSE is 3818022. Let's look at $R^2$:


```r
tss <- mean((College[!train, ]$Outstate - mean(College[!train, ]$Outstate))^2)
1 - mse / tss
```

```
## Error in eval(expr, envir, enclos): object 'mse' not found
```

The $R^2$ is approximately 0.79, so the GAM model explains approximately 79% of the variation in the data.

## Part (d)
For which variables, if any, is there evidence of a non-linear relationship with the response?


```r
summary(gam.fit)
```

```
## Error in eval(expr, envir, enclos): object 'gam.fit' not found
```

Strong relationship between response and quintic `Expend` spline, and a significant non-linear relationship between `Grad.Rate` and the response.



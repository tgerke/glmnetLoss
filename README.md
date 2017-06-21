# glmnetLoss
a hack of the glmnet package that allows the user to customize the objective function in cross validation.

Current usage simply replaces the cv.lognet function (hence, only works with logistic regression), but would generalize easily to the other cv... functions. 
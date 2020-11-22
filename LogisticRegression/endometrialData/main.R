
#' Example: estimation of binomial model with extreme data

library(brglm2nc)

ML <- glm(HG ~ NV + PI + EH, data = endometrial, family = binomial("probit"))
BR_mean <- update(ML, method = "brglmFit", type = "AS_mean")
BC <- update(ML, method = "brglmFit", type = "correction")
BR_median <- update(ML, method = "brglmFit", type = "AS_median")
BR_nc <- update(ML, method = "brglmFit", type = "AS_nc")

summary(ML)
summary(BC)
summary(BR_mean)
summary(BR_median)
summary(BR_nc)

require(stargazer)
require(dplyr)

stargazer(PE=ML, BR=BR_mean, MU=BR_median, NC=BR_nc, column.labels = c("PE", "BR", "MU", "NC"), title="Logistic regression on endometrial data", align=TRUE, label="tab:endometrial", out = "output/endometrial.tex")

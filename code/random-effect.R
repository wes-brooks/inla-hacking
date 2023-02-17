library(ggplot2)
library(INLA)

data(pulp, package="faraway")
summary(pulp)


imod <- inla(bright ~ f(operator, model="iid"),
             family="gaussian",
             data=pulp)
summary(result)



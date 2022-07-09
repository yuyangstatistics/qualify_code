# References: https://towardsdatascience.com/insurance-risk-pricing-tweedie-approach-1d71207268fc

library(insuranceData)
data(dataCar)

library(tweedie)
library(statmod)
est_p <-tweedie.profile(claimcst0 ~ veh_value+veh_body+veh_age+ 
                          gender+area+agecat,data=dataCar,
                        link.power = 0,do.smooth = TRUE, 
                        do.plot = TRUE)
tweedie_model <- glm(claimcst0 ~ veh_value+veh_body+veh_age+
                       gender+ area+agecat,data=dataCar, 
                     family = tweedie(var.power=1.6, link.power=0),
                     offset=log(exposure))
summary(tweedie_model)



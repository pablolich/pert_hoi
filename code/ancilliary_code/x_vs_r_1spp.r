library(tidyverse)

safe_sqrt <- function(x){
  result = c()
  for (i in x){
    if (i<0){
      res_i = 0
    }
    else{
      res_i = sqrt(i)
    }
    result = c(result, res_i)
  }
  return(result)
}

get_x <- function(alpha, r){
  xplus = 1/(2*alpha)*(alpha - 1 + safe_sqrt((1-alpha)^2 + 4*alpha*r))
  xminus = 1/(2*alpha)*(alpha - 1 - safe_sqrt((1-alpha)^2 + 4*alpha*r))
  return(data.frame(alpha = alpha, r = r, xplus = xplus, xminus = xminus))
}

alpha = seq(0.001, 1, length.out = 5)
r = seq(-1, 2, length.out = 1000)
data = get_x(alpha, r)

lapply(data, is.nan)

to_plot = data %>% pivot_longer(c(xplus, xminus), 
                                names_to = "type",
                                values_to = "sol") %>% 
  filter(sol > -100)

ggplot(to_plot, aes(x = r, y = sol))+
  geom_line(aes(group = interaction(as.factor(alpha), as.factor(type)),
                col = as.factor(alpha)))+
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")

            
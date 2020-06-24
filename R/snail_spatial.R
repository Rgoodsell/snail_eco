# Analysis of snail ecological data
library(tidyverse)
library(Rtsne)
library(cluster)
library(logisticPCA)
library(DHARMa)
library(ggeffects)
library(sp)
library(glmmfields)
library(gstat)




# Is there within site spatial autocorrelation?
# Quick & dirty way to look is to split data by location

# Functions --------------------------------- 
#  Fit mods & test spatacf 

fitMods <- function(modDat, modForm){
   modDat %>% 
    glm(modForm , data=. , family = binomial) %>% 
    simulateResiduals() %>% 
    testSpatialAutocorrelation(.,modDat$x , modDat$y)
}

# Function to plot site specific variograms 
plotVariograms <- function(modDat,modForm){
  mod <- modDat %>% glm(modForm , data=. , family = binomial) 
  dat <- tibble(lon = modDat$x,lat = modDat$y,resids=resid(mod))
        coordinates(dat)<-c('lon','lat')
        variogram(resids~1,data=dat,alpha=c(0,45,90,135)) %>% plot()
}


modDat <- readRDS("snailModDat.rds") %>% ungroup() %>% drop_na()
str(modDat)

table(modDat$year,modDat$location)

modDat_list <- modDat %>% rowwise() %>% 
  mutate(x = x + rnorm(1 , 0 , 0.00001), # add teeny bit of random noise to xy coords (dharma doesn't like identical spatcoords)
         y = y + rnorm(1 , 0 , 0.00001)) %>% 
  split(.$location)


#  Refit model (without z*location interaction & year terms) to each location
modForm <- formula(species ~ month + y  + z   + aspect_norm + Incline + LPC1 + LPC2)

# Only spatial autocor in England... 
modDat_list %>% 
      map(~fitMods(.,modForm))

  

# Check variograms ----- 
# Suggests spatial patterns in France.. 
modDat_list %>%
  map(~plotVariograms(.,modForm)) 


# Fit spatially explicit models
spatDat <- modDat  %>% mutate(species = as.numeric(species)-1) %>% drop_na()

spatMod <- glmmfields(species ~ month  + z * location + year + aspect_norm + Incline + LPC1+LPC2, 
                      data = spatDat,
                      family = binomial(link="logit"),
                      lat = "y", lon = "x",iter = 500, chains = 1,
                      seed = 1 
)

plot(spatMod$model)
summary(spatMod$model)

plot(spatMod, type = "spatial-residual", link = FALSE) + geom_point(size = 3)
plot(spatMod, type = "residual-vs-fitted" , link = FALSE)
plot(spatMod, type = "prediction", link = FALSE) +
  viridis::scale_colour_viridis() +
  geom_point(size = 3)



pred_grid <- expand.grid(x = seq(min(spatDat$x), max(spatDat$x), length.out = 30) ,
                         y = seq(min(spatDat$y), max(spatDat$y), length.out = 30) , 
                         z = 1.9, 
                         month = levels(spatDat$month),
                         location = c("E","F","W"), 
                         LPC1 = 1 ,LPC2 = 1 , year = levels(spatDat$year) ,
                         aspect_norm = mean(spatDat$aspect_norm), Incline = 0) 




str(pred_grid)
str(spatDat)

pred_grid$prediction <- predict(spatMod, newdata = pred_grid, 
                                type = "response",estimate_method = "median")$estimate


dev.new()
pred_grid %>%
    ggplot(aes(x,y))+
    geom_tile(aes(fill = prediction)) +
  scale_fill_viridis_c()+
    facet_grid(location~month)

  
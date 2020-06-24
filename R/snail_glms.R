# Analysis of snail ecological data
library(tidyverse)
library(Rtsne)
library(cluster)
library(logisticPCA)
library(DHARMa)
library(ggeffects)
library(sp)
library(ggpubr)
library(mvabund)


# Read in data 
eco_dat <- read.csv("data/data_wave_only_arc_sax_only_Rob.csv") %>% distinct() %>% 
            mutate(Aspect  = (droplevels(na_if(Aspect , y="?"))), # tidy NAs 
                    Aspect = as.numeric(levels(Aspect))[Aspect]) %>%  # coerce to numeric
            group_by(location) %>% drop_na(z) %>% 
            mutate(z = z - min(z)) %>% ungroup() %>% mutate(year=factor(year)) # relative elevation to station


# Transform aspect
eco_dat<- eco_dat %>% mutate(aspect_norm = case_when(
                                    location == "E" & Aspect >= 220 ~ Aspect -220,
                                    location == "E" & Aspect <= 220 ~ 360 - Aspect , 
                                    location == "F" & Aspect >= 180 ~ Aspect -180,
                                    location == "F" & Aspect <= 180 ~ 360 - Aspect, 
                                    location == "W" & Aspect >= 50  ~ Aspect - 50, 
                                    location == "W" & Aspect <= 50  ~ 360 - Aspect)) %>% # Normalise so 0 degrees = facing sea
                        mutate(aspect_norm =  case_when(aspect_norm > 180 ~ abs(aspect_norm-360) , # flip measurements when > 180 degrees 
                                                                         TRUE ~ aspect_norm)) %>% 
  rownames_to_column()



# TSNE  ---------------------------------- 
# Visualise differences in communities 
commDatR <- eco_dat %>% select(18:30) %>% rownames_to_column()
cl_dist <- daisy(commDatR , metric="gower")
out     <- Rtsne(commDatR, perplexity=30)
tRES    <- cbind(eco_dat , X= out$Y[,1] ,  Y=out$Y[,2])

  tRES %>% 
    ggplot(aes(X , Y)) +
        geom_point(aes(fill = season),size=1.2,pch=21)+
        scale_fill_brewer(type = "div",palette = 5)+
        coord_equal()+ 
        theme_minimal()+
        facet_grid(location~species) 
  

# Mulitvariate models  ---------------------------------- 

envDat  <- eco_dat %>%select(rowname,location , z , Incline , Aspect , season , month, year) %>% drop_na()
commDatA <- eco_dat %>% filter(rowname %in% envDat$rowname) %>% select(18:30) %>% mvabund()
commMod <- manyglm(commDatA ~ z + Incline + location + Aspect + year + month , data=envDat, family="binomial")

# Community composition varies by height,location,and month (sig dif in july and october). 
multinomTab <- summary(commMod)
saveRDS(multinomTab,"data/multinomTab.rds")

# ---------------------------------- 

# Summarise community data as LPCCS 
commDat <- eco_dat %>% select(18:30)
logpca_cv  <-  cv.lpca(commDat, ks = 2, ms = 1:10)
plot(logpca_cv)
logpca_model <-  logisticPCA(commDat, k = 2, m = which.min(logpca_cv))
plot(logpca_model, type = "trace")

# inspect PCS
species <- eco_dat$species 
plot(logpca_model, type = "scores") + geom_point(aes(colour = species)) + 
  ggtitle("Logistic PCA") + scale_colour_manual(values = c("blue", "red"))


# Set up data
# select vars and incorporate community structure scores
modDat <- eco_dat %>% select(month,season,species, x ,y ,z , location , year , aspect_norm , Incline , Crevice ,  18:30)
modDat$LPC1 <- logpca_model$PCs[,1]
modDat$LPC2 <- logpca_model$PCs[,2]

toscale <- c("aspect_norm" , "Incline" , "LPC1" , "LPC2")
scaleAttr <- scale(modDat[toscale] , scale = TRUE , center = TRUE)
modDat[toscale] <- scale(modDat[toscale] , scale = TRUE , center = TRUE)
modDat <- modDat %>% group_by(location) %>% mutate(x = x-median(x) , y=y-median(y)) %>% drop_na()

# Save
saveRDS(modDat , "data/snailModDat.rds") 

modDat %>% 
    ggplot(aes(z,LPC1))+
    geom_point() +
    facet_wrap(~location, scales = "free_x") +
    theme_minimal()

# quick & dirty model selection
mod   <-  glm(species ~ month + y + z + location + year + aspect_norm + Incline + LPC1+LPC2, data = modDat , family = binomial)
mod1  <-  glm(species ~ month + year + y + z*location   + aspect_norm + Incline + LPC1+LPC2, data = modDat , family = binomial)
mod2  <-  glm(species ~ x + y + z + location + year + month , modDat, family="binomial")
mod3  <-  glm(species ~ x + y + z + location + year + month + aspect_norm + Incline, data = modDat , family = "binomial")
mod4  <-  glm(species ~ x + y + z + location + year + month + LPC1*LPC2, data = modDat , family = "binomial")

AIC(mod,mod1, mod2 , mod3 , mod4) # saturated model - x coords best fit
simRES <- DHARMa::simulateResiduals(mod1) ; plot(simRES) 
modSum<- summary(mod1) # effect sizes

saveRDS(mod1 , "data/snail_mod.rds")
modS <- mod1


# -----------  Simulate marginal effects

# Function to rescale scaled vars
rescale <- function(modDat , scaleAttr){
  t(apply(modDat[toscale], 1,
          function(r)r*attr(scaleAttr,'scaled:scale') + attr(scaleAttr, 'scaled:center')))}

# Effect of aspect ratio
newData <- expand.grid(z = 1.9, 
                       location = levels(modDat$location), 
                       month = "july",
                       x = 1 , y = 0 ,  
                       LPC1 = 1 ,LPC2 = 1 , year = "2018",
                       aspect_norm = seq(-2,2,length.out = 100), Incline = 0 , Crevice = "1")

preds <- predict(modS , newdata = newData , type = "response", se.fit = TRUE)
predDat <- cbind(newData , preds = preds$fit , SE = preds$se.fit)
predDat[toscale] <- rescale(predDat , scaleAttr)

P1 <- predDat %>% 
  ggplot(aes(aspect_norm , preds, group = location))+
  scale_y_continuous(limits = c(0,1))+
  geom_line(aes(lty=location))+
  scale_fill_brewer(type = "qual",palette = 6)+
  geom_ribbon(aes(x= aspect_norm, ymin=preds-SE,ymax=preds+SE, fill = location),alpha=0.2)+
  theme_classic()+
  labs(y="Pr(Saxatalis)", x="Aspect")
  

# Effect of elevation
newData <- expand.grid(z = seq(0,4,length.out = 100), 
                       location = levels(modDat$location), 
                       month = "july",
                       x = 1 , y = 0 ,  
                       LPC1 = 1 ,LPC2 = 1 , year = "2018",
                       aspect_norm = 0, Incline = 0 , Crevice = "1")





preds <- predict(modS , newdata = newData ,type = "response", se.fit = TRUE)
predDat <- cbind(newData , preds = preds$fit , SE = preds$se.fit)

P2 <- predDat %>% 
  ggplot(aes(z , preds, group = location))+
  scale_y_continuous(limits = c(0,1))+
  geom_line(aes(lty=location))+
  geom_ribbon(aes(x= z, ymin=preds-SE,ymax=preds+SE, fill = location),alpha=0.2)+
  scale_fill_brewer(type = "qual",palette = 6)+
  theme_classic()+
  theme(axis.title.y = element_blank() , 
        axis.text.y = element_blank())+
  labs(y="Pr(Saxatalis)", x="Height")




# Effect of incline
newData <- expand.grid(z = 1.9, 
                       location = levels(modDat$location), 
                       x = 1 , y = 0 ,  
                       month = "july", year = "2018",
                       LPC1 = seq(-2,2,length.out = 100) ,LPC2 = 1 , 
                       aspect_norm = 0, Incline = 0 , Crevice = "1")

preds <- predict(modS , newdata = newData ,type = "response", se.fit = TRUE)
predDat <- cbind(newData , preds = preds$fit , SE = preds$se.fit)
predDat[toscale] <- rescale(newData , scaleAttr )


P3 <- predDat %>% 
  ggplot(aes(LPC1, preds, group = location))+
  scale_y_continuous(limits = c(0,1))+
  geom_line(aes(lty=location))+
  geom_ribbon(aes(x= LPC1, ymin=preds-SE,ymax=preds+SE, fill = location),alpha=0.2)+
  scale_fill_brewer(type = "qual",palette = 6)+
  theme_classic()+
  theme(axis.title.y = element_blank() , 
        axis.text.y = element_blank())+
  labs(y="Pr(Saxatalis)", x="LPC1")




ggarrange(P1,P3,P2 , common.legend = TRUE, nrow = 1)




library(plotly)
# x y z plots 
filter(eco_dat , location == "F") %>% 
 plot_ly(x = .$x, 
                y = .$y, 
                z = .$z ,
                color = .$species
                  ) %>% add_markers(size = 0.5)

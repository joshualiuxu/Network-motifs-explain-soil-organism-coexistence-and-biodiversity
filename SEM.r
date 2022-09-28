#SEM

library(piecewiseSEM)

metanew <- read.csv("metanew.csv", row.names=1)
data=metanew
data[,5:44]=scale(data[,5:44],center=T,scale=T)
data=as.data.frame(data)

#SEM model for cycfac motif 
sem.model1 <- psem (
lm(biodiversity ~ fragility+ cycfac_total + modularity +Distance_from_equator+Elevation+TSEA+MAP+pH+Texture+Soil_C,  data = data),
lm(fragility~ cycfac_total + modularity +Distance_from_equator+Elevation+TSEA+MAP+pH+Texture+Soil_C, data = data),
lm(modularity~  cycfac_total +Distance_from_equator+Elevation+TSEA +MAP  +pH+Texture+Soil_C ,data = data),
lm(cycfac_total~  Distance_from_equator+Elevation+TSEA+MAP+pH+Texture+Soil_C, data = data),
lm(pH~ Distance_from_equator+Elevation+TSEA+MAP , data = data),
lm(Soil_C~ Distance_from_equator+Elevation+TSEA+MAP+pH +Texture  , data = data),
lm(Texture~ Distance_from_equator+Elevation+TSEA+MAP , data = data),
lm(TSEA~ Distance_from_equator+Elevation, data = data),
lm(MAP~Distance_from_equator+Elevation +TSEA  , data = data)
)
#export result
summary(sem.model1)



#SEM model for facmcom motif
sem.model2 <- psem (
lm(biodiversity ~ fragility+ facmcom_total + modularity +Distance_from_equator+Elevation+TSEA+MAP+pH+Texture+Soil_C,  data = data),
lm(fragility~ facmcom_total + modularity +Distance_from_equator+Elevation+TSEA+MAP+pH+Texture+Soil_C, data = data),
lm(modularity~  facmcom_total +Distance_from_equator+Elevation+TSEA +MAP  +pH+Texture+Soil_C ,data = data),
lm(facmcom_total~  Distance_from_equator+Elevation+TSEA+MAP+pH+Texture+Soil_C, data = data),
lm(pH~ Distance_from_equator+Elevation+TSEA+MAP , data = data),
lm(Soil_C~ Distance_from_equator+Elevation+TSEA+MAP+pH +Texture  , data = data),
lm(Texture~ Distance_from_equator+Elevation+TSEA+MAP , data = data),
lm(TSEA~ Distance_from_equator+Elevation, data = data),
lm(MAP~Distance_from_equator+Elevation +TSEA  , data = data)
)
#export result
summary(sem.model2)


#SEM model for trancom motif
sem.model3 <- psem (
lm(biodiversity ~ fragility+ trancom_total + modularity +Distance_from_equator+Elevation+TSEA+MAP+pH+Texture+Soil_C,  data = data),
lm(fragility~ trancom_total + modularity +Distance_from_equator+Elevation+TSEA+MAP+pH+Texture+Soil_C, data = data),
lm(modularity~  trancom_total +Distance_from_equator+Elevation+TSEA +MAP  +pH+Texture+Soil_C ,data = data),
lm(trancom_total~  Distance_from_equator+Elevation+TSEA+MAP+pH+Texture+Soil_C, data = data),
lm(pH~ Distance_from_equator+Elevation+TSEA+MAP , data = data),
lm(Soil_C~ Distance_from_equator+Elevation+TSEA+MAP+pH +Texture  , data = data),
lm(Texture~ Distance_from_equator+Elevation+TSEA+MAP , data = data),
lm(TSEA~ Distance_from_equator+Elevation, data = data),
lm(MAP~Distance_from_equator+Elevation +TSEA  , data = data)
)
#export result
summary(sem.model3)


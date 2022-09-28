#Driving forces for network motifs 
#Mantel test

library(linkET)
library(dplyr)
library(ggplot2)

top <- read.csv("top.csv", row.names=1)

spec=top[,20:24]
chem=top[,10:19]
mantel_top <- mantel_test(spec, chem,
spec_dist =  dist_func(.FUN = "vegdist", method = "euclidean"),
env_dist = dist_func(.FUN = "vegdist", method = "euclidean"),
spec_select = list(cycfac = 1:1,
facmcom = 2:2,
tranfac = 3:3,
trancom = 4:4,
trancomfac = 5:5
)) %>%
mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf),
labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

#output
set_corrplot_style()
p_top= qcorrplot(correlate(chem,  method = "spearman"), type = "upper", diag = FALSE) +
geom_square() +
geom_couple(aes(colour = pd, size = rd),
data = mantel_top, drop = TRUE, nudge_x = 0.5,
curvature = nice_curvature()) +
scale_size_manual(values = c(0.5, 1, 2)) +
scale_colour_manual(values = c("#D95F02" ,"#1B9E77" ,"#CCCCCC99")) +
guides(size = guide_legend(title = "Mantel's r",
override.aes = list(colour = "grey35"),
order = 2),
colour = guide_legend(title = "Mantel's P",
override.aes = list(size = 3),
order = 1),
fill = guide_colorbar(title = "Spearman's rho", order = 3))
p_top
ggsave("p_top.pdf", p_top,height = 6,width=8)

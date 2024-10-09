mtcars.pca <- prcomp(mtcars[,c(1:7,10,11)], center = TRUE, scale. = TRUE)
summary(mtcars.pca)
str(mtcars.pca)

ggbiplot(mtcars.pca, labels=rownames(mtcars))

ggbiplot(mtcars.pca,ellipse = TRUE,labels = rownames(mtcars), groups = mtcars.country)
ggbiplot(mtcars.pca,ellipse = TRUE, choices = c(3,4),labels = rownames(mtcars), groups = mtcars.country)

ggbiplot(mtcars.pca,ellipse=TRUE,circle=TRUE, labels=rownames(mtcars), groups=mtcars.country)
ggbiplot(mtcars.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1, labels=rownames(mtcars),groups=mtcars.country)

ggbiplot(mtcars.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE,labels=rownames(mtcars), groups=mtcars.country)

ggbiplot(mtcars.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1, labels=rownames(mtcars),
         groups=mtcars.country) +
  scale_colour_manual(name="Origin", values= c("forest green", "red3", "dark blue"))+
  ggtitle("PCA of mtcars dataset")+
  theme_minimal()+
  theme(legend.position = "bottom")




# --------------kidney disease-------------
kidneyData <- read.csv("kidney_disease.csv", header = TRUE)
str(kidneyData)

kidneyData$pcv <- sapply(kidneyData$pcv, as.numeric)
kidneyData$wc <- sapply(kidneyData$wc, as.numeric)
kidneyData$rc <- sapply(kidneyData$rc, as.numeric)

kidneyData <- kidneyData %>%mutate(age = na.approx(age)) %>% mutate(bp = na.approx(bp)) %>% mutate(sg = na.approx(sg)) %>% mutate(al = na.approx(al)) %>% mutate(su = na.approx(su)) %>% mutate(bgr = na.approx(bgr))%>% mutate(bu = na.approx(bu))%>% mutate(sc = na.approx(sc)) %>% mutate(hemo = na.approx(hemo))%>% mutate(pcv = na.approx(pcv))%>% mutate(wc = na.approx(wc))%>% mutate(rc = na.approx(rc))
kidneyData$sod[is.na(kidneyData$sod)]<-mean(kidneyData$sod,na.rm=TRUE)
kidneyData$pot[is.na(kidneyData$pot)]<-mean(kidneyData$pot,na.rm=TRUE)

missingPercentage <- colMeans(is.na(kidneyData))
data <- kidneyData[,c(2:3,11:19)]
pca <- prcomp(data, center = TRUE, scale. = TRUE)
ggbiplot(pca,ellipse = TRUE,obs.scale = 2, var.scale = 1, choices = c(1,2), groups = kidneyData$classification)

fviz_pca_ind(pca, col.ind="cos2", geom = "point") +
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=0.6)+ theme_minimal()


fviz_pca_biplot(pca, habillage=as.factor(kidneyData$classification), addEllipses=TRUE,
                label = "var", col.var = "red", col.ind = "#696969",  repel = TRUE) +
  theme_minimal() + theme_gray(base_size =12) + labs(title="", x ="PC1 (34.5% explained var.)", y = "PC2 (12.0% explained var.)") +
  coord_fixed()

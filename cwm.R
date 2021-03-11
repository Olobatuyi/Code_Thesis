library(flexCWM)
library(parallel)
library(OpenMPController)
library(data.table)
library(tsne)
library(Rtsne)
library(factoextra)
library(FactoMineR)
library(xtable)
library(MBCbook); library(data.table)
library(caret); library(clues);library(HDclassif)


#### Abalon data
abal <- read.csv("C:\\Users\\Olobatuyi\\Downloads\\ONGOING WRITEUP\\Data\\abalone-data.txt", header = F)

colnames(abal) <- c("Sex","Length","Diameter","Height","Whole_weight","Shucked_weight",
                    "Viscera_weight","Shell_weight","Rings")

#### DIM reduction
# data.frame as input
head(abal)
abal$Sex <- log(as.numeric(abal$Sex)+0.5) / max(log(as.numeric(abal$Sex)+0.5))
x <- abal[,-c(9)]
plot(abal[,c(3,4)], col = sex)

sex <- as.numeric(aba$Sex)
for(i in 1:15){

  tsne = Rtsne(x, dims = 2, perplexity=30, verbose=TRUE, max_iter = 1500, pca=T)

  colors = rainbow(length(unique(aba$Sex)))
  names(colors) = unique(aba$Sex)
  plot(tsne$Y, t='n', main="tsne")
  text(tsne$Y, labels=aba$Sex, col=colors[aba$Sex])

  readline(prompt="Press [enter] to continue")
}

colors = rainbow(length(unique(sex)))
tsne_out <- Rtsne(x, dims = 3, pca = FALSE, theta=0.0, check_duplicates = F)
plot(tsne_out$Y, col = colors, asp=0, pch = "+", bty = "n", xaxt = "n",ann = FALSE, yaxt = "n")


abal <- cbind(abal, tsne$Y)
colnames(abal) <- c("Sex", "Length", "Diameter", "Height", "Whole_weight", "Shucked_weight",
                    "Viscera_weight", "Shell_weight", "Rings",  "tsne.1","tsne.2")
head(abal)
attach(abal)

summary(abal)


#log(as.numeric(Sex)+0.5)

fit_EII <- cwm(Rings ~ Sex + Length+Diameter+Height+Whole_weight+Shucked_weight+
               Viscera_weight+Shell_weight, familyY = gaussian, k = 3, initialization = "random.soft",
            Xnorm = cbind(Length,Diameter,Height,Whole_weight,Shucked_weight,
                          Viscera_weight,Shell_weight), modelXnorm = "VEE", data = abal)

set.seed(1)

fit_EII <- cwm(scale(Rings) ~ tsne.1+tsne.2, familyY = gaussian, k = 1:5, initialization = "random.soft", 
               Xnorm = cbind(tsne.1,tsne.2), modelXnorm = "VVV")


getIC(fit_EII, "BIC")
clus <- getCluster(fit_EII, "BIC")
adjustedRand(Rings, clus)


d <- data.frame(EII,VII,EEI,VEI,EVI,VVI,EEE,VEE,EVE,EEV,VVE,VEV,EVV,VVV)
Comp <- 1:5

ggplot(d, aes(x = Comp, BIC))+
  geom_line(aes(y = EII, colour = "EII"))+
  geom_line(aes(y = VII, colour = "VII"))+
  geom_line(aes(y = EEI, colour = "EEI"))+
  geom_line(aes(y = VEI, colour = "VEI"))+
  geom_line(aes(y = EVI, colour = "EVI"))+
  geom_line(aes(y = VVI, colour = "VVI"))+
  geom_line(aes(y = EEE, colour = "EEE"))+
  geom_line(aes(y = VEE, colour = "VEE"))+
  geom_line(aes(y = EVE, colour = "EVE"))+
  geom_line(aes(y = EEV, colour = "EEV"))+
  geom_line(aes(y = VVE, colour = "VVE"))+
  geom_line(aes(y = VEV, colour = "VEV"))+
  geom_line(aes(y = EVV, colour = "EVV"))+
  geom_line(aes(y = VVV, colour = "VVV"))



# fit_kmeans1 <- cwm(Rings ~ tsne.1,familyY = gaussian, k = 1:5, initialization = "random.soft", Xnorm = Length,
#                   modelXnorm = "E", data = abal)

#c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV")

table(abal$Rings)
table(wee)
unique(aba$Rings)
clus <- getCluster(fit_VVI, "BIC")
plot(abal[-c(9,10,11)], col = Rings)
text(abal[,c(2,5)], labels=aba$Sex, col=clus)
xtable(table(clus, abal$Sex))

for (i in 1:length(Rings)) {
  
  if(Rings[i] ==1 || Rings[i] ==2|| Rings[i] ==3 || Rings[i] ==4|| Rings[i] ==5 || Rings[i] ==6 || Rings[i] ==7|| Rings[i] ==8) Rings[i] <- 1
  else if(Rings[i] ==9 || Rings[i] ==10) Rings[i] <- 2
  else Rings[i] <- 3
}

table(clus)
# colors = rainbow(length(unique(sex)))
# tsne_out <- Rtsne(abal[,c(10,11, 12)], dims = 3, pca = FALSE, theta=0.0, check_duplicates = F)
# plot(tsne_out$Y, col = colors, asp=0, pch = "+", bty = "n", xaxt = "n",ann = FALSE, yaxt = "n")

levels(wee) <- 1:3
confusionMatrix(table(Rings,clus))
library(caret); library(clues);library(HDclassif)

getIC(fit_VVI)
adjustedRand(Rings, clus)


r <- PCA(abal[,-c(9,10,11)], ncp = 5, graph = FALSE)

fviz_pca_ind(r, geom = "point", axes = c(2,5),habillage = as.factor(clus),
             palette = c("red", "blue", "gold","magenta", "black", "cyan","orange", 
                         "grey", "green","brown"),addEllipse = F, repel = F,
             ggtheme = theme_minimal(), title = "")

# Extractor
clus <- getCluster(fit_kmeans)
table(clus, abal$Sex)




EII <- c(-87696.1,-2109.2,-85062.3,-84460.1,-32541.1)
VII <- c(-87696.1,-2230.8,-2024.8,-917.4,-341.9)
EEI <- c(-87688.0,-2115.9,-2101.1,-84412.6,NA)
VEI <- c(-87688.0,-2184.1,176940.8,-838.7,NA)
EVI <- c(-87688.0,-2006.1,-84632.3,-1225.5,NA)
VVI <- c(-87688.0,-2085.3,177254.5,NA,-32220.5)
EEE <- c(-87687.8,-2122.8,1793.2,NA,-32554.9)
VEE <- c(-87687.8,-2146.5,176932.7,-32459.0,NA)
EVE <- c(-87712.7,-1999.1,-84415.3,-31949.4,-32016.0)
EEV <- c(-87687.8,-2086.2,2651.1,-32021.3,NA)
VVE <- c(-87734.3,-85777.6,-1219.1,-118.8,-32002.6)
VEV <- c(-87687.8,-2186.5,-1337.2,182260.9,-32072.1)
EVV <- c(-87687.8,-85606.0,-1196.0,NA,-1063.1)
VVV <- c(-87687.8,-2108.3,177462.9,715,-32028.74)

we <- cbind(EII, VII, EEI, VEI, EVI, VVI, EEE, VEE, EVE, EEV, VVE, VEV, EVV, VVV)
matplot(we, type = "l", xlab = "Components", ylab = "BIC", col = 1:14)
legend("topleft", legend = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", 
                                  "EEV", "VVE", "VEV", "EVV", 'VVV'), col = 1:14, lty = 1:14)

colnames(tsne_out$Y) <- c("a","b")
redu <- data.frame(tsne_out$Y, abal$Rings)
colnames(redu) <- c("a","b","R")
attach(redu)
head(redu)

fit_kmeans <- cwm(R ~ a+b, familyY = gaussian, k = 1:5, initialization = "kmeans", 
                  Xnorm = cbind(a,b), modelXnorm = "VVI", data = redu)

#c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV")
summary(fit_kmeans)
plot(fit_kmeans)



#### Protein data_mclust
prot <- read.csv("C:\\Users\\Olobatuyi\\Downloads\\ONGOING WRITEUP\\Data\\Protein\\ecoli-data.txt", 
                 header = F, sep = "")

dim(prot)
str(prot)

colnames(prot) <- c("Seq_name", "mcg", "gvh", "lip", "chg", "aac", "alm1", "alm2", 
                    "class", "tsne.1","tsne.2","tsne.3")

set.seed(1)

for(i in 1:15){
  
  tsne = Rtsne(prot[,-c(1,9)], dims = 3, perplexity=i, verbose=TRUE, max_iter = 1000, pca=T)
  
  colors = rainbow(length(unique(prot$class)))
  names(colors) = unique(prot$class)
  plot(tsne$Y, t='n', main="tsne")
  text(tsne$Y, labels=prot$class, col=colors[prot$class])
  
  readline(prompt="Press [enter] to continue")
}

dim(tsne$Y)
colors = rainbow(length(unique(clus)))
names(colors) = unique(clus)
tsne_out <- Rtsne(prot[,-c(1,9)], dims = 3, pca = FALSE, theta=0.0, check_duplicates = F)
plot(tsne_out$Y, col = colors, asp=0, pch = "+", bty = "n", xaxt = "n",ann = FALSE, yaxt = "n")
prot <- cbind(prot, tsne$Y)

attach(prot)
str(prot)

fit2 <- cwm(mcg ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:8, initialization = "random.soft", 
            Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VVE", data = prot)

summary(fit2)

# Extractor
clus <- getCluster(fit2, "AIC")
plot(prot[,c(2,10,11,12)], col = clus)
xtable(table(clus, prot$class))
getIC(fit2)

wee <- prot$class
levels(wee) <- 1:8
confusionMatrix(table(sort(wee), sort(clus)))
library(caret); library(clues);library(HDclassif)
adjustedRand(as.vector(wee), clus)



EII <- c(-8291.1,-7270.3,-7217.3,-6433.0,-6479.2,-6309.6,-6937.0,-6815.8)
VII <- c(-7960.9,-7223.6,-6565.9,-6474.5,-6336.2,-6394.1,-6337.5,-6183.0)
EEI <- c(-7960.9,-7422.0,-6577.6,-6437.5,-6389.8,-6453.8,-6210.9,-6213.5)
VEI <- c(-7960.9,-7086.6,-7262.7,-6411.8,-6404.0,-6304.2,-6213.3,-6346.9)
EVI <- c(-7960.9,-7196.4,-6566.5,-6644.1,-6440.1,-6361.3,-6214.6,-6230.0)
VVI <- c(-7960.9,-7223.6,-6565.9,-6474.5,-6336.2,-6394.1,-6337.5,-6183.0)
EEE <- c(-7412.5,-7005.6,-6460.8,-6327.7,-6284.3,-6201.5,-6295.4,-6214.2)
VEE <- c(-7412.5,-7143.0,-6457.2,-6413.9,-6219.8,-6222.3,-6215.4,-6155.1)
EVE <- c(-8182.5,-6997.2,-6452.2,-6192.0,-6147.5,-6168.1,-6186.7,-6227.5)
EEV <- c(-7412.5,-6745.1,-6429.6,-6431.5,-6180.8,-6182.5,-6150.6,-6247.5)
VVE <- c(-8189.6,-6994.1,-6954.8,-6176.8,-6132.6,-6222.1,-6185.7,-6198.5)
VEV <- c(-7412.5,-6941.4,-6503.3,-6199.0,-6174.0,-6174.7,-6324.7,-6129.2)
EVV <- c(-7412.5,-6733.9,-6473.3,-6229.9,-6242.5,-6241.3,-6184.3,-6299.9)
VVV <- c(-7412.5,-6822.0,-6412.8,-6363.1,-6154.2,-6183.6,-6241.4,-6192.4)


we <- cbind(EII, VII)
matplot(we, type = "l", xlab = "Components", ylab = "BIC", col = 1:14)
legend("topleft", legend = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", 
                             "EEV", "VVE", "VEV", "EVV", 'VVV'), col = 1:14, lty = 1:14)
d <- data.frame(EII,VII,EEI,VEI,EVI,VVI,EEE,VEE,EVE,EEV,VVE,VEV,EVV,VVV)
Comp <- 1:8

ggplot(d, aes(x = Comp, BIC))+
  geom_line(aes(y = EII, colour = "EII"))+
  geom_line(aes(y = VII, colour = "VII"))+
  geom_line(aes(y = EEI, colour = "EEI"))+
  geom_line(aes(y = VEI, colour = "VEI"))+
  geom_line(aes(y = EVI, colour = "EVI"))+
  geom_line(aes(y = VVI, colour = "VVI"))+
  geom_line(aes(y = EEE, colour = "EEE"))+
  geom_line(aes(y = VEE, colour = "VEE"))+
  geom_line(aes(y = EVE, colour = "EVE"))+
  geom_line(aes(y = EEV, colour = "EEV"))+
  geom_line(aes(y = VVE, colour = "VVE"))+
  geom_line(aes(y = VEV, colour = "VEV"))+
  geom_line(aes(y = EVV, colour = "EVV"))+
  geom_line(aes(y = VVV, colour = "VVV"))
  
  
################## SEIZURE

seizure <- read.csv("C:\\Users\\Olobatuyi\\Downloads\\ONGOING WRITEUP\\Data\\Epileptic Seizure Recognition\\Epileptic Seizure Recognition.csv", 
                    header = T)
head(seizure)
dim(seizure)
seizure <- seizure[,-1]

for(i in 9:15){
  
  tsne = Rtsne(seizure[,-179], dims = 2, perplexity=250, verbose=TRUE, max_iter = 10000, pca=T, 
               check_duplicates = F)
  
  colors = rainbow(length(unique(seizure$y)))
  names(colors) = unique(seizure$y)
  plot(tsne$Y, t='n', main="tsne")
  text(tsne$Y, labels=seizure$y, col=colors[seizure$y])
  
  readline(prompt="Press [enter] to continue")
}

?optimSA
seizure <- cbind(seizure, tsne$Y)

# colors = rainbow(length(unique(seizure$y)))
# names(colors) = unique(seizure$y)
# tsne_out <- Rtsne(seizure[,-179], dims = 3, pca = F, theta=0.2, check_duplicates = F)
# plot(tsne_out$Y, col = colors, asp=0, pch = "+", xaxt = "n",ann = FALSE, yaxt = "n")
# seizure <- cbind(seizure, tsne_out$Y)

dim(seizure)

nam <- colnames(seizure)
nam[c(180, 181)] <- c("tsne.1","tsne.2")
colnames(seizure) <- nam
seizure$y <- log(seizure$y)
attach(seizure)


set.seed(1)

fit_EIIs <- cwm(scale(y) ~., familyY = gaussian, k = 1:3, initialization = "random.soft", 
               Xnorm = cbind(X1,X2,X3), modelXnorm = "EII", data = seizure)


fit_EI <- cwm(scale(y) ~ tsne.1+tsne.2, familyY = gaussian, k = 1:3, initialization = "random.soft", 
                Xnorm = cbind(tsne.1,tsne.2), modelXnorm = "EEI")

fit_EI <- cwm(y1 ~ tsne.1+tsne.2, familyY = binomial, k = 1:3, initialization = "random.soft", 
              Xnorm = cbind(tsne.1,tsne.2), modelXnorm = "VVE")

getIC(fit_EI, "ICL")
clus <- getCluster(fit_EI, "BIC")
adjustedRand(y1, clus)
confusionMatrix(table(y1,clus-1))
table(clus)



EII_BIC <- c(-169881,-166546,-163231);EII_ICL <- c(-169881,-169186,-164750)
VII_BIC <- c(-169881,-166556,-162870);VII_ICL <- c(-169881,-169185,-164866)
EEI_BIC <- c(-168667,-164041,-162360);EEI_ICL <- c(-168667,-165176,-165164)
VEI_BIC <- c(-168667,-164050,-162591);VEI_ICL <- c(-168667,-165176,-164680)
EVI_BIC <- c(-168667,-164050,-160881);EVI_ICL <- c(-168667,-165185,-162118)
VVI_BIC <- c(-168667,-164059,-159935);VVI_ICL <- c(-168667,-165183,-161259)
EEE_BIC <- c(-168593,-163889,-163367);EEE_ICL <- c(-168593,-165081,-166885)
VEE_BIC <- c(-168593,-163898,-162569);VEE_ICL <- c(-168593,-165090,-164475)
EVE_BIC <- c(-169086,-164816,-161938);EVE_ICL <- c(-169086,-166081,-163565)
EEV_BIC <- c(-168593,-163898,-162271);EEV_ICL <- c(-168593,-165090,-163699)
VVE_BIC <- c(-169109,-163927,-159226);VVE_ICL <- c(-169109,-165145,-160531)
VEV_BIC <- c(-168593,-163908,-162034);VEV_ICL <- c(-168593,-165099,-163612)
EVV_BIC <- c(-168593,-163908,-161937);EVV_ICL <- c(-168593,-165103,-163613)
VVV_BIC <- c(-168593,-163149,-159204);VVV_ICL <- c(-168593,-164215,-160499)


r <- PCA(seizure[,c(180,181)], ncp = 5, graph = FALSE)

fviz_pca_ind(r, geom = "point", axes = c(1,2),habillage = as.factor(clus),
             palette = c("red", "black", "blue","cyan","orange", 
                         "grey", "green","brown"),addEllipse = F, repel = T,
             ggtheme = theme_minimal(), title = "")

y1 <- ifelse(y==1, 0, 1)
table(y1, clus)

library(devtools)

install_github("fawda123/ggord")

library(ggord)
library(klaR)

ldaw <- lda(y ~., data = seizure)
ggord::ggord(ldaw, as.factor(clus-1))
klaR::partimat(as.factor(clus-1)~ X1+X2+X3, method = "lda")
ldaw$prior
ldaw$call
ldaw$counts
ldaw
ldaw$terms[1]
table(y, )

fit_eii <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
            Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EII", data = seizure)
fit_vii <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VII", data = seizure)
fit_eei <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EEI", data = seizure)
fit_vei <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VEI", data = seizure)
fit_evi <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EVI", data = seizure)
fit_vvi <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VVI", data = seizure)
fit_eee <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EEE", data = seizure)
fit_vee <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VEE", data = seizure)
fit_eve <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EVE", data = seizure)
fit_vve <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VVE", data = seizure)
fit_eev <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EEV", data = seizure)
fit_evv <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EVV", data = seizure)
fit_vvv <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:6, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VVV", data = seizure)
fit_vev <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VEV", data = seizure)

summary(fit2)

# Extractor
clus <- getCluster(fit_vev, "AIC")
adjustedRand(as.vector(wee), clus)
plot(seizure[,c(180,181,182)], col = clus)

xtable(table(clus, seiz$y))

getIC(fit_vev)
table(seiz$y)
wee <- seiz$y
confusionMatrix(table(sort(wee), sort(clus)))
adjustedRand(as.vector(wee), clus)

r <- PCA(fas[,-1], ncp = 5, graph = FALSE)

fviz_pca_ind(r, geom = "point", axes = c(1,3), 
             habillage = as.factor(fas$label),
             palette = c("red", "blue", "gold","magenta", "black", "cyan","orange", 
                         "grey", "green","brown"),
             addEllipse = F, repel = TRUE,
             ggtheme = theme_minimal(), title = "")

############# GENE EXPRESSION 

gene <- fread("C:\\Users\\Olobatuyi\\Downloads\\ONGOING WRITEUP\\Data\\TCGA-PANCAN-HiSeq-801x20531\\data.csv", 
              header = T)
gene_label <- fread("C:\\Users\\Olobatuyi\\Downloads\\ONGOING WRITEUP\\Data\\TCGA-PANCAN-HiSeq-801x20531\\labels.csv",
              header = T)

gene <- cbind(gene, gene_label$Class)
gene <- gene[,-1]
dim(gene)
head(gene)
colnames(gene)[20532] <- "Class"

gene <- as.matrix(gene)

#### Tsne for gene

colors = rainbow(length(unique(gene$Class)))
names(colors) = unique(gene$Class)
tsne_out <- Rtsne(gene[,-20532], dims = 3, pca = F, theta=0.5, check_duplicates = F)
plot(tsne_out$Y, col = colors, asp=0, pch = "+", xaxt = "n",ann = FALSE, yaxt = "n")

dim(tsne_out$Y)
gene <- cbind(gene, tsne_out$Y)

dim(seizure)

nam <- colnames(seizure)
nam[c(180, 181,182)] <- c("tsne.1","tsne.2","tsne.3")
colnames(seizure) <- nam
seizure$y <- log(seizure$y)
attach(seizure)
levels(gene$Class)

gene <- as.matrix(gene)
class(gene)

for(i in 1:15){
  
  tsne = Rtsne(gene[,-20532], dims = 2, perplexity=1000, verbose=TRUE, max_iter = 4000, pca=T, 
               check_duplicates = F)
  
  colors = rainbow(length(unique(gene$Class)))
  names(colors) = unique(gene$Class)
  plot(tsne$Y, t='n', main="tsne")
  text(tsne$Y, labels=gene$Class, col=colors[gene$Class])
  
  readline(prompt="Press [enter] to continue")
}


############## Fashion Image data

fash <- fread("D:\\traindcm\\fashionmnist\\fashion-mnist_train.csv", header = T)
head(fash)
dim(fash)

for(i in 1:15){
  
  tsne = Rtsne(fash[,-1], dims = 2, verbose=TRUE, max_iter = 10000, pca=T, perplexity = 500,
               check_duplicates = F)
  
  colors = rainbow(length(unique(fash$label)))
  names(colors) = unique(fash$label)
  plot(tsne$Y, t='n', main="tsne")
  text(tsne$Y, labels=fash$label, col=colors[fash$label])
  
  readline(prompt="Press [enter] to continue")
}

dim(tsne$Y)
fas <- cbind(fash$label, tsne$Y)
colnames(fas) <- c("label", "tsne.1", "tsne.2", "tsne.3")
fas <- as.data.frame(fas)
fas$label <- log(fas$label + 0.5)
attach(fas)
fas <- as.data.frame(fas)
fit_eii <- cwm(label ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:10, initialization = "mclust", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EII", data = fas)

fit_vii <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "mclust", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VII", data = seizure)
fit_eei <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "mclust", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EEI", data = seizure)
fit_vei <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VEI", data = seizure)
fit_evi <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EVI", data = seizure)
fit_vvi <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VVI", data = seizure)
fit_eee <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EEE", data = seizure)
fit_vee <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VEE", data = seizure)
fit_eve <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EVE", data = seizure)
fit_vve <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VVE", data = seizure)
fit_eev <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EEV", data = seizure)
fit_evv <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "EVV", data = seizure)
fit_vvv <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VVV", data = seizure)
fit_vev <- cwm(y ~ tsne.1+tsne.2+tsne.3, familyY = gaussian, k = 1:5, initialization = "random.hard", 
               Xnorm = cbind(tsne.1,tsne.2,tsne.3), modelXnorm = "VEV", data = seizure)

summary(fit2)

# Extractor
clus <- getCluster(fit_eev, "AIC")
adjustedRand(as.vector(wee), clus)
plot(seizure[,c(180,181,182)], col = clus)

xtable(table(clus, seiz$y))

getIC(fit_vvv)

wee <- seiz$y
confusionMatrix(table(sort(wee), sort(clus)))
adjustedRand(as.vector(wee), clus)

EII <- c(-312028,-306019,1351293,1353885,1354440)
VII <- c(-312028,-306005,1353698,1356425,1359043)
EEI <- c(-311558,-305527,1351286,1353888,1354422)
VEI <- c(-311558,-305016,1354012,1492195,1625142)
EVI <- c(-311558,-305453,1352182,1355063,1623220)
VVI <- c(-311558,1348003,1354857,1358142,1495020)
EEE <- c(-308496,-303896,1351608,1354223,1623402)
VEE <- c(-308496,-303905,1354625,1493409,1357974)
EVE <- c(-311212,-303802,1492072,1625703,NA)

EEV <- c(-7412.5,-6745.1,-6429.6,-6431.5,-6180.8,-6182.5,-6150.6,-6247.5)

VVE <- c(-310951,-303811,1355108,-298905,1359483)
VEV <- c(-308496,1353297,1356458,13557391,1358841)
EVV <- c(-308496,-303364,-298330,1626054,NA)
VVV <- c(-308496,1353558,1357871,1626007,NA)

d <- data.frame(EII_ICL,VII_ICL,EEI_ICL,VEI_ICL,EVI_ICL,VVI_ICL,EEE_ICL,
                VEE_ICL,EVE_ICL,VVE_ICL,VEV_ICL,EVV_ICL,VVV_ICL)#,EEV,VVE,VEV,EVV,VVV)
Comp <- 1:3

ggplot(d, aes(x = Comp, ICL))+
  geom_line(aes(y = EII_ICL, colour = "EII"))+
  geom_line(aes(y = VII_ICL, colour = "VII"))+
  geom_line(aes(y = EEI_ICL, colour = "EEI"))+
  geom_line(aes(y = VEI_ICL, colour = "VEI"))+
  geom_line(aes(y = EVI_ICL, colour = "EVI"))+
  geom_line(aes(y = VVI_ICL, colour = "VVI"))+
  geom_line(aes(y = EEE_ICL, colour = "EEE"))+
  geom_line(aes(y = VEE_ICL, colour = "VEE"))+
  geom_line(aes(y = EVE_ICL, colour = "EVE"))+
  geom_line(aes(y = EEV_ICL, colour = "EEV"))+
  geom_line(aes(y = VVE_ICL, colour = "VVE"))+
  geom_line(aes(y = VEV_ICL, colour = "VEV"))+
  geom_line(aes(y = EVV_ICL, colour = "EVV"))+
  geom_line(aes(y = VVV_ICL, colour = "VVV"))

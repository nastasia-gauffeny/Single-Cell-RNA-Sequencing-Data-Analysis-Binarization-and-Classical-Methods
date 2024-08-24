
## Nastasia Gauffeny - nastasia.gauffeny@live.com 

### Pipeline scRNA-seq : binarisation VS méthodes classiques

# 8. Fin ####

library(rio)
library(ggplot2)
library(scales)
library(dplyr)
library(ggpubr)

data <- import_list("./Data/BP_human_mm.xlsx")

## Visualisations des critères de choix de résolution et distances #### 

## Nombre de clusters par résolution

df1 <- data[[1]]
df1$Résolution <- as.factor(df1$Résolution)
df1$`Matrice de comptage` <- as.factor(df1$`Matrice de comptage`)
df1$Intégration <- as.factor(df1$Intégration)
df1$Intégration <- factor(df1$Intégration, levels = c("Avant", "Après"))

# Plot

g1 <- ggplot(df1, aes(x = `Matrice de comptage`, y = `Nb clusters`, color = Résolution)) +
  geom_jitter(size = 3, width = 0.05, height = 0.01) +
  scale_y_continuous(name = "Nombre de clusters", breaks = c(seq(4, 17, 1))) +
  scale_x_discrete(name = "Données") +
  labs(title = "Nombre de clusters par résolution\nselon le format des données et l'intégration") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot(g1)

p1 <- ggplot(df1, aes(x = Résolution, y = `Nb clusters`, color = `Matrice de comptage`)) +
  geom_jitter(size = 3, width = 0.05, height = 0.01) +
  scale_y_continuous(name = "Nombre de clusters", breaks = c(seq(4, 17, 1))) +
  scale_x_discrete(name = "Résolution") +
  labs(title = "Nombre de clusters par résolution\nselon le format des données et l'intégration") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(colour = "Données")
plot(p1)

## ARI par résolution

df2 <- data[[2]]
df2$Résolution <- as.factor(df2$Résolution)
df2$`Adjusted Rand Index` <- as.numeric(df2$`Adjusted Rand Index`)
df2$Comparaison <- as.factor(df2$Comparaison)
df2$Comparaison <- factor(df2$Comparaison, levels = c("Comptage-Binaire", "Binaire-Label", "Comptage-Label"))
df2 <- df2 %>%
  mutate(Comparaison = recode(Comparaison, 
                              "Comptage-Binaire" = "Comptage\nBinaire",
                              "Binaire-Label" = "Binaire\nLabel",
                              "Comptage-Label" = "Comptage\nLabel"))
df2$Intégration <- as.factor(df2$Intégration)
df2$Intégration <- factor(df2$Intégration, levels = c("Avant", "Après"))

# Plot

g2 <- ggplot(df2, aes(x = Comparaison, y = `Adjusted Rand Index`, color = Résolution)) +
  geom_jitter(size = 3, width = 0.1, height = 0) +
  scale_y_continuous(name = "ARI", breaks = breaks_pretty()) +
  scale_x_discrete(name ="Clusterings comparés") +
  labs(title = "Index de Rand Ajusté (ARI) par résolution\nselon le format des données et l'intégration") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot(g2)

p2 <- ggplot(df2, aes(x = Résolution, y = `Adjusted Rand Index`, color = Comparaison)) +
  geom_jitter(size = 3, width = 0.1, height = 0) +
  scale_y_continuous(name = "ARI", breaks = breaks_pretty()) +
  scale_x_discrete(name = "Résolution") +
  labs(title = "Index de Rand Ajusté (ARI) par résolution\nselon le format des données et l'intégration") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        legend.key.spacing.y = unit(2, 'mm'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(colour = "Clusterings\ncomparés")
plot(p2)

## Largeur moyenne de silhouette par résolution

df3 <- data[[3]]
df3$Résolution <- as.factor(df3$Résolution)
df3$`Largeur moyenne de silhouette` <- as.numeric(df3$`Largeur moyenne de silhouette`)
df3$`Matrice de comptage` <- as.factor(df3$`Matrice de comptage`)
df3$Intégration <- as.factor(df3$Intégration)
df3$Intégration <- factor(df3$Intégration, levels = c("Avant", "Après"))

# Plot

g3 <- ggplot(df3, aes(x = `Matrice de comptage` , y = `Largeur moyenne de silhouette`, color = Résolution)) +
  geom_jitter(size = 3, width = 0.1, height = 0) + 
  ylim(0, 1) +
  scale_x_discrete(name ="Données") +
  labs(title = "Largeur moyenne de silhouette par résolution\nselon le format des données et l'intégration") +
  xlab("Largeur moyenne de silhouette") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot(g3)

p3 <- ggplot(df3, aes(x = Résolution , y = `Largeur moyenne de silhouette`, color = `Matrice de comptage`)) +
  geom_jitter(size = 3, width = 0.1, height = 0) + 
  ylim(0, 1) +
  scale_x_discrete(name = "Résolution") +
  labs(title = "Largeur moyenne de silhouette par résolution\nselon le format des données et l'intégration") +
  xlab("Largeur moyenne de silhouette") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(colour = "Données")
plot(p3)

## Score de spécificité / AUC (scMiko) par résolution

df4 <- data[[4]]
df4$Résolution <- as.factor(df4$Résolution)
df4$`Score de spécificité` <- as.numeric(df4$`Score de spécificité`)
df4$`Matrice de comptage` <- as.factor(df4$`Matrice de comptage`)
df4$Intégration <- as.factor(df4$Intégration)
df4$Intégration <- factor(df4$Intégration, levels = c("Avant", "Après"))

# Plot

g4 <- ggplot(df4, aes(x = `Matrice de comptage` , y = `Score de spécificité`, color = Résolution)) +
  geom_jitter(size = 3, width = 0.1, height = 0) +
  ylim(0, 1) +
  scale_x_discrete(name ="Données") +
  labs(title = "Score de spécificité par résolution\nselon le format des données et l'intégration") +
  xlab("Score de spécificité") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot(g4)

p4 <- ggplot(df4, aes(x = Résolution , y = `Score de spécificité`, color = `Matrice de comptage`)) +
  geom_jitter(size = 3, width = 0.1, height = 0) +
  ylim(0, 1) +
  scale_x_discrete(name = "Résolution") +
  labs(title = "Score de spécificité par résolution\nselon le format des données et l'intégration") +
  xlab("Score de spécificité") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(colour = "Données")
plot(p4)

## GGarrange

plot1 <- ggarrange(g1, g2, g3, g4, legend = "right",
                   common.legend = TRUE)
print(plot1)
ggsave("./Output/Plots/08-ggplots_critères.png", plot1,
       height = 10, width = 12, bg = "white")

plot2 <- ggarrange(p1, p2, p3, p4, 
                   common.legend = FALSE)
print(plot2)
ggsave("./Output/Plots/08-ggplots_critères_2.png", plot2,
       height = 10, width = 12, bg = "white")

rm(data)
rm(df1)
rm(df2)
rm(df3)
rm(df4)
rm(g1)
rm(g2)
rm(g3)
rm(g4)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(plot1)
rm(plot2)
gc()

#
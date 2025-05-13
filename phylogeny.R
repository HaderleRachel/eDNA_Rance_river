library("U.PhyloMaker")
library(ape)
library(picante)
library(phangorn)
library(vegan)
library(hillR)
library(dplyr)


sp.list<-read.csv("sp.list.csv")

# Créer un nouveau dataframe avec les colonnes species, genus, family, species.relative, et genus.relative
sp.list <- sp.list %>%
  dplyr::select(species, genus, family) %>%  # Sélectionner les colonnes species, genus, et family du dataframe existant
  dplyr::mutate(
    species.relative = NA,  # Créer une nouvelle colonne species.relative, initialisée à NA ou à la valeur souhaitée
    genus.relative = NA     # Créer une nouvelle colonne genus.relative, initialisée à NA ou à la valeur souhaitée
  )

megatree<-read.tree("fish_megatree.tre")
gen.list<-read.csv("fish_genus_list.csv")
result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)
result
phylo <- result$phylo
phylo

#########Vertébrés
PAtaxa=read.csv2('fiche code 2.csv', row.names = NULL, stringsAsFactors = FALSE)

test_taxa_list = tibble::tibble(
  species = PAtaxa$Taxonomie,
  genus = PAtaxa$Genus,
  family = PAtaxa$Family,
  order = PAtaxa$Order
)
sp.list <- test_taxa_list %>%
  dplyr::select(species, genus, family) %>%  # Sélectionner les colonnes species, genus, et family du dataframe existant
  dplyr::mutate(
    species.relative = NA,  # Créer une nouvelle colonne species.relative, initialisée à NA ou à la valeur souhaitée
    genus.relative = NA     # Créer une nouvelle colonne genus.relative, initialisée à NA ou à la valeur souhaitée
  )

megatree_fish<-read.tree("~/Documents/R/ADNeRance/AdneRance/phylogeny/fish_megatree.tre")
megatree_mm<-read.tree("~/Documents/R/ADNeRance/AdneRance/phylogeny/mammal_megatree.tre")
megatree_oi<-read.tree("~/Documents/R/ADNeRance/AdneRance/phylogeny/bird_megatree.tre")
megatree_amp<-read.tree("~/Documents/R/ADNeRance/AdneRance/phylogeny/amphibian_megatree.tre")
megatree_rep<-read.tree("~/Documents/R/ADNeRance/AdneRance/phylogeny/reptile_megatree.tre")


#taxons communs ?
# Obtenir les labels des tips pour chaque arbre
tips_fish <- megatree_fish$tip.label
tips_mm <- megatree_mm$tip.label
tips_oi <- megatree_oi$tip.label
tips_amp <- megatree_amp$tip.label
tips_rep <- megatree_rep$tip.label

# Vérifier les taxons communs
common_fish_mm <- intersect(tips_fish, tips_mm)
common_fish_oi <- intersect(tips_fish, tips_oi)
common_fish_amp <- intersect(tips_fish, tips_amp)
common_fish_rep <- intersect(tips_fish, tips_rep)

# Afficher les résultats
cat("Taxons communs entre poissons et mammifères:", common_fish_mm, "\n")
cat("Taxons communs entre poissons et oiseaux:", common_fish_oi, "\n")
cat("Taxons communs entre poissons et amphibies:", common_fish_amp, "\n")
cat("Taxons communs entre poissons et reptiles:", common_fish_rep, "\n")

# Combiner les arbres
combined_tree <- bind.tree(megatree_fish, megatree_mm)
combined_tree <- bind.tree(combined_tree, megatree_oi)
combined_tree <- bind.tree(combined_tree, megatree_amp)

# Notez que si le dernier arbre (reptiles) n'a pas de taxon commun avec les poissons,
# il peut être utile d'évaluer son intégration à part.
# Pour le combiner, vous devrez peut-être l'intégrer différemment ou le laisser séparé.

# Afficher l'arbre combiné
plot(combined_tree, show.tip.label = FALSE)



gen.list_fish<-read.csv("~/Documents/R/ADNeRance/AdneRance/phylogeny/fish_genus_list.csv")
gen.list_mm<-read.csv("~/Documents/R/ADNeRance/AdneRance/phylogeny/mammal_genus_list.csv")
gen.list_oi<-read.csv("~/Documents/R/ADNeRance/AdneRance/phylogeny/bird_genus_list.csv")
gen.list_amp<-read.csv("~/Documents/R/ADNeRance/AdneRance/phylogeny/amphibian_genus_list.csv")

gen.list <- rbind(gen.list_fish, gen.list_mm, gen.list_oi, gen.list_amp)

result <- phylo.maker(sp.list, combined_tree, gen.list, nodes.type = 1, scenario = 3)
result
phylo <- result$phylo
phylo

jaccard<-read.csv("fiche code 3.csv", header = TRUE, row.names = 1)
jaccard[,] <- lapply(jaccard[,], as.numeric)
# Remplacer les points par des tirets bas dans les noms des colonnes de jaccard
colnames(jaccard) <- gsub("\\.", "_", colnames(jaccard))
# Définir l'ordre souhaité des noms de lignes
ordre_stations <- c("STATION 1", "STATION 2", "STATION 3", "STATION 4", "STATION 5")

# Réorganiser les lignes du dataframe en utilisant l'ordre défini
jaccard <- jaccard[ordre_stations, ]

# Afficher le dataframe réorganisé
print(jaccard)

comm <- jaccard


prunedphy <- prune.sample(comm, phylo)
prunedphy


# Réorganiser les colonnes de comm pour qu'elles correspondent à l'ordre des tip.label dans prunedphy
comm <- comm[, prunedphy$tip.label]
print(colnames(comm))


par(mfrow=c(2,3))
# Boucler sur chaque station pour créer un arbre avec des points pour les espèces présentes

stations <- c("STATION 1", "STATION 2", "STATION 3", "STATION 4", "STATION 5")

# Boucle pour chaque station dans l'ordre défini
for (station in stations) {   
  # Créer un vecteur binaire pour la présence d'espèces dans la station
  presence <- comm[station, ] == 1
  
  # Tracer l'arbre
  plot(prunedphy, show.tip.label = FALSE, main = paste(station), direction = "rightwards")
  
  # Ajouter des points bleus aux extrémités pour les espèces présentes
  tiplabels(pch = 16, col = "blue", tip = which(presence))
}


pd.result <- pd(comm, phylo, include.root = FALSE)

# Vérifiez les résultats
print(pd.result)



# Vérifiez si l'arbre est enraciné
if (!is.rooted(phylo)) {
  # Convertir en arbre bifilaire si nécessaire
  phylo_bifurcated <- multi2di(phylo, random = TRUE)
  
  # Enraciner l'arbre à une espèce spécifique (par exemple, le premier tip.label)
  phylo_rooted <- root(phylo_bifurcated, outgroup = phylo_bifurcated$tip.label[1], resolve.root = TRUE)
  
  # Maintenant, votre arbre est enraciné
  phylo <- phylo_rooted
}


# Calculer la diversité phylogénétique
pd.result <- pd(comm, phylo, include.root = TRUE)
print(pd.result)
print((pd.result$PD)/(pd.result$SR))


#####PD, MPD, MNTD avec SES
ses.pd.result <- ses.pd(comm, phylo, null.model = "taxa.labels", run = 1000)
ses.pd.result

phydist <- cophenetic(phylo)

ses.mpd.result <- ses.mpd(comm, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result

ses.mntd.result <- ses.mntd(comm, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result


####Beta div
comdist.result <- comdist(comm, phydist)
comdist.result
comdist.clusters <- hclust(comdist.result)
plot(comdist.clusters)

comdistnt.result <- comdistnt(comm, phydist)
comdistnt.result
comdistnt.clusters <- hclust(comdistnt.result)
plot(comdistnt.clusters)


####PD globale tt vert

# Initialiser une matrice de présence pour toutes les espèces de l'arbre phylo
comm_full <- matrix(1, nrow = 1, ncol = length(combined_tree$tip.label)) # Une seule ligne, toutes les espèces présentes
colnames(comm_full) <- combined_tree$tip.label
rownames(comm_full) <- "All_species"

# Calculer la diversité phylogénétique maximale
pd_result <- pd(comm_full, combined_tree, include.root = FALSE)
pd_result
pd_result$PD/pd_result$SR


######Unifrac - jaccard

unifrac <- unifrac(comm, phylo)
unifrac


#####Hill number

hill_numbers <-hillR:: hill_phylo(comm, tree = prunedphy, q = 0)  # q = 0 pour la richesse
hill_numbers
hill_parti <-hillR:: hill_phylo_parti(comm, tree = prunedphy, q = 0)  # q = 0 pour la richesse
hill_parti
hill_parti_pariwise <-hillR:: hill_phylo_parti_pairwise(comm, tree = prunedphy, q = 0)  # q = 0 pour la richesse
hill_parti_pariwise









############
# Calculer la matrice des distances phylogénétiques
phydist <- cophenetic(phylo)

# Initialiser des listes pour stocker les résultats VPD, SES et p-values
vpd_results <- list()
ses_results <- list()
p_values <- list()

# Nombre de répétitions pour le modèle nul
num_reps <- 1000

# Boucle pour chaque station
for (station in stations) {
  presence <- comm[station, , drop = FALSE] # Sous-matrice pour la station
  present_tips <- which(presence == 1) # Indices des espèces présentes
  
  if (length(present_tips) > 1) { # Vérifier qu'il y a au moins deux espèces
    # Extraire les distances pairwise pour les espèces présentes
    distances <- phydist[present_tips, present_tips]
    # Calculer la variance des distances
    vpd_observed <- var(distances[upper.tri(distances)]) # Variance des distances pairwise
    vpd_results[[station]] <- vpd_observed
    
    # Simuler la VPD
    null_vpd <- numeric(num_reps)
    for (i in 1:num_reps) {
      # Échantillonnage aléatoire d'espèces
      random_sample <- sample(colnames(comm), length(present_tips))
      random_tips <- which(colnames(comm) %in% random_sample)
      random_distances <- phydist[random_tips, random_tips]
      null_vpd[i] <- var(random_distances[upper.tri(random_distances)])
    }
    
    # Calculer la moyenne et l'écart-type des valeurs nulles
    null_mean <- mean(null_vpd)
    null_sd <- sd(null_vpd)
    
    # Calculer la SES
    ses <- (vpd_observed - null_mean) / null_sd
    ses_results[[station]] <- ses
    
    # Calculer la p-value
    p_value <- mean(null_vpd >= vpd_observed) # Test d'une différence supérieure
    p_values[[station]] <- p_value
  } else {
    vpd_results[[station]] <- NA # Pas assez d'espèces présentes
    ses_results[[station]] <- NA
    p_values[[station]] <- NA
  }
}

# Afficher les résultats VPD, SES et p-values
vpd_results
ses_results
p_values


library(ggplot2)
library(dplyr)
library(tidyr)

# Créer un DataFrame avec les valeurs SES
ses_pd <- ses.pd.result$pd.obs.z
ses_mpd <- ses.mpd.result$mpd.obs.z
ses_vpd  <- unlist(ses_results)


# Combiner les résultats dans un seul DataFrame
ses_df <- data.frame(
  Station = c("STATION 1", "STATION 2", "STATION 3", "STATION 4", "STATION 5"),
  SES = c(ses_pd, ses_mpd, ses_vpd),
  Biodiv_Index = rep(c("PD", "MPD", "VPD"), each = length(ses_pd))
)

# Convertir Station en facteur pour l'ordre
ses_df$Station <- factor(ses_df$Station, levels = c("STATION 1", "STATION 2", "STATION 3", "STATION 4", "STATION 5"))
ses_df
# Créer le bargraph
ggplot(ses_df, aes(x = Station, y = SES, fill = Biodiv_Index)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Valeurs SES pour PD, MPD et VPD", y = "Valeur SES", x = "Stations") +
  theme_minimal()


library(ape)
library(PhyloMeasures)

compute_metrics <- function(tree, comm, stations, num_reps = 1000) {
  phydist <- cophenetic(tree)  # Compute the phylogenetic distance matrix
  results <- list()
  
  for (station in stations) {
    presence <- comm[station, , drop = FALSE]
    present_tips <- which(presence == 1)
    
    if (length(present_tips) > 1) {
      # Prune the tree
      pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, colnames(comm)[which(presence == 0)]))
      
      # Compute metrics
      NTI <- mntd.query(tree = pruned_tree, matrix = presence, standardize = TRUE, reps = num_reps)
      NRI <- mpd.query(tree = pruned_tree, matrix = presence, standardize = TRUE, reps = num_reps)
      PD <- pd.query(tree = pruned_tree, matrix = presence, reps = num_reps)
      
      distances <- phydist[present_tips, present_tips]
      VPD <- var(distances[upper.tri(distances)])
      
      # Store results
      results[[station]] <- data.frame(SR = sum(presence), PD = PD, NTI = NTI, NRI = NRI, VPD = VPD)
    } else {
      results[[station]] <- NA  # Not enough species present
    }
  }
  return(results)
}

# Usage example:
# results <- compute_metrics(phylo, comm, stations)




#####Package  daijiang/lirrr
##Calcul de VPD

mvpd <- function(samp, dis, abundance.weighted = FALSE){
  N <- dim(samp)[1]
  mpd = vpd = numeric(N)
  samp_stad = vegan::decostand(samp, method = 'total', MARGIN = 1)
  
  for (i in 1:N) {
    # cat("row ", i)
    spp <- names(samp[i, samp[i, ] > 0])
    if (length(spp) > 1) {
      sample.dis <- dis[spp, spp]
      if (abundance.weighted) {
        sample.weights <- crossprod(as.matrix(samp_stad[i, spp, drop = FALSE]))
        wm = weighted.mean(sample.dis[lower.tri(sample.dis)], sample.weights[lower.tri(sample.weights)])
        mpd[i] <- wm
        vpd[i] = weighted.mean((sample.dis[lower.tri(sample.dis)] - wm)^2, 
                               sample.weights[lower.tri(sample.weights)]) * 
          (sum(lower.tri(sample.dis))/(sum(lower.tri(sample.dis)) - 1)) # n-1 instead of n
      }
      else {
        mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
        vpd[i] <- var(sample.dis[lower.tri(sample.dis)])
      }
    } else {
      mpd[i] = NA
      vpd[i] = NA
    }
  }
  data.frame(site = row.names(samp), mpd = mpd, vpd = vpd, stringsAsFactors = FALSE)
}


mvpd.result <- mvpd(samp = comm, dis = phydist, abundance.weighted = FALSE)
mvpd.result

#on check avec les valeurs de mpd calculées avec picante
ses.mpd.result




###On vole le code de la fonction ses.mpd de picante pour calculer les SES de VPD : 

ses.vpd <- function(samp, dis, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", 
                                              "phylogeny.pool", "independentswap", "trialswap"), 
                    abundance.weighted = FALSE, runs = 999, iterations = 1000) {
  dis <- as.matrix(dis)
  
  # Calculer les valeurs observées de VPD en extrayant la colonne `vpd` des résultats de `mvpd`
  vpd.obs <- mvpd(samp, dis, abundance.weighted = abundance.weighted)$vpd
  null.model <- match.arg(null.model)
  
  # Calculer les valeurs VPD pour chaque modèle nul
  vpd.rand <- switch(null.model,
                     taxa.labels = t(replicate(runs, mvpd(samp, taxaShuffle(dis), abundance.weighted=abundance.weighted)$vpd)),
                     richness = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="richness"), dis, abundance.weighted)$vpd)),
                     frequency = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="frequency"), dis, abundance.weighted)$vpd)),
                     sample.pool = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="richness"), dis, abundance.weighted)$vpd)),
                     phylogeny.pool = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="richness"), 
                                                             taxaShuffle(dis), abundance.weighted)$vpd)),
                     independentswap = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="independentswap", iterations), dis, abundance.weighted)$vpd)),
                     trialswap = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="trialswap", iterations), dis, abundance.weighted)$vpd))
  )
  
  # Calculer la moyenne et l'écart-type des valeurs simulées
  vpd.rand.mean <- apply(X = vpd.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
  vpd.rand.sd <- apply(X = vpd.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
  
  # Calculer les valeurs SES pour VPD
  vpd.obs.z <- (vpd.obs - vpd.rand.mean) / vpd.rand.sd
  vpd.obs.rank <- apply(X = rbind(vpd.obs, vpd.rand), MARGIN = 2, FUN = rank)[1, ]
  vpd.obs.rank <- ifelse(is.na(vpd.rand.mean), NA, vpd.obs.rank)
  
  data.frame(ntaxa = specnumber(samp), vpd.obs, vpd.rand.mean, vpd.rand.sd, vpd.obs.rank, 
             vpd.obs.z, vpd.obs.p = vpd.obs.rank / (runs + 1), runs = runs, row.names = row.names(samp))
}

ses.vpd.result <- ses.vpd(comm, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.vpd.result


#Pour vérifier que ces valeurs de ses sont bien calculées on va faire le test avec mpd :
ses.mvpd <- function(samp, dis, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", 
                                              "phylogeny.pool", "independentswap", "trialswap"), 
                    abundance.weighted = FALSE, runs = 999, iterations = 1000) {
  dis <- as.matrix(dis)
  
  # Calculer les valeurs observées de VPD en extrayant la colonne `vpd` des résultats de `mvpd`
  vpd.obs <- mvpd(samp, dis, abundance.weighted = abundance.weighted)$mpd
  null.model <- match.arg(null.model)
  
  # Calculer les valeurs VPD pour chaque modèle nul
  vpd.rand <- switch(null.model,
                     taxa.labels = t(replicate(runs, mvpd(samp, taxaShuffle(dis), abundance.weighted=abundance.weighted)$mpd)),
                     richness = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="richness"), dis, abundance.weighted)$mpd)),
                     frequency = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="frequency"), dis, abundance.weighted)$mpd)),
                     sample.pool = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="richness"), dis, abundance.weighted)$mpd)),
                     phylogeny.pool = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="richness"), 
                                                             taxaShuffle(dis), abundance.weighted)$mpd)),
                     independentswap = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="independentswap", iterations), dis, abundance.weighted)$mpd)),
                     trialswap = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="trialswap", iterations), dis, abundance.weighted)$mpd))
  )
  
  # Calculer la moyenne et l'écart-type des valeurs simulées
  vpd.rand.mean <- apply(X = vpd.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
  vpd.rand.sd <- apply(X = vpd.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
  
  # Calculer les valeurs SES pour VPD
  vpd.obs.z <- (vpd.obs - vpd.rand.mean) / vpd.rand.sd
  vpd.obs.rank <- apply(X = rbind(vpd.obs, vpd.rand), MARGIN = 2, FUN = rank)[1, ]
  vpd.obs.rank <- ifelse(is.na(vpd.rand.mean), NA, vpd.obs.rank)
  
  data.frame(ntaxa = specnumber(samp), vpd.obs, vpd.rand.mean, vpd.rand.sd, vpd.obs.rank, 
             vpd.obs.z, vpd.obs.p = vpd.obs.rank / (runs + 1), runs = runs, row.names = row.names(samp))
}


ses.mvpd.result <- ses.mvpd(comm, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mvpd.result

ses.mvd.result <- ses.mpd(comm, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result
####ça marche !!!!!!

#############################

library(ggplot2)
library(dplyr)
library(tidyr)

# Créer un DataFrame avec les valeurs SES
ses_pd <- ses.pd.result$pd.obs.z
ses_mpd <- ses.mpd.result$mpd.obs.z


# Combiner les résultats dans un seul DataFrame
ses_df <- data.frame(
  Station = c("STATION 1", "STATION 2", "STATION 3", "STATION 4", "STATION 5"),
  SES = c(ses_pd, ses_mpd, ses_vpd),
  Biodiv_Index = rep(c("PD", "MPD", "VPD"), each = length(ses_pd))
)

# Convertir Station en facteur pour l'ordre
ses_df$Station <- factor(ses_df$Station, levels = c("STATION 1", "STATION 2", "STATION 3", "STATION 4", "STATION 5"))

# Créer le bargraph
ggplot(ses_df, aes(x = Station, y = SES, fill = Biodiv_Index)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Valeurs SES pour PD, MPD et VPD", y = "Valeur SES", x = "Stations") +
  theme_minimal()


library(ggplot2)
library(dplyr)
library(factoextra)
library(ggrepel)

# Prepare data for PCA
pca_data <- ses_df %>%
  pivot_wider(names_from = Biodiv_Index, values_from = SES) %>%
  select(-Station) %>%
  scale() # Normalize data

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PCA scores for plotting
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Station <- ses_df$Station[1:length(pca_scores$PC1)]

# Plot PCA
library(factoextra)

# Use fviz_pca_biplot to plot PCA results
phylo.PCA <- fviz_pca_biplot(pca_result, 
                             geom = "point", 
                             label = "var", 
                             habillage = pca_scores$Station, 
                             labelsize = 4, 
                             repel = TRUE, 
                             legend.title = "Station", 
                             title = NULL) +
  geom_text_repel(aes(label = pca_scores$Station), size = 3) + # Add station labels
  theme(text = element_text(size = 10))

# Display plot
phylo.PCA


library(betapart)
# Example with a phylogenetic tree and community data
beta_uni <- phylo.beta.pair(comm, phylo, index.family = "jaccard")
beta_uni
# UniFrac Turnover
unifrac_turn <- beta_uni$phylo.beta.jtu
unifrac_turn
# UniFrac Phylogenetic Diversity
unifrac_pd <- beta_uni$phylo.beta.jne
unifrac_pd





#####################
##Functionnal diversity

traits_updated <- read.csv("traits_updated.csv")
#On garde les valeurs de Gf calculées précédemment
#Le truc c'est que ça ne considère que les poissons car limites Fishbase ...
#On va déjà essayé avec les poissons et on verra ce que ça vaut après


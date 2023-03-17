# Illustrate GRAMPA dataset in ternary plots

## Working directory
getwd()
setwd('~/.../')

## packages
install.packages("ggplot2")
install.packages('ggtern')
library(ggplot2)
library(ggtern)

## Load datasets
grampa_db <- read.csv('./Data/GRAMPA/grampa_pep2d.csv', header=TRUE)
grampa_db

pep2D_db <- read.csv("./Data/GRAMPA/pep2d_one.csv", header=TRUE)
pep2D_db

pdbe_db <- read.csv("./Data/GRAMPA/oneid_pdbe.csv", header=TRUE)
pdbe_db

#Subsetting
grampa_sub <- grampa_db[,4:6]
colnames(grampa_sub) <- c('helix_H', 'sheet_E', 'coil_C')

# Change null values from 0.0 to 0.5
grampa_sub[grampa_sub == 0.0] <- 0.5
grampa_sub

# Ternary plots
TernDens <- ggtern::ggtern(data = grampa_sub,
                           aes(x = sheet_E,
                               y = helix_H, 
                               z = coil_C),
                           aes(x,y,z))  + 
  #geom_point(size=0.4, colour='darkgrey', alpha=0.8) +
  stat_density_tern(geom='polygon', 
                    #color='grey',
                    #n=300,
                    bins=50,
                    expand = 1, h=0.1,
                    base='identity',
                    aes(fill   = ..level.., alpha = 0.5),
                    na.rm = TRUE) +
  #geom_point(size=0.4, colour='darkgrey', alpha=0.8) +
  #scale_fill_distiller(palette = 'RdYlBu') +
  scale_fill_continuous(low = 'darkolivegreen1', high='red', name = "density") +
  #scale_fill_viridis(option= "B", direction = -1) +
  theme_bw() +
  theme_showarrows() +
  theme_anticlockwise() +
  ggtitle('Density Map GRAMPA N=5980')

TernDens
#ggsave('./Figures/TernaryPlots/GRAMPA_density.png')

TernDens2 <- ggtern::ggtern(data = grampa_sub,
                           aes(x = sheet_E,
                               y = helix_H, 
                               z = coil_C),
                           aes(x,y,z))  + 
  #geom_point(size=0.4, colour='darkgrey', alpha=0.8) +
  stat_density_tern(geom='polygon', 
                    #color='black',
                    #n=300, 
                    bins=50,
                    expand = 1, h=0.15,
                    base='identity',
                    aes(fill   = ..level.., alpha = 0.5),
                    na.rm = TRUE) +
  geom_point(size=0.8, colour='darkgrey', alpha=0.8) +
  #scale_fill_distiller(palette = 'RdYlBu') +
  scale_fill_continuous(low = 'darkolivegreen1', high='red', name = "density", limits=c(0,25)) +
  #scale_fill_viridis(option= "B", direction = -1) +
  theme_bw() +
  theme_showarrows() +
  theme_anticlockwise() +
  ggtitle('Density Map E_coli (Navy) N=3367')

TernDens2
#ggsave('./Figures/TernaryPlots/GRAMPA_density_2.png')

#Overlapping ternary plots between PEP2D and PDBe geom points
Overlay <- ggtern(NULL,aes(x = sheet_E,y = helix_H,z = coil_C)) +  
  geom_point(data = pdbe_sub, size=1.4, colour='dodgerblue3', alpha=0.7) +
  geom_point(data = pep2D_sub, size=1.4, colour='darkslategray3', alpha=0.7) +
  theme_bw() +
  theme_showarrows() +
  theme_anticlockwise() +
  ggtitle('Distributions PDBe and PEP2D (N=261)')

Overlay
#ggsave('./Figures//TernaryPlots/PDBe&Pep2D_points.png')

#Create a common dataframe to PDBe and Pep2D dataframes
##Make sure both dataframes share a column (pdb_id)
pep2D_sub$pdb_id <- pep2D_db$peptide_ID
pdbe_sub$pdb_id <- pdbe_db$pdb_id

## Merge their subsets into a new dataframe, rename columns
struc_db <- merge(x=pdbe_sub,y=pep2D_sub, by = 'pdb_id')
colnames(struc_db) <- c('pdb_id', 'pdbe_helix_H', 'pdbe_sheet_E', 'pdbe_coil_C', 'pep2d_helix_H', 'pep2d_sheet_E', 'pep2d_coil_C' )
#write.csv(struc_db, file='./Data/GRAMPA/structures_pdbe_pep2d_HEC_values.csv')


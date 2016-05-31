
library(QFASA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(rCharts)


# Fatty Acid Names
fa.names = read.csv(file=system.file("exdata", "FAset.csv", package="QFASA"), as.is=TRUE)

# Predators
predators = read.csv(file=system.file("exdata", "predatorFAs.csv", package="QFASA")) %>%
    dplyr::select_(.dots = fa.names$FA) %>% 
    dplyr::mutate(id=row_number()) %>%
    tidyr::gather(FattyAcid, percentage, -id, factor_key=TRUE) %>%
    dplyr::group_by(id) %>% 
    dplyr::mutate(percentage=percentage/sum(percentage)) %>%
    tidyr::spread(FattyAcid, percentage) %>% 
    dplyr::ungroup() %>%
    dplyr::select(-id)

# Prey
prey <- read.csv(file=system.file("exdata", "preyFAs.csv", package="QFASA")) %>%
    dplyr::select_('Species', .dots = fa.names$FA) %>% # filter fatty acids
    dplyr::mutate(id=row_number()) %>%  # add a row id for sample grouping in long format
    tidyr::gather(FattyAcid, proportion, -id, -Species, factor_key=TRUE) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(proportion=proportion/sum(proportion)) %>% # normalize sample proportions 
    dplyr::group_by(Species, FattyAcid) %>% 
    dplyr::summarize(proportion=mean(proportion)) %>% # average sample proportions
    tidyr::spread(FattyAcid, proportion)

# Coefficents are the same for all predators in same group (species)
cal = read.csv(file=system.file("exdata", "CC.csv", package="QFASA"), as.is=TRUE) %>%
    dplyr::filter(FA %in% fa.names$FA) %>%
    transmute(FattyAcid = as.factor(FA), CalCoeff = CC)
cal.mat = replicate(dim(predators)[1], cal$CalCoeff)

# Account for fat content of prey species
FC <- read.csv(file=system.file("exdata", "preyFAs.csv", package="QFASA")) %>%
    dplyr::select(Species, lipid) %>%
    dplyr::group_by(Species) %>%
    dplyr::summarize(lipid=mean(lipid))


########################################
# KL Distance Measure
########################################
QFASA.KL = p.QFASA(as.matrix(predators), 
                as.matrix(prey %>% dplyr::ungroup() %>% dplyr::select(-Species)), 
                cal.mat, 
                dist.meas=1,# KL
                gamma=1, 
                FC$lipid, 
                start.val=rep(1, dim(prey)[1]), 
                fa.names$FA)

DietEst = QFASA.KL$'Diet Estimates'
colnames(DietEst) = prey$Species
DietEst.KL <- data.frame(DietEst, dist='KL')

########################################
# AIT Distance Measure
########################################
QFASA.AIT = p.QFASA(as.matrix(predators), 
                as.matrix(prey %>% dplyr::ungroup() %>% dplyr::select(-Species)), 
                cal.mat, 
                dist.meas=2, # AIT
                gamma=1, 
                FC$lipid, 
                start.val=rep(1, dim(prey)[1]), 
                fa.names$FA)

DietEst = QFASA.AIT$'Diet Estimates'
colnames(DietEst) = prey$Species
DietEst.AIT <- data.frame(DietEst, dist='AIT')


########################################
# CS Distance Measure
########################################
QFASA.CS = p.QFASA(as.matrix(predators), 
                as.matrix(prey %>% dplyr::ungroup() %>% dplyr::select(-Species)), 
                cal.mat, 
                dist.meas=3,# CS
                gamma=1, 
                FC$lipid, 
                start.val=rep(1, dim(prey)[1]), 
                fa.names$FA)

DietEst = QFASA.CS$'Diet Estimates'
colnames(DietEst) = prey$Species
DietEst.CS <- data.frame(DietEst, dist='CS')


########################################
# Combine results from models using different distance measures
########################################
# Diet Estimates
DietEst <- dplyr::rbind_list(DietEst.KL, DietEst.AIT, DietEst.CS) %>%
    tidyr::gather(prey, proportion, -dist) %>%
    dplyr::mutate(dist=as.factor(dist), prey=as.factor(prey)) %>%
    dplyr::group_by(dist, prey) %>%
    dplyr::summarize(proportion=mean(proportion)) %>%
    dplyr::ungroup()
    
ggplot(data=DietEst, aes(x=prey, y=proportion, fill=dist)) +
    geom_bar(stat='identity', position='dodge') +
    theme(text=element_text(size=8), axis.text.x  = element_text(angle=90, hjust = 1))

nPlot(data=DietEst, proportion ~ prey, group='dist', type = 'multiBarHorizontalChart')
nPlot(data=DietEst, proportion ~ dist, group='prey', type = 'multiBarHorizontalChart')



# Additional Measures: ModFAS
Add.meas.KL = plyr::ldply(QFASA.KL$'Additional Measures', data.frame) %>% 
    dplyr::select(contains('ModFAS')) %>%
    dplyr::mutate(dist='KL')
    
Add.meas.AIT = plyr::ldply(QFASA.AIT$'Additional Measures', data.frame) %>% 
    dplyr::select(contains('ModFAS')) %>%
    dplyr::mutate(dist='AIT')

Add.meas.CS = plyr::ldply(QFASA.CS$'Additional Measures', data.frame) %>% 
    dplyr::select(contains('ModFAS')) %>%
    dplyr::mutate(dist='CS')

Add.meas <- dplyr::rbind_list(Add.meas.KL, Add.meas.AIT, Add.meas.CS) %>%
    dplyr::mutate(dist=as.factor(dist)) %>%
    tidyr::gather(ModFAS, proportion, -dist, factor_key=TRUE) %>%
    dplyr::group_by(dist, ModFAS) %>%
    dplyr::summarize(proportion=mean(proportion)) %>%
    dplyr::ungroup()

ggplot(data=Add.meas, aes(x=ModFAS, y=proportion, fill=dist)) + 
    geom_bar(stat='identity', position='dodge') +
    theme(text=element_text(size=8), axis.text.x  = element_text(angle=90, hjust = 1))

nPlot(data=Add.meas, proportion ~ ModFAS, group='dist', type = 'multiBarHorizontalChart')
nPlot(data=Add.meas, proportion ~ dist, group='ModFAS', type = 'multiBarHorizontalChart')




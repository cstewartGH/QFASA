## server.R ##
library(shiny)
library(QFASA)
library(dplyr)
library(tidyr)
library(rCharts)


########################################
## Server
########################################
shinyServer(function(input, output) {

########################################
## Data
########################################
    getmodels <- eventReactive(input$goButton, {
        print("getmodels()")

        # Progress indicator
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Reading data...", value = 0)
        
        # Fatty Acid Names
        fa.names <- data.frame(FA=input$fattyacid, stringsAsFactors = FALSE)

        ## Predators
        data(predatorFAs)
        predators = predatorFAs %>%
            dplyr::select_(.dots = fa.names$FA) %>% 
            dplyr::mutate(id=row_number()) %>%
            tidyr::gather(FattyAcid, percentage, -id, factor_key=TRUE) %>%
            dplyr::group_by(id) %>% 
            dplyr::mutate(percentage=percentage/sum(percentage)) %>%
            tidyr::spread(FattyAcid, percentage) %>% 
            dplyr::ungroup() %>%
            dplyr::select(-id)

        ## Prey
        data(preyFAs)
        prey <- preyFAs %>%
            dplyr::select_('Species', .dots = fa.names$FA) %>% # filter fatty acids
            dplyr::mutate(id=row_number()) %>%  # add a row id for sample grouping in long format
            tidyr::gather(FattyAcid, proportion, -id, -Species, factor_key=TRUE) %>% 
            dplyr::group_by(id) %>% 
            dplyr::mutate(proportion=proportion/sum(proportion)) %>% # normalize sample proportions 
            dplyr::group_by(Species, FattyAcid) %>% 
            dplyr::summarize(proportion=mean(proportion)) %>% # average sample proportions
            tidyr::spread(FattyAcid, proportion) %>%
            dplyr::filter(Species %in% input$prey)

        ## Coefficents are the same for all predators in same group (species)
        data(CC)
        cal = CC %>%
            dplyr::filter(FA %in% fa.names$FA) %>%
            transmute(FattyAcid = as.factor(FA), CalCoeff = CC)
        cal.mat = replicate(dim(predators)[1], cal$CalCoeff)
        
        # Account for fat content of prey species
        FC <- preyFAs %>%
            dplyr::select(Species, lipid) %>%
            dplyr::group_by(Species) %>%
            dplyr::summarize(lipid=mean(lipid)) %>%
            dplyr::filter(Species %in% input$prey)


        ########################################
        # KL Distance Measure
        ########################################
        progress$set(message = "Reading running KL...", value = 1)
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
        progress$set(message = "Reading running AIT...", value = 2)
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
        progress$set(message = "Reading running CS...", value = 3)
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
        progress$set(message = "Combining data...", value = 4)
        DietEst <- dplyr::rbind_list(DietEst.KL, DietEst.AIT, DietEst.CS) %>%
            tidyr::gather(prey, proportion, -dist) %>%
            dplyr::mutate(dist=as.factor(dist), prey=as.factor(prey)) %>%
            dplyr::group_by(dist, prey) %>%
            dplyr::summarize(proportion=mean(proportion)) %>%
            dplyr::ungroup()

        # Additional Measures: PropDistCont
        Add.meas.KL = plyr::ldply(QFASA.KL$'Additional Measures', data.frame) %>% 
            dplyr::select(contains('PropDistCont')) %>%
            dplyr::mutate(dist='KL')
    
        Add.meas.AIT = plyr::ldply(QFASA.AIT$'Additional Measures', data.frame) %>% 
            dplyr::select(contains('PropDistCont')) %>%
            dplyr::mutate(dist='AIT')

        Add.meas.CS = plyr::ldply(QFASA.CS$'Additional Measures', data.frame) %>% 
            dplyr::select(contains('PropDistCont')) %>%
            dplyr::mutate(dist='CS')

        Add.meas <- dplyr::rbind_list(Add.meas.KL, Add.meas.AIT, Add.meas.CS) %>%
            dplyr::mutate(dist=as.factor(dist)) %>%
            tidyr::gather(PropDistCont, proportion, -dist, factor_key=TRUE) %>%
            dplyr::group_by(dist, PropDistCont) %>%
            dplyr::summarize(proportion=mean(proportion)) %>%
            dplyr::ungroup()
        
        return(list(DietEst=DietEst, AddMeas=Add.meas))
    })


########################################
## Visualization
########################################

    output$dietestbyprey <- renderChart({
        print("renderChart(dietestbyprey)")
        qfasa <- getmodels()
        p <- nPlot(data=qfasa$DietEst,
                   proportion ~ dist,
                   group='prey',
                   type = 'multiBarChart')
        p$chart(stacked='true')
        p$set(dom='dietestbyprey')
        return(p)
    })


    output$addmeasbyprey <- renderChart({
        print("renderChart(addmeasbyprey)")
        qfasa <- getmodels()
        p <- nPlot(data=qfasa$AddMeas,
                   proportion ~ dist,
                   group='PropDistCont',
                   type = 'multiBarChart')

        p$chart(stacked='true')
        p$set(dom='addmeasbyprey')
        p$addParams(title = 'Prop Dist Cont')
        return(p)
    })
    
})



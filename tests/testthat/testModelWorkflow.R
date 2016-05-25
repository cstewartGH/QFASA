test_that("Model Workflow Regression", {

    dist.meas = 1

    # Fatty Acids
    fa.set = as.vector(unlist(read.csv(file=system.file("exdata",
                                       "FAset.csv",
                                       package="QFASA"), as.is=TRUE)))
    print("fa.set")
    str(fa.set)

    # Predators
    predators = read.csv(file=system.file("exdata",
                                          "predatorFAs.csv",
                                          package="QFASA"), as.is=TRUE)
    tombstone.info = predators[,1:4]
    predator.matrix = predators[,5:(ncol(predators))]
    npredators = nrow(predator.matrix)
    print("predator.matrix")
    str(predator.matrix)

    # Prey
    prey = read.csv(file=system.file("exdata", "preyFAs.csv", package="QFASA"), as.is=TRUE)
    prey.sub = (prey[,4:(ncol(prey))])[fa.set]
    prey.sub = prey.sub/apply(prey.sub,1,sum) 
    group = as.vector(prey$Species) 
    prey.matrix = MEANmeth(cbind(group,prey.sub))
    print("prey.matrix")
    str(prey.matrix)
    
    FC = prey[, c(2,3)] 
    FC = as.vector(tapply(FC$lipid, FC$Species, mean, na.rm=TRUE))
    print("FC")
    str(FC)

    # Calibration Coefficients
    cal = read.csv(file=system.file("exdata", "CC.csv", package="QFASA"), as.is=TRUE)
    cal.vec = cal[,2]
    cal.mat = replicate(npredators, cal.vec)
    print("cal.mat")
    str(cal.mat)

    # Run QFASA
    Q = p.QFASA(predator.matrix,
                prey.matrix,
                cal.mat,
                dist.meas,
                gamma=1,
                FC,
                start.val = rep(1,nrow(prey.matrix)),
                fa.set)

    # Diet Estimates
    DietEst = Q$'Diet Estimates'
    DietEst = round(DietEst*100,digits=2)
    colnames(DietEst) = (as.vector(rownames(prey.matrix)))
    DietEst = cbind(tombstone.info,DietEst)
    print("DietEst")
    print(DietEst)

    # Check diet estimates
    DietEstCheck = read.csv(file=system.file("exdata", "DietEst.csv", package="QFASA"),
                            as.is=TRUE)
    expect_equal(DietEst, DietEstCheck)

    # Additional Measures
    AdditionalMeasures = plyr::ldply(Q$'Additional Measures', data.frame)
    print("AdditionalMeasures")
    print(AdditionalMeasures)
    
})

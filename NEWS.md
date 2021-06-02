---
output:
  html_document: default
  pdf_document: default
---
# News

## QFASA 1.1.0

* Made changes to DESCRIPTION including updating date, version, adding "cre" to Connie Stewart's role, removing "cre" from Justing Kamerman's role, adding Holly Steeves as an author, adding new imports and packages to LinkingTo. 
* Cleaned up p.QFASA: added reference, added to example, etc...
* Cleaned up pseudo.pred: changed description of parameters, added more explanation to "Value" and example.
* * Vignettes: 
  + Changed name of modelling workflow vingette to "QFASA_Workflow_Example" and updated it.  
  + Added a comment to the Parallel_Execution_for_Confidence_Intervals vignette.  
  + Added MUFASA_Worflow_Example vignette.
* beta.meths.CI is deprecated and replaced by conf.meth. Hid this function and bias.all from index. Made several corrections to conf.meth code.
* prey.cluster was missing code.  Added example and default values.
* Added functions to compute repeatability.
* Added functions to compute diet estimates and calibration coefficients simultaneously.
* Added functions to compute diet estimates via MUFASA (i.e. ML estimation)


Added the following functions:

## QFASA 1.0.3
* Added functions to test a difference between two independent samples of compositional data.
* Added functions to calculate confidence intervals on diet estimates.
* Added functions to perform hierarchical clustering on prey database.
* Added functions to calculate bias in diet estimates.
* Fixed: Add species names to diet estimates that are returned in p.QFASA
* Fixed: Set gamma=1 in p.QFASA
* Fixed: Add "Value" to p.QFASA function
* Fixed: Fix incorrect spelling of Shelley Lang and Chris Field
* Fixed: Package citation incorrect

## QFASA 1.0.2
* Fix bug whereby using data frames masks dimension incompatibilities of inputs.

## QFASA 1.0.1
* Fix documentation based on feddback from Connie Stewart. 
* Added separate tests for each distance measure.
* Fixed author list
* Added bibliographic references
* Add Shiny sample application


## QFASA 1.0.0
Initial release of code created in support of article

> Iverson, S. J., Field, C., Don Bowen, W. and Blanchard, W. (2004)
> QUANTITATIVE FATTY ACID SIGNATURE ANALYSIS: A NEW METHOD OF
> ESTIMATING PREDATOR DIETS Ecological Monographs 74: 211–235.
> doi:10.1890/02-4105

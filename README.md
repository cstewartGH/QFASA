---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



[![Travis-CI Build Status](https://travis-ci.org/justinkamerman/QFASA.svg?branch=master)](https://travis-ci.org/justinkamerman/QFASA)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/QFASA)](http://cran.r-project.org/package=QFASA)


# Overview
Accurate estimates of the diets of predators are required in many
areas of ecology, but for many species current methods are imprecise,
limited to the last meal, and often biased. The diversity of fatty
acids and their patterns in organisms, coupled with the narrow
limitations on their biosynthesis, properties of digestion in
monogastric animals, and the prevalence of large storage reservoirs of
lipid in many predators, led us to propose the use of quantitative
fatty acid signature analysis (QFASA) to study predator diets.


# Installing
## Via GitHub

```r
devtools::install_github('justinkamerman/QFASA')
```

## Via CRAN

```r
install.packages('QFASA')
```

# Load Package


```r
library(QFASA)
```

# Modeling Inputs
Prior to starting make sure that:

* Fatty acid names in all files are the same (contain the exact same
  numbers/characters and punctuation)
* There are no fatty acids in the prey file that do not appear in the
  predator file and visa versa


## Distance Measure
Choose from one of three distance measures:  
1=KL (Kullback-Leibler) || 2=AIT (Aitchison) || 3=CSD (Chi-Squared) 


```r
dist.meas=1
```

## Fatty Acid Set
* This is the list of FAs to be used in the modelling.
* The simplest alternative is to load a .csv file which contains a
  single column with a header row and the names of the fatty acids
  listed below (see example file __"FAset.csv"__).
* A more complicated alternative is to load a .csv file with the full
  set of FAs and then add code to subset the FAs you wish to use from
  that set --> this alternative is useful if you are planning to test
  multiple FA sets.
* Regardless of how you load the FA set it must be converted to a
  vector.


```r
data(FAset)
fa.set = as.vector(unlist(FAset))
```

## Matrix of Predator FA signatures
* The FA signatures in the originating .csv file should be in
  percentages not proportions (i.e. adding to ~100 not 1).
  
* Each predator signature is a row with the FAs in columns (see
  example file __"predatorFAs.csv"__).
  
* the FA signatures are subsetted for the chosen FA set (created
  above) and renormalized during the modelling so there is no need to
  subset and/or renormalize prior to loading the .csv file or running
  p.QFASA BUT make sure that the the same FAs appear in the predator
  and prey files (if a FA appears in one but not the other the code
  will give you an error).

* Unlike the original QFASApack code the predator FA .csv file can
  contain as much tombstone data in columns as you wish but the
  predator FA signatures must be extracted as a separate input in
  order to run in p.QFASA. For example: in the code below the predator
  .csv file ("predatorFAs.csv") has 4 tombstone columns (SampleCode,
  AnimalCode, SampleGroup, Biopsy). Prior to running QFASA the
  tombstone (columns 1-4) and FA data (columns 5 onward) are each
  extracted from the original data frame. The FA data becaomes the the
  predator.matrix (which is passed to p.QFASA) and the tombstone data
  is retained so that it can be recombined with the model output later
  on.
  

```r
data(predatorFAs)
tombstone.info = predatorFAs[,1:4]
predator.matrix = predatorFAs[,5:(ncol(predatorFAs))]

# number of predator FA signatures this is used to create the matrix of CC values (see section 6 below)
npredators = nrow(predator.matrix)
```


## Matrix of Prey FA signatures
* The FA signatures in the originating .csv file should be in
  percentages not proportions (i.e. adding to ~100 not 1).
  
* The prey file should contain all of the individual FA signatures of
  the prey and their lipid contents (where appropriate) - a matrix of
  the mean values for the FAs (prey.matrix) by the designated prey
  modelling group is then calculated using the MEANmeth function
  loaded above.
  
* Like the predator .csv file you can have as many tombstone data
  columns as required but there must be at least one column that
  identifies the modelling group to be used (in the example file used
  below __"preyFAs.csv"__ it is the "Species" column).
  
* Unlike the predator data, the prey data is not subsetted and
  renomalized during the modelling so the prey file needs to be
  subsetted for the desired FA set (created above) and renormalized to
  sum to 1 prior to calculating the mean values (see code
  below). Example: in the code below the "preyFAs.csv" file has 3
  tombstone columns. The full FA set is extracted from the data frame
  (columns 4 onward), subsetted for the FA set in use and then
  renormalized over 1. The modelling group names (the "Species" column
  in this case) is then added back to the subsetted and renormalized
  data (as the first column) and the average values calculated using
  the MEANmeth function. Note that for the MEANmeth function to work
  the modelling group name must be in the first column.
    

```r
#full file
data(preyFAs)

#extract prey FA only from data frame and subset them for the FA set designated above
prey.sub=(preyFAs[,4:(ncol(preyFAs))])[fa.set]

#renormalize over 1
prey.sub=prey.sub/apply(prey.sub,1,sum) 

#extract the modelling group names from the full file
group=as.vector(preyFAs$Species)

#add modelling group names to the subsetted and renormalized FAs
prey.matrix=cbind(group,prey.sub)

#create an average value for the FA signature for each designated modelling group
prey.matrix=MEANmeth(prey.matrix) 
```

## Prey Lipid Content
* mean lipid content by modelling group is calculated from the full
  prey file using the modelling group as a summary variable (see code
  below).
* **Note:** if no lipid content correction is going to be applied then
  a vector of '1's of length equal to the number of modelling groups
  is used as the vector instead i.e. FC=rep(1,nrow(prey.matrix))


```r
#numbers are the column which identifies the modelling group, and the column which contains the lipid contents
FC = preyFAs[,c(2,3)] 
FC = as.vector(tapply(FC$lipid,FC$Species,mean,na.rm=TRUE))
```


## Calibration Coefficients
* Originating .csv file should contain 2 columns (with headers). The
  first contains the FA names, the second the value of the CC for each
  FA (see example file __"CC.csv"__).
* __IMPORTANT:__ the FAs in the CC.csv file __MUST__ be exactly the
  same as the FAs in the originating predator.csv file __AND__ they
  __MUST__ BE IN THE __*EXACT*__ SAME ORDER.
  

```r
data(CC)
cal.vec = CC[,2]
cal.mat = replicate(npredators, cal.vec)
```


# Run QFASA


```r
Q = p.QFASA(predator.matrix, prey.matrix, cal.mat, dist.meas, gamma=1, FC, start.val=rep(1,nrow(prey.matrix)), fa.set)
```

## p.QFASA Output
The QFASA output is a list with 2 components:

* Diet Estimates
* Additional Measures

### Diet Estimates
This is a matrix of the diet estimate for each predator (by rows, in
the same order as the input file) by the modelling groups (by column,
in the same order as the prey.matrix file). The estimates are
expressed as a proportion (they will sum to 1). In the code below the
Diet Estimate matrix is extracted from the QFASA output and the
modelling group identities and predator tombstone data (created above)
are added to the matrix:


```r
DietEst = Q$'Diet Estimates'

#estimates changed from proportions to percentages
DietEst = round(DietEst*100,digits=2)
colnames(DietEst) = (as.vector(rownames(prey.matrix)))
DietEst = cbind(tombstone.info,DietEst)
knitr::kable(DietEst)
```



|SampleCode |AnimalCode |SampleGroup |Biopsy | capelin| coho| eulachon| herring| mackerel| pilchard| pollock| sandlance| squid| surfsmelt_lg| surfsmelt_s|
|:----------|:----------|:-----------|:------|-------:|----:|--------:|-------:|--------:|--------:|-------:|---------:|-----:|------------:|-----------:|
|3-01A      |P031       |T           |A      |   35.30|    0|        0|   44.90|     9.25|     2.36|       0|      8.19|     0|         0.00|           0|
|3-01B      |P031       |T1          |B      |   43.06|    0|        0|   44.79|     3.72|     4.25|       0|      4.18|     0|         0.00|           0|
|3-01C      |P031       |T1          |C      |   50.15|    0|        0|   34.34|     6.14|     3.97|       0|      5.40|     0|         0.00|           0|
|3-02A      |P032       |T           |A      |   37.80|    0|        0|   47.14|     1.42|     5.04|       0|      8.60|     0|         0.00|           0|
|3-02B      |P032       |T1          |B      |   39.86|    0|        0|   45.50|     3.69|     5.77|       0|      5.18|     0|         0.00|           0|
|3-02C      |P032       |T1          |C      |   47.99|    0|        0|   35.59|     4.47|     5.39|       0|      6.55|     0|         0.00|           0|
|3-04A      |P034       |T           |A      |   41.29|    0|        0|   46.33|     8.24|     2.65|       0|      1.49|     0|         0.00|           0|
|3-04B      |P034       |T2          |B      |   31.96|    0|        0|   28.00|     2.95|     2.40|       0|     34.70|     0|         0.00|           0|
|3-04C      |P034       |T2          |C      |   21.07|    0|        0|   14.45|     0.27|     2.00|       0|     37.03|     0|        25.17|           0|
|3-05A      |P035       |T           |A      |   28.91|    0|        0|   63.50|     0.00|     3.93|       0|      3.66|     0|         0.00|           0|

### Additional Measures
This is a list of lists where each list (one per predator) is itself a list of four outputs:

* __ModFAS__: the value of the modelled FA (i.e. after CCs have been
  applied and the FA subsetted and renormalised over the designated FA
  set). These are expressed as proportions (they will sum to 1).
  
* __DistCont__: the contribution of each FA to the final minimized distance.

* __PropDistCont__: the contribution of each FA to the final minimized
  distance as a proportion of the total.
  
* __MinDist__: the final minimized distance in the code below the
  'ldply' function from the plyr package is used to compile the lists
  within 'Additional Measures' into a data frame with one row per
  predator (in the same order as the input predator matrix) and the
  values for each of the 4 lists arranged into columns. The 'ldply'
  function automatically names the columns of the data frame with a
  concatenation of the originating list name and the FA name so that
  the 4 sets of outputs can be easily identified within the data
  frame.
  

```r
# plyr package
library(plyr)
Add.meas = ldply(Q$'Additional Measures', data.frame)
knitr::kable(Add.meas)
```



| ModFAS.c14.0| ModFAS.c16.0| ModFAS.c16.1w7| ModFAS.c16.2w6| ModFAS.c16.2w4| ModFAS.c16.3w6| ModFAS.c17.0| ModFAS.c16.3w4| ModFAS.c16.4w3| ModFAS.c16.4w1| ModFAS.c18.0| ModFAS.c18.1w9| ModFAS.c18.1w7| ModFAS.c18.2w6| ModFAS.c18.2w4| ModFAS.c18.3w6| ModFAS.c18.3w4| ModFAS.c18.3w3| ModFAS.c18.3w1| ModFAS.c18.4w3| ModFAS.c18.4w1| ModFAS.c20.1w11| ModFAS.c20.1w9| ModFAS.c20.1w7| ModFAS.c20.2w6| ModFAS.c20.3w6| ModFAS.c20.4w6| ModFAS.c20.3w3| ModFAS.c20.4w3| ModFAS.c20.5w3| ModFAS.c22.1w11| ModFAS.c22.1w9| ModFAS.c22.1w7| ModFAS.c21.5w3| ModFAS.c22.4w6| ModFAS.c22.5w6| ModFAS.c22.4w3| ModFAS.c22.5w3| ModFAS.c22.6w3| DistCont.c14.0| DistCont.c16.0| DistCont.c16.1w7| DistCont.c16.2w6| DistCont.c16.2w4| DistCont.c16.3w6| DistCont.c17.0| DistCont.c16.3w4| DistCont.c16.4w3| DistCont.c16.4w1| DistCont.c18.0| DistCont.c18.1w9| DistCont.c18.1w7| DistCont.c18.2w6| DistCont.c18.2w4| DistCont.c18.3w6| DistCont.c18.3w4| DistCont.c18.3w3| DistCont.c18.3w1| DistCont.c18.4w3| DistCont.c18.4w1| DistCont.c20.1w11| DistCont.c20.1w9| DistCont.c20.1w7| DistCont.c20.2w6| DistCont.c20.3w6| DistCont.c20.4w6| DistCont.c20.3w3| DistCont.c20.4w3| DistCont.c20.5w3| DistCont.c22.1w11| DistCont.c22.1w9| DistCont.c22.1w7| DistCont.c21.5w3| DistCont.c22.4w6| DistCont.c22.5w6| DistCont.c22.4w3| DistCont.c22.5w3| DistCont.c22.6w3| PropDistCont.c14.0| PropDistCont.c16.0| PropDistCont.c16.1w7| PropDistCont.c16.2w6| PropDistCont.c16.2w4| PropDistCont.c16.3w6| PropDistCont.c17.0| PropDistCont.c16.3w4| PropDistCont.c16.4w3| PropDistCont.c16.4w1| PropDistCont.c18.0| PropDistCont.c18.1w9| PropDistCont.c18.1w7| PropDistCont.c18.2w6| PropDistCont.c18.2w4| PropDistCont.c18.3w6| PropDistCont.c18.3w4| PropDistCont.c18.3w3| PropDistCont.c18.3w1| PropDistCont.c18.4w3| PropDistCont.c18.4w1| PropDistCont.c20.1w11| PropDistCont.c20.1w9| PropDistCont.c20.1w7| PropDistCont.c20.2w6| PropDistCont.c20.3w6| PropDistCont.c20.4w6| PropDistCont.c20.3w3| PropDistCont.c20.4w3| PropDistCont.c20.5w3| PropDistCont.c22.1w11| PropDistCont.c22.1w9| PropDistCont.c22.1w7| PropDistCont.c21.5w3| PropDistCont.c22.4w6| PropDistCont.c22.5w6| PropDistCont.c22.4w3| PropDistCont.c22.5w3| PropDistCont.c22.6w3|   MinDist|
|------------:|------------:|--------------:|--------------:|--------------:|--------------:|------------:|--------------:|--------------:|--------------:|------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|---------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|---------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|----------------:|----------------:|----------------:|----------------:|--------------:|----------------:|----------------:|----------------:|--------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|-----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|-----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|------------------:|------------------:|--------------------:|--------------------:|--------------------:|--------------------:|------------------:|--------------------:|--------------------:|--------------------:|------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|---------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|---------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|--------------------:|---------:|
|    0.0547664|    0.1853966|      0.0701887|      0.0009476|      0.0017854|      0.0069045|    0.0056682|      0.0052232|      0.0013693|      0.0093276|    0.0232666|      0.1480562|      0.0370854|      0.0096465|      0.0015922|      0.0012739|      0.0008316|      0.0047609|      0.0007809|      0.0141435|      0.0014503|       0.0178365|      0.0650470|      0.0040928|      0.0017300|      0.0007539|      0.0061731|      0.0004931|      0.0038550|      0.1147319|       0.0820286|      0.0074502|      0.0022416|      0.0040179|      0.0003608|      0.0011030|      0.0004209|      0.0110993|      0.0920990|      0.0088791|      0.0150925|        0.0001173|         6.05e-05|        0.0002396|        0.0000650|      0.0003814|        0.0002666|        0.0000305|        0.0005486|      0.0001031|        0.0183438|        0.0000862|        0.0000401|         2.07e-05|          2.0e-07|        0.0001591|        0.0000213|        0.0000010|        0.0013139|         7.90e-06|         0.0031149|        0.0155393|        0.0003308|        0.0001363|         4.91e-05|        0.0001160|        0.0001174|        0.0000440|        0.0006008|         0.0212495|        0.0000001|        0.0035884|        0.0002199|         1.26e-05|        0.0001506|        0.0001160|        0.0057975|        0.0147840|          0.0794585|          0.1350617|            0.0010495|            0.0005416|            0.0021443|            0.0005815|          0.0034132|            0.0023854|            0.0002733|            0.0049093|          0.0009225|            0.1641565|            0.0007714|            0.0003593|            0.0001856|             1.60e-06|            0.0014239|            0.0001909|            0.0000086|            0.0117581|            0.0000703|             0.0278749|            0.1390601|            0.0029605|            0.0012198|            0.0004395|            0.0010377|            0.0010502|            0.0003934|            0.0053766|             0.1901595|            0.0000006|            0.0321122|            0.0019676|            0.0001128|            0.0013478|            0.0010383|            0.0518812|            0.1323003| 0.1117455|
|    0.0550460|    0.1841167|      0.0720054|      0.0009930|      0.0016304|      0.0070830|    0.0053360|      0.0056322|      0.0012572|      0.0104616|    0.0233679|      0.1478895|      0.0367106|      0.0093263|      0.0016973|      0.0013089|      0.0008664|      0.0043575|      0.0007758|      0.0138342|      0.0015588|       0.0142726|      0.0682818|      0.0042321|      0.0016245|      0.0007779|      0.0061654|      0.0004480|      0.0038504|      0.1163810|       0.0814998|      0.0076830|      0.0022185|      0.0041182|      0.0003539|      0.0010787|      0.0004058|      0.0119736|      0.0893801|      0.0063032|      0.0171896|        0.0036939|         2.20e-06|        0.0003432|        0.0002371|      0.0004968|        0.0003633|        0.0000041|        0.0010276|      0.0000137|        0.0178350|        0.0001426|        0.0005043|         1.04e-05|          2.2e-06|        0.0002163|        0.0001160|        0.0000737|        0.0011233|         9.90e-06|         0.0023613|        0.0139278|        0.0003984|        0.0000013|         7.77e-05|        0.0003517|        0.0000505|        0.0002181|        0.0000005|         0.0244185|        0.0000336|        0.0032871|        0.0001371|         5.53e-05|        0.0001587|        0.0001128|        0.0064767|        0.0121005|          0.0553517|          0.1509500|            0.0324377|            0.0000191|            0.0030137|            0.0020822|          0.0043629|            0.0031900|            0.0000363|            0.0090237|          0.0001206|            0.1566182|            0.0012518|            0.0044287|            0.0000910|             1.94e-05|            0.0018996|            0.0010187|            0.0006473|            0.0098644|            0.0000866|             0.0207353|            0.1223068|            0.0034983|            0.0000115|            0.0006822|            0.0030885|            0.0004436|            0.0019155|            0.0000041|             0.2144305|            0.0002949|            0.0288652|            0.0012042|            0.0004860|            0.0013933|            0.0009902|            0.0568754|            0.1062608| 0.1138759|
|    0.0560279|    0.1789890|      0.0723739|      0.0009591|      0.0017444|      0.0067273|    0.0049163|      0.0053720|      0.0012702|      0.0101304|    0.0228679|      0.1378388|      0.0355206|      0.0096423|      0.0016793|      0.0012881|      0.0008372|      0.0043167|      0.0008086|      0.0137151|      0.0015843|       0.0152122|      0.0732648|      0.0046766|      0.0016630|      0.0007612|      0.0059978|      0.0004516|      0.0039308|      0.1163864|       0.0861903|      0.0086955|      0.0022833|      0.0040713|      0.0003448|      0.0010633|      0.0004023|      0.0123792|      0.0936161|      0.0017016|      0.0107110|        0.0061674|         6.90e-06|        0.0004523|        0.0002586|      0.0004484|        0.0004357|        0.0000136|        0.0008326|      0.0000424|        0.0126327|        0.0000073|        0.0009051|         1.39e-05|          1.7e-06|        0.0000459|        0.0001008|        0.0000944|        0.0012928|         3.40e-06|         0.0029103|        0.0144433|        0.0007050|        0.0000070|         1.93e-05|        0.0000384|        0.0000488|        0.0002088|        0.0004433|         0.0285698|        0.0000262|        0.0046912|        0.0000356|         1.12e-05|        0.0000781|        0.0000473|        0.0052183|        0.0080927|          0.0167216|          0.1052544|            0.0606055|            0.0000674|            0.0044448|            0.0025410|          0.0044065|            0.0042818|            0.0001335|            0.0081813|          0.0004165|            0.1241387|            0.0000720|            0.0088944|            0.0001370|             1.68e-05|            0.0004512|            0.0009907|            0.0009272|            0.0127039|            0.0000333|             0.0285992|            0.1419305|            0.0069274|            0.0000687|            0.0001898|            0.0003778|            0.0004791|            0.0020517|            0.0043557|             0.2807480|            0.0002577|            0.0460989|            0.0003498|            0.0001098|            0.0007675|            0.0004645|            0.0512792|            0.0795249| 0.1017630|
|    0.0549324|    0.1886819|      0.0713946|      0.0010226|      0.0017054|      0.0071749|    0.0053501|      0.0058160|      0.0012653|      0.0109474|    0.0243862|      0.1504301|      0.0368751|      0.0092908|      0.0017885|      0.0013530|      0.0008909|      0.0046221|      0.0008388|      0.0148194|      0.0016100|       0.0128313|      0.0632252|      0.0039940|      0.0016646|      0.0007920|      0.0063065|      0.0004683|      0.0040066|      0.1192530|       0.0750986|      0.0069816|      0.0021633|      0.0042342|      0.0003746|      0.0011081|      0.0004206|      0.0121530|      0.0897284|      0.0130936|      0.0188271|        0.0015858|         4.30e-05|        0.0002297|        0.0000408|      0.0005260|        0.0001712|        0.0000203|        0.0002750|      0.0005548|        0.0239421|        0.0000178|        0.0000801|         1.40e-06|          9.1e-06|        0.0001440|        0.0002830|        0.0000192|        0.0001850|         1.30e-05|         0.0016961|        0.0158307|        0.0005330|        0.0001075|         0.00e+00|        0.0000122|        0.0000066|        0.0001355|        0.0012861|         0.0139420|        0.0000800|        0.0025324|        0.0002661|         1.00e-07|        0.0001324|        0.0001049|        0.0080911|        0.0131539|          0.1109885|          0.1595887|            0.0134418|            0.0003643|            0.0019472|            0.0003454|          0.0044584|            0.0014508|            0.0001722|            0.0023313|          0.0047031|            0.2029459|            0.0001513|            0.0006794|            0.0000121|             7.71e-05|            0.0012205|            0.0023990|            0.0001630|            0.0015682|            0.0001103|             0.0143772|            0.1341896|            0.0045183|            0.0009115|            0.0000000|            0.0001030|            0.0000559|            0.0011484|            0.0109018|             0.1181795|            0.0006784|            0.0214662|            0.0022555|            0.0000010|            0.0011227|            0.0008889|            0.0685845|            0.1114991| 0.1179728|
|    0.0551564|    0.1860906|      0.0718531|      0.0010309|      0.0016327|      0.0072122|    0.0052826|      0.0059086|      0.0012421|      0.0113276|    0.0241941|      0.1489066|      0.0369172|      0.0092367|      0.0018067|      0.0013635|      0.0009050|      0.0044098|      0.0008202|      0.0144234|      0.0016661|       0.0139787|      0.0644588|      0.0040255|      0.0016226|      0.0008094|      0.0063120|      0.0004537|      0.0040279|      0.1191223|       0.0771774|      0.0072574|      0.0021436|      0.0042723|      0.0003667|      0.0010828|      0.0004133|      0.0126588|      0.0884308|      0.0075153|      0.0150745|        0.0012235|         2.24e-05|        0.0003131|        0.0002711|      0.0006123|        0.0003867|        0.0000006|        0.0008265|      0.0002306|        0.0167278|        0.0000755|        0.0004254|         3.30e-06|          1.1e-06|        0.0002367|        0.0001534|        0.0000901|        0.0008256|         2.25e-05|         0.0018837|        0.0120798|        0.0005337|        0.0000055|         3.60e-06|        0.0002554|        0.0000568|        0.0001010|        0.0000000|         0.0186703|        0.0000472|        0.0022597|        0.0002001|         1.65e-05|        0.0001432|        0.0000487|        0.0077364|        0.0103954|          0.0755503|          0.1515409|            0.0122994|            0.0002253|            0.0031475|            0.0027254|          0.0061557|            0.0038871|            0.0000064|            0.0083088|          0.0023185|            0.1681610|            0.0007593|            0.0042765|            0.0000331|             1.09e-05|            0.0023794|            0.0015418|            0.0009061|            0.0082998|            0.0002260|             0.0189360|            0.1214356|            0.0053647|            0.0000554|            0.0000359|            0.0025671|            0.0005709|            0.0010153|            0.0000000|             0.1876884|            0.0004741|            0.0227160|            0.0020120|            0.0001655|            0.0014394|            0.0004898|            0.0777722|            0.1045025| 0.0994748|
|    0.0560505|    0.1817995|      0.0723815|      0.0010029|      0.0017353|      0.0068926|    0.0048537|      0.0056954|      0.0012426|      0.0110586|    0.0237688|      0.1398926|      0.0357112|      0.0094975|      0.0017979|      0.0013444|      0.0008776|      0.0043657|      0.0008541|      0.0143584|      0.0016866|       0.0138687|      0.0694806|      0.0044664|      0.0016535|      0.0007911|      0.0061452|      0.0004535|      0.0040879|      0.1193432|       0.0810038|      0.0081787|      0.0022098|      0.0042271|      0.0003589|      0.0010720|      0.0004093|      0.0130023|      0.0923806|      0.0023871|      0.0087871|        0.0021140|         5.80e-06|        0.0003958|        0.0003019|      0.0005416|        0.0003837|        0.0000211|        0.0005778|      0.0004160|        0.0114337|        0.0000748|        0.0007350|         6.00e-06|          3.7e-06|        0.0000521|        0.0001235|        0.0002417|        0.0008650|         2.36e-05|         0.0018387|        0.0112123|        0.0005756|        0.0000040|         2.06e-05|        0.0000082|        0.0000061|        0.0000562|        0.0005961|         0.0189813|        0.0000657|        0.0036601|        0.0000206|         5.07e-05|        0.0000651|        0.0000114|        0.0064457|        0.0064989|          0.0299852|          0.1103792|            0.0265544|            0.0000726|            0.0049721|            0.0037927|          0.0068028|            0.0048198|            0.0002655|            0.0072581|          0.0052258|            0.1436248|            0.0009394|            0.0092332|            0.0000759|             4.70e-05|            0.0006548|            0.0015511|            0.0030362|            0.0108653|            0.0002962|             0.0230963|            0.1408430|            0.0072303|            0.0000503|            0.0002586|            0.0001032|            0.0000768|            0.0007066|            0.0074877|             0.2384340|            0.0008247|            0.0459761|            0.0002588|            0.0006369|            0.0008178|            0.0001431|            0.0809679|            0.0816362| 0.0796084|
|    0.0547001|    0.1817555|      0.0715331|      0.0009510|      0.0015960|      0.0070167|    0.0056379|      0.0053222|      0.0013040|      0.0095083|    0.0224540|      0.1486681|      0.0371563|      0.0094265|      0.0015574|      0.0012542|      0.0008343|      0.0043098|      0.0007035|      0.0130305|      0.0014459|       0.0174640|      0.0699567|      0.0042452|      0.0016257|      0.0007582|      0.0060870|      0.0004500|      0.0036742|      0.1126239|       0.0867604|      0.0079380|      0.0022620|      0.0039682|      0.0003395|      0.0010656|      0.0003999|      0.0112330|      0.0889833|      0.0047755|      0.0174695|        0.0036603|         6.00e-06|        0.0001590|        0.0001289|      0.0007356|        0.0004164|        0.0000281|        0.0012425|      0.0003970|        0.0201024|        0.0002863|        0.0000133|         3.82e-05|          6.6e-06|        0.0001921|        0.0000508|        0.0000030|        0.0010064|         5.00e-07|         0.0032654|        0.0156070|        0.0005452|        0.0000089|         2.14e-05|        0.0001790|        0.0000451|        0.0000147|        0.0000077|         0.0291126|        0.0000320|        0.0045975|        0.0001798|         4.68e-05|        0.0001565|        0.0001103|        0.0093666|        0.0164309|          0.0366090|          0.1339216|            0.0280598|            0.0000456|            0.0012191|            0.0009878|          0.0056393|            0.0031924|            0.0002151|            0.0095252|          0.0030435|            0.1541055|            0.0021948|            0.0001017|            0.0002928|             5.09e-05|            0.0014730|            0.0003898|            0.0000228|            0.0077150|            0.0000035|             0.0250325|            0.1196438|            0.0041797|            0.0000683|            0.0001644|            0.0013720|            0.0003457|            0.0001125|            0.0000588|             0.2231776|            0.0002455|            0.0352442|            0.0013786|            0.0003587|            0.0012000|            0.0008454|            0.0718042|            0.1259598| 0.1304458|
|    0.0563824|    0.1968077|      0.0671719|      0.0009391|      0.0026376|      0.0061370|    0.0047962|      0.0048937|      0.0015236|      0.0090623|    0.0264545|      0.1343605|      0.0343480|      0.0105995|      0.0018357|      0.0013559|      0.0008015|      0.0063743|      0.0011786|      0.0189172|      0.0015480|       0.0120723|      0.0547957|      0.0042769|      0.0021356|      0.0007097|      0.0062515|      0.0006337|      0.0046436|      0.1262907|       0.0638349|      0.0066688|      0.0022629|      0.0042235|      0.0004413|      0.0012529|      0.0004936|      0.0113655|      0.1095215|      0.0036647|      0.0099124|        0.0014762|         3.04e-05|        0.0003230|        0.0000780|      0.0000663|        0.0001142|        0.0000072|        0.0000934|      0.0000010|        0.0132954|        0.0001011|        0.0004996|         3.46e-05|          5.0e-07|        0.0001643|        0.0010005|        0.0000480|        0.0002908|         1.76e-05|         0.0015343|        0.0150237|        0.0002299|        0.0000085|         2.14e-05|        0.0000647|        0.0000328|        0.0003402|        0.0000457|         0.0094349|        0.0005935|        0.0038285|        0.0002098|         4.90e-06|        0.0000490|        0.0001888|        0.0051552|        0.0066288|          0.0491150|          0.1328495|            0.0197849|            0.0004070|            0.0043294|            0.0010451|          0.0008882|            0.0015309|            0.0000971|            0.0012514|          0.0000134|            0.1781888|            0.0013554|            0.0066959|            0.0004638|             7.10e-06|            0.0022024|            0.0134090|            0.0006430|            0.0038974|            0.0002360|             0.0205638|            0.2013518|            0.0030816|            0.0001140|            0.0002874|            0.0008667|            0.0004402|            0.0045593|            0.0006120|             0.1264499|            0.0079545|            0.0513105|            0.0028117|            0.0000660|            0.0006569|            0.0025302|            0.0690916|            0.0888414| 0.0746140|
|    0.0560089|    0.1993966|      0.0621763|      0.0009045|      0.0036597|      0.0054941|    0.0043480|      0.0044963|      0.0017555|      0.0084072|    0.0311406|      0.1224158|      0.0334652|      0.0108072|      0.0020314|      0.0015495|      0.0007844|      0.0073389|      0.0016374|      0.0207108|      0.0015147|       0.0085509|      0.0418152|      0.0056803|      0.0024959|      0.0007186|      0.0086272|      0.0007766|      0.0050396|      0.1375296|       0.0454824|      0.0066172|      0.0024761|      0.0043525|      0.0006707|      0.0021301|      0.0005243|      0.0122009|      0.1342687|      0.0019381|      0.0035051|        0.0027773|         1.35e-05|        0.0004209|        0.0000144|      0.0008635|        0.0000194|        0.0001671|        0.0000288|      0.0003219|        0.0082416|        0.0010711|        0.0015401|         3.05e-05|          1.5e-06|        0.0000059|        0.0015270|        0.0001205|        0.0012847|         4.43e-05|         0.0010478|        0.0105549|        0.0004779|        0.0000206|         1.00e-07|        0.0001562|        0.0000012|        0.0005087|        0.0002161|         0.0046874|        0.0002321|        0.0035418|        0.0000906|         9.70e-06|        0.0001440|        0.0001202|        0.0026136|        0.0019608|          0.0385154|          0.0696546|            0.0551910|            0.0002689|            0.0083636|            0.0002852|          0.0171597|            0.0003853|            0.0033203|            0.0005718|          0.0063978|            0.1637800|            0.0212859|            0.0306061|            0.0006069|             2.94e-05|            0.0001174|            0.0303453|            0.0023954|            0.0255310|            0.0008809|             0.0208233|            0.2097509|            0.0094974|            0.0004101|            0.0000025|            0.0031035|            0.0000240|            0.0101087|            0.0042948|             0.0931494|            0.0046122|            0.0703841|            0.0018012|            0.0001924|            0.0028619|            0.0023879|            0.0519379|            0.0389659| 0.0503210|
|    0.0536044|    0.1916365|      0.0711010|      0.0010240|      0.0015210|      0.0075065|    0.0059829|      0.0058547|      0.0012665|      0.0105092|    0.0239668|      0.1617185|      0.0383445|      0.0089801|      0.0016851|      0.0013215|      0.0008915|      0.0045276|      0.0007313|      0.0140701|      0.0014783|       0.0133808|      0.0611832|      0.0036268|      0.0016014|      0.0007879|      0.0063729|      0.0004538|      0.0037150|      0.1157722|       0.0753060|      0.0062637|      0.0021495|      0.0041354|      0.0003677|      0.0011093|      0.0004140|      0.0110648|      0.0845737|      0.0112392|      0.0103927|        0.0061986|         7.60e-05|        0.0002106|        0.0000047|      0.0007563|        0.0000227|        0.0000918|        0.0001002|      0.0001074|        0.0177674|        0.0001146|        0.0000145|         4.05e-05|          1.6e-06|        0.0001552|        0.0001715|        0.0002282|        0.0002335|         4.50e-06|         0.0014146|        0.0124359|        0.0001925|        0.0014653|         3.50e-06|        0.0001069|        0.0000219|        0.0000435|        0.0022707|         0.0075854|        0.0000000|        0.0030467|        0.0001723|         1.50e-06|        0.0000405|        0.0000800|        0.0059879|        0.0092741|          0.1220660|          0.1128721|            0.0673211|            0.0008259|            0.0022876|            0.0000508|          0.0082139|            0.0002463|            0.0009975|            0.0010880|          0.0011666|            0.1929665|            0.0012449|            0.0001579|            0.0004397|             1.74e-05|            0.0016857|            0.0018630|            0.0024789|            0.0025355|            0.0000494|             0.0153636|            0.1350625|            0.0020904|            0.0159142|            0.0000380|            0.0011607|            0.0002380|            0.0004725|            0.0246620|             0.0823833|            0.0000000|            0.0330888|            0.0018716|            0.0000159|            0.0004396|            0.0008686|            0.0650327|            0.1007230| 0.0920750|



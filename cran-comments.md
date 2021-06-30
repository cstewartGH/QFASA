
# Submit Version 1.1.1

## Test Environments
### Rhub

```
Build ID:	QFASA_1.1.1.tar.gz-23fe7bf608b59f95970bd71e758f9a52
Platform:	Windows Server 2008 R2 SP1, R-release, 32/64 bit
Submitted:	9 minutes 12.9 seconds ago
Build time:	8 minutes 55.9 seconds
```

```
Build ID:	QFASA_1.1.1.tar.gz-0c35ee472464bc4a0d72a22e103eddca
Platform:	macOS 10.13.6 High Sierra, R-release, CRAN's setup
Submitted:	10 minutes 49.4 seconds ago
Build time:	10 minutes 47.8 seconds

```

## R CMD check results
There was 1 warning:

```
 WARNING
  'qpdf' is needed for checks on size reduction of PDFs
```

There were 2 notes:

```
> checking installed package size ... NOTE
    installed size is 56.3Mb
    sub-directories of 1Mb or more:
      doc    1.5Mb
      libs  54.3Mb
  NB: this package is only installed for sub-architecture 'x64'

> checking compiled code ... NOTE
  Note: information on .o files for x64 is not available
  File 'H:/QFASA R Package/QFASA.Rcheck/QFASA/libs/x64/CommonDiet.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)
  File 'H:/QFASA R Package/QFASA.Rcheck/QFASA/libs/x64/ErrorModelSimpleEquant.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)
  File 'QFASA/libs/x64/CommonDiet.dll':
    Found non-API call to R: 'R_RunExitFinalizers'
  File 'QFASA/libs/x64/ErrorModelSimpleEquant.dll':
    Found non-API call to R: 'R_RunExitFinalizers'
  
  Compiled code should not call entry points which might terminate R nor
  write to stdout/stderr instead of to the console, nor use Fortran I/O
  nor system RNGs. The detected symbols are linked into the code but
  might come from libraries and not actually be called.
  Compiled code should not call non-API entry points in R.
```

# Submit Version 1.1.0
## Test Environments
* local OS X build, check, install, R version 4.0.3. 
* win-builder

## R CMD check results
There was 1 warning:

```
 WARNING
  'qpdf' is needed for checks on size reduction of PDFs
```

There were 2 notes:

```
checking installed package size ... NOTE
    installed size is 21.8Mb
    sub-directories of 1Mb or more:
      doc    1.5Mb
      libs  19.9Mb
      

> checking compiled code ... NOTE
  Note: information on .o files for i386 is not available
  Note: information on .o files for x64 is not available
  File 'H:/QFASA R Package/QFASA.Rcheck/QFASA/libs/i386/QFASA.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)
  File 'H:/QFASA R Package/QFASA.Rcheck/QFASA/libs/x64/QFASA.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)
  File 'QFASA/libs/i386/QFASA.dll':
    Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
  File 'QFASA/libs/x64/QFASA.dll':
    Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
  
  Compiled code should not call entry points which might terminate R nor
  write to stdout/stderr instead of to the console, nor use Fortran I/O
  nor system RNGs. The detected symbols are linked into the code but
  might come from libraries and not actually be called.
  It is good practice to register native routines and to disable symbol
  search.
```

# Submit Version 1.0.3
## Test Environments
* local OS X build, check, install, R version 3.5.1
* travis-ci Ubuntu Linux 16.04.6 LTS build, check, install, R version 3.6.0
* win-builder

## R CMD check results
There was 1 NOTE:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Justin Kamerman <justin@kaleco.net>'

New maintainer:
  Justin Kamerman <justin@kaleco.net>
Old maintainer(s):
  Justin Kamerman <justin.kamerman@unb.ca>
```


# Submit Version 1.0.2
## Test Environments
* local OS X build, check, install, R 3.3.0
* local AWS Linux build, check, install, R 3.2.2
* win-builder

## R CMD check results
There was 1 NOTE:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Justin Kamerman <justin.kamerman@unb.ca>’

License components with restrictions and base license permitting such:
  MIT + file LICENSE
File 'LICENSE':
  YEAR: 2016
  COPYRIGHT HOLDER: Connie Stewart, Justin Kamerman
```


# Submit Version 1.0.1
## Test Environments
* local OS X build, check, install, R 3.3.0
* local AWS Linux build, check, install, R 3.2.2
* win-builder

## R CMD check results
There was 1 NOTE:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Justin Kamerman <justin.kamerman@unb.ca>’

License components with restrictions and base license permitting such:
  MIT + file LICENSE
File 'LICENSE':
  YEAR: 2016
  COPYRIGHT HOLDER: Connie Stewart, Justin Kamerman
``

# Submit Version 1.0.0 
## Resubmission
* Added @examples code to main function p.QFASA() detailing use of package:

```
* checking files in ‘vignettes’ ... OK
* checking examples ... OK
* checking for unstated dependencies in ‘tests’ ... OK
* checking tests ...
  Running ‘testthat.R’
 OK
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in ‘inst/doc’ ... OK
* checking re-building of vignette outputs ... OK
* DONE

Status: OK

R CMD check results
0 errors | 0 warnings | 0 notes
```

## Resubmission
* Incorporate late authorship change request in DESCRIPTION
* Updated LICENSE to comply with CRAN template

## Test Environments
* local OS X build, check, install, R 3.2.4
* local OS X build, check, install, R 3.3.0
* local AWS Linux build, check, install, R 3.2.2
* win-builder

## R CMD check results
There were no NOTES, ERRORS or WARNINGS.

Submitted from Linux because of LaTeX errors on OSX:

```
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
LaTeX errors found:
! Font T1/phv/m/n/14.4=phvr8t at 14.4pt not loadable: Metric (TFM) file not fou
nd.
<to be read again> 
                   relax 
l.45 \Rdcontents{\R{} topics documented:}
* checking PDF version of manual without hyperrefs or index ... ERROR
* DONE

Status: 1 ERROR, 1 WARNING
See
  ‘/private/var/folders/g6/7xdk55_s7d5b4q17985h2b1r0000gp/T/Rtmpqktvy0/QFASA.Rcheck/00check.log’
for details.

R CMD check results
1 error  | 1 warning  | 0 notes
checking PDF version of manual without hyperrefs or index ... ERROR

checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
LaTeX errors found:
! Font T1/phv/m/n/14.4=phvr8t at 14.4pt not loadable: Metric (TFM) file not fou
nd.
<to be read again> 
                   relax 
l.45 \Rdcontents{\R{} topics documented:}
```



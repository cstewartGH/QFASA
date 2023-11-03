
# Submit Version 1.2.0

## Test Environments
### Rhub

> check_for_cran()

  ‘connie.stewart@unb.ca’ is not validated, or does not match the package maintainer's email. To
  validate it now, please enter the email address below. Note that R-hub will send a token to this address.
  If the address does not belong to you, quit now by pressing ENTER .

  Email address: connie.stewart@unb.ca

Please check your emails for the R-hub access token.
Token: 9b25182aed00494587c75a7821e577ea
Token added for ‘connie.stewart@unb.ca’

For info the token(s) and email(s) are stored at C:\Users\cstewart\AppData\Local/rhub/rhub/validated_emails.csv

─  Building package
─  Uploading package
─  Preparing build, see status at
   https://builder.r-hub.io/status/QFASA_1.2.0.tar.gz-f855c924a1ae4fc18eedc2f0d6eaa554
   https://builder.r-hub.io/status/QFASA_1.2.0.tar.gz-3db280cee5c64994ba62b753cef10bcb
   https://builder.r-hub.io/status/QFASA_1.2.0.tar.gz-23c427e30bec44ca89cec72603dc2e2b
   https://builder.r-hub.io/status/QFASA_1.2.0.tar.gz-00685c9085c04c62ac64b3982220985b
─  Build started
─  Creating new user
─  Downloading and unpacking package file
─  Querying package dependencies
─  Installing package dependencies
─  Running R CMD check
   setting _R_CHECK_FORCE_SUGGESTS_ to false
   setting R_COMPILE_AND_INSTALL_PACKAGES to never
   setting R_REMOTES_STANDALONE to true
   setting R_REMOTES_NO_ERRORS_FROM_WARNINGS to true
   setting _R_CHECK_FORCE_SUGGESTS_ to true
   setting _R_CHECK_CRAN_INCOMING_USE_ASPELL_ to true
   'getOption("repos")' replaces Bioconductor standard repositories, see
   'help("repositories", package = "BiocManager")' for details.
   Replacement repositories:
       CRAN: https://cloud.r-project.org
─  using log directory 'C:/Users/USERfwXEvjmzNU/QFASA.Rcheck'
─  using R Under development (unstable) (2023-10-14 r85331 ucrt)
─  using platform: x86_64-w64-mingw32
─  R was compiled by (788ms)
       gcc.exe (GCC) 12.2.0
       GNU Fortran (GCC) 12.2.0
─  running under: Windows Server 2022 x64 (build 20348)
─  using session charset: UTF-8 (766ms)
─  using option '--as-cran'
✔  checking for file 'QFASA/DESCRIPTION'
─  this is package 'QFASA' version '1.2.0'
─  package encoding: UTF-8 (770ms)
─  checking CRAN incoming feasibility ... [10s] Note_to_CRAN_maintainers (2.3s)
   Maintainer: 'Connie Stewart <connie.stewart@unb.ca>'
✔  checking package namespace information
✔  checking package dependencies (3.1s)
✔  checking if this is a source package
✔  checking if there is a namespace
✔  checking for executable files (1.6s)
✔  checking for hidden files and directories
✔  checking for portable file names
─  checking whether package 'QFASA' can be installed ... [107s] OK (1m 50.4s)
─  used C++ compiler: 'G__~1.EXE (GCC) 12.2.0'
N  checking installed package size
     installed size is 26.2Mb
     sub-directories of 1Mb or more:
       libs  25.5Mb
✔  checking package directory
✔  checking for future file timestamps
✔  checking 'build' directory (1.5s)
✔  checking DESCRIPTION meta-information
✔  checking top-level files
✔  checking for left-over files
✔  checking index information (789ms)
✔  checking package subdirectories
✔  checking R files for non-ASCII characters
✔  checking R files for syntax errors
✔  checking whether the package can be loaded (774ms)
✔  checking whether the package can be loaded with stated dependencies
✔  checking whether the package can be unloaded cleanly
✔  checking whether the namespace can be loaded with stated dependencies
✔  checking whether the namespace can be unloaded cleanly (779ms)
✔  checking loading without being on the library search path (1.5s)
✔  checking use of S3 registration (3.1s)
✔  checking dependencies in R code (3.8s)
✔  checking S3 generic/method consistency (1.5s)
✔  checking replacement functions
✔  checking foreign function calls (2.3s)
─  checking R code for possible problems ... [13s] OK (12.8s)
✔  checking Rd files (1.6s)
✔  checking Rd metadata
✔  checking Rd line widths
✔  checking Rd cross-references (768ms)
✔  checking for missing documentation entries (783ms)
✔  checking for code/documentation mismatches (2.3s)
✔  checking Rd \usage sections (1.5s)
✔  checking Rd contents
✔  checking for unstated dependencies in examples (1.5s)
✔  checking contents of 'data' directory
✔  checking data for non-ASCII characters
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves (1.5s)
✔  checking line endings in C/C++/Fortran sources/headers
✔  checking line endings in Makefiles
✔  checking for GNU extensions in Makefiles
✔  checking for portable use of $(BLAS_LIBS) and $(LAPACK_LIBS)
✔  checking use of PKG_*FLAGS in Makefiles
✔  checking use of SHLIB_OPENMP_*FLAGS in Makefiles
✔  checking include directives in Makefiles
✔  checking pragmas in C/C++ headers and code
✔  checking compilation flags used
✔  checking compiled code (33.7s)
✔  checking installed files from 'inst/doc'
✔  checking files in 'vignettes'
─  checking examples ... [349s] NOTE (5m 49.4s)
   Examples with CPU (user + system) or elapsed time > 5s
                          user system elapsed
   forward.selection    174.70   8.05  182.67
   backward.elimination 150.77   8.11  159.02
✔  checking for unstated dependencies in 'tests'
─  checking tests
✔  Running 'testthat.R' (3.8s)
✔  checking for unstated dependencies in vignettes
✔  checking package vignettes in 'inst/doc'
─  checking re-building of vignette outputs ... [10s] OK (11.5s)
─  checking PDF version of manual ... [13s] OK (13.1s)
─  checking HTML version of manual ... [11s] OK (10.7s)
N  checking for non-standard things in the check directory
   Found the following files/directories:
     ''NULL''
N  checking for detritus in the temp directory
   Found the following files/directories:
     'lastMiKTeXException'
   
─  Done with R CMD check
─  Cleaning up files and user

# Submit Version 1.1.2

## Test Environments
### Rhub

rhub::check_for_cran()

* checking installed package size ... NOTE
  installed size is 26.2Mb
  sub-directories of 1Mb or more:
    libs  25.6Mb
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'


── R CMD check results ───────────── QFASA 1.1.2 ────
Duration: 2m 57.6s

❯ checking installed package size ... NOTE
    installed size is 26.2Mb
    sub-directories of 1Mb or more:
      libs  25.6Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

R CMD check succeeded


# Submit Version 1.1.1




```
Build ID:	QFASA_1.1.1.tar.gz-e4407e595ad7101d0dce4c5a080e691e
Platform:	Windows Server 2008 R2 SP1, R-release, 32/64 bit
Submitted:	8 minutes 29.7 seconds ago
Build time:	7 minutes 13.5 seconds

```
```
Build ID:	QFASA_1.1.1.tar.gz-9ecb3094393c126230287d7406453bea
Platform:	macOS 10.13.6 High Sierra, R-release, CRAN's setup
Submitted:	20 minutes 32 seconds ago
Build time:	20 minutes 12.6 seconds



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
    installed size is 29.3Mb
    sub-directories of 1Mb or more:
      libs  28.7Mb
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
  
  Compiled code should not call entry points which might terminate R nor
  write to stdout/stderr instead of to the console, nor use Fortran I/O
  nor system RNGs. The detected symbols are linked into the code but
  might come from libraries and not actually be called.
  
  
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



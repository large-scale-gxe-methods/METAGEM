# METAGEM


METAGEM (META-analysis of GEM summary statistics) is a software program for meta-analysis of large-scale gene-environment interaction testing results, including multi-exposure interactions, joint (main effect and interactions) tests, and marginal tests. It uses results directly from [GEM](https://github.com/large-scale-gxe-methods/GEM) output.


Current version: 1.0

## Contents 
- [Dependencies](#dependencies)
- [Quick Installation](#quick-installation)
- [Usage](#usage)
- [Contact](#contact)
- [License](#license)

<br />  

## Dependencies

- Compiler with C++11 support
- [Boost C++ Libraries](https://www.boost.org/) (Versions 1.70.0 - 1.79.0)
- Intel Math Kernal Library (MKL)

<br />

## Quick Installation

To install METAGEM, run the following lines of code:
 ```
git clone https://github.com/large-scale-gxe-methods/METAGEM
cd METAGEM
cd src
make
 ```
 
<br />

## Usage

### Running METAGEM

- [Command Line Options](#command-line-options)
- [Input Files](#input-files)
- [Output File](#output-file)
- [Example](#example)

<br />

### Command Line Options

Once METAGEM is installed, the executable ```./METAGEM``` can be used to run the program.
For a list of options, use ```./METAGEM --help```.

<details>
     <summary> <b>List of Options</b> </summary>

```
General Options:

   --help 
     Prints available options and exits.


Input/Output File Options:

   --input-files         
     Output files from GEM 'meta' or 'full' option separated by space. At least two files are required.
     
   --input-file-list     
     A no header text file containing a single file name per line. This file should contain at least two file names.
     
   --exposure-names
     The names of the exposure(s) to be included in the meta-analysis.
     
   --out                 
     Full path and extension to where METAGEM output results.
     Default: metagem.out
   
   --meta-option         
     Integer value indicating which summary statistics should be used for meta-analysis.
                        
     0: Both model-based and robust summary statistics.                      
     1: model-based summary statistics.                       
     2: robust summary statistics.                       
     Default: 0

   --additional-test
     The variable name(s) (may or may not include the genetic main effect G at first) and the full path of the output file for one additional test.
   
```
</details>

<br /> 

### Input Files

METAGEM accepts output files from GEM (v1.4.1 or later) with '--output-style' set to 'meta' or 'full'. Multiple GEM output files can be specified using the '--input-files' flag, separated by spaces. Alternatively, the '--input-file-list' option can be used to specify a text file without headers, where each line contains a single input file name.

<br />

### Output File

METAGEM will write results to the output file specified with the --out parameter (or 'metagem.out' if no output file is specified).
Below are details of the possible column headers in the output file.

```diff 
SNPID              - The SNP identifier as retrieved from the input files.
CHR                - The chromosome of the SNP.
POS                - The physical position of the SNP. 
Non_Effect_Allele  - The reference allele in association testing.  
Effect_Allele      - The coding allele in association testing.  
N_Samples          - The combined sample size from all studies in the meta-analysis.
AF                 - The summary effect allele frequency in all studies combined in the meta-analysis.  

Beta_Marginal           - The summary marginal genetic effect estimate (i.e., from a model with no interaction terms) from univariate meta-analysis using model-based results from each study.
SE_Beta_Marginal        - SE for the summary marginal genetic effect estimate from univariate meta-analysis using model-based results from each study.  
robust_Beta_Marginal    - The summary marginal genetic effect estimate (i.e., from a model with no interaction terms) from univariate meta-analysis using robust results from each study.
robust_SE_Beta_Marginal - SE for the summary marginal genetic effect estimate from univariate meta-analysis using robust results from each study.

Beta_G             - The summary genetic main effect (G) estimate from joint (main effect and interactions) meta-analysis using model-based results from each study.
Beta_G-*           - The summary GxE interaction effect estimate(s) from joint (main effect and interactions) meta-analysis using model-based results from each study.
SE_Beta_G          - SE for the summary genetic main effect (G) estimate from joint meta-analysis using model-based results from each study.  
SE_Beta_G-*        - SE for the summary GxE interaction effect estimate(s) from joint meta-analysis using model-based results from each study.
Cov_Beta_G_G-*     - Covariance(s) between the summary genetic main effect (G) estimate and the summary GxE interaction effect estimate(s) from joint meta-analysis using model-based results from each study.  
Cov_Beta_G-*_G-*   - Covariance(s) between the summary GxE interaction effect estimate(s) from joint meta-analysis using model-based results from each study.
robust_Beta_G      - The summary genetic main effect (G) estimate from joint (main effect and interactions) meta-analysis using robust results from each study.
robust_Beta_G-*    - The summary GxE interaction effect estimate(s) from joint (main effect and interactions) meta-analysis using robust results from each study.
robust_SE_Beta_G   - SE for the summary genetic main effect (G) estimate from joint meta-analysis using robust results from each study.  
robust_SE_Beta_G-* - SE for the summary GxE interaction effect estimate(s) from joint meta-analysis using robust results from each study.
robust_Cov_Beta_G_G-*   - Covariance(s) between the summary genetic main effect (G) estimate and the summary GxE interaction effect estimate(s) from joint meta-analysis using robust results from each study.
robust_Cov_Beta_G-*_G-* - Covariance(s) between the summary GxE interaction effect estimate(s) from joint meta-analysis using robust results from each study.

P_Value_Marginal           - The summary marginal genetic effect test p-value from univariate meta-analysis using model-based results from each study.
P_Value_Interaction        - The summary GxE interaction effect test p-value (K degrees of freedom test) from joint (main effect and interactions) meta-analysis using model-based results from each study. (K is the number of GxE interaction terms)
P_Value_Joint              - Joint (main effect and interactions) test p-value (K+1 degrees of freedom test) from joint meta-analysis using model-based results from each study.
robust_P_Value_Marginal    - The summary marginal genetic effect test p-value from univariate meta-analysis using robust results from each study.
robust_P_Value_Interaction - The summary GxE interaction effect test p-value (K degrees of freedom test) from joint (main effect and interactions) meta-analysis using robust results from each study. (K is the number of GxE interaction terms)
robust_P_Value_Joint       - Joint (main effect and interactions) test p-value (K+1 degrees of freedom test) from joint meta-analysis using robust results from each study.
```

<br />

The '--meta-option' flag can be used to specify which columns should be included in the output file:

* 0 - Meta-analyses will be performed on both model-based and robust results from each study, and all column headers listed above will be available in the output file.
* 1 - Only model-based results from each study will be used in the meta-analysis, and columns above containing the 'robust_' prefix will be excluded from the output file.
* 2 - Only robust results from each study will be used in the meta-analysis, and summary statistics columns above without the 'robust_' prefix will be excluded from the output file.
* Default: 0 
 
<br />

The '--additional-test' flag can be used to define an additional test to run. The genetic main effect 'G' may or may not be specified at first, indicating whether the additional test is a joint test or an interaction-only test. Then the name(s) of interaction(s) which will be included in the additional test should be listed. At last, the full path of the output file of the additonal test should be specified.

For an example situation which has a total of 2 covariates (cov1 and cov2) in the full test, the '--additional-test' flag may be used as:
```unix
--additional-test G cov1 metagem2.out
```
This is to define an additional joint test with only cov1 as interaction. The flag may also be used as:
```unix
--additional-test cov1 metagem2.out
```
This is to define an additional interaction-only test with only cov1 as interaction.

<br />

### Example
```unix
./METAGEM --input-files file1.out file2.out file3.out --exposure-names cov1 --out metagem.out
```
<br />
<br />

## Contact 
For comments, suggestions, bug reports and questions, please contact Han Chen (Han.Chen.2@uth.tmc.edu), Alisa Manning (AKMANNING@mgh.harvard.edu), Kenneth Westerman (KEWESTERMAN@mgh.harvard.edu) or Cong Pan (Cong.Pan@uth.tmc.edu). For bug reports, please include an example to reproduce the problem without having to access your confidential data.

<br />
<br />

## License 

 ```
 METAGEM: META-analysis of GEM summary statistics
 Copyright (C) 2021-2023 Duy T. Pham and Han Chen
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ```

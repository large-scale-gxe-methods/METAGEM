# METAGEM


METAGEM (META-analysis of GEM summary statistics) is a software program for meta-analysis of large-scale gene-environment testing results, including multi-exposure interaction, joint, and marginal tests. It uses results directly from [GEM](https://github.com/large-scale-gxe-methods/GEM) output.


Current version: 1.0

## Contents

- [Quick Installation](#quick-installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Contact](#contact)
- [License](#license)

## Quick Installation

To install METAGEM, run the following lines of code:
 ```
git clone https://github.com/large-scale-gxe-methods/METAGEM
cd METAGEM
cd src
make
 ```
## Dependencies
- boost-1.70.0 or up, 
- GCC-5.2 or up, 
- Intel Math Kernel Library: tested with mkl-2019.3.199

## Usage

### Running METAGEM

1. [Command Line Options](#command-line-options)

### Command Line Options

Once METAGEM is installed, the executable ```./METAGEM``` can be used to run the program.
For a list of options, use ```./METAGEM -help```.

<details>
     <summary> <b>List of Options</b> </summary>

```
General Options:

   --help 
   Prints available options and exits.
   
   --version 
   Prints the version of METAGEM and exits.


Input/Output File Options:
   --input-files         
     Output files from GEM 'meta' or 'full' option.
     
   --input-file-list     
     A no header text file containing a single file name per line.
     
   --exposure-names
     
   --out                 
   Full path and extension to where METAGEM output results.
   Default: metagem.out
   
   --meta-option         
     Integer value indicating which summary statistics should be used for meta-analysis.
                        
     0: Both model-based and robust summary statistics.                      
     1: model-based summary statistics.                       
     2: robust summary statistics.                       
     Default: 0
   
```
</details>

<br /> 

### Input Files

METAGEM uses output files from GEM (v1.4.1 and up, with use of --output-style: 'meta' or 'full').

### Input File List

A no header txt file containing a single input file name per line.

### Output File Format

METAGEM will write results to the output file specified with the --out parameter (or 'metagem.out' if no output file is specified).
Below are details of the possible column headers in the output file.

```diff 
SNPID              - The SNP identifier as retrieved from the genotype file.
CHR                - The chromosome of the SNP.
POS                - The physical position of the SNP. 
Non_Effect_Allele  - The allele not counted in association testing.  
Effect_Allele      - The allele that is counted in association testing.  
N_Samples          - The number of samples without missing genotypes.
AF                 - The allele frequency of the effect allele.  

Beta_Marginal           - The coefficient estimate for the marginal genetic effect (i.e., from a model with no interaction terms).
SE_Beta_Marginal        - The model-based SE associated with the marginal genetic effect estimate.  
robust_Beta_Marginal  
robust_SE_Beta_Marginal - The robust SE associated with the marginal genetic effect estimate.

Beta_G             - The coefficient estimate for the genetic main effect (G).
Beta_G-*           - The coefficient estimate for the interaction or interaction covariate terms.
SE_Beta_G          - Model-based SE associated with the the genetic main effect (G).  
SE_Beta_G-*        - Model-based SE associated with any GxE or interaction covariate terms.
Cov_Beta_G_G-*          - Model-based covariance between the genetic main effect (G) and any GxE or interaction covariate terms.  
robust_Beta_G  
robust_Beta_G-*    
robust_SE_Beta_G   - Robust SE associated with the the genetic main effect (G).  
robust_SE_Beta_G-* - Robust SE associated with any GxE or interaction covariate terms.
robust_Cov_Beta_G_G-*   - Robust covariance between the genetic main effect (G) and any GxE or interaction covariate terms.   

P_Value_Marginal           - Marginal genetic effect p-value from model-based SE.
P_Value_Interaction        - Interaction effect p-value (K degrees of freedom test of interaction effect) from model-based SE. (K is number of major exposures)
P_Value_Joint              - Joint test p-value (K+1 degrees of freedom test of genetic and interaction effect) from model-based SE.
robust_P_Value_Marginal    - Marginal genetic effect p-value from robust SE.
robust_P_Value_Interaction - Interaction effect p-value from robust SE.
robust_P_Value_Joint       - Joint test p-value (K+1 degrees of freedom test of genetic and interaction effect) from robust SE.
```

<br />

The --meta-option flag can be used to specify which columns should be included in the output file:

* 0: Both model-based and robust summary statistics.
* 1: model-based summary statistics.
* 2: robust summary statistics.
* Default: 0 

### Examples
<br />

```unix
./METAGEM --input-files file1.out file2.out file3.out --exposure-names cov1 --out metagem.out
```
<br />
<br />

## Contact 
For comments, suggestions, bug reports and questions, please contact Han Chen (Han.Chen.2@uth.tmc.edu), Alisa Manning (AKMANNING@mgh.harvard.edu), Kenny Westerman (KEWESTERMAN@mgh.harvard.edu) or Cong Pan (Cong.Pan@uth.tmc.edu). For bug reports, please include an example to reproduce the problem without having to access your confidential data.

<br />
<br />

## License 

 ```
 METAGEM: META-analysis of GEM summary statistics
 Copyright (C) 2022 Duy T. Pham and Han Chen
 
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
 

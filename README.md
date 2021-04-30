# PremPS
## About
<font size=4> 
  
PremPS evaluates the effects of single mutations on protein stability by calculating the changes in unfolding Gibbs free energy. It can be applied to a large number of tasks, including finding functionally important variants, understanding their molecular mechanisms and protein design. 3D structure of a protein is required for this method.
  
</font>

## Scoring mutations with PremPS
<font size=4> 

We recommend that most users who just want to obtain PremPS predictions use [PremPS website](https://lilab.jysw.suda.edu.cn/research/PremPS/) to obtain scores.

</font>

## Source code releases
<font size=4> 
  
You can download [releases](https://github.com/minghuilab/PremPS/releases) on github.

</font>

## Installation

#### I. PREREQUISITES

<font size=4>
 
PremPS requires the following software and packages.

1. DSSP

   This is available at the DSSP website.

   https://swift.cmbi.umcn.nl/gv/dssp/

2. PROVEAN

   This is available at the PROVEAN website.

   http://provean.jcvi.org/index.php/

3. NCBI BLAST 2.4.0

   This is available at the NCBI ftp site.

   ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/

4. FoldX

   This is available at the FoldX website.

   http://foldxsuite.crg.eu/

5. VMD

   This is available at the VMD website.

   https://www.ks.uiuc.edu/Research/vmd/

6. Python packages: pandas and rpy2

   To install these packages you can use the following command:
</font>

<font size=4>

	$ conda install -c conda-forge pandas
	$ conda install -c r rpy2

</font> 

<font size=4>

7. R packages: randomForest, e1071, xgboost and stringr

</font>

<font size=4>

	$ install.packages('randomForest')
	$ install.packages('e1071')
	$ install.packages('xgboost')
	$ install.packages('stringr')

</font> 

#### II. INSTALLATION INSTRUCTIONS

<font size=4>

1. Download and/or install prerequisites described above.

2. Download and unpack the distribution:

</font>

<font size=4>

	$ wget https://github.com/minghuilab/PremPS/archive/v1.0.0.tar.gz
	$ tar -zxvf v1.0.0.tar.gz

</font> 

<font size=4>

3. Change to the source directory:

</font>

<font size=4>

	$ cd PremPS-1.0.0

</font> 

<font size=4>

4. Change the path parameters in PremPS.py (line 15-20):

</font>

<font size=4>

	workdir = Your working directory
	pathvmd = path for running VMD software  # /usr/local/bin/vmd
	pathmkdssp = path for running DSSP software  # /usr/local/bin/mkdssp
	pathpsiblast = path for running PSI-BLAST software  # /usr/local/bin/blast/psiblast
	pathblastdb = path for blastdb  # /usr/local/bin/blastdb/nr
	pathrscript = path for running Rscript  # /usr/local/bin/Rscript
	
</font>

&nbsp; &nbsp; The FoldX software needs to be installed in the working directory.

#### III. RUNNING PremPS

<font size=4>

	$ python PremPS.py -i 2020100417132606935696574

</font> 

## Platform

<font size=4>

PremPS is only intended to run on *linux* operating systems.

</font>

## Issues

<font size=4>

You will need to have Python 2 (or 3) and R 3.4.0 (or higher) installed.

</font>

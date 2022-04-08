# PlasmaCNA

PlasmaCNA provides a toolkit for analyzing tumor-derived copy number changes using low pass whole-genome sequencing of plasma cell-free DNA.
 - Ability to generate and analyze both fixed size windows (ichorCNA, tMAD) and variable sized windows (PlasmaSeq)
 - A set of functions for window read depth normalization
 - A julia implementation of the PlasmaSeq algorithm

## Citations
If PlasmaCNA is used in published research, please cite:  
>Favaro, Patricia F., et al. Feasibility of circulating tumor DNA analysis in dogs with naturally-occurring malignant and benign splenic lesions. (2022) Scientific Reports doi:10.1038/s41598-022-09716-6

Publications using the PlasmaSeq implementation provided in PlasmaCNA should also cite:
>Farris, C., Trimarchi, J.M. Plasma-seq: a novel strategy for metastatic prostate cancer analysis. Genome Med 5, 35 (2013). https://doi.org/10.1186/gm439

## Dependencies
Julia package dependencies are handled by the julia package manager, except QuickArgParse  
https://github.com/brmcdonald/QuickArgParse  
The R package DNACopy is used by the PlasmaSeq implementation for segmentation  
GenMap is required for calculating mappability:  
https://github.com/cpockrandt/genmap  

## Usage

The primary workflows are provided in the /bin directory of the package, 
with inputs specified in their usage statements.  

**fixbin_generation.jl**  
Generates a table of fixed size windows across the genome, 
annotated with GC content and average mappability. Used by ichorCNA and other methods.

**flexbin_generation.jl**  
Generates a table of windows across the genome such that each has approximately the 
same number of mappable bases. Windows are annotated with GC content and average mappability. 
Used by PlasmaSeq.  

**plasmaseq.jl**  
An implementation of the PlasmaSeq method for identifying large tumor-derived copy number changes.
The original PlasmaSeq implementation can be found here: https://github.com/PeterUlz/PlasmaSeq

## Generating a mappability wig file using GenMap
```Bash
genmap index -F /path/to/refgenome.fna -I genmap.index
genmap -K 50 -E 1 -w -T $NumCores -I genmap.index -O genmap.k50.e1.wig
```



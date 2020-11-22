# Thesis2020

This repository collects the codes used in simulations by Michele Lambardi di San Miniato for his Ph.D thesis.

The thesis will be available for consultation by 2020-2021 on http://paduaresearch.cab.unipd.it (stub).

The **main** branch contains the original codes, so they lack proper documentation. Some suitable **side** branches will be created to make the code more readable and useful for the public (stub).

The code can be run as is, but care should be taken in changing all paths to inputs and outputs for the program, as by default it **overwrites contents in** a folder named Temp in the user's home directory, which is "**~/Temp**" in Linux environments.

The code uses the R package pbapply for parallelization purposes, so a small adaptation is require to run, for instance, on **Windows** systems. One can replace all occurrences of pbapply functions with just lapply and it will do fine, it will just take more time to run. Consider running codes in a plain R console instead of more consuming software such as RStudio.

# Subject

The code is written in **R**. It uses some base packages, plus some **brglm2nc** package which is not available on CRAN but modified the cran package **brglm2**, in a way that is illustrated in the directory brglm2nc of this repository.

It was mentioned that the code is poorly documented and made available for reproducibility checks mainly. But it is also meant to provide running examples that implement some enhanced statistical methods that were either discussed or introduced in the author's Ph.D thesis.

Some codes deal with bias reduction techniques. All of the codes deal with the novel method discussed in the author's thesis, which approximately improves some statistical properties of so called profile score functions. For more details about these methods, one can refer to the material linked in the references.

# References

(stub)

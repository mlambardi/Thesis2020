# Readme

Instructions follows that allow to install R package brglm2 with a slight modification that implements the novel statistical method illustrated in the author's theory at least in the case of Generalized Linear Models (glm). The original package provides a function brglmFit() that can be fed to R package stats' glm() function that implements:

- maximum likelihood estimation, with type="ML",
- mean bias reduction, with type="AS_mean",
- median bias reduction, with type="AS_median",

among others. In general, as new versions of the package are developed, more statistical methods are implemented that enhance maximum likelihood estimation in some sense. Our modification to that package implements

- nuisance control, with type="AS_nc"

You can install the original R package brglm2 in R with:

        install.packages("brglm2")

Or you can install a modified version as follows, or both. The following version has a different name, it is the R package **brglm2nc**, so it won't conflict with brglm2 if you want to install both. Consider that brglm2nc is based on brglm2 v0.6, so you may install brglm2 too, so to catch up with the any other additional method that has been implemented in the meantime.

The package brglm2nc is not available on CRAN, but you may install it manually as follows.

# Instructions

Commands are provided for Linux-like command lines.

If some commands do not work, please, open a discussion on the issues, as the uploader run the codes on Ubuntu 20.* which may come with some pre-installed relevant utilities.

1. Download R package **brglm2** in its version 0.6 as a zip file, so the file at:

        https://github.com/ikosmidis/brglm2/archive/v0.6.zip

   Save it to some local download folder I'll denote >>download<<. In Linux:

        wget 'https://github.com/ikosmidis/brglm2/archive/v0.6.zip' -O brglm2.zip

2. Unzip the archive into a directory brglm2, so:

        unzip -o brglm2.zip -d brglm2

3. Enter the newly create directory **brglm2**, then into the subfolder **brglm2-0.6**. You'll be at **>>download<</brglm2/brglm2-0.6**.

4. Replace files:

    - **DESCRIPTION**
    - **R/brglmControl.R**
    - **R/brglmFit.R**
    
   with the counterparts in the same subfolder of this repository as this installation instructions. This action only expands the functionalities of the package, it shouldn't hamper its pre-existing features.

5. Consider using **diff** to compare the two versions. The actual modifications are just few.

6. In the file explorer, leave the folder brglm2-0.6 so that you are in >>download<</brglm2 and you see the folder brglm2-0.6 it contains.

7. In >>download<</brglm2, open a terminal, then R. If you run the command

        list.files()

   you should see the brglm2-0.6 folder of the R package you are about to install.

8. In R, run:

        install.packages(pkgs = "brglm2-0.6", repos = NULL, type = "source")
    
   to install the package.

9. If the installation failed, maybe you have to install additional R packages. If so, do it, then retry. Please write on the issues of this repo if you find more obstacles.

10. In R, try:

        require(brglm2nc)
        glm(NV ~ HG + PI + EH, data = endometrial, method = "brglmFit", type="AS_nc")

    If it works, the installation is done.

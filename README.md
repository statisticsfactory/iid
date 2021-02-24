## What does it do?

* The code performs Canonical Correlation Analysis (CCA) using the CCA R package
* It also includes the regularized CCA from the CCA R package
* It gives a rich set of plots and tables for interpreting the CCA results, including legends for everything
* It simplifies your work by saving all results in well organized folders
* It allows for reproducibility of your CCA by creating a folder ready to be sent for publication by scientific journals

## Files and folder description

* *demo-data*: folder containing example data for testing the package
* *CCApackage.R*: file containing the functions defining the package
* *main_CCA.R*: file containing code for running the CCA

# How to run a CCA using the CCApackage

1. Put the files *CCApackage.R* and *main_CCA.R* in your current R work directory
2. Open R version >= 4.0.1
3. Open the file named *main_CCA.R*
4. Read the comments and choose to run one of the described ```runCCA()``` function modes
5. Go to the folder you used as your R work directory and look for the results inside the newly created folder named *output_CCA* (or *R_CCA*)

## The ```runCCA()``` function and its modes of action

The ```runCCA()``` function performs the routines of checking for data consistence, estimating the CCA parameters and plotting the results.
It works accepts two arguments:
1. *list_param*, which is a list of parameters and their values for controlling aspects of plots, data transformation and type of CCA
2. *publish*, which is a logical parameter controlling the type of output generated by the function
The file *main_CCA.R* shows and describes all parameters that can be controlled through *list_param*.
The modes of action of the ```runCCA()``` function are determined by *publish*:
* *publish = FALSE* serves for routine analysis, e.i. for trying out different sets of values for *list_param* until you find good results
* *publish = TRUE* serves for the final and definitive analysis, after you have found the right values for *list_param*
With *publish = TRUE* a folder named *R_CCA* will be created in your R work directory. This folder has all anyone needs for reproducing your CCA analysis,
including your datasets and R codes with the list of parameters you chose to fine tune your CCA. A *README* file inside *R_CCA* explaines how to reproduce
the CCA analysis you did. All figures, tables, legends and suggested references are also included in *R_CCA*. With *publish = FALSE* a folder named
*output_CCA* is created in your R work directory, it is analogous to the folder *R_CCA*, although it does not include all one needs to reproduce your CCA
in another computer.

## Changing the plots without needing to re-run ```runCCA()```

The outputs of ```runCCA()``` are automaticly save in the folders *output_CCA* or *R_CCA*, depending on the mode you choose for the *publish* parameter.
But if the results need to be hold in an R object for further analysis or modifications, just do ```results <- runCCA()```. The *results* object is an R list
containing all plots, tables, estimates and legends generates by the ```runCCA()```.

All plots cam be accessed through ```results$plots```, which will have 13 plots save as *p1, p2, ..., p13*. All plots are generated using the functionalities of 
the ggplot2 R package. It means, new layers of graphics can be added to enhance any plot. For example, suppose you wat to enlarge the x axis text of the plot *p1*,
you could do this:
~~~
new_plot1 <- results$plots$p1 + theme(axis.text.x = element_text(size = 24))
~~~
Now suppose you want to change the legend position on plot *p13*:
~~~
new_plot1 <- results$plots$p1 + theme(legend.position = "bottom")
~~~
Any other suitable ggplot2 layer could be added to *p1*. For a complete list of features you can change using the theme() layer, look-uk the 
<https://ggplot2.tidyverse.org/reference/theme.html>. For any other changes, make a google search for **how to do something with ggplot2**.

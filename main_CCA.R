# READ-ME: before anything, set your R working directory as the CCApackage folder (it is assumed you know how to do it!)
# This script has 3 parts:
# parte 1) EXTERNAL SOURCE: here you will call functions necessary for running CCA
# parte 2) PARAMETERS: here you will load your datasets (2 datasets) and set parameters for running your CCA
# parte 3) CCA: here you will actually run your CCA using your datasets and previously defined parameters
# The main function runCCA is in parte 3. It can be run in several modes:
# mode 1) Default mode: runCCA(publish = FALSE). It does not use any parameters by the use, instead it uses default parameters (not adequate for all analysis)
# mode 2) User defined mode: runCCA(list_param, publish = FALSE). It uses a list of parameters defined by the user in part 2 of this script
# mode 3) Publishing mode: runCCA(list_param, publish = TRUE). It creats a set of folders and files ready for publication (only use it after finding an awesome set of parameters for your CCA!)
# After runCCA, go look inside your working directory. A folder named output_CCA was created there. It includes all folders and files generated for your CCA.

# EXTERNAL SOURCE ---------------------------------------------------------------------------------------

# this will load all required functions for CCA
source("CCApackage.R")

# PARAMETERS --------------------------------------------------------------------------------------------

# use read.csv() to upload your datasets or use NA. If you choose NA, you will be directed to dynamically choose an Excel file containing 2 spreadsheets corresponding to data_X and data_Y
# if you use data_X <- NA or data_Y <- NA, after running the function runCCA a pop-up window will appear for you to choose the appropriate Excel file in your computer. 
# this will be automatically put as part of the DATA PARAMETERS in the object 'list_param' below.
data_X <- read.table("D:/statistical-analysis/Otavio/Gabriel/ClockGenes_COVID_nonICU_WBL.txt", header = TRUE, dec = ",")[, 1:10]
data_Y <- read.table("D:/statistical-analysis/Otavio/Gabriel/Citocinas_COVID_nonICU_WBL.txt", header = TRUE, dec = ",")[, 1:10]

# this list of parameters will determine the behavior of the function runCCA (i.e. aspect of figures and numerical results of the CCA)
list_param <- list(
    # FIXED PARAMETERS
    empty = FALSE,                                                           # keep this parameter as FALSE
    # GENERIC PARAMETERS
    force_install_all_packages = FALSE,                                      # if packages do not load correctly, set this parameter to TRUE
    print_descriptive_statistics = TRUE,                                     # print out in the R console descriptive measures for each variable in 'data_X' and 'data_Y'  
    # PLOT: HEATMAP OF PEARSON CORRELATIONS
    heatmap_clust_method = "ward.D2",                                        # heatmap: method for clustering the left lower block of correlation matrix (i.e. only X against Y) 
    heatmap_color_ramp = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"),   # heatmap: gradient of colors for filling correlation matrix
    heatmap_border_color = grey(0.4),                                        # heatmap: color for borders of the correlation matrix
    heatmap_fontsize = 5,                                                    # heatmap: font size for the variables' names
    heatmap_width = 13,                                                      # figure dimension: width
    heatmap_height = 7,                                                      # figure dimension: height
    # PLOT: UNDIRECTED GRAPH OF PEARSON CORRELATIONS
    plotCorCCA_n_canonical_variates = 2,                                     # Correlation and CCA plots: number of canonical variates to plot (if greater than 2, it will override 'plots_canonical_variate' for the correlation graphs)
    plotCorCCA_corr_min = 0.7,                                               # Correlation and CCA plots: threshold determining the minimum correlation to be considered between variables
    plotCorCCA_seed = 10,                                                    # Correlation and CCA plots: choose any positive number whatsoever
    plotCorCCA_node_size = 5,                                                # Correlation and CCA plots: size of nodes drawn in the graph
    plotCorCCA_label_size = 3,                                               # Correlation and CCA plots: label size to be shown over the nodes
    plotCorCCA_net_method = "kamadakawai",                                   # Correlation and CCA plots: method for determining the placement of nodes in the graph (choose between: "kamadakawai" and "circle") 
    plotCorCCA_width = 13,                                                   # figure dimension: width
    plotCorCCA_height = 7,                                                   # figure dimension: height
    # PLOT: CORRELATION OF VARIABLES WITH CANONICAL VARIATES
    plotCCA_corr_min = 0.7,                                                  # CCA plot: threshold determining the minimum correlation to be considered between variables and canonical variates 
    plotCCA_width = 13,                                                      # CCA plot: figure dimension: width
    plotCCA_height = 7,                                                      # CCA plot: figure dimension: height
    plotCCA_max_overlaps = 1000,                                             # CCA plot: maximum number of labels that are aloud to overlap. reduce it for avoiding labels that overlap too much 
    # PLOT: ESTIMATES OF CANONICAL VARIATES (COEFFICIENTS AND CORRELATIONS)
    plotEst_standardize = FALSE,                                             # Should estimated coefficients for the canonical variates be standardized?: TRUE, FALSE
    plotEst_width = 20,                                                      # figure dimension: width
    plotEst_height = 15,                                                     # figure dimension: height
    # PLOT: CUMULATIVE VARIANCE OF CANONICAL VARIATES
    plotCumVar_width = 20,                                                   # figure dimension: width
    plotCumVar_height = 10,                                                  # figure dimension: height
    #PLOT: CANONICAL CORRELATIONS
    plotCC_width = 20,                                                       # figure dimension: width
    plotCC_height = 10,                                                      # figure dimension: height
    # GENERIC PLOT PARAMETERS
    plots_canonical_variate = c(1, 2),                                       # The two main canonical variates chosen for plotting
    plots_x_title = "Dataset X",                                             # title describing the first dataset (data_X)
    plots_y_title = "Dataset Y",                                             # title describing the second dataset (data_Y)
    plots_short_x_title = "x",                                               # short label for the first dataset (data_X)
    plots_short_y_title = "y",                                               # short label for the second dataset (data_Y)
    plots_x_color = "green",                                                 # color for variables of the first dataset (data_X)
    plots_y_color = "pink",                                                  # color for variables of the second dataset (data_Y)
    plots_cv_color = "gold",                                                 # color for canonical variates (and their associated short titles)
    plots_neg_color = "red",                                                 # color for highlighting negative correlations or negative coefficient estimates 
    plots_pos_color = "blue",                                                # color for highlighting positive correlations or positive coefficient estimates 
    plots_file_extension = "pdf",                                            # file extension for figures: "pdf" or "png"
    # DATA SPECIFICATIONS
    data_log_trans = FALSE,                                                  # log transformation for data: TRUE or FALSE
    data_log_base = 2,                                                       # log transformation for data: if data_log_trans = TRUE, choose a base for the log 
    data_sqrt_trans = FALSE,                                                 # square root transformation for data: TRUE or FALSE
    data_scale_trans = TRUE,                                                 # scale transformation for data: TRUE or FALSE (subtract mean and divide by standard deviation)
    data_X = data_X,                                                         # your first dataset (data_X): choose between an object holding your data or NA
    data_Y = data_Y,                                                         # your second dataset (data_X): choose between an object holding your data or NA                                              
    data_sheet1 = 1,                                                         # if data_X = NA, this parameter will point to the number of a spreadsheet (corresponding to data_X) inside an Excel file that you will be asked to choose (a pop-up will show up)
    data_sheet2 = 2                                                          # if data_X = NA, this parameter will point to the number of a spreadsheet (corresponding to data_Y) inside an Excel file that you will be asked to choose (a pop-up will show up)
  )

# CCA ---------------------------------------------------------------------------------------------------
## uncomment only one of the runCCA below and run it to start a CCA analysis

#cca <- runCCA(publish = FALSE) # run CCA using default values for parameters (list_param is not considered here)

#cca <- runCCA(list_param, publish = TRUE) # run CCA using custom values for parameters (using list_param)

#choose publish = TRUE if you have already found a good set of values for list_param and achieved a CCA you consider is great for publication.
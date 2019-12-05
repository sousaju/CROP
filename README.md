# CROP
R-script for paper **CROP: Correlation-based reduction of feature multiplicities in untargeted metabolomic data** (Š. Kouřil, J. de Sousa, J. Václavík, D. Friedecký and T. Adam; *in second review*)
***

**CROP** (**C**orrelation-based **R**emoval **O**f multi**P**licities) is a visual post-processing tool that removes redundant features from untargeted metabolomic data sets. It is based on a grouping of highly correlated features within a defined retention time window avoiding the condition of specific m/z difference making it a second-tier strategy for multiplicities reduction.
Graphical representation of correlation network for better understanding of the clusters composition and parameter tuning is provided.

![CROPped example data - correlation network](Example_correlation_network.PNG)
![CROPped example data - correlation network](example_data_CROPped_ccth_0.75_rtw+-0.02_correlation_network.pdf)

After CROPping your data set, you can directly continue with statistical pre-processing and analysis using our package [Metabol](https://github.com/AlzbetaG/Metabol).
***

### Input data
The user can freely use both mzTab and csv formatted data to be CROPped. 

- **mzTab**: For mzTab input, the CROP function returns data with modified SML table where all multiplicities are assigned to their cluster representatives in `SMF_ID_REFS` column and deleted from the data set. Data from `theoretical_neutral_mass` column are used as labels for the correlation network. Missing values in `abundance_assay` columns of SML table should be imputed before CROPping, otherwise they get treated as zeros when computing correlations. Missing values in `SML_ID`, `SMF_ID_REFS`, `SMF_ID`, and `retention_time_in_seconds/minutes` of SML and SMF table, respectively, result in errors since they prevent computing correlations in a given retention time window and/or assigning features to their clusters.

- **csv**: For csv input, the CROP function returns two csv files. The `output_table` stores a modified data set where all multiplicities are assigned to their cluster representant and deleted from the data set. The `list_of_clusters` provides contents of all identified clusters and their chosen representatives. The input table needs to have names/codes of features in the first column, retention time of features in minutes in second column, and peak areas of the samples (can include QCs) after imputation of missing values in the rest of the columns. Column header should contain names/codes of samples. Data from the first column are used as labels for the correlation network.

Please note that in both cases CROP is only usable if the majority of the abundances is filtered out by another preprocessing software first (otherwise many features belonging to different compounds would be connected into the “stretched clusters” due to similar retention times). 

***

### Usage
`CROP(mscsv="example_data.csv", name="project1", ccth=0.75, rtw=0.02, maxrtw = 0.04,  mcs=100, rtunit="min", funit="MW")`

### Parameters
* __`ccth`__

            threshold for correlation coefficient values
            default: ccth = 0.75

* __`funit`__
            
            units in which your features are assigned ("MW", "m/z" or anything else)
            only for column headers in the outputs
            default: funit = "MW"

* __`maxrtw`__

            maximal allowed retention time window to color stretched clusters phenomenon
            recommended to start with 2*rtw (higher number equals milder condition)
            default: maxrtw = NULL

* __`mcs`__

            maximal allowed cluster size
            should not be set high; if default is not enough, rather consider setting ccth bigger 
            and/or rtw smaller than changing mcs 
            default: mcs = 100

* __`mscsv`__

           filepath to an input mzTab or csv table of your MS data 
           if using csv, make sure you have:
              names of features in the first column
              retention time of features in minutes in second column
              samples in the rest of columns (including QCs)
              column header as names of samples
           if using mzTab, make sure you have:
              values in SML table after imputation of zeros
              no missing values in columns with SML_ID and SMF_ID_REFS of SML table
              no missing values in columns with SMF_ID and retention_time_in_seconds/minutes of SMF table

* __`name`__

            a note which will be in names of all files, i.e. a project name

* __`rtunit`__  

            units of retention time you are using in your data ("min" or "s")
            default: rtunit = "min"

* __`rtw`__     

            retention time window where +-rtw will be considered
            default: rtw = 0.02

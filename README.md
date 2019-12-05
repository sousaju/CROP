# CROP
R-script for paper **CROP: Correlation-based reduction of feature multiplicities in untargeted metabolomic data** (S. Kouril, J. Rendlova, J. Vaclavik, D. Friedecky and T. Adam; *submitted*)
***

CROP (**C**orrelation-based **R**emoval **O**f multi**P**licities) is a visual post-processing tool that removes redundant features from untargeted metabolomic data sets. It is based on a grouping of highly correlated features within a defined retention time window. Graphical representation of correlation network for better understanding of the clusters composition and parameter tuning is provided.
![CROPped example data - correlation network](example_data_CROPped_ccth_0.75_rtw+-0.02_correlation_network.pdf)
![CROPped example data - correlation network](Example_correlation_network.png)

After CROPping your data set, you can directly continue with statistical pre-processing and analysis using our package [Metabol](https://github.com/AlzbetaG/Metabol).
***

***

### Usage
```CROP(mscsv="example_data.csv", name="project1", ccth=0.75, rtw=0.02, maxrtw = 0.04,  mcs=100, rtunit="min", funit="MW")```

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
            recommended to set as max 2*rtw
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

# CROP
R-script for paper **CROP: Correlation-based reduction of feature multiplicities in untargeted metabolomic data** (S. Kouril, J. Rendlova, J. Vaclavik, D. Friedecky and T. Adam; *submitted*)
***

CROP (**C**orrelation-based **R**emoval **O**f multi**P**licities) is a visual post-processing tool that removes redundant features from untargeted metabolomic data sets. It is based on a grouping of highly correlated features within a defined retention time (RT) window. Graphical representation of correlation network for better understanding and parameter tuning is provided.
![CROPped example data - correlation network](example_data_CROPped_ccth_0.75_rtw+-0.02_correlation_network.pdf)
***

### Usage
```CROP(data="example_data.csv", name="project_name", ccth=0.75, rtw=0.02, maxrtw = NULL,  mcs=100, rtunit="min", funit="MW")```

### Parameters
* `ccth`

            threshold for correlation coefficient values
            default: `ccth = 0.75`

* `funit`
            
            units in which your features are measured ("MW", "m/z" or anything else)
            only for column headers in the outputs
            default: `funit = "MW"`

* `maxrtw`

            maximal allowed RT window to color stretched clusters phenomenon
            recommended to set as max `2\*rtw`
            default: `maxrtw = NULL`

* `mcs`

            maximal allowed cluster size
            should not be set high; if default is not enough, rather consider setting `ccth` bigger and/or `rtw` smaller than changing `mcs` 
            default: `mcs = 100`

* `data`

            filepath to an input csv table of your MS data with:
            names of features in the first column            
            retention time of features in minutes in second column          
            samples in the rest of columns (including QCs)           
            column header as names of samples

* `name`

            a note which will be in names of all files, i.e. a project name

* `rtunit`  

            units of retention time you ar using in your data ("min" or "s")
            default: `rtunit = "min"`

* `rtw`     

            retention time window where +-rtw will be considered
            default: `rtw = 0.02`

# Tuning Manifold Charting: the Official Code Base for the Paper 


## Dataset Structure
Our formatted reduced dataset is stored in `.mat` format as a structure array, which can be accessed from Matlab or Python. 

Each experimental session has one struct corresponding to Evolution, and one struct corresponding to Manifold experiment. 



## Code Structure
Our analysis code is written in Matlab and Python. Matlab for most statistical analysis, Python esp. pytorch for the *in silico* experiment and modelling.

### Statistical Analysis 
* Evolution trajectory, successfulness, best generation. 
* Tuning Map, basic stats (ANOVA, F).
* Kent fitting of tuning map, extract Kent parameters, population analysis of Kent parameter. 
* Radial Tuning curve, AUC for tuning curve. 
* Non-parametric statistics
    * Volume under the Surface (VUS). 

### Plotting 
* -> Figure 2B
* -> Figure 2C
* 


### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/Animadversio/Tuning-Manifold-Charting/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.

# Millennial

This model develops upon the conceptual framework presented in the Biogeochemistry paper Abramoff et al. 2018 (http://dx.doi.org/10.1007/s10533-017-0409-7). 
This code is publicly available but should be used for training purposes only until the Millennial Version 2 description and evaluation paper is published. Once published, the citation and licensing information will be added here. Until then, please let Rose Abramoff know if you have any activities planned using this model version.

The model is organized in folders, first by programming language and then by model version.
### Fortran:
#### MillennialV1 Files:
&nbsp;&nbsp;&nbsp; <b>main.F90</b> - original model code\
&nbsp;&nbsp;&nbsp; <b>simulation</b> - folder with model input, output, and run scripts\
&nbsp;&nbsp;&nbsp; <b>Table_3.xlsx, Table_3_with_params.xlsx</b> - some documentation

#### MillennialV2 Files:
&nbsp;&nbsp;&nbsp; <b>millennialv2.F90</b> - Rose's development version (subject to change)\
&nbsp;&nbsp;&nbsp; <b>simulationv2</b> - folder with model input, output, and run scripts

### R (under development):
#### models:
&nbsp;&nbsp;&nbsp; <b>derivs_Century.R</b> - Century model equations\
&nbsp;&nbsp;&nbsp; <b>derivs_V1.R</b> - Millennial V1 model equations\
&nbsp;&nbsp;&nbsp; <b>derivs_V2_ECA.R</b> - Millennial V2 model equations with ECA kinetics\
&nbsp;&nbsp;&nbsp; <b>derivs_V2_MM.R</b> - Millennial V2 model equations with Michaelis-Menten kinetics\
&nbsp;&nbsp;&nbsp; <b>derivs_V2_LIN.R</b> - Millennial V2 model equations with linear kinetics

#### simulation:
&nbsp;&nbsp;&nbsp; <b>run_functions.R</b> - run script

This repository branches from the original repository (https://github.com/email-clm/Millennial), containing the first version of the Millennial model framework written by Xiaofeng Xu. There are some differences between the equations presented in the Appendix of Abramoff et al. 2018 and the repository. They are noted as issues here (https://github.com/PNNL-TES/millenial/issues) where the model is being translated into R. This repository will further modify the original equations as part of ongoing model development.

## Original Readme

March 14, 2017
This is a repository for the newly developed Millennial model (version 1.0)

The Millennial model is a product of the third group of participants for the International Decade of Soils (IDOS) workshop held in Boulder, CO during Mar 14-16. The model was conceived by the entire group including Rose Abramoff, Sarah O'Brien, Eric Davidson, Wenting Feng, Adrien Finzi, Melannie Hartman, Daryl Moorhead, Josh Schimel, Melanie Mayes, and Xiaofeng Xu. The code is created in May - June 2016, by Xiaofeng Xu (xxu@sdsu.edu). The model has been updated with several modifications on processes and parameters during May 2016 - January 2017. 

Contact Xiaofeng Xu for any technical questions or if you want to test this model.

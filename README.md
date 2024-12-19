# Millennial
Nov 11, 2021\
Millennial Version 2 develops upon the conceptual framework presented in the Biogeochemistry paper Abramoff et al. (2018) (http://dx.doi.org/10.1007/s10533-017-0409-7). 
The Millennial Version 2 description and evaluation paper is published open access here: https://doi.org/10.1016/j.soilbio.2021.108466. You can cite the paper as:

Abramoff, Rose Z., Bertrand Guenet, Haicheng Zhang, Katerina Georgiou, Xiaofeng Xu, Raphael A. Viscarra Rossel, Wenping Yuan, and Philippe Ciais. "Improved global-scale predictions of soil carbon stocks with Millennial Version 2." Soil Biology and Biochemistry (2021) 164: 108466.

Please contact Rose Abramoff (rose.abramoff@igdore.org) with any questions, issues, comments, inquiries.

The model is organized in folders, first by programming language and then by model version.
### Fortran:
#### MillennialV1 Files:
&nbsp;&nbsp;&nbsp; <b>main.F90</b> - original model code\
&nbsp;&nbsp;&nbsp; <b>simulation</b> - folder with model input, output, and run scripts\
&nbsp;&nbsp;&nbsp; <b>Table_3.xlsx, Table_3_with_params.xlsx</b> - some documentation

#### MillennialV2 Files:
&nbsp;&nbsp;&nbsp; <b>millennialv2.F90</b> - main model script\
&nbsp;&nbsp;&nbsp; <b>simulationv2</b> - folder with model input, output, and run scripts

Similar to Millennial V1, Millennial V2's code is fully contained within the fortran script <b>millennialv2.F90</b>. The model can be compiled by calling from the folder <b>Fortran/MillennialV2/</b>:
```
gfortran -o name_of_executable millennialv2.F90
```
The model can be run by calling from the folder <b>Fortran/MillennialV2/simulationv2</b>:
```
../name_of_executable<run_control.txt
```
The executable takes one argument: <b>run_control.txt</b>\
This is an example of an input file that can be modified by the user.\
The arguments contained in the input file provide the information prompted by <b>millennialv2.F90</b> lines 429-443 and lines 520-537.\
Specifically,
1. Number of total simulation steps
2. Name of the parameter file
3. Whether or not to save output (1 = Yes, 0 = No)
4. Annual or daily output (1 = Annual, 0 = Daily)
5. Name of the file for saving model output
6. Name of the file for initializing the model pools
7. Input data containing the soil temperature in degrees Celsius, volumetric soil moisture in m3/m3, and NPP in gC/m2/day.

Examples of all files and their formats are provided in the folder <b>Fortran/MillennialV2/simulationv2</b>.

### R:
#### models:
&nbsp;&nbsp;&nbsp; <b>derivs_Century.R</b> - Century model equations\
&nbsp;&nbsp;&nbsp; <b>derivs_V1.R</b> - Millennial V1 model equations\
&nbsp;&nbsp;&nbsp; <b>derivs_V2_ECA.R</b> - Millennial V2 model equations with ECA kinetics\
&nbsp;&nbsp;&nbsp; <b>derivs_V2_MM.R</b> - Millennial V2 model equations with Michaelis-Menten kinetics\
&nbsp;&nbsp;&nbsp; <b>derivs_V2_LIN.R</b> - Millennial V2 model equations with linear kinetics

#### simulation:
&nbsp;&nbsp;&nbsp; <b>run_functions.R</b> - run script\
&nbsp;&nbsp;&nbsp; <b>model_tutorial.Rmd</b> - tutorial with explanation of how to run the model as well as examples and exercises

#### analysis:
&nbsp;&nbsp;&nbsp; <b>Site_Data_for_Millennial.Rmd</b> - Prepares input data needed to run the Millennial model as well as other data used for analysis\
&nbsp;&nbsp;&nbsp; <b>Millennial_Version_Analysis.Rmd</b> - Annotated analysis scripts\
&nbsp;&nbsp;&nbsp; <b>X.Rdata</b> - Multiple Rdata files called by <b>Site_Data_for_Millennial.Rmd</b>

This repository branches from the original repository (https://github.com/email-clm/Millennial), containing the first version of the Millennial model framework written by Xiaofeng Xu. There are some differences between the equations presented in the Appendix of Abramoff et al. (2018) and the repository. They are noted as issues here (https://github.com/PNNL-TES/millenial/issues) where the model is being translated into R. This repository may further modify the original equations as part of ongoing model development.

## Millennial V1 Archived Readme

March 14, 2017
This is a repository for the newly developed Millennial model (version 1.0)

The Millennial model is a product of the third group of participants for the International Decade of Soils (IDOS) workshop held in Boulder, CO during Mar 14-16. The model was conceived by the entire group including Rose Abramoff, Sarah O'Brien, Eric Davidson, Wenting Feng, Adrien Finzi, Melannie Hartman, Daryl Moorhead, Josh Schimel, Melanie Mayes, and Xiaofeng Xu. The code is created in May - June 2016, by Xiaofeng Xu (xxu@sdsu.edu). The model has been updated with several modifications on processes and parameters during May 2016 - January 2017. 

Contact Xiaofeng Xu for any technical questions or if you want to test this model.

## Notes
Dec 19, 2024
Kdes units used in Equations 10 and 12 should be gC/m2/d, and not mgC/L/d.

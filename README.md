# CO2BLOCKSEISM
An analytical tool for screening subsurface CO<sub>2</sub> storage resources constrained by induced seismicity

CO2BLOCKSEISM extends the tool [CO2BLOCK](https://github.com/co2block/CO2BLOCK), which estimates and optimizes CO<sub>2</sub> storage capacity subject to reservoir injectivity limitations. To learn more about the theory behind CO2BLOCKSEISM, please refer to the following paper:
- Kivi, I.R., De Simone, S. and Krevor, S., 2024. A simplified physics model for estimating subsurface CO<sub>2</sub> storage resources constrained ‎by fault slip potential. EarthArXiv. https://doi.org/10.31223/X5KM65

The CO2BLOCKSEISM tool is provided here in the context of two demonstration studies: [Oklahoma seismicity](https://github.com/imanrahimzadeh/CO2BLOCKSEISM/tree/main/Oklahoma%20seismicity) and [Utsira storage capacity](https://github.com/imanrahimzadeh/CO2BLOCKSEISM/tree/main/Utsira%20storage%20capacity). A general version of the code that can be applied to any study will be provided soon in this repository. 

Should you have any question/comment/suggestion, please contact Iman R. Kivi through the email i.rahimzadeh-kivi@imperial.ac.uk

CO2BLOCKSEISM is an open-source software. If you use it for academic purposes, please cite the reference paper above. 


**Code structure**

The tool comprises of a number scripts and functions as detailed below:
- "CO2BLOCKSEISM.m": this is the main script, where all input data for the study are specified. This script calls several functions to perform the required calculations,
- "read_data.m": this function reads the input parameters from the input excel file,
- "hydrogeology.m": this function involves analytical solutions of the pressure response of saline aquifers to multi-site CO<sub>2</sub> injection at time-varying rate,
- "Geomech_prob.m": this function develops a Monte Carlo-type probabilistic model for evaluating the probability of fault slip under inherent uncertainties of the geomechanical parameters of the subsurface,
- "Analytical_solution.m" and "FD_Nor_2zones.m": these functions calculate pressure changes in a saline aquifer in response to two-phase flow of CO<sub>2</sub> and water around a single injection site based on analytical solution by Nordbotten et al. (2005). A simplified solution for single phase flow of water is used for the demonstration study of induced seismicity in Oklahoma,
- "stress_projection.m": this function projects shear and normal stress components onto the fault planes for different realizations of the Monte Carlo simulation,
- "eos.m": this function returns brine viscosity and CO<sub>2</sub> density and viscosity using appropriate equations of state,
- "plot.m": this script visualizes the output data through appropriate plots. Scientific colormaps developed by Fabio Crameri (2018) is used in these plots. The functions, database and instructions of using the Crameri's colormaps are also provided.     

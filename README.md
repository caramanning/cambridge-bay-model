# cambridge-bay-model
model for methane cycling in Cambridge Bay, Nunavut

This model code is provided as a supplement to the following paper: 
Manning, C. C., Preston, V. L., Jones, S. F., Michel, A. P. M., Nicholson, D. P., Duke, P. J., et al. (2020). River inflow
dominates methane emissions in an Arctic coastal system. Geophysical Research Letters, 47, e2020GL087669. https://doi.org/
10.1029/2020GL087669

To run the code:
1) download this repository and add to your MATLAB path. Also download the Gas Toolbox (https://github.com/dnicholson/gas_toolbox) and add it to your path. 
2) run the m file CBfm.m.  The value of k_ox (oxidation rate constant) is set to 0 by default (matching Figure 4 in the paper) but can be adjusted to 0.01 d^-1 or 0.1 d^-1 to reproduce Figure S9 in the paper.

Troubleshooting: 
In certain versions of MATLAB you may need to convert CBw.DateTime to a char format. To do this, on line 119 of CBfm_importdata.m, change datenum(CBw.DateTime) to datenum(char(CBw.DateTime)).


The gas data used in the model was generated by the authors and is available on PANGAEA: 
Manning, C. C.; Preston, V. L., Jones, S. F., Michel, A. P. M., Nicholson, D. P., Duke, P. J., Ahmed, M. M. M., Manganini, K., Else, B. G. T., & Tortell, P. T. (2019): Dissolved methane, nitrous oxide, carbon dioxide, water isotope, salinity and temperature data from Cambridge Bay, Nunavut, Canada (2017-2018). PANGAEA, https://doi.org/10.1594/PANGAEA.907159

Other data sources are provided below. The preliminary 2018 discharge data for Freshwater Creek was provided by A. Pippy of Environment and Climate Change Canada (personal communication).

Environment and Climate Change Canada Historical Hydrometric Data website, River discharge data for station 10TF001, Freshwater Creek near Cambridge Bay https://wateroffice.ec.gc.ca/mainmenu/historical_data_index_e.html. Water Survey of Canada, Environment and Climate Change Canada. Downloaded on January 10, 2019.

Environment and Climate Change Canada Historical Database (2019) Hourly weather data for Cambridge Bay, Nunavut from Jan 1, 2018 to Jan 1, 2019. Station Climate ID: 2400602, WMO ID: 71288 TC ID: XCM  http://climate.weather.gc.ca. Downloaded on February 25, 2019.

Ocean Networks Canada Data Archive (2019) Corrected ice draft data from Cambridge Bay, Nunavut, from Jan 1, 2018 to Jan 1, 2019. http://www.oceannetworks.ca, Ocean Networks Canada, University of Victoria, Canada. Downloaded on May 25, 2019.


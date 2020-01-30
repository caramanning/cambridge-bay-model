% Matlab functions for modeling CH4 in Cambridge Bay estuary
% Author: Cara Manning, University of British Columbia, cmanning@eoas.ubc.ca
% Date updated: 30 Jan 2020

% To run this function you will need to install:
% Gibbs-SeaWater Oceanographic Toolbox http://www.teos-10.org/software.htm
% gas_toolbox https://github.com/dnicholson/gas_toolbox

% set variables
dt = 1/24; % model time increment in days
ydmod =datenum(2018,1,1):dt:datenum(2019,1,1); %model time steps
nd = length(ydmod);

k_ox = 0; % microbial oxidation first-order rate constant in day^-1
ox1 = datenum(2018,6,15); % date when oxidation should start (set to 0 prior to significant river inflow)
                           
nb = 10; % number of boxes in model
bsa = 5e5; % box surface area in m2
bd = 2; % box depth in meters
sasa = bsa.*nb; % study area surface area in m2
bv = bsa.*bd; % box volume

CH4_CB = nan(nd,nb); % initialize with NaN values in all box
CH4_CB(1,:) = 6e-6; % day 0 give concentration 6e-6 mol/m3
CH4_ice = 6e-6; % CH4 concentration in ice

% ice off dates for regions 1, 2 and 3
io1 = datenum(2018,6,25);
io2 = datenum(2018,7,8);
io3 = datenum(2018,7,20);

% ice on date at end of year for all regions
ft = datenum(2018,10,18);

% boxes included in each region for ice off
iob1 = 1:2;
iob2 = 3:10;
iob3 = 1000; % for now, no boxes use the third date

run('CBfm_importdata.m'); % import data from various files
run('CBfm_runmodel.m'); % run model for CH4 in estuary
run('CBfm_plot_top2panels.m'); % plot model and data

% now repeat the model incorporating the local ocean
nb = 130; % number of boxes in model, now includes 45 km2 coast oce + 5 km2 estuary
bsa = 5e5; % box surface area in m2
bd = 2; % box depth in meters
sasa = bsa.*nb; % study area surface area in m2
bv = bsa.*bd; % box volume

CH4_CB = nan(nd,nb); % initialize with NaN values in all box
CH4_CB(1,:) = 6e-6; % day 0 give concentration 6e-6 mol/m3
CH4_ice = 6e-6; % CH4 concentration in ice

io1 = datenum(2018,6,25);
io2 = datenum(2018,7,8);
io3 = datenum(2018,7,20);
ft = datenum(2018,10,18);

% boxes included in each region for ice off
iob1 = 1:2;
iob2 = 3:10;
iob3 = 11:nb; 

run('CBfm_importdata.m'); % import data from various files
run('CBfm_runmodel.m'); % run model for CH4 in estuary

% plot model results for domain including coastal ocean
if k_ox > 0
run('CBfm_plot_bottom2panels_withox.m'); % legend with microbial oxidation
else
run('CBfm_plot_bottom2panels.m'); % legend without microbial oxidation
    
end

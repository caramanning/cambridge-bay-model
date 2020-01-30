% script to import data

% load in gas concentration data
load CBg;

% load in Freshwater Creek discharge data 1970-2016
% the end of the spreadsheet had water level values for 2002-2016 so these
% have been removed
load freshwater_creek.mat;

% load in 2018 Freshwater Creek Discharge data
% see FWC_2017_2018.m
load FWCdischarge18.mat FWC18;

% load in weather data
load CBweather_WMO71288;

% load in ice data and interpolated values
load ONCice.mat;

% interpolate ice data
% date vector for ice thickness
dvi = datevec(ONCinterpice.time);
dni = datenum(dvi(:,1),dvi(:,2),dvi(:,3),dvi(:,4),dvi(:,5),dvi(:,6));

% modeled ice thickness
ice_mod = interp1(dni',ONCinterpice.data,ydmod);
ice_mod(isnan(ice_mod)) = 0; % set values at end to 0;

% calculate volume of water (melted ice) added to the box per time step
% equal to -1*(change in ice thickness).*0.9 as ice is ~10% less dense than
% water
dicedt = -1.*[0,diff(ice_mod)].*0.9; 

% -------------------------------------------
% move on to the time-series data
% find data for station B1
station_num = nan(size(CBg18.yyyy));
%B1=ismember(CBg18.station,'B1');
%station_num(B1) = 1;

ind_B1 = find(CBg18.station_num==1);
i_2 = find(CBg18.depth<=2);
B1_2 = intersect(ind_B1,i_2);

% salinity and temperature for station B1 at 2 m
B1_S = CBg18.ctd_S(B1_2);
B1_T = CBg18.ctd_T(B1_2);

% depth for B1 at 2 m
B1_2_D = 2.*ones(size(B1_S));

% equilibrium CH4 concentration in mol/m3
CH4_eq = CH4sol(B1_S,B1_T)'.*sw_dens(B1_S,B1_T,B1_2_D)./1e6;
CH4_eq(13:14) = CH4_eq(12);  % manually fill where T/S data is missing
CH4_CB_eq = interp1(CBg18.sampledn(B1_2),CH4_eq,ydmod)';

% since B1 data doesn't span the entire year we need to fill in model
% set values at beginning and end to the first or last non-NaN value
first_non_NaN_index= find(~isnan(CH4_CB_eq), 1);
CH4_CB_eq(1:first_non_NaN_index-1) = CH4_CB_eq(first_non_NaN_index);
CH4_CB_eq_lr = flipud(CH4_CB_eq);
first_non_NaN_index=find(~isnan(CH4_CB_eq_lr), 1);
CH4_CB_eq_lr(1:first_non_NaN_index-1) = CH4_CB_eq_lr(first_non_NaN_index);
CH4_CB_eq = flipud(CH4_CB_eq_lr);
CH4_CB_eq = repmat(CH4_CB_eq,1,nb);

% find data for station FWC/FWCB(same station)
FWC=find(CBg18.station_num==2|CBg18.station==3);
q = FWC;
c_CH4_FWC_18 = CBg18.c_CH4(q);
dnFWC_18 = [CBg18.yyyy(q),CBg18.mm(q),CBg18.dd(q),CBg18.HH(q),CBg18.MM(q),0.*CBg18.dd(q)];
dnFWC_18(isnan(dnFWC_18)) = 0;
datenumFWC_18 = datenum(dnFWC_18);
ydFWC_18 = datenumFWC_18;

% model FWC data by interpolation
FWCmod = interp1(ydFWC_18,c_CH4_FWC_18,ydmod);
FWC_molm3 = FWCmod ./1e9 * 1e3; % nmo/L * 10^3 L/m3 * (1 mol/1e9 nmol)

% since FWC data doesn't span the entire year we need to fill in model
% set values at beginning and end to the first or last non-NaN value
first_non_NaN_index= find(~isnan(FWC_molm3), 1);
FWC_molm3(1:first_non_NaN_index-1) = FWC_molm3(first_non_NaN_index);
FWC_molm3_lr = fliplr(FWC_molm3);
first_non_NaN_index= find(~isnan(FWC_molm3_lr), 1);
FWC_molm3_lr(1:first_non_NaN_index-1) = FWC_molm3_lr(first_non_NaN_index);
FWC_molm3 = fliplr(FWC_molm3_lr);

% -----------------------------------
% process river discharge data
doyCB = 1:366; % day 60 is Feb 29
dv = datevec(CB.daten);
CB.doy = datenum(0*dv(:,1),dv(:,2),dv(:,3));
for ii = 1:length(doyCB)
    d = CB.doy == ii & CB.PARAM == 1 & ~isnan(CB.Value);
    m(ii) = nanmean(CB.Value(d));
    s(ii) = nanstd(CB.Value(d));
    qh(ii) = quantile(CB.Value(d),0.75);
    ql(ii) = quantile(CB.Value(d),0.25);
end
m(m<0) = 0;
m(60) = []; % remove February 29
md = m; % mean discharge in m3/s

s(s<0) = 0;
s(60) = []; % remove February 29
stdd = s; % standard deviation of discharge in m3/s

dist = datenum(2018,1,1,12,0,0):datenum(2018,12,31,12,0,0);
% discharge per time step dt
disdt = interp1(dist,md,ydmod).*dt*86400; % discharge/second * time step (day) * (seconds/day)
disdt = naninterp(disdt);
disdt(isnan(disdt)) = 0;
% note: day 152 is the first one with measurable discharge

CBw = CBweather_WMO71288;

CBw_dn = [datenum(CBw.DateTime) 
    datenum(2019,1,1)]; % add on final date
CBw_u10 = CBw.WindSpdkmh.*1000/3600; % (1000 m/1 km) * (1 hr / 3600s)
CBw_u10 = [CBw_u10
    CBw_u10(end)]; % add in wind speed for 2019-01-01

CBw_yd = CBw_dn; 


u10mod = interp1(CBw_yd',CBw_u10',ydmod);
u10mod = naninterp(u10mod); 

% calculate schmidt number using fixed T and S
[~,Sc_CH4] = gasmoldiff(0,2,'CH4');
Sc_CH4 = Sc_CH4.*ones(size(u10mod));
kCH4 = kgas(u10mod,Sc_CH4,'W14').*86400; % m/d

% initial concentration in each box is 6 nmol/L = 6e-6 mol/m3
nd = length(u10mod);
CH4_CB = nan.*CH4_CB_eq;
CH4_CB(1,:) = 6e-6; % day 0 give concentration 6e-6 mol/m3

% set ice-corrected k values
kiceCH4 = repmat(kCH4,nb,1)';
kiceCH4(ydmod<io1,iob1) = 0; % box group 1 ice-covered Jan 1-June 25
kiceCH4(ydmod<io2,iob2) = 0; % box group 2 ice-covered Jan 1 through Jul 8

if nb>iob2
    kiceCH4(ydmod<io3,iob3) = 0; % box group 3 ice-covered Jan 1 through Jul 24
end

kiceCH4(ydmod>ft,:) = 0; % all boxes ice covered Oct 18 - Dec 31

% calculate dicedt for box closer to river that melts first
dicedt_box1 = interp1([ydmod-(io2-io1) max(ydmod)],[dicedt 0],ydmod);
dicedt_box1 = dicedt_box1';
% now model ice thickness
dicedt_m = repmat(dicedt,nb,1)';

for i =iob1
dicedt_m(:,i) = dicedt_box1;
end

% volume of ice added to each box per time step
vicedt_m = dicedt_m.*bsa;

% river input of CH4 per time step
% mol (hr)^-1 as the time step is hourly
CH4_in = disdt.*FWC_molm3;


% ChemYak mean and standard deviation in upper 1 m each day, taken from
% chemyak data files
cy_ch4 =[407 19
    362 75
 268 48 
198 23
151 68] ;

cy_dn = datenum(2018,6,28,12,0,0):1:datenum(2018,7,2,12,0,0); % approximate sampling time


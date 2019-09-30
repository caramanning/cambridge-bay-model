
% At t = 1, add river derived water and ice melt to box 1 and remove CH4 by
% gas exchange.
% The new CH4 concentration in box 1 is determined by mass conservation
%(assuming some water is lost from the mixed layer to balance the ice melt
% being added.
% I = ice melt volume (size T x N)  m3
%C(1,1) = [D(0) * Cr(0) + I(0,1) * Ci + (B - D(0) - I(0,1)) * C(0,1)]/B - Fc(0,1);
%CH4_CB(i+1,1) = [CH4_in(i)]; 


% fluxes in mol hr^-1 or mol (timestep)^-1
% positive flux adds CH4 to box
% negative flux removes CH4 from box

F_riv = zeros(size(CH4_CB)); % river flux in (will be NaN for all except box 1)
F_in = zeros(size(CH4_CB)); % flux in to each box from the adjacent box (due to river and ice input)
F_out_riv = zeros(size(CH4_CB)); %flux out to the adjacent box (due to river input)
F_out_ice = zeros(size(CH4_CB)); %flux out to adjacent box(due to ice melt)
F_ge = zeros(size(CH4_CB)); % flux due to gas exchange
F_ice =zeros(size(CH4_CB)); % flux in due to ice melt


for n = 1:nb
if n==1
for i = 1:nd-1
F_riv(i,n) = CH4_in(i);
F_out_riv(i,n) = -1*(disdt(i)).*CH4_CB(i,n);
F_out_ice(i,n) = -1*vicedt_m(i,n).*CH4_CB(i,n);
F_in(i,n) = 0;
F_ge(i,n) = kiceCH4(i,n).*(CH4_CB_eq(i,n)-CH4_CB(i,n)).*dt.*bsa;
F_ice(i,n) = vicedt_m(i,n).*CH4_ice;

CH4_CB(i+1,n) = (CH4_CB(i,n).*bv + F_out_riv(i,n) + F_out_ice(i,n)+F_riv(i,n)+F_in(i,n)+F_ge(i,n)+F_ice(i,n))./bv;

end

else 
    CH4_CB(n,n:nb) = CH4_CB(1,n:nb);
    for i = n:nd-1
        F_in(i,n) = (disdt(i)+vicedt_m(i,n-1)).*CH4_CB(i-1,n-1);        
        F_out_riv(i,n) = -1*(disdt(i)).*CH4_CB(i,n);
        F_out_ice(i,n) = -1*vicedt_m(i,n).*CH4_CB(i,n);
        F_ge(i,n) = kiceCH4(i,n).*(CH4_CB_eq(i,n)-CH4_CB(i,n)).*dt.*bsa;
        F_ice(i,n) = vicedt_m(i,n).*CH4_ice;
        CH4_CB(i+1,n) = (CH4_CB(i,n).*bv + F_out_riv(i,n) + F_out_ice(i,n)+F_in(i,n)+F_ge(i,n)+F_ice(i,n))./bv;
    end
end
end

% calculate fluxes per day
F_in_ps = F_in./dt;
F_out_riv_ps = F_out_riv./dt;
F_out_ice_ps = F_out_ice./dt;
F_ge_ps = F_ge./dt;
F_ice_ps = F_ice./dt;
F_riv_ps = F_riv./dt;

dn = nd/24;

% calculate max and min at each time step, for plotting
minCH4_CB = nan.*ydmod;
maxCH4_CB = nan.*ydmod;
for i = 1:length(CH4_CB)
    minCH4_CB(i) = min(CH4_CB(i,:));
    maxCH4_CB(i) = max(CH4_CB(i,:));
end


% calculate net fluxes over entire model domain
F_riv_sum = sum(F_riv(:,1));
F_ge_sum = sum(sum(F_ge));
F_out_riv_sum = sum(F_out_riv(:,nb));
F_out_ice_sum = sum(sum(F_out_ice));
F_ice_sum = sum(sum(F_ice));

% save out cumulative fluxes from river in, gas ex, river out, ice in, ice
% out
fluxes=[F_riv_sum F_ge_sum F_out_riv_sum F_ice_sum F_out_ice_sum];
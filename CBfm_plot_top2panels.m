xl1 = [150 250];
xl2 = [datenum(2018,5,29) datenum(2018,9,6)];
figure(4)
clf; 
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 8]);

subplot(2,2,1)
hold on; box on;
set(gca,'fontsize',12);
set(gca,'tickdir','out');
s1 = scatter(datenum(0*dv(:,1),dv(:,2),dv(:,3),dv(:,4),dv(:,5),dv(:,6)),CB.Value,20,'MarkerFaceColor',[0.5 0.5 0.5]);
s1.MarkerFaceAlpha = 0.2;
s1.MarkerEdgeColor = 'none';

h1=plot(1:365,md,'k','linewidth',1);
datetick;

set(gca,'ylim',[-1 120]);
set(gca,'xlim',xl1);

%plot(FWC17.doy,FWC17.discharge,'.');
h2=plot(FWC18.doy,FWC18.discharge,'-r');
ylabel('Discharge [m^3 s^{-1}]');

legend([h1,h2],'mean 1970-2016','2018');

set(gca,'xtick',[datenum(0,6,1) datenum(0,7,1) datenum(0,8,1) datenum(0,9,1)]);
%set(gca,'xticklabel',['Jun';'Jul';'Aug';'Sep']);
set(gca,'xticklabel',[]);



subplot(2,2,2)
hold on;
set(gca,'fontsize',12);
set(gca,'tickdir','out');

% plot polygon with min and max values for each model time step
h3 = plot([datenum(2018,1,1) datenum(2018,6,4)],[6 6],'color',[0.8 0.8 0.8],'linewidth',1.5);

for i =1:length(ydmod)
rectangle('Position',[ydmod(i)-0.5*dt,minCH4_CB(i).*1e6,dt,(maxCH4_CB(i)-minCH4_CB(i)).*1e6],'facecolor',[0.8 0.8 0.8],'edgecolor','none');
end

% plot range of values from chemyak
%for i = 1:5
%s1=rectangle('Position',[cy_dn(i)-0.5-datenum(2018,0,0) cy_ch4(i,1)-cy_ch4(i,2) 1 2.*cy_ch4(i,2)],'facecolor',[0 0 1],'edgecolor','none');
%end

h2=errorbar(cy_dn,cy_ch4(:,1),cy_ch4(:,2),'^','color',[0.2 0.6 1],'markerfacecolor',[0.2 0.6 1]);
plot(cy_dn,cy_ch4(:,1),'-','color',[0.2 0.6 1],'linewidth',1); % data, river

% plot equilibrium
plot(ydmod,CH4_CB_eq(:,1).*1e6,':','color',[0.5 0.5 0.5]); % equilibrium, estuary
h4=plot(ydmod,FWC_molm3.*1e6,'--k','linewidth',1); % data, river
h5=plot(CBg18.sampledn(B1_2)+1,CBg18.c_CH4(B1_2),'or','markerfacecolor','r'); % data, estuary

h6=plot(datenumFWC_18,c_CH4_FWC_18,'ok'); % data, river

box on;
datetick;
set(gca, 'YScale', 'log')
%set(gca,'xlim',[150 250]);

legend([h3,h4,h5,h2,h6],'model, estuary','model, river','data, estuary','data, ChemYak','data, river');


%legend('river','B1','model, box1','model, box2','model, box3','model, box4');set(gca,'xlim',[160 200]); % day 151 is June 1, day 243 is Aug 1

datetick;
ylabel('CH_4 [nM]');
box on;
datetick;
set(gca, 'YScale', 'log')
set(gca,'xlim',xl2); % day 151 is June 1, day 243 is Aug 1
set(gca,'xtick',[datenum(2018,6,1) datenum(2018,7,1) datenum(2018,8,1) datenum(2018,9,1)]);

set(gca,'xticklabel',[]);
ylim([1 10000]);
yticks([10^0 10^1 10^2 10^3 10^4]);
yticklabels({'10^0','10^1','10^2','10^3','10^4'});

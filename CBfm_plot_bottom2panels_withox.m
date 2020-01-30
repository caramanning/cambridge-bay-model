figure(4)
subplot(2,2,3)
hold on; box on;set(gca,'fontsize',12);
set(gca,'tickdir','out');
plot(ydmod,F_riv_ps(:,1).*16./1000,'color','k');
plot(ydmod,sum(F_ice_ps(:,1:max(iob2)),2).*16./1000,'color',[0 158 115]./255,'linewidth',1.5);
plot(ydmod,sum(F_ge_ps(:,1:max(iob2)),2).*16./1000,'color',[230,159,0]./255,'linewidth',1)
plot(ydmod,sum(F_ge_ps(:,iob3),2).*16./1000,'color',[0 114 178]./255,'linewidth',1)
plot(ydmod,sum(F_ox_ps(:,1:max(iob2)),2).*dt.*16./1000,'--','color',[0.5 0.5 0.5],'linewidth',1.5)
plot(ydmod,sum(F_ox_ps(:,iob3),2).*dt.*16./1000,':','color',[0.5 0.5 0.5],'linewidth',1.5)


ylabel('Modeled daily delivery [kg CH_4 d^{-1}]');
legend('river inflow','ice melt','gas exch, estuary','gas exch, coastal ocean','oxidation, estuary','oxidation, coastal ocean','location','southwest');

set(gca,'xlim',xl2); 
set(gca,'xtick',[datenum(2018,6,1) datenum(2018,7,1) datenum(2018,8,1) datenum(2018,9,1)]);
set(gca,'xticklabel',['Jun';'Jul';'Aug';'Sep']);
%%set(gca,'ylim',[-250 250]);
set(gca,'ylim',[-1000 100]);


subplot(2,2,4)
hold on; box on;set(gca,'fontsize',12);
set(gca,'tickdir','out');
plot(ydmod,cumsum(F_riv_ps(:,1)).*dt.*16./1000,'color','k');
plot(ydmod,cumsum(sum(F_ice_ps(:,1:max(iob2)),2)).*dt.*16./1000,'color',[0 158 115]./255,'linewidth',1.5);
plot(ydmod,cumsum(sum(F_ge_ps(:,1:max(iob2)),2)).*dt.*16./1000,'color',[230,159,0]./255,'linewidth',1.5)
plot(ydmod,cumsum(sum(F_ge_ps(:,iob3),2)).*dt.*16./1000,'color',[0 114 178]./255,'linewidth',1.5)
plot(ydmod,cumsum(sum(F_ox_ps(:,1:max(iob2)),2)).*dt.*16./1000,'--','color',[0.5 0.5 0.5],'linewidth',1.5)
plot(ydmod,cumsum(sum(F_ox_ps(:,iob3),2)).*dt.*16./1000,':','color',[0.5 0.5 0.5],'linewidth',1.5)

ylabel('Modeled cumulative delivery [kg CH_4]');

set(gca,'xlim',xl2); 
set(gca,'xtick',[datenum(2018,6,1) datenum(2018,7,1) datenum(2018,8,1) datenum(2018,9,1)]);
set(gca,'xticklabel',['Jun';'Jul';'Aug';'Sep']);
set(gcf,'renderer','Painters')
set(gca,'ylim',[-1000 1000]);
wysiwyg;

%print -depsc -r300 modeloutput_withox_v3.eps;
%print -dpdf -r300 modeloutput_withox_v3.pdf;

%print(gcf,'-depsc','-painters','out.eps');
%epsclean('out.eps'); % cleans and overwrites the input file
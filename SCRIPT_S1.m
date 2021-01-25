clear, clc

%Import Cac OD600 Data
CACOD = xlsread('TABLE S1.xlsx','Cac OD600','B2:B7');
CACtime = xlsread('TABLE S1.xlsx','Cac OD600','A2:A7');

%Import Clj OD600 Data
CLJOD = xlsread('TABLE S1.xlsx','Clj OD600','B2:B8');
CLJtime = xlsread('TABLE S1.xlsx','Clj OD600','A2:A8');

%Import Co-culture OD600 Data
COCOD = xlsread('TABLE S1.xlsx','Co-culture OD600','B2:B7');
COCtime = xlsread('TABLE S1.xlsx','Co-culture OD600','A2:A7');

%Import non-hybrid Clj relative abundnace from flow cytometry experiments
[CLJ_plot,~,~] = xlsread('TABLE S1.xlsx','RNA Exchange','B13:B19');
%Import non-hybrid Cac relative abundance from flow cytometry experiments
[CAC_plot,~,~] = xlsread('TABLE S1.xlsx','RNA Exchange','C13:C19');
%Import hybrid cell relative abundance from flow cytometry experiments
[MIXED_plot,~,~] = xlsread('TABLE S1.xlsx','RNA Exchange','D13:D19');
[t_RNA,~,~] =  xlsread('TABLE S1.xlsx','RNA Exchange','A13:A19');

%Import Cac/Clj ratio calculated from genome copy number
[GC_ratio,~,~] = xlsread('TABLE S1.xlsx','Genome Copy Number','G2:G6');
%Import Cac/Clj ratio standard deviation
[GC_stDEV,~,~] = xlsread('TABLE S1.xlsx','Genome Copy Number','H2:H6');
[GC_time,~,~] = xlsread('TABLE S1.xlsx','Genome Copy Number','A2:A6');
%Import genome copy number for fusion parameter conersion to cell/L
[GC_Cac,~,~] = xlsread('TABLE S1.xlsx','Genome Copy Number','B2:B2');
[GC_Clj,~,~] = xlsread('TABLE S1.xlsx','Genome Copy Number','D2:D2');

%Calculate relative abundances and standard deviations  from genome copy 
%number data
Cac_GC = (GC_ratio )./(GC_ratio + 1);
GC_stDEV = GC_stDEV./(GC_ratio + 1);

%interpolate OD600 values
tq = [0:0.1:100];
yq1 = interp1(CACtime,CACOD,tq,'pchip');
yq2 = interp1(CLJtime,CLJOD,tq,'pchip');
yq3 = interp1(COCtime,COCOD,tq,'pchip');


%Calculate time-dependent C. acetobutylicum specific growth rate
for i = 2:length(yq1)
    mu1(i,1) = log(yq1(i)/yq1(i-1))/(tq(i)-tq(i-1));
end
mu1 = min(mu1,3);
k = 2;
for i = 1:201
    muset1_init(i) = min(mu1(k),3.0);
    k = k + 1;
end
%Calculate time-dependent C. ljungdahlii specific growth rate
for i = 2:length(yq2)
    mu2(i,1) = log(yq2(i)/yq2(i-1))/(tq(i)-tq(i-1));
end
k = 2;
for i = 1:201
    muset2_init(i) = min(mu2(k),3.0);
    k = k + 1;
end
%Calculate time-dependent co-culture specific growth rate
for i = 2:length(yq2)
    mu3(i,1) = log(yq3(i)/yq3(i-1))/(tq(i)-tq(i-1));
end
k = 2;
for i = 1:201
    muset3_init(i) = min(mu3(k),3.0);
    k = k + 1;
end

%Initialize storage variable
fullx = zeros(0,301);

%Define fusion parameter
f = 0.00564;
%Convert fusion parameter to cell/L, assume genome copy number
%concentration is equivalent to cell concentration
muL_conversion = 1e-6;
f_converted = f*100/((GC_Cac/muL_conversion)+(GC_Clj/muL_conversion));

for i = 1:501
    if i == 1
        %initial cell abundance conditions, growth rates
        x0 = [10; 90; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
        tspan = [0 i*0.1];
        mu = [muset1_init(i) muset2_init(i) muset3_init(i)];
    else
        %iteratively updated cell abundances, growth rates, integration 
        %over 0.1 hr intervals
        tspan = [(i-1)*0.1 i*0.1];
        x0 = xa(end,:);
        mu = [mu1(i) mu2(i) mu3(i)];
    end
    
    %Cell death is not included in growth model
    mu = max(mu,0);
    
    %System of ODEs representing hybrid/non-hybrid cell growth and cell
    %fusion
    g=@(t,x) [mu(1)*max(x(1),0) - f*x(1)*max(x(2),0) + mu(1)*max(x(11),0);mu(2)*max(x(2),0) - f*x(1)*max(x(2),0) + mu(2)*x(12);...
    f*x(1)*max(x(2),0) - max(mu(1),0)*max(x(3),0);f*x(1)*max(x(2),0) - mu(2)*x(4);mu(1)*(x(3)-max(x(5),0));mu(2)*(x(4)-x(6));...
    mu(1)*(max(x(5),0)-x(7));mu(2)*(x(6)-x(8));mu(1)*(x(7)-max(x(9),0));mu(2)*(x(8)-x(10));mu(1)*(max(x(9),0)-x(11));mu(2)*(x(10)-x(12))];
    
    %ODE integration
    [t,xa] = ode45(g,tspan,x0);
    
    %Store solution
    fullx = [fullx;xa];
    
    %store ODE timesteps
    if i == 1
        fullt = t;
    else
        fullt = [fullt;t];
    end
end

%Calculate total Cac abundance
plotCAC = fullx(:,1) + fullx(:,3) + fullx(:,5) + fullx(:,7) + fullx(:,9) + fullx(:,11);
%Calculate total Clj abundance
plotCLJ = fullx(:,2) + fullx(:,4) + fullx(:,6) + fullx(:,8) + fullx(:,10) + fullx(:,12);
%Calculate total cell abundance
total = plotCAC + plotCLJ;
%Calculate total non-hybrid cell abundance
pure = fullx(:,1) + fullx(:,2);
%Calculate total hybrid cell abundance
mixed = fullx(:,3) + fullx(:,4) + fullx(:,5) + fullx(:,6) + fullx(:,7) + fullx(:,8) + fullx(:,9) + fullx(:,10) + fullx(:,11) + fullx(:,12);

%Calculate Cac relative abundance
Cac_frac = (plotCAC)./total;
%Calculate Clj relative abundance
Clj_frac = (plotCLJ)./total;
%Calculate total hybrid cell relative abundance
plotMIXED = mixed./total;

%Pull Cac, Clj relative abundances to compare to experimental data
Cac_genome = [Cac_frac(528);Cac_frac(938);Cac_frac(1266);Cac_frac(2008)];
Clj_genome = [Clj_frac(528);Clj_frac(938);Clj_frac(1266);Clj_frac(2008)];    


plotT = fullt;

plotpureCAC = fullx(:,1);
plotpureCLJ = fullx(:,2);

%Calculate non-hybrid Cac, non-hybrid Clj relative abundance
pureCACfrac = plotpureCAC./(plotCAC + plotCLJ);
pureCLJfrac = plotpureCLJ./(plotCAC + plotCLJ);



%Calculate relative abundance of each organism in hybrid and non-hybrid
%states
plot1 = fullx(:,1)./plotCAC;
plot2 = fullx(:,2)./plotCLJ;
plot3 = fullx(:,3)./plotCAC;
plot4 = fullx(:,4)./plotCLJ;
plot5 = fullx(:,5)./plotCAC;
plot6 = fullx(:,6)./plotCLJ;
plot7 = fullx(:,7)./plotCAC;
plot8 = fullx(:,8)./plotCLJ;
plot9 = fullx(:,9)./plotCAC;
plot10 = fullx(:,10)./plotCLJ;
plot11 = fullx(:,11)./plotCAC;
plot12 = fullx(:,12)./plotCLJ;


%Calculate weighted squared residual
errCac = ((Cac_GC(2) - Cac_frac(4940))/GC_stDEV(2))^2 + ((Cac_GC(3) - Cac_frac(9040))/GC_stDEV(3))^2 + ((Cac_GC(4) - Cac_frac(12320))/GC_stDEV(4))^2 + ((Cac_GC(5) - Cac_frac(19700))/GC_stDEV(5))^2;
errClj = ((1-Cac_GC(2) - Clj_frac(4940))/GC_stDEV(1))^2 + ((1-Cac_GC(3) - Clj_frac(9040))/GC_stDEV(3))^2 + ((1 - Cac_GC(4) - Clj_frac(12320))/GC_stDEV(4))^2 + ((1 - Cac_GC(5) - Clj_frac(19700))/GC_stDEV(5))^2;
err_total = errCac + errClj;



%Plot main manuscript Figure 3
figure;
plot(plotT,plot1,':','Color',[0/255,176/255,80/255],'LineWidth',3)
hold on
plot(plotT,plot3,'-','Color',[255/255,255/255,0/255],'LineWidth',3)
hold on
plot(plotT,plot5,'-','Color',[204/255,239/255,16/255],'LineWidth',3)
hold on
plot(plotT,plot7,'-','Color',[153/255,223/255,32/255],'LineWidth',3)
hold on
plot(plotT,plot9,'-','Color',[101/255,207/255,48/255],'LineWidth',3)
hold on
plot(plotT,plot11,'-','Color',[50/255,191/255,64/255],'LineWidth',3)
xlim([0 50])
ylim([0 1])
xlabel('\bfTime (hr)')
ylabel('\bfFractional abundance')
title('\bfa. {\itCac} in non-hybrid and hybrid states')
legend('Non-hybrid {\itCac}','{\itCac} Hybrid state 1','{\itCac} Hybrid state 2','{\itCac} Hybrid state 3','{\itCac} Hybrid state 4','{\itCac} Hybrid state 5')
set(gca,'FontSize',24)

figure;
plot(plotT,plot2,':','Color',[255/255 0/255 0/255],'LineWidth',3)
hold on
plot(plotT,plot4,'-','Color',[255/255,255/255,0/255],'LineWidth',3)
hold on
plot(plotT,plot6,'-','Color',[255/255,204/255,0/255],'LineWidth',3)
hold on
plot(plotT,plot8,'-','Color',[255/255,153/255,0/255],'LineWidth',3)
hold on
plot(plotT,plot10,'-','Color',[255/255,101/255,0/255],'LineWidth',3)
hold on
plot(plotT,plot12,'-','Color',[255/255,50/255,0/255],'LineWidth',3)
xlim([0 50])
ylim([0 1])
xlabel('\bfTime (hr)')
ylabel('\bfFractional abundance')
title('\bfb. {\itClj} in non-hybrid and hybrid states')
legend('Non-hybrid {\itClj}','{\itClj} Hybrid state 1','{\itClj} Hybrid state 2','{\itClj} Hybrid state 3','{\itClj} Hybrid state 4','{\itClj} Hybrid state 5')
%h = suptitle({'\bfPure and Mixed Cell Fractional Abundances vs. Time';'(Predicted Dilution Pool Abundances)'});
%set(h,'FontSize',20)
set(gca,'FontSize',24)


%Plot main manuscript figure 1
figure;
h(1) = plot(plotT,Cac_frac,'-','Color',[0/255,176/255,80/255],'LineWidth',3);
hold on
h(2) = plot(plotT,Clj_frac,'-','Color',[255/255 0/255 0/255],'LineWidth',3);
hold on
k = 4;
for i = 1:5
    h(k) = plot([GC_time(i), GC_time(i)],[Cac_GC(i)-GC_stDEV(i),Cac_GC(i)+GC_stDEV(i)],'-','Color',[0/255,176/255,80/255],'LineWidth',1);
    hold on
    k = k + 1;
    h(k) = plot([GC_time(i) - 1, GC_time(i) + 1],[Cac_GC(i)-GC_stDEV(i),Cac_GC(i)-GC_stDEV(i)],'-','Color',[0/255,176/255,80/255],'LineWidth',1);
    hold on
    k = k + 1;
    h(k) = plot([GC_time(i) - 1, GC_time(i) + 1],[Cac_GC(i)+GC_stDEV(i),Cac_GC(i)+GC_stDEV(i)],'-','Color',[0/255,176/255,80/255],'LineWidth',1);
    hold on
    k = k + 1;
    h(k) = plot([GC_time(i), GC_time(i)],[Cac_GC(i),Cac_GC(i)],'*','Color',[0/255,176/255,80/255],'LineWidth',3);
    hold on
    k = k + 1;
    h(k) = plot([GC_time(i), GC_time(i)],[1-Cac_GC(i)-GC_stDEV(i),1-Cac_GC(i)+GC_stDEV(i)],'-','Color',[255/255 0/255 0/255],'LineWidth',1);
    hold on
    k = k + 1;
    h(k) = plot([GC_time(i) - 1, GC_time(i) + 1],[1-Cac_GC(i)-GC_stDEV(i),1-Cac_GC(i)-GC_stDEV(i)],'-','Color',[255/255 0/255 0/255],'LineWidth',1);
    hold on
    k = k + 1;
    h(k) = plot([GC_time(i) - 1, GC_time(i) + 1],[1-Cac_GC(i)+GC_stDEV(i),1-Cac_GC(i)+GC_stDEV(i)],'-','Color',[255/255 0/255 0/255],'LineWidth',1);
    hold on
    k = k + 1;
    h(k) = plot([GC_time(i), GC_time(i)],[1-Cac_GC(i),1-Cac_GC(i)],'*','Color',[255/255 0/255 0/255],'LineWidth',3);
    hold on
    k = k + 1;
end


ylim([0 1])
xlabel('\bfTime (hr)')
ylabel('\bfFractional Abundance')
xlim([-1,50]);
legend(h([1 2]),'{\itCac} genome','{\itClj} genome','TextColor','k','EdgeColor','k','Color','None');
title({'\bf\color{black}Predicted vs. Experimental Genome Relative Abundance'})
set(gca,'FontSize',24)
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gca,'Color','None')


%Plot main manuscript Figure 2
figure;
p(1) = plot(plotT,pureCACfrac,'Color',[0/255 176/255 80/255],'LineWidth',3);
hold on
p(2) = plot(plotT,pureCLJfrac,'Color',[255/255 0/255 0/255],'LineWidth',3);
hold on

p(3) = plot(plotT,plotMIXED,'Color',[255/255,255/255,0/255],'LineWidth',3);
hold on

p(4) = plot(t_RNA,CLJ_plot,'*','Color',[255/255 0/255 0/255],'LineWidth',15);
hold on
p(5) = plot(t_RNA,CAC_plot,'*','Color',[0/255 176/255 80/255],'LineWidth',15);
hold on
p(6) = plot(t_RNA,MIXED_plot,'*','Color',[255/255,255/255,0/255],'LineWidth',15);
hold on
xlim([0 50])
ylim([0 1])
ylabel('\bfFractional Abundance')
xlabel('\bfTime (hr)')
title('\bfCell Type Abundance')
legend([p(1) p(5) p(2) p(4) p(3) p(6)],'Predicted non-hybrid {\itCac}','Experimental non-hybrid {\itCac}','Predicted non-hybrid {\itClj}','Experimental non-hybrid {\itClj}','Predicted hybrid','Experimental hybrid')
set(gca,'FontSize',24)


%Calculate squared residual between model prediction and experimental
%relative abundances measured using flow cytometry
err_CAC = abs(pureCACfrac(836) - CAC_plot(2))^2 + abs(pureCACfrac(1656) - CAC_plot(3))^2 + abs(pureCACfrac(3706) - CAC_plot(4))^2 + abs(pureCACfrac(8216) - CAC_plot(5))^2 + abs(pureCACfrac(11086) - CAC_plot(6))^2 + abs(pureCACfrac(16414) - CAC_plot(7))^2
err_CLJ = abs(pureCLJfrac(836) - CLJ_plot(2))^2 + abs(pureCLJfrac(1656) - CLJ_plot(3))^2 + abs(pureCLJfrac(3706) - CLJ_plot(4))^2 + abs(pureCLJfrac(8216) - CLJ_plot(5))^2 + abs(pureCLJfrac(11086) - CLJ_plot(6))^2 + abs(pureCLJfrac(16414) - CLJ_plot(7))^2
err_MIX = abs(plotMIXED(836) - MIXED_plot(2))^2 + abs(plotMIXED(1656) - MIXED_plot(3))^2 + abs(plotMIXED(3706) - MIXED_plot(4))^2 + abs(plotMIXED(8216) - MIXED_plot(5))^2 + abs(plotMIXED(11086) - MIXED_plot(6))^2 + abs(plotMIXED(16414) - MIXED_plot(7))^2
error1 = err_CAC + err_CLJ + err_MIX;




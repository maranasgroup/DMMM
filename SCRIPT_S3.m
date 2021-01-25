clear, clc

%Add path to cobratoolbox
addpath /path/to/cobratoolbox/
addpath /path/to/CPLEX/

%Initialize cobratoolbox
initCobraToolbox


changeCobraSolver ('gurobi', 'all');
global CBTDIR


%Import concentration data
[concentration,~,~] = xlsread('TABLE S1.xlsx','Fermentation Profiles','A2:L8');
%Import genome copy number for fusion parameter conersion to cell/L
[GC_Cac,~,~] = xlsread('TABLE S1.xlsx','Genome Copy Number','B2:B2');
[GC_Clj,~,~] = xlsread('TABLE S1.xlsx','Genome Copy Number','D2:D2');
%Interpolate concentration data
tq = [0:0.1:100];
yq1 = interp1(concentration(:,1),concentration(:,2),tq,'pchip');
yq2 = interp1(concentration(:,1),concentration(:,3),tq,'pchip');
yq3 = interp1(concentration(:,1),concentration(:,4),tq,'pchip');
yq4 = interp1(concentration(:,1),concentration(:,5),tq,'pchip');
yq5 = interp1(concentration(:,1),concentration(:,6),tq,'pchip');
yq6 = interp1(concentration(:,1),concentration(:,7),tq,'pchip');
yq7 = interp1(concentration(:,1),concentration(:,8),tq,'pchip');
yq8 = interp1(concentration(:,1),concentration(:,9),tq,'pchip');
yq9 = interp1(concentration(:,1),concentration(:,10),tq,'pchip');
yq10 = interp1(concentration(:,1),concentration(:,11),tq,'pchip');
yq11 = interp1(concentration(:,1),concentration(:,12),tq,'pchip');
pp_pchip_gluc = pchip(tq,yq1);
pp_pchip_fruc = pchip(tq,yq2);
pp_pchip_lac = pchip(tq,yq3);
pp_pchip_ac = pchip(tq,yq4);
pp_pchip_bdo = pchip(tq,yq5);
pp_pchip_acetoin = pchip(tq,yq6);
pp_pchip_etoh = pchip(tq,yq7);
pp_pchip_isoprop = pchip(tq,yq8);
pp_pchip_but = pchip(tq,yq9);
pp_pchip_acetone = pchip(tq,yq10);
pp_pchip_btoh = pchip(tq,yq11);
data(:,1) = tq(1:1000);
data(:,2) = pp_pchip_gluc.coefs(:,3);
data(:,3) = pp_pchip_fruc.coefs(:,3);
data(:,4) = pp_pchip_lac.coefs(:,3);
data(:,5) = pp_pchip_ac.coefs(:,3);
data(:,6) = pp_pchip_isoprop.coefs(:,3);
data(:,7) = pp_pchip_acetoin.coefs(:,3);
data(:,8) = pp_pchip_bdo.coefs(:,3);
data(:,9) = pp_pchip_btoh.coefs(:,3);
data(:,10) = pp_pchip_but.coefs(:,3);
data(:,11) = pp_pchip_etoh.coefs(:,3);
data(:,12) = pp_pchip_acetone.coefs(:,3);

%Indecies used to partition glucose (1) and fructose (2) uptake according
%to each stoichiometric model in DMMM framework
Cacpureind = [1 2];
Cljpureind = [2];
Cacfusedind = [1 2];
Cljfusedind = [1 2];

%List of cell types in simulation
species = {'Cacpure';'Cljpure'};

%0.1 hour timesteps
dt = 0.1;

%Define initial concentrations of each extracellular carbon source/sink
C_gluc0 = concentration(1,2);
C_fruc0 = concentration(1,3);
C_lac0 = concentration(1,4);
C_ac0 = concentration(1,5);
C_bdo0 = concentration(1,6);
C_acetoin0 = concentration(1,7);
C_etoh0 = concentration(1,8);
C_isoprop0 = concentration(1,9);
C_but0 = concentration(1,10);
C_acetone0 = concentration(1,11);
C_btoh0 = concentration(1,12);

%Define reaction list, index, and molecular weight for FVA, maximum titer
%calculations
rxnsList = {'ExCom_lac__L[u]';'ExCom_ac[u]';'ExCom_2ppoh[u]';'ExCom_actn__R[u]';'ExCom_btd_RR[u]';'ExCom_btoh[u]';'ExCom_but[u]';'ExCom_etoh[u]';'ExCom_acetone[u]'};
rxn_index = [5072; 4987; 4985; 4991; 5009; 5011; 5012; 5034; 4988];
MW = [0.09008; 0.060052; 0.0601; 0.08811; 0.090121; 0.074121; 0.08811; 0.04607; 0.05808];
C0 = [C_lac0; C_ac0; C_isoprop0; C_acetoin0; C_bdo0; C_btoh0; C_but0; C_etoh0; C_acetone0];
%Define initial glucose, fructose concentrations
Cgluc(1,1) = C_gluc0;
Cfruc(1,1) = C_fruc0;

%Define intial cell type total abundnace and fractional abundance
XX(1,:) = [10 90];
X(1,:) = [0.1 0.9];

%Define fusion parameter
f = 0.00564;
%Convert fusion parameter to cell/L, assume genome copy number
%concentration is equivalent to cell concentration
muL_conversion = 1e-6;
f_converted = f*100/((GC_Cac/muL_conversion)+(GC_Clj/muL_conversion));

%Define growth efficiency 
bm_eff = 1;

for k = 1:331
    %Define total glucose, fructose uptake at each iteration of loop
    ut(1,1) = data(k,2);
    ut(2,1) = data(k,3);
    
    %Exclude non-hybrid Clj from consideration when partitioning glucose
    %consumption between cell type/organisms
    frac(1,1) = [X(k,1)];
    
    %Partition fructose consumption between all cell type/organisms
    frac(2,1) = sum(X(k,:));
    
    %Initialize summed NGAM by hybrid Cac, Clj, to be used when predicting
    %fermentation yields
    CACfused_ATPuse = 0;
    CLJfused_ATPuse = 0;
    
    %Maximize growth for each organism/cell type included in the simulation
    %(Non-hybrid Cac, Non-hybrid Clj, 5 hybrid metabolism models for both
    %Cac and Clj
    for i = 1:length(species)
        
        %Flux balance analysis for Non-hybrid Cac
        if strcmp(species{i},'Cacpure') == 1
            %load model
            load model.mat
            %Assign substrate utilization based on the fractional abundance
            %of the organism/cell type in the culture
            for j = 1:length(eval(strcat(species{i},'ind')))
                if Cacpureind(j) == 1
                    if X(k,i) == 0
                        modelJoint.lb(4639) = 0; %Cacpure glucose uptake
                        modelJoint.lb(5046) = 0; %Community glucose uptake
                    else
                        modelJoint.lb(4639) = ut(Cacpureind(j))*(X(k,i)/frac(1,1)); %Cacpure glucose uptake
                        modelJoint.lb(5046) = ut(Cacpureind(j))*(X(k,i)/frac(1,1)); %Community glucose uptake
                    end
                else
                    if X(k,i) == 0
                        modelJoint.lb(4557) = 0; %Cacpure fructose uptake
                        modelJoint.lb(5041) = 0; %Community fructose uptake
                    else
                        modelJoint.lb(4557) = ut(Cacpureind(j))*X(k,i); %Cacpure fructose uptake
                        modelJoint.lb(5041) = ut(Cacpureind(j))*X(k,i); %Community fructose uptake
                    end
                end
            end
            
            %Block network activity for all other organisms/cell types in
            %the simulation
            for r = 1476:4550
                modelJoint.lb(r) = 0; %set all other species to zero
                modelJoint.ub(r) = 0; %set all other species to zero
            end
            
            %Assign NGAM based on the fractional abundance the 
            %organism/cell type occupies in the co-culuture
            modelJoint.lb(1100) = 8.39*X(k,i); %Cacpure ATPM
            %Store ATPM for use in fermentation yield FVA
            CACpure_ATPstore = modelJoint.lb(1100);%Cacpure ATPM

            %Perfrom flux balance analysis, store biomass flux
            modelJoint = changeObjective (modelJoint, 'CacpureEx_cpd11416_norm2[u]');
            FBA = optimizeCbModel (modelJoint, 'max');
            mu(k,i) = FBA.obj;
            modelJoint.lb(4668) = bm_eff*FBA.obj;
            
            %Identify and store feasible ranges for H2, CO2 export by 
            %non-hybrid Cac
            modelJoint = changeObjective (modelJoint, 'CacpureEx_h2[u]');
            h2_store = optimizeCbModel (modelJoint, 'max');
            h2_store_var = h2_store.obj;
            modelJoint = changeObjective (modelJoint, 'CacpureEx_co2[u]');
            co2_store = optimizeCbModel (modelJoint, 'max');
            co2_store_var = co2_store.v(4653); %Cac pure CO2 evolution-->Change to [e] ->[u]
            
            %Store non-hybrid Cac glucose and fructose uptake for use in
            %fermentation product FVA
            gluc_flux(k,i) = FBA.full(5046); %Community glucose uptake
            fruc_flux(k,i) = FBA.full(5041); %Community fructose uptake
            
            %X(k,i) = X(k,i)/frac(2,1);
        %Flux balance analysis for Non-hybrid Clj
        elseif strcmp(species{i},'Cljpure') == 1
            load model.mat
            %Assign substrate utilization based on the fractional abundance
            %of the organism/cell type in the culture
            for j = 1:length(eval(strcat(species{i},'ind')))
                if X(k,i) == 0
                    modelJoint.lb(4695) = 0; %Cljpure fructose uptake
                    modelJoint.lb(5041) = 0; %Community fructose uptake
                else
                    modelJoint.lb(4695) = ut(Cljpureind(j))*X(k,i); %Cljpure fructose uptake
                    modelJoint.lb(5041) = ut(Cljpureind(j))*X(k,i); %Community fructose uptake
                end
                modelJoint.lb(4699) = 0; %Cljpure glucose uptake
                modelJoint.lb(5046) = 0; %Community fructose uptake
            end
            %Block non-hybrid Cac flux
            for r = 1:1475
                modelJoint.lb(r) = 0; %set Cac pure bounds to zero
                modelJoint.ub(r) = 0; %set Cac pure bounds to zero
            end
            %Block Cac, Clj hybrid model flux
            for r = 2270:4550
                modelJoint.lb(r) = 0; %set fused model bounds to zero
                modelJoint.ub(r) = 0; %set fused model bounds to zero
            end

            %Place lower bound on Clj H2, CO2 uptake equal to the summed
            %lower bound of H2, CO2 export from all other organsism/cell
            %types in the simulation
            modelJoint.lb(5063) = -h2_store_var; %Cljpure H2 uptake
            modelJoint.lb(4708) = -h2_store_var; %Cljpure H2 uptake
            modelJoint.lb(4686) = -co2_store_var;
            modelJoint.lb(5020) = -co2_store_var;

            %Assign NGAM based on the fractional abundance the 
            %organism/cell type occupies in the co-culuture
            modelJoint.lb(1702) = 0.45*X(k,i); %Cljpure ATPM
            %Store ATPM for use in fermentation yield FVA
            CLJpure_ATPstore = modelJoint.lb(1702); %Cljpure ATPM
            
            %Non-hybrid Clj Flux balance analysis
            modelJoint = changeObjective (modelJoint, 'CljpureBIOMASS_Cl_DSM_WT_46p666M1');
            FBA = optimizeCbModel (modelJoint, 'max');
            
            %Store non-hybrid Clj glucose, fructore uptake, growth rate
            if isempty(FBA.obj)
                fruc_flux(k,i) = 0;
                gluc_flux(k,i) = 0;
                mu(k,i) = 0;
            else
                gluc_flux(k,i) = FBA.full(5046); %Community glucose uptake
                fruc_flux(k,i) = FBA.full(5041); %Community fructose uptake
                mu(k,i) = FBA.obj;
            end
        end
    end
    load model.mat    
    
    %Set lower bounds for community glucose and fructose uptake rates
    modelJoint.lb(5046) = sum(gluc_flux(k,:));
    modelJoint.lb(5041) = sum(fruc_flux(k,:));
    
    %Set lower bounds for organism/cell type glucose, fructose uptake rates
    modelJoint.lb(4639) = gluc_flux(k,1);
    modelJoint.lb(4557) = fruc_flux(k,1); 
    modelJoint.lb(4695) = fruc_flux(k,2);
    
    %Place lower bound on organism/cell type growth rate (99% of maximum
    %growth rate)
    modelJoint.lb(4668) = bm_eff*mu(k,1)-bm_eff*mu(k,1)*0.01;
    modelJoint.lb(1476) = bm_eff*mu(k,2)-bm_eff*mu(k,2)*0.01;

    %Block flux through hybrid states with zero relative abundance
    for r = 2270:4550
        modelJoint.lb(r) = 0;
        modelJoint.ub(r) = 0;
    end
    %Loop over all co-culture fermentationa products
    for rr = 1: length(rxnsList)
        %Define fermentation product export flux as FBA objective
        modelJoint = changeObjective (modelJoint, rxnsList{rr});
        
        %Place lower bounds on Clj H2 and CO2 uptake rate
        modelJoint.lb(4708) = -h2_store_var;
        modelJoint.lb(4686) = -co2_store_var; 
        %Store Clj CO2, H2 uptake rates
        co2_vals{k,1} = co2_store_var;
        h2_vals{k,1} = h2_store_var;
        
        %Place NGAM lower bounds for each organism/cell type
        modelJoint.lb(1100) = CACpure_ATPstore;
        modelJoint.lb(1702) = CLJpure_ATPstore;
        modelJoint.lb(3369) = CACfused_ATPuse;
        modelJoint.lb(3978) = CLJfused_ATPuse;
        
        %Maximize fermentation yield
        FBA = optimizeCbModel (modelJoint, 'max');
        %Store solution
        if ~isempty(FBA.obj)
            store_obj(k,rr) = FBA.obj;
            if k == 1
                dStore_obj(k,rr) = FBA.obj*dt;
                Cstore_obj(k,rr) = C0(rr) + dStore_obj(k,rr);
                Cstore_GL(k,rr) = Cstore_obj(k,rr) * MW(rr);
            else
                dStore_obj(k,rr) = FBA.obj*dt;
                Cstore_obj(k,rr) = Cstore_obj(k-1,rr) + dStore_obj(k,rr);
                Cstore_GL(k,rr) = Cstore_obj(k,rr) * MW(rr);
            end
                
        else
            store_obj(k,rr) = 0;
            dStore_obj(k,rr) = 0;
            Cstore_obj(k,rr) = 0;
            Cstore_GL(k,rr) = 0;
        end
    end
    %Euler Integration
    dCgluc = sum(gluc_flux(k,:))*dt;
    Cgluc(k+1,1) = Cgluc(k,1) + dCgluc;
    dCfruc = sum(fruc_flux(k,:))*dt;
    Cfruc(k+1,1) = Cfruc(k,1) + dCfruc;
    dX(k,:) = dt*(bm_eff*mu(k,:).*XX(k,:));
    XX(k+1,:) = max(XX(k,:)+dX(k,:),0);
    X(k+1,:) = XX(k+1,:)./(sum(XX(k+1,:)));
    XX(k+1,:) = max(XX(k,:) + dX(k,:),0);
    X(k+1,:) = XX(k+1,:)./(sum(XX(k+1,:)));
end


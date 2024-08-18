% to generalize to multiple species


%% Figure 4B
% Generate patterns and compute biomass of multispecies colonies

clear; 
num_species = 5;
initRatios = [1 0 0 0 0;
    0 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1;
    % 2 species
    1 1 0 0 0;
    1 0 1 0 0; 
    1 0 0 1 0; 
    1 0 0 0 1; 
    1 0 0 1 0; 
    % 3 species(more out there)
    1 1 1 0 0; 
    1 1 0 1 0;
    1 1 0 0 1;
    1 0 1 1 0;
    1 0 1 0 1;
    1 0 0 1 1;
    % 4 species
    1 1 1 1 0;
    1 1 1 0 1;
    1 1 0 1 1;
    1 0 1 1 1;
    0 1 1 1 1;
    % 5 species
    1 1 1 1 1
    ]; % initial ratio of WT:cheater:hyperswarmer


load 'parameters.mat'
Output_Biomass = zeros(size(initRatios, 1), num_species); % the final biomass of each species 

for iter = 22 %length(initRatios)
    
    fprintf('iter = %d\n', iter)
    initialRatio = initRatios(iter, :);
    initialFract = initialRatio / sum(initialRatio); % initial fraction of each species
    SetParameters
    BranchingColonyMultispecies_Core

    filename = ['results\Fig4B_' num2str(iter,'%02d')];    
    MakeFigure_patterns
    Output_Biomass(iter, :) = Biomass;
    
end

save('results\Fig4B_biomass.mat', 'initRatios', 'Output_Biomass')


%% Figure S9 (& Figure 5A)
% Competition between species on solid or in liquid
% compute final ratios when initiate with varying initial ratios

clear; 
ratios = 10.^(linspace(-3,3,61))';
initRatios = zeros(length(ratios), 3);
% prefix = 'CTvsWT'; ra = 2; rb = 1; % Cheater : WT
% prefix = 'CTvsHS'; ra = 2; rb = 3; % Cheater : hyperswarmer
prefix = 'HSvsWT'; ra = 3; rb = 1; % Hyperswarmer : WT
initRatios(:, ra) = ratios; initRatios(:, rb) = ones(length(ratios), 1);

load 'parameters.mat'
Output_Biomass_solid = zeros(size(initRatios, 1), 3); % the final biomass of each species
Output_Biomass_liquid = zeros(size(initRatios, 1), 3); 

for iter = 1 : length(initRatios)
    
    fprintf('iter = %d\n', iter)
    initialRatio = initRatios(iter, :);
    initialFract = initialRatio / sum(initialRatio); % initial fraction of each species
    SetParameters

    BranchingColonyMultispecies_Core
    Output_Biomass_solid(iter, :) = Biomass;
    
    LiquidCulture
    Output_Biomass_liquid(iter, :) = Biomass;

end

finalRatios_solid = Output_Biomass_solid(:, ra) ./ Output_Biomass_solid(:, rb);
finalRatios_liquid = Output_Biomass_liquid(:, ra) ./ Output_Biomass_liquid(:, rb);

save(['results\FigS9_' prefix '_biomass.mat'], 'initRatios', 'Output_Biomass_solid', 'finalRatios_solid', 'Output_Biomass_liquid', 'finalRatios_liquid')

%% Figure 5BC
% Show regions of cooperation on solid medium (Fig 5C)

clear; 
% initialRatio = [1 0 0]; prefix = 'WT'; % Wildtype alone
initialRatio = [1 1 0]; prefix = 'WT+CT'; % Wildtype + cheater
load 'parameters.mat'

initialFract = initialRatio / sum(initialRatio); % initial fraction of each species
SetParameters
saveintermediateimages = true;
BranchingColonyMultispecies_Core

% Period of cooperation during liquid growth (Fig 5B)
LiquidCulture
t = (0 : dt : totalt)';
cooperPeriod = aCsV(:,1) < max(aCs_act); % cooperation period
cooperPeriod(1) = 0;
Results = [t, C, N, double(cooperPeriod)];
save(['results\Fig5B_' prefix '_liquidgrowth.mat'], 'C', 'N', 't', 'cooperPeriod', 'Results')

%% Figure 6A
% Evolution of population structure over multiple cycles

clear;
ndays = 7;
N0V = [0.4 0.8 1.2 1.5]; 

for iN = 1 : length(N0V)

    fprintf('N0 = %f\n', N0V(iN))
    load 'parameters.mat'
    Parameters(3)  = Parameters(3) * N0V(iN);  % change N0
    Parameters(10) = Parameters(10) / N0V(iN); % keep the same N_upper
    Parameters(11) = Parameters(11) / N0V(iN); % keep the same N_lower
    
    initialRatio = [1 1e-7 1e-7];
    initialFract = initialRatio / sum(initialRatio); % initial fraction of each species
    SetParameters
    BiomassDynamics = []; % record of the biomass of each species over multiple days
    
    for iday = 1 : ndays
        BranchingColonyMultispecies_Core
        BiomassDynamics = [BiomassDynamics; BiomassV(1 : end - 1,:)]; 
        initialRatio = BiomassV(end, :);
        initialFract = initialRatio / sum(initialRatio);
        fprintf('Day %d\n', iday)
    end
    BiomassDynamics = [BiomassDynamics; BiomassV(end, :)];
    
    filename = ['results\Fig6A_N0=' num2str(N0V(iN)) '_' num2str(ndays) 'days'];
    save([filename '.mat'], 'BiomassDynamics','initialRatio','ndays')
    MakeFigure_evolution

end



num_species = 5;  

% assign cell names
speciesName = cell(1, num_species);  
for i = 1:num_species
    speciesName{i} = num2str(i);  
end

% other vectors will follow the same order

%% Tunable parameters
bN = Parameters(1);  % nutrient consumption rate
DN = Parameters(2);  % nutrient diffusivity
N0 = Parameters(3);  % initial nutrient conc.

% aCs_act = [1, Parameters(7), Parameters(12)] * Parameters(4);  % cell growth rate of each species
% gs      = [1, 1, Parameters(8)] * Parameters(5);   % autonomous motility
% hs_act  = [1, 0, Parameters(13)] * Parameters(6);  % cooperative motility coefficients
% rs      = [1, Parameters(14), Parameters(15)];     % response coefficients

aCs_act = [1, Parameters(7), Parameters(12), 1, 1.2] * Parameters(4);  % cell growth rate of each species
gs      = [1, 1, Parameters(8), 0.9, 1.2] * Parameters(5);   % autonomous motility
hs_act  = [1, 0, Parameters(13), 0.5, 1] * Parameters(6);  % cooperative motility coefficients
rs      = [1, Parameters(14), Parameters(15), 1, 0.9];     % response coefficients


N_upper = Parameters(10) * N0;           % upper bound of nutrient for cooperative motility
N_lower = N_upper - Parameters(11) * N0; % lower bound of nutrient for cooperative motility

%% Fixed parameters
L      = 90;    % domain size
totalt = 16;    % total time
dt     = 0.02;  % time step
nx     = 1001; ny = nx; % number of nodes

r0  = 5;   % initial colony radius
C0  = 8;   % initial cell density

% branch density & width of single-species colonies
Densities = [0.14, 0.14, 0.25, 0.2, 0.2];
Widths    = [5, 5, 20, 10, 15];

%% Others
noiseamp = 0;  % noise amplitude of branch direction
dt_updatebranch = 10 * dt;  % time step for updating branch locations
saveintermediateimages = false;
plotgap = 4;  % time interval for making plots while computing (same unit as totalt)
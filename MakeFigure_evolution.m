% modified for multiple species

% making figures for Fig 6A
figure(2); clf
set(gcf, 'position', [320.3333  345.6667  185.3333   82.6667])

Frat = BiomassDynamics ./ sum(BiomassDynamics, 2); 
area(linspace(0, ndays, (nt + 1) * ndays + 1), Frat)

% 12 color map, take the first num_species to label each strain
color_maps = [[66, 206, 227],
    [31, 120, 180],
    [178, 223, 138],
    [51, 160, 44],
    [251, 154, 153],
    [227, 26, 28],
    [253, 191, 111],
    [255, 127, 0],
    [202, 178, 214],
    [106, 61, 154],
    [255, 255, 153],
    [177, 89, 40]]

colororder = color_maps(1:num_species)/255;


ylim([0 1]); xlim([0 ndays])
xlabel 'Days'
ylabel 'Fraction'
set(gca, 'TickLength', [0,0])
set(gca, 'FontSize', 8)
set(gca, 'LabelFontSizeMultiplier', 10/8)
set(gca, 'xtick', 0 : ndays)
set(gca, 'ytick', [0 1])
set(gca, 'TitleFontWeight', 'bold')
set(gca, 'XGrid', 'on')
set(gca, 'Layer', 'top')
set(gca, 'GridAlpha', 0.5)

fig = gcf; fig.PaperPositionMode = 'auto';
print([filename '.svg'],'-dsvg')
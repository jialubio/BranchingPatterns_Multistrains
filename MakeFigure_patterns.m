% making figures for Fig 4B

multipliers = ones(num_species, 1);
multipliers(1) = 5;

Ctotal = 0; % Initialize total concentration
for i = 1:num_species
    Ctotal = Ctotal + multipliers(i) * C{i}; % Sum up the weighted concentrations of all species
end

all_p = cell(num_species, 1); % Initialize a cell array to hold the proportions for each species
for i = 1:num_species
    all_p{i} = multipliers(i) * C{i} ./ Ctotal; % Calculate the proportion for each species
    all_p{i}(isnan(all_p{i})) = 0; % Replace NaN values with 0
end
p_list = all_p(1: num_species);

ind = 1 : 2 : nx;


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
    [177, 89, 40]];
color_list = color_maps(1:num_species, :);
ColorMap = MarkMixing_3color(color_list, p_list);

figure(2); clf; set(gcf,'position',[360   227   391   391])
hold off; surf(xx(ind, ind), yy(ind, ind), ones(size(xx(ind, ind))), ColorMap(ind, ind, :))
view([0, 0, 1]); shading interp; axis equal; box on
axis([-L/2 L/2 -L/2 L/2]);
set(gca,'YTick',[], 'XTick',[], 'Visible', 'off')

print([filename '.tif'],'-dtiff','-r800')
figure(1)

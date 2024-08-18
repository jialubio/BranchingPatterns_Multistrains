
figure(1);
set(gcf,'position',[213  73  966  538]);
iz = 1 : 2 : nx; % plot less cols & rows to save time

% Plot each species
subplot(2, num_species, j)
    hold off; pcolor(xx(iz, iz), yy(iz, iz), C{j}(iz, iz));
    shading interp; axis equal; caxis([-max(C{j}(:)) max(C{j}(:))])
    axis([-L/2 L/2 -L/2 L/2]); hold on
    set(gca,'YTick',[], 'XTick',[])
    plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)
    title(speciesName{j})
    drawnow

% Plot all species
if j == find(initialRatio,1,'last')

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

Ctotal = 0;
for i = 1:num_species
    Ctotal = Ctotal + C{i};  % Add weights here if needed, e.g., weights(i)*C{i}
end

p_list = cell(num_species, 1);

% Calculate the proportion for each species
for i = 1:num_species
    p_list{i} = C{i} ./ Ctotal;
    p_list{i}(isnan(p_list{i})) = 0;  % Replace NaN with 0
end


subplot(2, num_species, num_species+1) % total cell density
    hold off; pcolor(xx(iz, iz), yy(iz, iz), Ctotal(iz, iz));
    shading interp; axis equal; caxis([-max(Ctotal(:)) max(Ctotal(:))])
    axis([-L/2 L/2 -L/2 L/2]); colormap('gray'); hold on
    set(gca,'YTick',[], 'XTick',[])
    plot(Tipx{j}(:,ib), Tipy{j}(:,ib), '.', 'markersize', 5)
    title('Total cell density')

% subplot(2, num_species, num_species+2) % show each species by color
%     ColorMap = MarkMixing_N_color(color_list, p_list);
%     hold off; surf(xx(iz, iz), yy(iz, iz), ones(size(xx(iz, iz))), ColorMap(iz, iz, :))
%     view([0, 0, 1]); shading interp; axis equal; box on
%     axis([-L/2 L/2 -L/2 L/2]);
%     set(gca,'YTick',[], 'XTick',[])
%     title('All')

subplot(2, num_species, num_species+3) % line graph of cell densities
    yyaxis left; hold off
    mid = (nx + 1) / 2;
    plot(x(mid:end), C{1}(mid:end,mid), '-', 'color', color_list(1, :)/255, 'linewidth', 2, 'DisplayName', speciesName{1}); hold on
    plot(x(mid:end), C{2}(mid:end,mid), '-', 'color', color_list(2, :)/255, 'linewidth', 2, 'DisplayName', speciesName{2});
    plot(x(mid:end), C{3}(mid:end,mid), '-', 'color', color_list(3, :)/255, 'linewidth', 2, 'DisplayName', speciesName{3});
    plot(x(mid:end), C{4}(mid:end,mid), '-', 'color', color_list(4, :)/255, 'linewidth', 2, 'DisplayName', speciesName{4});
    plot(x(mid:end), C{5}(mid:end,mid), '-', 'color', color_list(5, :)/255, 'linewidth', 2, 'DisplayName', speciesName{5});
    plot(x(mid:end), Ctotal(mid:end,mid), 'k-', 'linewidth', 2, 'DisplayName', 'Total C')
    ylabel 'Cell densities';
    yyaxis right; hold off
    plot(x(mid:end), N(mid:end,mid), '-', 'color', [0.7,0.7,0.7], 'linewidth', 2, 'DisplayName', 'nutrient'); 
    ylim([0 N0])
    ylabel 'Nutrient'
    xlabel 'Distance from center'
    legend

drawnow

end
    
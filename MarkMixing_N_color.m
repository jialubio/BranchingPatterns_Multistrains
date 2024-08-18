function colors = MarkMixing_N_color(color_list, p_list)
% Generating a smooth color map using multiple colors
% Based on Mark's method for generating a color gradient
% Inputs:
%   color_list - a cell array of colors, where each color is an RGB vector
%   p_list - a list of proportions for the colors

gamma = .43;
num_species = size(color_list, 1);
lin_colors = zeros(num_species, 3);
brightness = zeros(num_species, 1);

% Convert each color from sRGB and calculate its brightness
for i = 1:num_species
    lin_colors(i, :) = from_sRGB(color_list(i, :));
    brightness(i) = sum(lin_colors(i, :))^gamma;
end

% Initialize the output color matrix
colors = zeros(size(p_list{1}, 1), size(p_list{1}, 2), num_species);
intensity = zeros(size(p_list{1}));

% Calculate overall intensity based on brightness and proportions
for i = 1:num_species
    intensity = intensity + brightness(i) * p_list{i};
end
intensity = intensity .^ (1/gamma);

% Calculate color contributions
for dim = 1:3 % For each color dimension (R, G, B)
    color_dim = zeros(size(p_list{1}));
    for i = 1:num_species
        color_dim = color_dim + lin_colors(i, dim) * p_list{i};
    end
    colors(:, :, dim) = color_dim;
end

% Adjust each color channel based on the overall intensity
sumcolor = sum(colors, 3);
sumcolor(sumcolor == 0) = 1; % Avoid division by 0
for dim = 1:3
    colors(:, :, dim) = colors(:, :, dim) .* intensity ./ sumcolor;
    colors(:, :, dim) = to_sRGB_f(colors(:, :, dim)); % Convert back to sRGB
end

end

function f = to_sRGB_f(x)
%     ''' Returns a sRGB value in the range [0,1]
%         for linear input in [0,1].
f = 12.92*x;
f(x > 0.0031308) = (1.055 * (x(x > 0.0031308) .^ (1/2.4))) - 0.055;
end

function f = to_sRGB(x)
%     ''' Returns a sRGB value in the range [0,255]
%         for linear input in [0,1]
f = round(255.9999 * to_sRGB_f(x));
end


function y = from_sRGB(x)
%     ''' Returns a linear value in the range [0,1]
%         for sRGB input in [0,255].
x = x / 255.0;
y = x / 12.92;
y(x > 0.04045) = ((x(x > 0.04045) + 0.055) / 1.055) .^ 2.4;
end

function f = lerp(colors, proportions)
f = 0; % Initialize the interpolated color
for i = 1:length(colors)
    f = f + colors{i} * proportions{i};
end

end 
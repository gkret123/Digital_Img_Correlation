% ME360 Engineering Experimentation
% Team A2 - Digital Image Correlation

% Large Box Beam

clear all;

% filepath for 1000N:
file_path = "C:\Users\maria\Downloads\DIC_LargeBox_1000N.xlsx";

% filepath for 4000N:
%file_path = "C:\Users\maria\Downloads\DIC_LargeBox_4000N.xlsx";

opts = detectImportOptions(file_path);
opts.VariableNamingRule = 'preserve';
opts.Sheet = 1; % Adjust if needed
data = readtable(file_path, opts);

% Ensure column names are valid MATLAB identifiers
data.Properties.VariableNames = matlab.lang.makeValidName(data.Properties.VariableNames);

% Define coordinates to remove, note: these values change for each beam
x_val = 35.4785;
y_val = 17.7537;
z_val = 483.516;

% Check if required columns exist to avoid errors
if all(ismember({'X', 'Y', 'Z'}, data.Properties.VariableNames))
    % Filter out exact match row
    rowsToRemove = (data.X == x_val) & (data.Y == y_val) & (data.Z == z_val);
    data(rowsToRemove, :) = [];

    % Display number of remaining rows
    disp(['Remaining rows: ', num2str(height(data))]);
else
    disp("Error: Required columns are missing from the dataset.");
end

% Scatter plot X vs Y
figure;
scatter(data.X, data.Y, 1, 'filled'); % Marker size set to 1
xlabel('X');
ylabel('Y');
title('Scatter Plot of X vs Y');
grid on;

% Display number of rows in the dataset
disp(['Number of rows: ', num2str(height(data))]);

% Heat Map of exx Strain Field
figure;
scatter(data.X, data.Y, 10, data.exx, 'filled'); % Marker size set to 10
colormap('parula'); % MATLAB does not have 'inferno', using 'parula' as an alternative
colorbar;
xlabel('X');
ylabel('Y');
title('Heat Map of Axial Strain Field');
grid on;

% Heat Map of eyy Strain Field
figure;
scatter(data.X, data.Y, 10, data.eyy, 'filled'); % Marker size set to 10
colormap('hot'); % MATLAB does not natively support 'plasma', using 'hot' as an alternative
colorbar;
xlabel('X');
ylabel('Y');
title('Heat Map of Lateral Strain Field');
grid on;

% Heat Map of exy Strain Field
figure;
scatter(data.X, data.Y, 10, data.exy, 'filled'); % Marker size set to 10
colormap('jet'); % MATLAB does not natively support 'viridis', using 'jet' as an alternative
colorbar;
xlabel('X');
ylabel('Y');
title('Heat Map of Shear Strain Field');
grid on;

% Clean column names
data.Properties.VariableNames = matlab.lang.makeValidName(data.Properties.VariableNames);

% Create a single figure with three vertically stacked subplots
figure;

% Interpolated Shear Strain (ε_xy)
subplot(3,1,3);
data_filtered = data(~isnan(data.exy) & data.exy ~= 0, :);
xi = linspace(min(data_filtered.X), max(data_filtered.X), 300);
yi = linspace(min(data_filtered.Y), max(data_filtered.Y), 300);
[X_grid, Y_grid] = meshgrid(xi, yi);
Z = griddata(data_filtered.X, data_filtered.Y, data_filtered.exy, X_grid, Y_grid, 'linear');
contourf(X_grid, Y_grid, Z, 100, 'LineStyle', 'none');
colormap('hsv'); 
colorbar;
xlabel('x (mm)');
ylabel('y (mm)');
title('Interpolated Shear Strain \epsilon_{xy}');
axis equal;
grid on;

% Interpolated Lateral Strain (ε_yy)
subplot(3,1,2);
data_filtered = data(~isnan(data.eyy) & data.eyy ~= 0, :);
xi = linspace(min(data_filtered.X), max(data_filtered.X), 300);
yi = linspace(min(data_filtered.Y), max(data_filtered.Y), 300);
[X_grid, Y_grid] = meshgrid(xi, yi);
Z = griddata(data_filtered.X, data_filtered.Y, data_filtered.eyy, X_grid, Y_grid, 'linear');
contourf(X_grid, Y_grid, Z, 100, 'LineStyle', 'none');
colormap('hsv'); 
colorbar;
xlabel('x (mm)');
ylabel('y (mm)');
title('Interpolated Lateral Strain \epsilon_{yy}');
axis equal;
grid on;

% Interpolated Axial Strain (ε_xx)
subplot(3,1,1);
data_filtered = data(~isnan(data.exx) & data.exx ~= 0, :);
xi = linspace(min(data_filtered.X), max(data_filtered.X), 300);
yi = linspace(min(data_filtered.Y), max(data_filtered.Y), 300);
[X_grid, Y_grid] = meshgrid(xi, yi);
Z = griddata(data_filtered.X, data_filtered.Y, data_filtered.exx, X_grid, Y_grid, 'linear');
contourf(X_grid, Y_grid, Z, 100, 'LineStyle', 'none');
colormap('hsv'); 
colorbar;
xlabel('x (mm)');
ylabel('y (mm)');
title('Interpolated Axial Strain \epsilon_{xx}');
axis equal;
grid on;
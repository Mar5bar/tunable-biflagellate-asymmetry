addpath(genpath('../common'))
% Loop through all the .csv files in the ./data directory. For each file, load in the data, and compute the curvature.

% Get a list of all the .csv files in the ./data directory whose name ends with 'filtered'.
files = dir('./data/*filtered.csv');

multipliers = [];
curvatures = [];

% Loop through the files, parsing the name by splitting at _ and removing the 'Trans' prefix.
figure;
for i = 1:length(files)
    % Load the data from the file.
    [x,y] = parse_datafile(['./data/' files(i).name]);
    % Compute the curvature.
    curvature = curvature_from_traj(x, y);
    % Parse the filename for the frequency multiplier.
    name_parts = split(files(i).name, '_');
    freq_multiplier = 100 / str2double(name_parts{1}(6:end));
    % Append the multiplier and curvature to the lists.
    multipliers = [multipliers freq_multiplier];
    curvatures = [curvatures curvature];
    % Plot the path with the curvature as the title.
    nexttile();
    plot(x, y);
    axis equal
    xlim([0,300])
    ylim([0,300])
    title(['Mul: ' num2str(freq_multiplier)]);
end

% Fit a curve of the form y = a*(x-1)/(x+1) to the data.
% Define the function to fit.
f = fittype('a*(x-1)/(x+1)', 'independent', 'x', 'dependent', 'y');
% Fit the curve.
fitresult = fit(multipliers', curvatures', f);
% Plot the curve.
figure
plot(fitresult, multipliers, curvatures);
xlabel('Frequency Multiplier');
ylabel('Curvature');
grid on;
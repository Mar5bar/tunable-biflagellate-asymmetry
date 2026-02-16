addpath(genpath('../common'))
clear all
close all
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
    if (freq_multiplier == 1 && curvature > 0.01)
        continue
    end
    multipliers = [multipliers freq_multiplier];
    curvatures = [curvatures curvature];
    % Plot the path with the curvature as the title.
    nexttile();
    hold on
    plot(x, y);
    par = circleFit(x,y);
    t = linspace(0,2*pi,1e4);
    plot(par(3)*cos(t)+par(1), par(3)*sin(t)+par(2))
    axis equal
    xlim([0,300])
    ylim([0,300])
    title(['Mul: ' num2str(freq_multiplier)]);
    set(gca,'FontSize',12)
end

% Fit a curve of the form y = a*(x-1)/(x+1) to the data.
% Define the function to fit.
f = fittype('a*(x-1)/(x+1)', 'independent', 'x', 'dependent', 'y');
% Fit the curve.
fitresult = fit(multipliers', curvatures', f);
fun = @(x) fitresult.a * (x-1)./(x+1);
% Plot the curve.
figure
font_size = 14;
set(gcf, 'Position',  [611,895,500,400])
set(gcf,'color','w');
hold on
box on
xs = linspace(1,max(multipliers),1e2);
plot(multipliers, curvatures,'+','Color','black','LineWidth',2)
plot(xs, fun(xs), 'k:','LineWidth',1);

% interval = confint(fitresult);
% fun = @(x) interval(1) * (x-1)./(x+1);
% plot(xs, fun(xs), 'k:','LineWidth',1);
% fun = @(x) interval(2) * (x-1)./(x+1);
% plot(xs, fun(xs), 'k:','LineWidth',1);

xlabel('$k_L / k_R$','FontSize',font_size,'Interpreter','latex');
ylabel('$\kappa$','FontSize',font_size,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');   
set(gca,'FontSize',font_size);
grid on;
legend({'Data','Scaling law'},'location','southeast')

exportgraphics(gcf, 'fit_to_robot_data.png', 'Resolution', 400)

% Compute the mean curvature at a given multiplier:

uniq_mults = unique(multipliers);
avgs = [];
for i = 1 : length(uniq_mults)
    mult = uniq_mults(i);
    total = 0;
    count = 0;
    for j = 1 : length(multipliers)
        if (mult == multipliers(j))
            total = total + curvatures(j);
            count = count + 1;
        end
    end
    avgs(i) = total / count;
end
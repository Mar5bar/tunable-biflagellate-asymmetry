addpath(genpath('../common'))
clear all
close all
% Load the given csv, which contains curvature values.

data = readmatrix('./data/curvature_values_0731-1.csv');
multipliers = data(:,3);
curvatures = data(:,2); % Data is given in units of 1/m.

% Delete the third and fourth entries as they are erroneous.
% multipliers(3:4) = [];
% curvatures(3:4) = [];

init_guess = curvatures(end)*(multipliers(end)+1)/(multipliers(end)-1);

% Fit a curve of the form y = a*(x-1)/(x+1) to the data.
% Define the function to fit.
f = fittype('a*(x-1)/(x+1)', 'independent', 'x', 'dependent', 'y');
% Fit the curve.
fitresult = fit(multipliers, curvatures, f, 'StartPoint', init_guess);
fun = @(x) fitresult.a * (x-1)./(x+1);
% Plot the curve.
figure
font_size = 14;
set(gcf, 'Position',  [611,895,500,400])
set(gcf,'color','w');
hold on
box on
xs = linspace(1,max(multipliers),1e2);
% Compute the mean and range for each unique multiplier
uniq_mults = unique(multipliers);
mean_curvatures = zeros(size(uniq_mults));
min_curvatures = zeros(size(uniq_mults));
max_curvatures = zeros(size(uniq_mults));

for i = 1:length(uniq_mults)
    mult = uniq_mults(i);
    relevant_curvatures = curvatures(multipliers == mult);
    mean_curvatures(i) = mean(relevant_curvatures);
    min_curvatures(i) = min(relevant_curvatures);
    max_curvatures(i) = max(relevant_curvatures);
end

% Plot mean curvatures with error bars
errorbar(uniq_mults, mean_curvatures, mean_curvatures - min_curvatures, max_curvatures - mean_curvatures, 'o', 'Color', 'black', 'LineWidth', 1);
plot(xs, fun(xs), 'k','LineWidth',1);

% interval = confint(fitresult);
% fun = @(x) interval(1) * (x-1)./(x+1);
% plot(xs, fun(xs), 'k:','LineWidth',1);
% fun = @(x) interval(2) * (x-1)./(x+1);
% plot(xs, fun(xs), 'k:','LineWidth',1);

xlabel('$k_L / k_R$','FontSize',font_size,'Interpreter','latex');
ylabel('$\kappa$ (1/m)','FontSize',font_size,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');   
set(gca,'FontSize',font_size);
grid on;
legend({'Data','Scaling law'},'location','southeast')

exportgraphics(gcf, 'fit_to_robot_data.png', 'Resolution', 400)
exportgraphics(gcf, 'fit_to_robot_data.pdf', 'ContentType','vector')

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
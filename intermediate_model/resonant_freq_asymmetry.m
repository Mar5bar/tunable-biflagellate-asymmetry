close all
clear all

figure('Position',[812 334 556 470])
hold on
fontsize = 24;

ts = linspace(0,30*pi,1e5);
gamma = 1;
delta = 0.05;
kL = 1;

kRs = linspace(1.1,3,1e3);
vals = 0*kRs;
approxVals = 0*kRs;

K = @(t,k) -cos(k*t);

for i = 1 : length(kRs)
	kR = kRs(i);
	vals(i) = trapz(ts, cos(gamma*(kR - kL)*ts + delta*(K(ts,kR) - K(ts,kL))));
	approxVals(i) = trapz(ts, cos(gamma*(kR - kL)*ts));
end

plot(kRs, vals, 'LineWidth', 2, 'Color', 'black')
plot(kRs, approxVals, 'LineWidth', 1, 'Color', 0.5*[1,1,1], 'LineStyle','-')
legend({'$\delta \ll \gamma$', '$\delta = 0$'}, 'Interpreter', 'latex')
xlim([min(kRs),max(kRs)])
ylim([-9,7])
axis square
box on
xlabel('$k_R/k_L$')
ylabel('$\int\cos{\theta(t)} \, \mathrm{d}t$')
set(gca,'FontSize', fontsize)

inset = axes('Position',[0.471223021582734 0.22 0.366906474820144 0.238297872340426]);
plot(kRs, vals - approxVals, 'Color', 'black')
set(inset, 'XTick', [2])
set(inset, 'YTick', [])
xlim([min(kRs),max(kRs)])
box on
ylabel('Error')
yline(0,'Color','black')
set(gca,'FontSize', fontsize)

exportgraphics(gcf, 'resonant_freq_asymmetry.pdf')
clearvars, clc %, close all

Sys0 = struct('g', [1.985, 1.98, 1.97], 'lwpp', 2, 'weight', 1);
Sys1 = struct('g', 2.02, 'lwpp', 2, 'weight', 1);

Exp = struct('Field', 340, 'mwRange', [9.3, 9.8], 'nPoints', 10000, 'Harmonic', 1);
B0 = linspace(Exp.mwRange(1), Exp.mwRange(2), Exp.nPoints);
% gs = linspace(0.75, 1.75, 5);

clf
for ii = 1
    % Sys1.g = gs(ii);
    ysim = pepper({Sys0, Sys1}, Exp);
    ysim = ysim/max(ysim);
    ydi = cumtrapz(ysim);
    ydi = ydi/max(ydi);
    plot(B0, ysim)
    hold on
    plot(B0, ydi)
end

% ylim([-1, 1])
clearvars, clc %, close all

Sys0 = struct('g', 2.05, 'lwpp', 2, 'weight', 1);
Sys1 = struct('g', 2.00, 'lwpp', 2, 'weight', 1);

Exp = struct('Range', [320, 1000], 'mwFreq', 9.4, 'nPoints', 10000);
B0 = linspace(Exp.Range(1), Exp.Range(2), Exp.nPoints);
gs = linspace(0.75, 1.75, 5);

clf
for ii = 1:5
    Sys1.g = gs(ii);
    ysim = pepper({Sys0, Sys1}, Exp);
    ysim = ysim/max(ysim);
    ydi = cumtrapz(cumtrapz(ysim));
    ydi = ydi/max(ydi);
    plot(B0, ysim)
    hold on
    % plot(B0, ydi)
end
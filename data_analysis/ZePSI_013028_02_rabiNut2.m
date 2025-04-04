clearvars
addpath(genpath('util/'))

%% IMPORT

load('../data/processed/ZePSI-E-013028-Nut2.mat')
nMeas = numel(x);

%% FIT cosine
expcos = @(xx, p) [ones(length(xx), 1), exp(-xx/p(1)).*cos(2*pi*p(2)*xx)];

fitOpt = optimoptions('lsqnonlin','Display','off');

figure()
clf
w0 = 47*1e-3;
clear("pci")
for ii = 1:nMeas
    p0 = [40, w0];
    fitmodel = @(p) expcos(x{ii}', p);
    ydata = real(y{ii}');
    ydata = ydata/max(ydata);

    [yfit{ii}, pfit{ii}, pci(ii, :)] = ...
        lsqnonlin2steps(ydata, fitmodel, p0, fitOpt);

    nexttile
    plot(x{ii}, ydata, 'o')
    hold on
    plot(x{ii}, yfit{ii})
    title(sprintf('%d: %.4f MHz, %.2f ns', ...
          ii, pfit{ii}(2)*1e3, 1/2/pfit{ii}(2)))
end

nexttile
for ii = 1:nMeas
    freq(ii) = pfit{ii}(2);
end
% errorbar(1:nMeas, freq, pci(:, 2, 1) - freq', pci(:, 2, 2) - freq', 'o')
errorbar(1:nMeas, freq, pci(:, 2), 'o')
xlim([0.9, nMeas + 0.1])
fprintf("Average Rabi frequency: %.2f\n", sum(freq*1e3./pci(:, 2)')/sum(1./pci(:, 2)))

clearvars
addpath(genpath('util/'))

%% IMPORT

load('../data/processed/ZePSI-E-013012-Nut2.mat')

%% RABI NUTATIONS EXPERIMENT PULSE 2

iMax = 196;
integWidth = 90;
% integWidth2 = 80;

nTau = size(y{1}, 1);
nMeas = numel(y);
figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    [integWindow{ii}, integWindowPlot{ii}] = deal(zeros(size(y{ii})));
    for itau = 1:nTau
        valMax = max(real(y{ii}(itau, :)));
        iInteg = (iMax - integWidth/2):(iMax + integWidth/2);
        integWindow{ii}(itau, iInteg) = ones(1, integWidth + 1);
    end
    y2{ii} = integWindow{ii}.*y{ii};  % Signal in the window

    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x{1}, x{2}, real(y{ii}));
    hold on
    plot(sax, x{1}, x{2}, imag(y{ii}));
    plot(sax, x{1}, x{2}, real(integWindow{ii})*0.2e4);
    ylim([-0.5e4, 0.5e4])
    plotTitle = strsplit(Param{ii}.TITL, '-');
    title(plotTitle{end})
end

% Integrate to get Rabi nutations
xrabi = x{2}*2 + 2;  % Each step is an increment of 2 ns of the pulse length
rabii = zeros(nMeas, nTau);
% Create window function
winham = windowhamming(nTau, 1);
for ii = 1:nMeas
    rabii(ii, :) = sum(y2{ii}(1:nTau, :), 2);
    % rabii2(ii, :) = sum(y22{ii}(1:nTau, :), 2);

    nexttile
    plot(xrabi, real(rabii(ii, :)), 'o-')
    hold on
    plot(xrabi, imag(rabii(ii, :)), 'o-')
    plot(xrabi, winham*max(real(rabii(ii, :))))
    plot(xrabi, winham.*real(rabii(ii, :)))
end

%% FIT cosine

rabii = rabii(:, 1:nTau);  % Recover rabii before zero filling
expcos = @(xx, A, tau, w) A*exp(-xx/tau).*cos(2*pi*w*xx);
fitmodel = @(p) expcos(xrabi', 1, p(1), p(2));  % A = 1

fitOpt = optimoptions('lsqnonlin','Display','off');

clf
w0 = 47*1e-3;
for ii = 1:nMeas
    p0 = [40, w0];
    ydata = real(rabii(ii, :))/max(real(rabii(ii, :)));

    [pfit{ii}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
        @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
    pci(ii, :, :) = nlparci(pfit{ii}, residual, 'jacobian', jacobian);

    [yfit{ii}, pfit{ii}(3), pfit{ii}(4)] = mldividefun(...
        fitmodel, ydata, pfit{ii});
    nexttile
    plot(xrabi, ydata, 'o')
    hold on
    plot(xrabi, yfit{ii})
    title(sprintf('%d: %.4f MHz, %.2f ns', ...
          ii, pfit{ii}(2)*1e3, 1/2/pfit{ii}(2)))
end

nexttile
for ii = 1:nMeas
    freq(ii) = pfit{ii}(2);
end
% errorbar(1:nMeas, freq, pci(:, 2, 1) - freq', pci(:, 2, 2) - freq', 'o')
errorbar(1:nMeas, freq, pci(:, 2, 1) - freq', ...
    pci(:, 2, 2) - freq', 'o')
xlim([0.9, nMeas + 0.1])
fprintf("Average Rabi frequency: %.2f\n", mean(freq*1e3))

%% FFT
%{
APPLY_WINDOW = 1;
tStep = xrabi(2) - xrabi(1);
fSampl = 1/tStep;
nzf = 1024;  % Zero filling
if APPLY_WINDOW
    winrabii = winham.*rabii;  % Apply window function
else
    winrabii = rabii;
end
if nzf ~= 0
    fxrabi = fSampl/nzf*(-nzf/2:nzf/2 - 1);  % Freq axis
    frabii = zeros(nMeas, nzf);  % Initialize fft arrays
    winrabii(:, nzf) = zeros(nMeas, 1);  % Zero filling
else
    fxrabi = fSampl/nTau*(-nTau/2:nTau/2 - 1);  % Freq axis
    frabii = zeros(nMeas, nTau);  % Initialize fft arrays
end

for ii = 1:nMeas
    frabii(ii, :) = fft(winrabii(ii, :));

    nexttile
    plot(fxrabi, abs(fftshift(frabii(ii, :))), 'o-')
    xlim([-0.1, 0.1])
end
%}

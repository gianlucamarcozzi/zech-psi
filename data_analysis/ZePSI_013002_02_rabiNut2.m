clearvars
addpath(genpath('util/'))

%% IMPORT

load('../data/processed/ZePSI-E-013002-Nut2.mat')

%% RABI NUTATIONS EXPERIMENT PULSE 1 (ONLY IN-PHASE-CHANNEL FOR NOW)

iMax = 214;
integWidth = 120;
% integWidth2 = 80;

nTau = size(y{1}, 1);
nMeas = numel(y);
figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    [integWindow{ii}, integWindowPlot{ii}] = deal(zeros(size(y{ii})));
    % [integWindow2{ii}, integWindowPlot2{ii}] = deal(zeros(size(y{ii})));
    for itau = 1:nTau
        valMax = max(real(y{ii}(itau, :)));
        iInteg = (iMax - integWidth/2):(iMax + integWidth/2);
        % iInteg2 = (iMax - integWidth2/2):(iMax + integWidth2/2);
        integWindow{ii}(itau, iInteg) = ones(1, integWidth + 1);
        % integWindow2{ii}(itau, iInteg2) = ones(1, integWidth2 + 1);
    end
    y2{ii} = integWindow{ii}.*y{ii};  % Signal in the window
    % y22{ii} = real(integWindow2{ii}.*y{ii});  % Signal in the window

    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x{1}, x{2}, real(y{ii}));
    hold on
    plot(sax, x{1}, x{2}, imag(y{ii}));
    plot(sax, x{1}, x{2}, real(integWindow{ii})*0.9e4);
    plot(sax, x{1}, x{2}, real(y2{ii}));
    % ylim(setaxlim([min(real(y{ii}), [], 'all'), max(real(y{ii}), ...
        % [], 'all')], 1.1))
    ylim([-1e4, 1e4])
    % xline([-1, 1]/2*integWidth + x{ii}{1}(iMax), '--')
    % fprintf("%d max at:\t%d\n", [ii, iMax])
end

% Integrate to get Rabi nutations
rabii = zeros(nMeas, nTau);
% rabii2 = zeros(nMeas, nTau);
figure(22)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
xrabi = x{2}*2 + 2;  % Each step is an increment of 2 ns of the pulse length
for ii = 1:nMeas
    rabii(ii, :) = sum(y2{ii}(1:nTau, :), 2);
    % rabii2(ii, :) = sum(y22{ii}(1:nTau, :), 2);

    nexttile
    % plot(xrabi, rabii(ii, :) - mean(rabii(ii, :)), 'o-')
    plot(xrabi, real(rabii(ii, :)), 'o-')
    hold on

    % plot(xrabi, (rabii2(ii, :) - mean(rabii2(ii, :)))*0.5, 'o-')
    plot(xrabi, imag(rabii(ii, :)), 'o-')
end

% Apply window function
winham = windowhamming(nTau, 3/4*nTau);

% FFT
tStep = xrabi(2) - xrabi(1);
fSampl = 1/tStep;
nzf = 0;  % Zero filling
if nzf ~= 0
    fxrabi = fSampl/nzf*(-nzf/2:nzf/2 - 1);
else
    fxrabi = fSampl/nTau*(-nTau/2:nTau/2 - 1);
end

figure(4)
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    if nzf ~= 0
        rabii(ii, nzf) = 0;  % Zero filling
        % rabii2(ii, nzf) = 0;  % Zero filling
    end
    % frabii = zeros()
    % frabii(ii, :) = fft(rabii(ii, :) - mean(rabii(ii, :)));
    frabii(ii, :) = fft(rabii(ii, :));
    % frabii2(ii, :) = fft(rabii2(ii, :) - mean(rabii(ii, :)));

    [~, ifMax(ii)] = max(abs(fftshift(frabii(ii, :))));
    % [~, ifMax2(ii)] = max(abs(fftshift(frabii2(ii, :))));

    nexttile
    plot(fxrabi, abs(fftshift(frabii(ii, :))), 'o-')
    % yyaxis right
    % plot(fxrabi, abs(fftshift(frabii2(ii, :))), 'o-')
    xlim([-0.1, 0.1])
end

% for ii = 1:nMeas
%     fprintf("%d pi pulse length in ns:\t%.3f\t%.3f\n", ...
%         [ii, 1/fxrabi(ifMax(ii))/2, 1/fxrabi(ifMax2(ii))/2])
% end

%% FIT cosine

% expcos = @(xx, A, tau, w, phi) A*exp(-xx/tau).*cos(2*pi*w*xx + phi);
expcos = @(xx, A, tau, w) A*exp(-xx/tau).*cos(2*pi*w*xx);
fitmodel = @(p) expcos(xrabi, 1, p(1), p(2));  % A = 1

fitOpt = optimoptions('lsqnonlin','Display','off');
figure(7)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
w0 = [40, 40]*1e-3;
for ii = 1:nMeas

    p0 = [500, w0(ii)];
    ydata = real(rabii(ii, :))/max(real(rabii(ii, :)));
    ydata = ydata';
    % plot(xrabi, ydata, 'o') 
    % hold on
    % plot(xrabi, fitmodel(p0))

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
xlim(setaxlim(1:nMeas, 1.1))


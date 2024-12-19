clearvars

%% IMPORT

generalFolder = "../data/raw/";
expName = "ZePSI-E-013002";  % ZePSI-E-001002
measFolder = generalFolder + expName;

[x0, y0, Param0] = loadfolderelexsys(measFolder);
nMeas = numel(Param0);

%% Scrollables

disp("Plot...")
% Traces real
figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x0{ii}{1}, x0{ii}{2}, real(y0{ii}));
    hold on
end

% Traces imag
figure(2)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    sax = ScrollableAxes();
    plot(sax, x0{ii}{1}, x0{ii}{2}, imag(y0{ii}));
    hold on
end

%% SEPARATE ESEEM FROM NUTATIONS

sEseem = size(y0{1});
iEseem = [];
iRabi1 = [];
allMeas = 1:nMeas;
for ii = allMeas
    if size(y0{ii}) == sEseem
        iEseem(end + 1) = ii;
        iRabi1(end + 1) = ii + 1;
    end
end
iRabi2 = allMeas(~ismember(allMeas, [iEseem, iRabi1]));

x1 = x0{1};
x2 = x0{2};
x3 = x0{3};

%% ESEEM RINGDOWN

disp("Subtract ringdown...")

inew = 0;
for ii = iEseem
    inew = inew + 1;
    y1{inew} = y0{ii} - y0{ii}(end, :);
end

figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 3:numel(iEseem)
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x1{1}, x1{2}, real(y1{ii}));
    hold on
    plot(sax, x1{1}, x1{2}, imag(y1{ii}));
    ylim([-1e4, max([real(y1{ii}), imag(y1{ii})], [], 'all')])
end

%% SAVE

disp("Save...")

% ESEEM
clear('x', 'y', 'Param')
pathToProcessed = "../data/processed/";
savePath = pathToProcessed + expName + "-ESEEM.mat";
x{1} = x1{1};
x{2} = x1{2}(1:end - 1);
% Cut out the last trace that is zero everywhere
for ii = 1:numel(iEseem)
    y{ii} = y1{ii}(1:end - 1, :);
end
Param = Param0(iEseem);
% save(savePath, 'x', 'y', 'Param')

% NUT1
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut1.mat";
x = x2;
y = y0(iRabi1);
Param = Param0(iRabi1);
% save(savePath, 'x', 'y', 'Param')

% NUT2
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut2.mat";
x = x3;
y = y0(iRabi2);
Param = Param0(iRabi2);
% save(savePath, 'x', 'y', 'Param')

%% PLOTS

iy = iEseem(14);
nSpectrum = 5;
i1 = round(linspace(1, size(y0{iy}, 1), nSpectrum));
i2 = 30:size(y0{iy}, 2);
% cmap = viridis(nSpectrum);
% figure()
clf
tL = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
yplot = real(y0{iy})*1e-4;
for ii = 1:nSpectrum
    plot(x1{1}(i2), yplot(i1(ii), i2),  'DisplayName', append(...
        char(hex2dec('03C4')), sprintf(' = %d ns', 150 + i1(ii) - 1)))
    hold on
end
xlim(setaxlim(x1{1}(i2), 1))
ylim(setaxlim(yplot(i1(2), i2), 1.05))
xticks(100:200:1000)
yticks((-5:2:5))
legend()
nexttile
yplot = imag(y0{iy})*1e-4;
for ii = 1:nSpectrum
    plot(x1{1}(i2), yplot(i1(ii), i2))  % + (ii - 1)*0.6e4
    hold on
end
xlim(setaxlim(x1{1}(i2), 1))
ylim(setaxlim(yplot(i1(2), i2), 1.05))
xticks(100:200:1000)
yticks((-5:0))
labelaxesfig(tL, "Time / ns",...
    ["In-phase / a.u.", "Out-of-phase / a.u."])
% nexttile(1); ylabel("In-phase channel / a.u.")
% nexttile(2); ylabel("Out-of-phase channel / a.u.")
savePath = "../images/" + "ZePSI-013002-01-eseemRaw.png";
% saveas(gcf, savePath)

%% ESEEM INTEGRATION (IGNORING OUT-OF-PHASE CHANNEL FOR NOW)
%{
% Find and plot integration window
ies = iEseem(1:5);
nTau = 280;
integWidth = 20;
iMax = 213;
figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = ies
    [integWindow{ii}, integWindowPlot{ii}] = deal(zeros(size(y1{ii})));
    for itau = 1:nTau
        valMax = max(real(y1{ii}(itau, :)));
        iInteg = (iMax - integWidth/2 + itau*2 - 2):(iMax + integWidth/2 + ...
            itau*2 - 2);
        integWindow{ii}(itau, iInteg) = ones(1, integWidth + 1);
        integWindowPlot{ii}(itau, iInteg) = ...
            valMax*integWindow{ii}(itau, iInteg);
    end
    y2{ii} = real(integWindow{ii}.*y1{ii});  % Signal in the window

    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x0{ii}{1}, x0{ii}{2}, real(y1{ii}));
    hold on
    plot(sax, x0{ii}{1}, x0{ii}{2}, real(integWindowPlot{ii}));
    plot(sax, x0{ii}{1}, x0{ii}{2}, y2{ii});
    ylim(setaxlim([-1000, max(real(y1{ii}), [], 'all')], 1.1))
    % xline([-1, 1]/2*integWidth + x0{ii}{1}(iMax), '--')
    % fprintf("%d max at:\t%d\n", [ii, iMax])
end

% Integrate to get Hahn-echo decay
hedecay = zeros(nMeas, nTau);
figure(2)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
xdecay = 0:2:(nTau - 1)*2;
for ii = ies
    hedecay(ii, :) = sum(y2{ii}(1:nTau, :), 2);

    nexttile
    plot(xdecay, hedecay(ii, :))
end

% Fit Hahn-echo decay
expdecay = @(p) p(1)*exp(-(150 + xdecay)/p(2));
iplot = 1;
for ii = ies
    p0 = [1e5, 200];
    [pfit{ii}, ~, residual, ~, ~, ~, jacobian] = ...
        lsqnonlin(@(p) hedecay(ii, :) - expdecay(p), p0);
    pci(ii, :, :) = nlparci(pfit{ii}, residual, 'jacobian', jacobian);

    nexttile(iplot)
    hold on
    plot(xdecay, expdecay(pfit{ii}))
    iplot = iplot + 1;
end

% Plot best fit results
iarr = 1;
for ii = ies
    ampls(iarr) = pfit{ii}(1);
    taus(iarr) = pfit{ii}(2);
    iarr = iarr + 1;
end

figure(3)
clf
errorbar(1:numel(ies), ampls, ...
    pci(ies, 1, 1) - ampls', pci(ies, 1, 2) - ampls', 'o-')
yyaxis right
errorbar(1:numel(ies), taus, ...
    pci(ies, 2, 1) - taus', pci(ies, 2, 2) - taus', 'o-')
xlim(setaxlim(1:numel(ies), 1.1))


%% ESEEM PHASE
iarr = 1;
for ii = ies
    ampls(iarr) = pfit{ii}(1);
    taus(iarr) = pfit{ii}(2);
    iarr = iarr + 1;
end

figure(3)
clf
errorbar(1:numel(ies), ampls, pci(ies, 1, 1), pci(ies, 1, 2), 'o-')
yyaxis right
plot(1:numel(ies), taus, pci(ies, 2, 1), pci(ies, 2, 2), 'o-')

disp("Phase correction...")

shiftphase = @(y, p) y*exp(1i*p*pi/180);

phis = 1:360;
for iphi = phis
    shiftedy{1}(iphi, :) = shiftphase(y1{end - 2}(103, :), iphi);
end

figure(3)
clf
sax = ScrollableAxes('Index', 1);
plot(sax, x0{1}{1}, phis, real(shiftedy{1}));
hold on
plot(sax, x0{1}{1}, phis, imag(shiftedy{1}));

figure(4)
clf
sax = ScrollableAxes('Index', 1);
for ii = iEseem
    plot(sax, x0{ii}{1}, x0{ii}{2}, real(y1{ii}));
    hold on
end
%}


%% RABI NUTATIONS EXPERIMENT PULSE 2 (ONLY IN-PHASE-CHANNEL FOR NOW)

iRabi2 = [3, 14];

iMax = 214;
integWidth = 20;
integWidth2 = 80;
nTau = size(y0{2}, 1);
figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = iRabi2
    [integWindow{ii}, integWindowPlot{ii}] = deal(zeros(size(y0{ii})));
    [integWindow2{ii}, integWindowPlot2{ii}] = deal(zeros(size(y0{ii})));
    for itau = 1:nTau
        valMax = max(real(y0{ii}(itau, :)));
        iInteg = (iMax - integWidth/2):(iMax + integWidth/2);
        iInteg2 = (iMax - integWidth2/2):(iMax + integWidth2/2);
        integWindow{ii}(itau, iInteg) = ones(1, integWidth + 1);
        integWindow2{ii}(itau, iInteg2) = ones(1, integWidth2 + 1);
        integWindowPlot{ii}(itau, iInteg) = ...
            valMax*integWindow{ii}(itau, iInteg);
        integWindowPlot2{ii}(itau, iInteg2) = ...
            valMax*integWindow2{ii}(itau, iInteg2);
    end
    y2{ii} = real(integWindow{ii}.*y0{ii});  % Signal in the window
    y22{ii} = real(integWindow2{ii}.*y0{ii});  % Signal in the window

    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x0{ii}{1}, x0{ii}{2}, real(y0{ii}));
    hold on
    plot(sax, x0{ii}{1}, x0{ii}{2}, real(integWindowPlot2{ii}));
    plot(sax, x0{ii}{1}, x0{ii}{2}, y22{ii});
    ylim(setaxlim([min(real(y0{ii}), [], 'all'), max(real(y0{ii}), ...
        [], 'all')], 1.1))
    % xline([-1, 1]/2*integWidth + x0{ii}{1}(iMax), '--')
    % fprintf("%d max at:\t%d\n", [ii, iMax])
end

% Integrate to get Rabi nutations
rabi = zeros(nMeas, nTau);
rabi2 = zeros(nMeas, nTau);
figure(2)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
xrabi = 2:2:nTau*2;
for ii = iRabi2
    rabi(ii, :) = sum(y2{ii}(1:nTau, :), 2);
    rabi2(ii, :) = sum(y22{ii}(1:nTau, :), 2);

    nexttile
    plot(xrabi, rabi(ii, :) - mean(rabi(ii, :)), 'o-')
    yyaxis right
    plot(xrabi, (rabi2(ii, :) - mean(rabi2(ii, :))), 'o-')
end

% FFT
nzf = 1024;  % Zero filling
tStep = xrabi(2) - xrabi(1);
fSampl = 1/tStep;
fxrabi = fSampl/nzf*(-nzf/2:nzf/2 - 1);
figure(4)
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = iRabi2
    rabi(ii, nzf) = 0;  % Zero filling
    rabi2(ii, nzf) = 0;  % Zero filling
    frabi(ii, :) = fft(rabi(ii, :) - mean(rabi(ii, :)));
    frabi2(ii, :) = fft(rabi2(ii, :) - mean(rabi(ii, :)));

    [~, ifMax(ii)] = max(abs(fftshift(frabi(ii, :))));
    [~, ifMax2(ii)] = max(abs(fftshift(frabi2(ii, :))));

    nexttile
    plot(fxrabi, abs(fftshift(frabi(ii, :))), 'o-')
    yyaxis right
    plot(fxrabi, abs(fftshift(frabi2(ii, :))), 'o-')
end

for ii = iRabi2
    fprintf("%d pi pulse length in ns:\t%.3f\t%.3f\n", ...
        [ii, 1/fxrabi(ifMax(ii))/2, 1/fxrabi(ifMax2(ii))/2])
end

%% RABI NUTATIONS EXPERIMENT PULSE 2 (ONLY IN-PHASE-CHANNEL FOR NOW)

iRabi1 = [2, 5, 7, 9, 11, 13, 16, 18];

iMax = 214;
integWidth = 20;
integWidth2 = 80;
nTau = size(y0{2}, 1);
figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = iRabi1
    [integWindow{ii}, integWindowPlot{ii}] = deal(zeros(size(y0{ii})));
    [integWindow2{ii}, integWindowPlot2{ii}] = deal(zeros(size(y0{ii})));
    for itau = 1:nTau
        valMax = max(real(y0{ii}(itau, :)));
        iInteg = (iMax - integWidth/2):(iMax + integWidth/2);
        iInteg2 = (iMax - integWidth2/2):(iMax + integWidth2/2);
        integWindow{ii}(itau, iInteg) = ones(1, integWidth + 1);
        integWindow2{ii}(itau, iInteg2) = ones(1, integWidth2 + 1);
        integWindowPlot{ii}(itau, iInteg) = ...
            valMax*integWindow{ii}(itau, iInteg);
        integWindowPlot2{ii}(itau, iInteg2) = ...
            valMax*integWindow2{ii}(itau, iInteg2);
    end
    y2{ii} = real(integWindow{ii}.*y0{ii});  % Signal in the window
    y22{ii} = real(integWindow2{ii}.*y0{ii});  % Signal in the window

    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x0{ii}{1}, x0{ii}{2}, real(y0{ii}));
    hold on
    plot(sax, x0{ii}{1}, x0{ii}{2}, real(integWindowPlot{ii}));
    plot(sax, x0{ii}{1}, x0{ii}{2}, y2{ii});
    ylim(setaxlim([min(real(y0{ii}), [], 'all'), max(real(y0{ii}), ...
        [], 'all')], 1.1))
    % xline([-1, 1]/2*integWidth + x0{ii}{1}(iMax), '--')
    % fprintf("%d max at:\t%d\n", [ii, iMax])
end

% Integrate to get Rabi nutations
rabii = zeros(nMeas, nTau);
rabii2 = zeros(nMeas, nTau);
figure(22)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
xrabi = 2:2:nTau*2;
for ii = iRabi1
    rabii(ii, :) = sum(y2{ii}(1:nTau, :), 2);
    rabii2(ii, :) = sum(y22{ii}(1:nTau, :), 2);

    nexttile
    % plot(xrabi, rabii(ii, :) - mean(rabii(ii, :)), 'o-')
    plot(xrabi, rabii(ii, :) - mean(rabii(ii, :)), 'o-')
    yyaxis right
    % plot(xrabi, (rabii2(ii, :) - mean(rabii2(ii, :)))*0.5, 'o-')
    plot(xrabi, (rabii2(ii, :)) - mean(rabii2(ii, :)), 'o-')
end

% FFT
nzf = 1024;  % Zero filling
tStep = xrabi(2) - xrabi(1);
fSampl = 1/tStep;
fxrabi = fSampl/nzf*(-nzf/2:nzf/2 - 1);
figure(4)
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = iRabi1
    rabii(ii, nzf) = 0;  % Zero filling
    rabii2(ii, nzf) = 0;  % Zero filling
    frabii(ii, :) = fft(rabii(ii, :) - mean(rabii(ii, :)));
    frabii2(ii, :) = fft(rabii2(ii, :) - mean(rabii(ii, :)));

    [~, ifMax(ii)] = max(abs(fftshift(frabii(ii, :))));
    [~, ifMax2(ii)] = max(abs(fftshift(frabii2(ii, :))));

    nexttile
    plot(fxrabi, abs(fftshift(frabii(ii, :))), 'o-')
    yyaxis right
    plot(fxrabi, abs(fftshift(frabii2(ii, :))), 'o-')
end

for ii = iRabi1
    fprintf("%d pi pulse length in ns:\t%.3f\t%.3f\n", ...
        [ii, 1/fxrabi(ifMax(ii))/2, 1/fxrabi(ifMax2(ii))/2])
end

%{
for ii = iEseem
    [~, bestPhase(ii)] = correctphase(y1{ii}(1, 150:300));
    for it = 1:size(y1{ii}, 1)
        y2{ii}(it, :) = shiftphase(y1{ii}(it, :), bestPhase(ii));
    end
end

% for ii = iEseem
%     for it = 1
%         [~, bestPhase(ii, it)] = correctphase(y1{ii}(it, 150:end));
%     end
% end

figure(3)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = iEseem
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x0{ii}{1}, x0{ii}{2}, real(y2{ii}));
end

figure(4)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = iEseem
    nexttile
    sax = ScrollableAxes();
    plot(sax, x0{ii}{1}, x0{ii}{2}, imag(y2{ii}));
end

% figure(5)
% clf
% tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
% for ii = iEseem
%     nexttile
%     plot(bestPhase(ii, :))
% end
figure(5)
clf
plot(bestPhase)

%% SAVE

% savePath = "../data/processed/" + expName + ".mat";
% save(savePath, 'x', 'y', 'Param')

%}


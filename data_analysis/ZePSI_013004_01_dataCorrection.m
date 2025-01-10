clearvars

%% IMPORT

generalFolder = "../data/raw/";
expName = "ZePSI-E-013004";
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
for ii = 1:numel(iEseem)
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x1{1}, x1{2}, real(y1{ii}));
    hold on
    plot(sax, x1{1}, x1{2}, imag(y1{ii}));
    ylim([-1e4, max([real(y1{ii}), imag(y1{ii})], [], 'all')])
end

%% ADJUST PHASE

shiftphase = @(y, p) y*exp(1i*p*pi/180);
nph = 360;
yp = zeros(nph, numel(x3{1}));
for ii = 1:nph
    yp(ii, :) = shiftphase(y0{iRabi2(2)}(1, :), ii);
end

% figure()
clf
sax = ScrollableAxes('Index', 1);
plot(sax, x3{1}, 1:nph, real(yp));
hold on
plot(sax, x3{1}, 1:nph, imag(yp));

bestPhase = 100;
for ii = 1:numel(iEseem)
    for itau = 1:size(y1{ii}, 1)
        y2{iEseem(ii)}(itau, :) = shiftphase(y1{ii}(itau, :), bestPhase);
    end
end

for ii = iRabi1
    for itau = 1:size(y0{ii}, 1)
        y2{ii}(itau, :) = shiftphase(y0{ii}(itau, :), bestPhase);
    end
end

for ii = iRabi2
    for itau = 1:size(y0{ii}, 1)
        inew = inew + 1;
        y2{ii}(itau, :) = shiftphase(y0{ii}(itau, :), bestPhase);
    end
end

%% SAVE

disp("Save...")
disp("The first two measurements have to be thrown away because of a phase shift or something else that caused a change in the background between the first and last scan (see ringdown-subtracted spectra).")

% ESEEM
clear('x', 'y', 'Param')
pathToProcessed = "../data/processed/";
savePath = pathToProcessed + expName + "-ESEEM.mat";
x{1} = x1{1};
x{2} = x1{2}(1:end - 1);
% Cut out the last trace that is zero everywhere
for ii = 3:numel(iEseem)
    y{ii} = y2{ii}(1:end - 1, :);
end
Param = Param0(iEseem(3:end));
% save(savePath, 'x', 'y', 'Param')

% NUT1
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut1.mat";
x = x2;
y = y2(iRabi1(3:end));
Param = Param0(iRabi1(3:end));
% save(savePath, 'x', 'y', 'Param')

% NUT2
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut2.mat";
x = x3;
y = y2(iRabi2(2));
Param = Param0(iRabi2(2));
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
savePath = "../images/" + "ZePSI-013004-01-eseemRaw.png";
% saveas(gcf, savePath)

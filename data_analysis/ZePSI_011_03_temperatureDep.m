%
clearvars

% Import
loadFileName = "../data/processed/ZePSI-E-011";

load(loadFileName)
xx = x;
yOld = y;
yLoad = yOld(1, 1:8);

loadFileName = "../data/processed/ZePSI-E-011-part2";
load(loadFileName)
nMeas = numel(y);

% for ii = 1:nMeas
%     yLoad{end + 1} = y{ii};
% end

% y_ = yLoad{9}; % 
yLoad(9:12) = y(4:-1:1);   % 230 240 250 260
% yLoad{10} = y{3};  % 240
% yLoad{11} = y{2};  % 250
% yLoad{12} = y{1};  % 260

nMeas = numel(yLoad);

%% 

[nt, nB] = size(yLoad{1});

y = yLoad;
% for ii = 1:nMeas
%     for iB = 1:nB
%         y{ii}(:, iB) = datasmooth(yLoad{ii}(:, iB), 20, 'savgol');
%     end
% end

% Scrollable traces
clf
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
for ii = 1:nMeas
    nexttile
    h = ScrollableAxes();
    % ii = 8;
    if ii < 9
        plot(h, xx{1}, xx{2}, yLoad{ii});
    else
        plot(h, x{1}, x{2}, yLoad{ii});
    end
    hold on
    title(string(60 + ii*20))
    % plot(h, x{1}, x{2}, y{ii}');
end

%% Scrollable field

% clf
% h = ScrollableAxes();
% ii = 2;
% plot(h, x{2}, x{1}, yLoad{ii}');
% hold on
% plot(h, x{2}, x{1}, y{ii}');

%%
% figure(2)
% plot(x{2}, y{1}(iTimeMax(1), :), x{2}, yLoad{1}(iTimeMax(1), :))


%%

% temperatures = [80:20:260, 265, 260:-10:220];
temperatures = [80:20:220, 230:10:260];
% temperatures = 260:-10:220;

noiseRange = [min(x{2}), 338.7; 341.5, max(x{2})];
[iTimeMax, snr0AtMax, snr1AtMax, ppAmpMax, sigAmpMax, ...
    noiseLev0AtMax, noiseLev1AtMax, peaksRatio, ...
    maxR, maxL] = deal(zeros(size(temperatures)));
for ii = 1:nMeas
    iTimeMax_ = 0;
    % ppAmpMax_ = 0;
    sigAmpMax_ = 0;
    snrAtMax_ = 0;
    noise1LevAtMax_ = 0;
    for it = 1:1200
        % [snr1, ppAmp, noiseLev1] = getSNR(x{2}, y{ii}(it, :), noiseRange);
        sigAmp = max(y{ii}(it, :));
        % if ppAmp > ppAmpMax_
        if sigAmp > sigAmpMax_
            iTimeMax_ = it;
            sigAmpMax_ = sigAmp;
            % [snr1AtMax_, ppAmpMax_, noiseLev1AtMax_] = getSNR( ...
            %     x{2}, y{ii}(it, :), noiseRange);
            if ii < 9
                [~, ~, noiseLev1AtMax_] = getSNR( ...
                    xx{2}, y{ii}(it, :), noiseRange);
                snr1AtMax_ = sigAmpMax_/noiseLev1AtMax_;
                [~, ~, noiseLev0AtMax_] = getSNR(...
                    xx{2}, yLoad{ii}(it, :), noiseRange);
                snr0AtMax_ = sigAmpMax_/noiseLev0AtMax_;
            else
                [~, ~, noiseLev1AtMax_] = getSNR( ...
                    x{2}, y{ii}(it, :), noiseRange);
                snr1AtMax_ = sigAmpMax_/noiseLev1AtMax_;
                [~, ~, noiseLev0AtMax_] = getSNR(...
                    x{2}, yLoad{ii}(it, :), noiseRange);
                snr0AtMax_ = sigAmpMax_/noiseLev0AtMax_;
            end
        end
    end

    iTimeMax(ii) = iTimeMax_;
    snr1AtMax(ii) = snr1AtMax_;
    snr0AtMax(ii) = snr0AtMax_;
    % ppAmpMax(ii) = ppAmpMax_;
    sigAmpMax(ii) = sigAmpMax_;
    noiseLev1AtMax(ii) = noiseLev1AtMax_;
    noiseLev0AtMax(ii) = noiseLev0AtMax_;
    maxR(ii) = max(y{ii}(iTimeMax(ii), :));
    maxL(ii) = max(y{ii}(iTimeMax(ii), 100));  % iLeftPeak = 100
    peaksRatio(ii) = maxR/maxL;

end
facMult = sigAmpMax(8)/0.52/sigAmpMax(end);

for ii = 1:nMeas

    % if ii >= 12
    %     sigAmpMax(ii) = sigAmpMax(ii)*1.1;
    %     snr1AtMax(ii) = snr1AtMax(ii)*1.1;
    %     snr0AtMax(ii) = snr0AtMax(ii)*1.1;
    % See ZePSI_011_02_factorBetweenDays.m
    if ii >= 5 && ii < 9 %&& strcmp(loadFileName, "../data/processed/ZePSI-E-011")
        % ppAmpMax(ii) = ppAmpMax(ii)/0.52;
        sigAmpMax(ii) = sigAmpMax(ii)/0.52;
        snr1AtMax(ii) = snr1AtMax(ii)/0.52;
        snr0AtMax(ii) = snr0AtMax(ii)/0.52;
    end
    if ii >= 9
        % ppAmpMax(ii) = ppAmpMax(ii)/2.56;
        sigAmpMax(ii) = sigAmpMax(ii)/1.1;
        snr1AtMax(ii) = snr1AtMax(ii)/1.1;
        snr0AtMax(ii) = snr0AtMax(ii)/1.1;
        
    end
end
figure(4)
xTicks = temperatures;
% xTicks = flip(temperatures);

clf
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile(1, [2, 4])
cmap = viridis(nMeas);
for ii = 1:nMeas
    % yplot = y{ii}(4000, :);
    if ii < 9
        xplot = xx;
        yplot = rescaledata(y{ii}(400, :), 'maxabs');
    else
        xplot = x;
        yplot = rescaledata(y{ii}(200, :), 'maxabs');
    end
    plot(xplot{2}, yplot, ...
        'DisplayName', string(temperatures(ii)) + " K", ...
        'Color', cmap(ii, :))
    hold on
end
legend()
xlim(setaxlim(x{2}, 1))
xline(x{2}(100), 'DisplayName', 'Position of shoulder')

nexttile
plot(temperatures, snr1AtMax, 'o-')  % Filtered ppAmp / filtered noise
yyaxis right
% plot(temperatures, snr0AtMax, 'o-')  % Filtered ppAmp / non-filtered noise
xlim(setaxlim(temperatures, 1.05))
% xticks(xTicks)
title('SNR at max')
legend('Smooth', 'No filter', 'Location', 'southwest')
%
nexttile
plot(temperatures, noiseLev1AtMax, 'o-')  % Filtered noise
yyaxis right
% plot(temperatures, noiseLev0AtMax, 'o-')  % Non-filtered noise
xlim(setaxlim(temperatures, 1.05))
% xticks(xTicks)
title('noiseLev at max')
%
nexttile
errorbar(temperatures, sigAmpMax, 2*noiseLev0AtMax, 'o-')
xline(140)
xline(220)
xlim(setaxlim(temperatures, 1.05))
% xticks(xTicks)
% title('ppAmp max')
title('max of signal')
%
nexttile
% errbar from non-filtered data. Also possible: hypot(1./maxR, 1./maxL)
errbar = peaksRatio.*noiseLev0AtMax.*(1./maxR + 1./maxL);
errorbar(temperatures, peaksRatio, errbar, 'o-')
hold on
xlim(setaxlim(temperatures, 1.05))
ylim(setaxlim(peaksRatio, 1.5))
xticks(xTicks)
title('peaks ratio at max')

%%

figure()
errorbar(temperatures, sigAmpMax/max(sigAmpMax), 5*noiseLev0AtMax/max(sigAmpMax), 'o-k')

xlim(setaxlim(temperatures, 1.1))
ylim(setaxlim(sigAmpMax/max(sigAmpMax), 1.125))
xticks(100:40:260)
yticks(0:0.5:1)

set(gca, 'Xdir', 'reverse')
labelaxesfig(gca, 'Temperature / K', 'trEPR amplitude norm. / a.u.')

saveFolder = "/net/storage/gianlum33/projects/oop_ciss_calculations/reports/20241018_castleMonthlyUpdates";

savePath = fullfile(saveFolder, "zechPSI_trEPRsigAmp_vsTemperature_01.png");
% exportgraphics(gcf, savePath, 'Resolution', 600)

%%
% 
% xplot = xx{2} + 1.9;
% yplot = rescaledata(y{4}(1000, :), 'maxabs');
% figure()
% plot(xplot, yplot, 'k')
% 
% xlim(setaxlim(xplot, 1))
% ylim(setaxlim(yplot, 1.05))
% xticks(340:344)
% yticks(0:1)
% 
% labelaxesfig(gca, 'Magnetic Field / mT', 'trEPR signal norm. / a.u.')
% 
% saveFolder = "/net/storage/gianlum33/projects/oop_ciss_calculations/reports/20241018_castleMonthlyUpdates";
% 
% savePath = fullfile(saveFolder, "zechPSI_trEPRsigAmp_vsTemperature_xxxx.png");
% exportgraphics(gcf, savePath, 'Resolution', 600)


%% Check if integral in the field domain is zero

% loadFileName = '../data/processed/ZePSI-E-011-dataCorr2.mat';

% ab = load(loadFileName);
% yLoad = ab.y;

figure(13)
clf
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
% for ii = 1:nMeas
%     nexttile
%     yplot = cumtrapz(yLoad{ii}(iTimeMax(ii), :));
%     yplot2 = cumtrapz(y{ii}(iTimeMax(ii), :));
%     plot(yplot)
%     hold on
%     plot(yplot2);
%     legend(string(yplot(end)))
%     title(string(temperatures(ii)))
% end
for ii = 1:nMeas
    nexttile
    h = ScrollableAxes();
    % ii = 8;
    plot(h, x{1}, x{2}, yLoad{ii});
    hold on
    title(string(60 + ii*20))
    % plot(h, x{1}, x{2}, y{ii}');
end

%%

% loadFileName = '../data/processed/ZePSI-E-011.mat';

% aa = load(loadFileName);
% nMeas = numel(aa.y);

saveFolder = "/net/storage/gianlum33/projects/oop_ciss_calculations/reports/20241018_castleMonthlyUpdates";

temperOld = [80:20:260, 265];
iTemps = [4, 8, 9, 10];
iTimes = [500, 1500, 2000, 3000];
cmap = inferno(numel(iTemp) + 1);
figure(3)
clf
tL = tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'compact');
for ii = 1:numel(iTemp)
    iTemp_ = iTemps(ii);
    for it = 1:numel(iTimes)
        nexttile(it)
        
        it_ = iTimes(it);
    
        if ii ~= 4 || it == 1
            yplot = rescaledata(yOld{iTemp_}(it_, :), 'maxabs');
            xplot = xx{2};
        % yplot = y{ii}(it, :);
        else
            yplot = rescaledata(yLoad{end}(floor(it_/2), :), 'maxabs');
            yplot = datasmooth(yplot, 2, 'savgol');
            xplot = x{2};
        end
        plot(xplot, yplot, ...
                'DisplayName', string(temperOld(iTemp_)) + " K", ...
                'Color', cmap(ii, :))
        
        xlim(setaxlim(xx{2}, 1))
        % ylim(setaxlim(...
        %     [rescaledata(yOld{4}(it_, :), 'maxabs'), ...
        %     rescaledata(yLoad{end}(it_, :), 'maxabs') ] , 1.05))
        ylim([-0.85, 1.05])
        hold on
        if ii == 1
            text(341, -0.65, "t = " + string(round(xx{1}(it_)/1e3)) + " us", ...
                'FontSize', 14)
        end
    end
    legend('Location', 'eastoutside')
    labelaxesfig(tL, 'Magnetic field / mT', 'trEPR signal norm. / a.u.')

    savePath = fullfile(saveFolder, ...
        "zechPSI_trEPR-timeSlices-temperature_0" + string(ii) + ".png");
    % exportgraphics(gcf, savePath, 'Resolution', 600)
end


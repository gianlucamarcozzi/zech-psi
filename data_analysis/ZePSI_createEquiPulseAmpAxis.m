clearvars

loadFolder = "../data/processed/ZePSI_*mpfuYch.mat";
loadDir = dir(loadFolder);
nData = numel(loadDir);
xx = linspace(78, 83.5, 1000);

for ii = 1:nData
    filepath = fullfile(loadDir(ii).folder, loadDir(ii).name);
    aa = load(filepath);
    sigfun{ii} = aa.sigfun;
    pfit{ii} = aa.pfitsig;
    ympfus(ii, :) = sigfun{ii}(xx, pfit{ii});
end

ympfu = mean(ympfus, 1);

figure(1)
clf
plot(xx, ympfu, 'DisplayName', 'Average', 'Color', 'black')
hold on
for ii = 1:nData
    plot(xx, sigfun{ii}(xx, pfit{ii}), "DisplayName", 'off')
end
legend('Location', 'northwest')

ny = 25;
ygoal = linspace(ympfu(1), ympfu(end), ny);
ixgoal = zeros(ny, 1);
for ii = 1:ny
    [~, ixgoal(ii)] = min(abs(ympfu - ygoal(ii)));
end
xgoal = xx(ixgoal);

plot(xgoal, ympfu(ixgoal), 'k', 'Marker', 'square', 'LineStyle', 'none')

% SAVE TO FILE
fileID = fopen('ZePSI_createEquiPulseAmpAxis.txt', 'w');
fprintf(fileID, '[');
for ii = 1:ny - 1
    fprintf(fileID, '%.3f, ', xgoal(ii));
end
fprintf(fileID, '%.3f', xgoal(end));
fprintf(fileID, ']');
fclose(fileID);


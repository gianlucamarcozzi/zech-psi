clearvars

%% IMPORT
generalFolder = "../data/raw/";
expName = "ZePSI-E-014003";
measFolder = generalFolder + expName;
filepathParam = append(generalFolder, expName, '/', expName, '-param.txt');

% Load separately each measurement type
[xel0, yel0, paramel0] = loadfolderelexsys(measFolder, "*1.DTA");
[xed0, yed0, paramed0] = loadfolderelexsys(measFolder, "*2.DTA");
[xst0, yst0, paramst0] = loadfolderelexsys(measFolder, "*3.DTA");
[xro0, yro0, paramro0] = loadfolderelexsys(measFolder, "*4.DTA");
[xrt0, yrt0, paramrt0] = loadfolderelexsys(measFolder, "*5.DTA");

% Store mpfu amplitude of first pulse in param structure
xAmp = getamplitudefirstpulse(filepathParam);
for ii = 1:numel(xAmp)
    paramel0{ii}.mpfuAmpPulse1 = xAmp(ii);
    paramro0{ii}.mpfuAmpPulse1 = xAmp(ii);
end

%% SAVE
pathToProcessed = "../data/processed/";

% ESEEM
clear('x', 'y', 'param')
savePath = pathToProcessed + expName + "-ESEEM-light.mat";
x = xel0{1};
y = yel0;
param = paramel0;
save(savePath, 'x', 'y', 'param')

% ESEEM-dark
clear('x', 'y', 'param')
savePath = pathToProcessed + expName + "-ESEEM-dark.mat";
x = xed0{1};
y = yed0;
param = paramed0;
save(savePath, 'x', 'y', 'param')

% Standing ESE
clear('x', 'y', 'param')
savePath = pathToProcessed + expName + "-Standing-ESE.mat";
x = xst0{1};
y = yst0;
param = paramst0;
save(savePath, 'x', 'y', 'param')

% NUT1
clear('x', 'y', 'param')
savePath = pathToProcessed + expName + "-Nut1.mat";
x = xro0;
y = yro0;
param = paramro0;
save(savePath, 'x', 'y', 'param')

% NUT2
clear('x', 'y', 'param')
savePath = pathToProcessed + expName + "-Nut2.mat";
x = xrt0;
y = yrt0;
param = paramrt0;
save(savePath, 'x', 'y', 'param')

%%
function xAmp = getamplitudefirstpulse(filepathParam)
    % Load the amplitude of the mpfu of the first pulse
    xAmp = table2array(readtable(filepathParam, "Range", "B1"));
    % Check for nans and cut the array as soon as one is found
    isnanxAmp = isnan(xAmp);
    if sum(isnanxAmp)
        for ii = 1:numel(xAmp)
            if isnanxAmp(ii) == 1
                xAmp = xAmp(1:ii - 1);
                break
            end
        end
    end
end
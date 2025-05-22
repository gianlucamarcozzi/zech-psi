clearvars

%% IMPORT
generalFolder = "../data/raw/";
pathToProcessed = "../data/processed/";
expName = "ZePSI-E-014-007";
measFolder = generalFolder + expName;
filepathParam = append(generalFolder, expName, '/', expName, '-param.txt');

% Store mpfu amplitude of first pulse in param structure
xAmp = getamplitudefirstpulse(filepathParam);

titles = ["-ESEEM-light-plus", "-ESEEM-light-minus", ...
 "-ESEEM-dark-plus", "-ESEEM-dark-minus", ...
 "-Standing-ESE-plus", "-Standing-ESE-minus", ...
 "-Nut1-plus", "-Nut1-minus", "-Nut2"];
    
for ii = 1:numel(titles)
    ending = sprintf("*%i.DTA", ii);
    [x0, y, Param] = loadfolderelexsys(measFolder, ending);
    for iamp = 1:numel(Param)
        Param{iamp}.pAmpPulse1 = xAmp(iamp);
    end
    if contains(titles(ii), "ESEEM") || contains(titles(ii), "Standing")
        x = x0{1};
    else
        x = x0;
    end
    savePath = pathToProcessed + expName + titles(ii) + ".mat";
    save(savePath, "x", "y", "Param");
end

% Load separately each measurement type
%{
[xelp0, yelp0, paramelp0] = loadfolderelexsys(measFolder, "*1.DTA");
[xelm0, yelm0, paramelm0] = loadfolderelexsys(measFolder, "*2.DTA");
[xedp0, yedp0, paramedp0] = loadfolderelexsys(measFolder, "*3.DTA");
[xedm0, yedm0, paramedm0] = loadfolderelexsys(measFolder, "*4.DTA");
[xstp0, ystp0, paramstp0] = loadfolderelexsys(measFolder, "*5.DTA");
[xstm0, ystm0, paramstm0] = loadfolderelexsys(measFolder, "*6.DTA");
[xrop0, yrop0, paramrop0] = loadfolderelexsys(measFolder, "*7.DTA");
[xrom0, yrom0, paramrom0] = loadfolderelexsys(measFolder, "*8.DTA");
[xrt0, yrt0, paramrt0] = loadfolderelexsys(measFolder, "*9.DTA");
%}

%% SAVE

% ESEEM
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-ESEEM-light.mat";
x = xel0{1};
y = yel0;
Param = paramel0;
save(savePath, 'x', 'y', 'Param')

% ESEEM-dark
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-ESEEM-dark.mat";
x = xed0{1};
y = yed0;
Param = paramed0;
save(savePath, 'x', 'y', 'Param')

% Standing ESE
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Standing-ESE.mat";
x = xst0{1};
y = yst0;
Param = paramst0;
save(savePath, 'x', 'y', 'Param')

% NUT1
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut1.mat";
x = xro0;
y = yro0;
Param = paramro0;
save(savePath, 'x', 'y', 'Param')

% NUT2
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut2.mat";
x = xrt0;
y = yrt0;
Param = paramrt0;
save(savePath, 'x', 'y', 'Param')

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
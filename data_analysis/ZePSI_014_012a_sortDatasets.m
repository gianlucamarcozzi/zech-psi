clearvars

%% IMPORT
generalFolder = "../data/raw/";
pathToProcessed = "../data/processed/";
expName = "ZePSI-E-014-012";
measFolder = generalFolder + expName;
filepathParam = append(generalFolder, expName, '/', expName, '-param.txt');

% Store mpfu amplitude of first pulse in param structure
xAmp = getamplitudefirstpulse(filepathParam);

titles = ["-ESEEM-light", "-ESEEM-dark", "-Standing-ESE", ...
 "-Nut1", "-Nut2"];

ifile = 0;
for ii = 1:numel(titles)
    if ~strcmp(titles(ii), "-Nut2")
        for jj = 1:2
            % BrX and MinBrX will be saved in the same structure
            ifile = ifile + 1;
            ending = sprintf("*%i.DTA", ifile);
            [x0, y{jj}, Param{jj}] = loadfolderelexsys(measFolder, ending);
        end
    else
        % No plus and minus for Nut2
        ifile = ifile + 1;
        ending = sprintf("*%i.DTA", ifile);
        [x0, y, Param] = loadfolderelexsys(measFolder, ending);        
    end
    % Save pulse amplitudes
    if numel(Param) == 2
        for iamp = 1:numel(Param{1})
            for jj = 1:2
                Param{jj}{iamp}.pAmpPulse1 = xAmp(iamp);
            end
        end
    else
        Param{iamp}.pAmpPulse1 = xAmp(iamp);
    end            
    % Save only one x-axis for ESEEM and Standing
    if contains(titles(ii), "ESEEM") || contains(titles(ii), "Standing")
        x = x0{1};
    else
        x = x0;
    end
    % Save
    savePath = pathToProcessed + expName + titles(ii) + ".mat";
    save(savePath, "x", "y", "Param");
end

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
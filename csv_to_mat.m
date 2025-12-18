%% =======================================================================
%   Convert all CSV files → signals.mat
%   - Skip first 44 rows
%   - Row 45 = header row
%   - Use column 1 (time) and column 2 (magnitude)
%   - Store each signal as a struct with fields: time, magnitude
% =======================================================================

clear; clc;
addpath('C:\Users\...'); %add path of raw CSV data

inputFolder = uigetdir(pwd, "Select folder containing CSV files");

% inputFolder = "raw_data";        % folder containing CSV files
outputFile  = "signals.mat"; % final MAT file

csvFiles = dir(fullfile(inputFolder, "*.csv"));
numFiles = length(csvFiles);

signals = cell(numFiles,1);  % cell array to hold all samples
ref_bp = cell(numFiles, 1);

fprintf("Found %d CSV files.\n", numFiles);

for k = 1:numFiles

    filePath = fullfile(inputFolder, csvFiles(k).name);
    fprintf("Reading %s...\n", csvFiles(k).name);

    % --------------------------------------------------------------
    % 1. Read table starting from row 45 (Header row)
    % --------------------------------------------------------------
    opts = detectImportOptions(filePath);
    opts.DataLines = [46 Inf];    % skip first 44 rows
    opts.SelectedVariableNames = opts.VariableNames(1:8); 
    % selects first 2 columns only (time, magnitude)

    T = readtable(filePath, opts);

    % --------------------------------------------------------------
    % 2. Extract relevant columns
    % --------------------------------------------------------------
    time = T{:,1};
    magnitude = T{:,2};

    % Ensure column vectors
    time = time(:);
    magnitude = magnitude(:);

    % Store into cell
    signals{k}.time = time;
    signals{k}.magnitude = magnitude;
    signals{k}.filename = csvFiles(k).name;

    % Repeat for ref cuff readings
    % opts.SelectedVariableNames = opts.VariableNames(5:8); 
    time_ref = T{:,5};
    HR_ref = T{:,6};
    SBP_ref = T{:,7};
    DBP_ref = T{:,8};

    time_ref = time_ref(:);
    HR_ref = HR_ref(:);
    SBP_ref = SBP_ref(:);
    DBP_ref = DBP_ref(:);

    ref_bp{k}.time = time_ref;
    ref_bp{k}.HR = HR_ref;
    ref_bp{k}.SBP = SBP_ref;
    ref_bp{k}.DBP = DBP_ref;
    ref_bp{k}.filename = csvFiles(k).name;


end

% --------------------------------------------------------------
% 3. Save everything
% --------------------------------------------------------------
save(outputFile, "signals");
save("ref_bp.mat", "ref_bp");

fprintf("\nConversion complete! Saved to %s\n", outputFile);
fprintf("\nSaved reference BP data → ref_bp.mat\n");


%% =======================================================================
%   Pair each signals{k} with ref_bp{k}
%   Output: paired_data.mat containing paired{k}.signal and paired{k}.ref
% =======================================================================

clear; clc;

load signals.mat      % signals{k}.time, signals{k}.magnitude
load ref_bp.mat       % ref_bp{k}.time, ref_bp{k}.HR, ref_bp{k}.SBP, ref_bp{k}.DBP

numSamples = length(signals);

paired = cell(numSamples,1);

for k = 1:numSamples
    
    paired{k}.filename = signals{k}.filename;  % same filename
    
    % ---- Signal data ----
    paired{k}.signal.time      = signals{k}.time;
    paired{k}.signal.magnitude = signals{k}.magnitude;
    
    % ---- Reference BP data ----
    paired{k}.ref.time = ref_bp{k}.time;
    paired{k}.ref.HR   = ref_bp{k}.HR;
    paired{k}.ref.SBP  = ref_bp{k}.SBP;
    paired{k}.ref.DBP  = ref_bp{k}.DBP;

end

save("paired_data.mat","paired");

fprintf("SUCCESS: paired_data.mat created.\n");

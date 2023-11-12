%% Compress all given aMALDIms files from the aMALDIms (hdf5 file) to the reduced aMALDIms format (mat file)

% List of files to convert
files = {   
    'C:\Users\ChrisTow\Downloads\Converted\e0165p01t03.hdf5',...
    };

% Path to save log
path.temp = 'C:\Users\ChrisTow\Downloads\Converted'; %'__temporary__'; 
if ~exist(path.temp,'dir')
    mkdir(path.temp);  
end

% Start log
diary(fullfile(path.temp,'compression_diary.txt'));

% Convert given aMALDIms files to aMALDImsRed files
amaldims = aMALDIms('');
[compressfiles, duration] = amaldims.compress_multiple_respectra_by_threshold(files,...
    'Algorithm','PeakThreshold',... % compression method
    'Threshold',30,...              % signal cutoff value
    'FindPeaks',true,...            % also save all spectra peaks
    'NormalizeIntensity',true,...   % intensity (false) or snr (true)
    'IncludeRawData',false);        % save raw time-series data

% End log
save(fullfile(path.temp,'compression_log.mat'),'compressfiles','duration')
diary off

%% Review a compressed file
layout = '1536'; % plate layout (if ignored will default to NxN where N is the square root of the number of positions)
query_mass = 442.267627; % mass to look for (mass lock: Lib1.F13: 442.267627 Da)
mass_tolerance = 0.01; % mass tolerance (window will be query_mass +/- mass_tolerance Da)
object = aMALDImsRed(compressfiles{1},layout); % load a compressed file
object.summary_plot(); % plots of measurement stats
object.summary_plateview(); % plate view of measurement stats
object.heatmap(); % heatmap of all reduced spectra
object.view_mass_window(query_mass,mass_tolerance); % mass search

%% Previously Compressed Files

% % 2019-12-31
% files = {   
%     ['C:\Users\ChrisTow\Desktop\Converted\0147_Ugi_Reaction_Timed_2_redo\e0147p01t01.hdf5'],...
%     ['C:\Users\ChrisTow\Desktop\Converted\0144_Ugi_Reaction_Timed_2\e0144p01t02.hdf5'],...
%     ['C:\Users\ChrisTow\Desktop\Converted\0144_Ugi_Reaction_Timed_2\e0144p01t03.hdf5'],...
%     ['C:\Users\ChrisTow\Desktop\Converted\0143_Ugi_Reaction_Timed\e0143p01t01.hdf5'],...
%     ['C:\Users\ChrisTow\Desktop\Converted\0143_Ugi_Reaction_Timed\e0143p01t02.hdf5'],...
%     ['C:\Users\ChrisTow\Desktop\Converted\0143_Ugi_Reaction_Timed\e0143p01t03.hdf5'],...
%     ['C:\Users\ChrisTow\Desktop\Converted\0143_Ugi_Reaction_Timed\e0143p01t04.hdf5'],...    
% };

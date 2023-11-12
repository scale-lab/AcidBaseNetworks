%{
@title: test_convert_spectra_read.m
@description: Example MATLAB script to test to reading from previously converted raw data (made by running "test_convert_raw.py")
@author: chrisarcadia 
@created: 2018/10/29
%}

%% Set file
filename = 'C:\Users\ChrisTow\Desktop\Converted\single\Mix_1_100_1.hdf5';
h5disp(filename,'/','simple') % display strucutre of file
hinfo = h5info(filename); % grab file info

%% Load data

% grab file info
positions = h5read(filename,'/spectra/position');
Npos = length(positions);
index = 1; 
position = positions{index};
sMZ = '/spectra/mass_to_charge';
sSignal = '/spectra/signal';
hMZ = h5info(filename,sMZ);
hSignal = h5info(filename,sSignal);
hShape = hSignal.Dataspace.Size;
Npts = hShape(1);
mz = h5read(filename,sMZ,[1,index],[Npts,index]);
signal = h5read(filename,sSignal,[1,index],[Npts,1]);

% grab scan log
TIC = h5read(filename,'/scan/TIC');
elapsed = h5read(filename,'/scan/elapsedMinutes');

% grab a few settings
laser_power = h5readatt(filename,'/method/laser','power');
laser_shots = h5readatt(filename,'/method/laser','shot_count');
laser_freq = h5readatt(filename,'/method/laser','frequency');
laser_focus = h5readatt(filename,'/method/laser','focus');
ion_polarity = h5readatt(filename,'/method/ion','polarity');
walk_distance = h5readatt(filename,'/method/walk','grid_width');

% grab all available method info
info = struct();
attributes = hinfo.Groups(1).Attributes;
for n = 1:length(hinfo.Groups)
    if strcmp(hinfo.Groups(n).Name,'/method')
        method_info = hinfo.Groups(n).Groups;
        for m = 1:length(method_info)
            category_info = method_info(m).Attributes;
            category = method_info(m).Name;
            category = strsplit(category,'/');
            category = category{end};
            for o = 1:length(category_info)
                info.(category).(category_info(o).Name) = category_info(o).Value;
            end
        end
        for m = 1:length(attributes)
            info.(attributes(n).Name) = attributes(n).Value;
        end
        
    end
end
disp(info)

%% Plot one of the spectra
figure('color','w','name','Mass Spectrum')
plot(mz,signal)
xlabel('m/z')
ylabel('Signal')
box on
grid on
axis tight
title(position)

%% Plot the total ion current 
figure('color','w','name','Total Ion Current')
plot(elapsed,TIC)
xlabel('Time Elapsed [min]')
ylabel('TIC')
box on
grid on
axis tight
title(position)






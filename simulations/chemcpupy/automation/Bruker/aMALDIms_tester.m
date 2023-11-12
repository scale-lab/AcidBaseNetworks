%{
@title: aMALDIms_tester.m
@description: Tester script for the aMALDIms class
@author: Chris Arcadia (christopher_arcadia@brown.edu)
@created: 2018/11/13
%}

%% Set file to read

filename = 'C:\Users\ChrisTow\Desktop\Converted\lib6\e0116p03t01.hdf5';
mass_lock = 442.267627; % UgiLib1_F13

%% Get file information
object = aMALDIms(filename);
%object.print(); % display file structure
%tree = object.tree(); % get handle to file structure
ion_mode = object.method.ion.polarity;

%% Plot the total ion current at all positions

if object.has.scan
    figure('color','w','name','Total Ion Current')
    plot(object.scan.elapsedMinutes,object.scan.TIC,'ob')
    xlabel('Time Elapsed [min]')
    ylabel('Total Ion Current (TIC)')
    box on
    grid on
    axis tight
    set(gca,'yscale','log')
end

%% PlateView demonstrations - Plot the total ion current and time elapsed at each location

if object.has.scan
    %zlimits = max(object.scan.TIC)*[0.1,1]; % only show spots with TIC values between 10-100% of the max
    zlimits = []; % show all spots with TIC values
    markersize = 30; % set spot size
    object.plateview('Total Ion Current (TIC)',object.scan.TIC,'markersize',markersize,'logscale',true,'zLimits',zlimits);
    object.plateview('Time Elapsed [min]',object.scan.elapsedMinutes,'markersize',markersize,'logscale',false,'zLimits',zlimits);
    object.plateview('Time per Spot [min]',[0,diff(object.scan.elapsedMinutes)'],'markersize',markersize,'logscale',false,'zLimits',zlimits);
    custom_var = [1:30]'; 
    custom_pos = [custom_var,custom_var];
    custom_pos_label = object.positionInt2Str(custom_pos);
    object.plateview('Custom Variable and Positions',custom_var,'positions',custom_pos,...
                    'markersize',markersize,'markeredgecolor','k','colormap','jet',...
                    'linewidth',0.25,'logscale',false,'zLimits',zlimits);    
end

%% Get raw (time) and spectra data for a single position

% set position to view
index = 1;
position = object.positions{index};
count = 1;

% get time domain
if object.has.raw
    [raw.signal, raw.time, found] = object.get_data_raw(index,count);    
end

% get spectrum
if object.has.spectra
    [spectra.signal, spectra.mz, found] = object.get_data_spectra(index,count);    
end

% get resampled spectrum
if object.has.respectra
    [respectra.signal, respectra.mz, found] = object.get_data_respectra(index,count);    
end

% Plot one of the spectra and time series
if object.has.raw || object.has.spectra || object.has.respectra
    figure('color','w','name',position)
    subplot(2,1,1) % time series
    if object.has.raw
        plot(raw.time,raw.signal,'-b')
    end
    xlabel('Time [s]')
    ylabel('Signal')
    box on
    grid on
    axis tight
    title('Raw Time Series')
    subplot(2,1,2) % spectra
    hold on
    if object.has.spectra
        plot(spectra.mz,spectra.signal,'-b')
    end
    if object.has.respectra
        plot(respectra.mz,respectra.signal,'-r')
    end
    if object.has.spectra && object.has.respectra
        legend('original','resampled')
    end
    hold off
    xlabel('m/z')
    ylabel('Signal')
    box on
    grid on
    axis tight
    title('Saved Mass Spectrum')
end


%% Plot the intensity of all samples at a certain mass

mass_to_query = mass_lock; % [Da/C]
mass_window_half_width = 0.000; % 0.005; % width of the window [Da/C]

if object.has.respectra
    
    % get window data
    [intensities, window] = object.mass_shift_for_intensity(mass_to_query, 'WindowHalfWidth', mass_window_half_width);
    
    % plot intensities
    figure('color','w','name','Intensities at Given Mass')
    plot(intensities,'ob')
    xlabel('Sample Number')
    ylabel('Intensity')
    box on
    grid on
    % axis tight
    set(gca,'yscale','log')    
    title(['mz: ' num2str(window.mass.center) ' [-' num2str(abs(window.mass_delta(1)))  ', +' num2str(abs(window.mass_delta(2))) ']']);
    
    % plot hit map
    object.plateview(['Intensity at mz: ' num2str(window.mass.center)],intensities,'markersize',30,'logscale',true,'zLimits',[]);    
end





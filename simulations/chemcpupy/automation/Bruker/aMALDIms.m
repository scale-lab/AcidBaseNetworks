
classdef aMALDIms < handle
    % @title: AMALDIMS.m
    % @description: Class for accessing the data in an aMALDIms file (hdf5 file)
    % @author: Chris Arcadia (christopher_arcadia@brown.edu)
    % @created: 2018/11/02
    
    properties % class variables        
        verbose = true; % turn on or off verbose warnings and messages      
        
        filename = ''; % filename of the aMALDIms file
        loaded = false; % boolean indicating if file was successfully loaded
        
        info = struct(); % info about the file
        method = struct(); % info about measurement methods
        scan = struct(); % the measurement scan log
        options = struct(); % conversion options used to generated the file
        raw = struct(); % info about the raw (time domain) data
        spectra = struct(); % info about the spectral (mass/charge domain) data
        respectra = struct(); % info about the resampled spectral (mass/charge domain) data
        
        positions = {}; % MALDI spot positions
        duplicates = []; % position index (sample number) of repeated measurements (if present)        
        plate = struct(); % MALDI plate layout info {'1536', '384', 'unknown'}
        
        has = struct(); % structure indicating what datasets and groups have been found and loaded
        num = struct(); % stucture providing the number of points in (length of) the datasets
    end
    
    properties (Constant) % class constants        
                
        versions = {'1.0'}; % compatible versions of aMALDIms files, those supported by this class
                
        location = struct(  'Root', '/',... % location of root/top-most level of the file
                         'Spectra', '/spectra',... % location of spectra
                       'SpectraMZ', '/spectra/mass_to_charge',... % location of spectra mass-to-charge vector(s)
                      'SpectraSig', '/spectra/signal',... % location of spectra signal vectors
                      'SpectraPos', '/spectra/position',... % location of spectra position vector
                      'SpectraPar', '/spectra/parameters',... % location of spectra parameters
                      'SpectraSta', '/spectra/statistics',... % location of spectra statistics                      
                      'SpectraBac', '/spectra/background',... % location of spectra background                                                                  
                         'ReSpectra', '/respectra',... % location of resampled spectra
                       'ReSpectraMZ', '/respectra/mass_to_charge',... % location of resampled spectra mass-to-charge vector(s)
                      'ReSpectraSig', '/respectra/signal',... % location of resampled spectra signal vectors
                      'ReSpectraPos', '/respectra/position',... % location of resampled spectra position vector
                      'ReSpectraPar', '/respectra/parameters',... % location of resampled spectra parameters     
                      'ReSpectraSta', '/respectra/statistics',... % location of resampled spectra statistics                      
                      'ReSpectraBac', '/respectra/background',... % location of resampled spectra background  
                      'ReSpectraSum', '/respectra/signal_sum',... % location of resampled spectra background                                                                  
                             'Raw', '/raw',... % location of raw         
                          'RawTim', '/raw/time',... % location of raw time vector
                          'RawSig', '/raw/signal',... % location of raw signal vectors 
                          'RawPos', '/raw/position',... % location of raw position vector
                          'Method', '/method',... % location of method info
                            'Scan', '/scan',... % location of scan log
                         'ScanPos', '/scan/position',... % location of scan log position vector
                         'ScanTIC', '/scan/TIC',... % location of scan total ion current vector
                         'ScanTim', '/scan/elapsedMinutes',... % location of scan elapsed time [min] vector                            
                         'ScanMax', '/scan/maxPeak',... % location of scan max peak vector
                         'ScanNum', '/scan/number',... % location of scan number vector
                          'Option', '/options',... % location of conversion options
                             'etc', ''); % location strings
    end
    
    methods % main
                
        function object = aMALDIms(filename) 
        % construct an instance of this class   
            % set filename
            object.filename = filename;
            if exist(filename, 'file') == 2
                object.load();
            else
                if isempty(filename)
                    % do nothing (using the class for its functions only)
                else
                    warning('File not found.')                
                end
            end
        end
        
        function load(object)
        % load the contents of a file (excluding datasets)
            % load info about the file and its contents (group info and presence)
            [object.info, object.has.info] = object.load_file_info();
            [object.options, object.has.options] = object.load_options();            
            [object.method, object.has.method] = object.load_method_info();
            [object.scan, object.has.scan] = object.load_scan_log();
            [object.raw, object.has.raw] = object.load_raw_info();
            [object.spectra, object.has.spectra] = object.load_spectra_info();
            [object.respectra, object.has.respectra] = object.load_respectra_info();            
            [object.positions, object.has.positions, object.duplicates] = object.load_positions(); % get positions while checking that the various versions of the position vector are equivalent
            object.has.duplicates = not(isempty(object.duplicates));
            object.plate = object.determine_plate(); 
            % get the size of each dataset
            object.num.positions = length(object.positions);    
            fieldsData = {'raw','spectra','respectra'}; % datasets that can be loaded laster, as needed (partially or fully)
            for n = 1:length(fieldsData)
                object.num.(fieldsData{n}) = 0; % initialize field
                if object.has.(fieldsData{n})
                    if isfield(object.(fieldsData{n}),'points')
                        object.num.(fieldsData{n}) = double(object.(fieldsData{n}).points);
                    end
                end
            end
            % ensure that certain fields are present 
            fieldsRequired = {'info','options','method','scan','positions'}; % a valid aMALDIms file must have these fields        
            object.require(fieldsRequired);
            object.loaded = true;
        end
                
        function print(object) 
        % display structure of the file            
            h5disp(object.filename,'/','simple') 
        end
        
        function result = tree(object) 
        % get file structure tree
            result = h5info(object.filename); 
        end 
        
    end
    
    methods % for loading data
                
        function [info, found] = get_info(object,location)
        % loads info about a location in the HDF5 file (a wrapper for the h5info function, that handles errors)
            % initialize output            
            info = struct();
            found = false;   
            % try to get 
            try
                info = h5info(object.filename,location); 
                found = true;
            catch err
                object.handle_error(err,location);
            end
        end
        
        function [data, found] = get_data(object,location,varargin)
        % loads a dataset (fully or partially) from a location in the HDF5 file (a wrapper for the h5read function, that handles errors)
            % initialize output            
            data = [];
            found = false;            
            % try to get 
            try            
                Narg = length(varargin);
                if  Narg == 2 % read in a portion of data
                    start = varargin{1};
                    count = varargin{2};
                    data = h5read(object.filename,location,start,count);
                elseif Narg == 3 % read in portion of data with stride
                    start = varargin{1};
                    count = varargin{2};
                    stride = varargin{3};                    
                    data = h5read(object.filename,location,start,count,stride);
                else % read all data
                    data = h5read(object.filename,location);
                end
                found = true;
            catch err
                object.handle_error(err,location);
            end
        end        
        
        function [data, found] = get_all_data(object,location) 
        % loads all datasets (in their entirety) under a group location in the HDF5 file 
            % initialize output
            data = struct(); 
            % get group info
            [info, found] = object.get_info(location); 
            if found
                % get scan datasets
                datasets = info.Datasets;
                for m = 1:length(datasets)
                    dataset_name = datasets(m).Name;
                    dataset_location = [location '/' dataset_name];
                    [dataset_values, dataset_found] = object.get_data(dataset_location);
                    if dataset_found
                        data.(dataset_name) = dataset_values;
                    end
                end
            end
        end
                
        function [result, found, info] = get_attributes(object,location) 
        % get all availabe attributes of a group specified by its location
            % initialize output
            result = struct();   
            % get group info
            [info, found] = object.get_info(location); 
            if found
                % get attributes                
                attributes = info.Attributes;
                for m = 1:length(attributes)
                    attribute_name = object.parse_name(attributes(m).Name);
                    attribute_value = object.parse_value(attributes(m).Value);                
                    result.(attribute_name) = attribute_value;                
                end
            end
        end        
        
        function [result, found] = load_file_info(object)
        % grab all available info about the file
            [result, found] = object.get_attributes(object.location.Root); 
        end
        
        function [result, found] = load_options(object) 
        % grab all availabe info about conversion options used to generated aMALDIms file
            [result, found] = object.get_attributes(object.location.Option);        
        end
        
        function [result, found] = load_raw_info(object) 
        % grab all availabe info about raw (time domain) measurement
            [result, found] = object.get_attributes(object.location.Raw);               
        end     
        
        function [result, found] = load_spectra_info(object) 
        % grab all availabe info about spectra (mass/charge domain) measurement
            [result, found, info] = object.get_attributes(object.location.Spectra); 
            if found
                result.parameters = orderfields(object.get_attributes(object.location.SpectraPar)); % sort fields by name                 
                infoMZ = h5info(object.filename, object.location.SpectraMZ); % 
                result.has.every_mz = length(infoMZ.Dataspace.Size)>1;
                [stats, result.has.statistics] = object.get_all_data(object.location.SpectraSta);  
                if result.has.statistics
                    result.statistics = stats;
                end
                [background, result.has.background] = object.get_all_data(object.location.SpectraBac);  
                if result.has.background
                    result.background = background;                    
                end
            end
        end     
        
        function [result, found] = load_respectra_info(object) 
        % grab all availabe info about spectra (mass/charge domain) measurement
            [result, found] = object.get_attributes(object.location.ReSpectra);   
            if found
                result.parameters = orderfields(object.get_attributes(object.location.ReSpectraPar)); % sort fields by name                           
                [stats, result.has.statistics] = object.get_all_data(object.location.ReSpectraSta);  
                if result.has.statistics
                    result.statistics = stats;
                end
                [background, result.has.background] = object.get_all_data(object.location.ReSpectraBac);  
                if result.has.background
                    result.background = background;                    
                end          
                datasets = {object.get_info(object.location.ReSpectra).Datasets.Name};
                result.has.signal_sum = any(contains(datasets,'signal_sum'));
                %[~,result.has.signal_sum] = object.get_info(object.location.ReSpectraSum)
            end
        end  
        
        function [result, found] = load_method_info(object) 
        % grab all available info about measurement methods
            % initialize output
            result = struct();                        
            % get group info
            [info, found] = object.get_info(object.location.Method); 
            if found
                % get attributes of method
                attributes = info.Attributes;            
                for m = 1:length(attributes)
                    result.(attributes(m).Name) = attributes(m).Value;
                end            
                % get attributes of each subgroup of methods
                subgroups = info.Groups;
                for m = 1:length(subgroups)
                    subgroup_attr = subgroups(m).Attributes;
                    subgroup_name = subgroups(m).Name;
                    subgroup_name = strsplit(subgroup_name,'/'); % parse and remove slash
                    subgroup_name = subgroup_name{end};
                    for n = 1:length(subgroup_attr)
                        subgroup_attr_name = subgroup_attr(n).Name;
                        subgroup_attr_value = object.parse_value(subgroup_attr(n).Value);
                        result.(subgroup_name).(subgroup_attr_name) = subgroup_attr_value;
                    end
                end
            end
        end
        
        function [result, found] = load_scan_log(object) 
        % grab all available info from the scan log
            [result, found] = get_all_data(object,object.location.Scan);        
        end
                
        
        function [positions, found, indDuplicates, equivalent] = load_positions(object) 
        % load position vector (while checking if the various positions vectors are equivalent)    
            % initialize outputs
            positions = {};   
            indDuplicates = []; % index of duplicate positions
            % get the position vectors from every possible location
            vfound = [];
            vlabel = {};
            vector = {};
            n = 1;
            if object.has.scan
                vlabel{n} = 'scan'; 
                [vector{n}, vfound(n)] = object.get_data(object.location.ScanPos); 
                n = n + 1; 
            end
            if object.has.spectra
                vlabel{n} = 'spectra'; 
                [vector{n}, vfound(n)] = object.get_data(object.location.SpectraPos); 
                n = n + 1; 
            end           
            if object.has.raw
                vlabel{n} = 'raw'; 
                [vector{n}, vfound(n)] = object.get_data(object.location.RawPos); 
                n = n + 1; 
            end        
            if object.has.respectra
                vlabel{n} = 'respectra'; 
                [vector{n}, vfound(n)] = object.get_data(object.location.ReSpectraPos); 
                n = n + 1; 
            end   
            % check that at least one of the position vectors was present in the file
            num_found = sum(vfound);
            found = num_found>0;
            index = find(vfound>0);   
            equivalent = true; % initially assume all vectors are equivalent
            if found
                % choose a default position vector
                positions = vector{index(1)};  
                for n=2:length(index) 
                    equivalent = equivalent && isequal(positions,vector{index(n)});
                end
                % if the vectors differ, then return all
                if not(equivalent)
                    positions = struct();
                    for n=1:length(index)
                        positions.(vlabel{n}) = vector{index(n)}; 
                    end
                    object.warn('Positions vectors not the same. Returned all vectors instead of single vector.')                
                end  
                % check that each position is unique 
                if equivalent
                    positionsUnique = {};
                    for n = 1:length(positions)
                        pos = positions{n};
                        if sum(ismember(positionsUnique,pos))>0
                            indDuplicates = [indDuplicates, n];
                            object.warn(['A repeat at position ' pos ' was found.'])                                        
                        else
                            positionsUnique = [positionsUnique, {pos}];
                        end
                    end                
                end                        
            end
        end
        
        function [plate] = determine_plate(object)
        % determine the type of plate measured (the layout of MALDI grid)
            if not(isempty(object.positions))
                letter_pattern = @(str) num2str(isstrprop(str,'alpha'));
                pos = object.positions{1};
                pos_pattern = letter_pattern(pos);
                pos_count = length(pos);
                switch pos_pattern
                    case letter_pattern('X00Y00')
                        plate.layout = '1536';
                        plate.dimensions = [32,48];                        
                    case letter_pattern('A00')
                        plate.layout = '384';
                        plate.dimensions = [16,24];
                    otherwise
                        plate.layout = 'unknown';
                        plate.dimensions = [1,1]*ceil(sqrt(pos_count));                        
                end
            end
        end     
        
        function [positionsInt] = positionStr2Int(object,positionsStr)
        % convert positions in string form ('A1' or 'X01Y01') to integer form ([1,1])
            if strcmp(object.plate.layout,'384')
            	positionsInt = object.positionAlphaStr2Int(positionsStr);
            elseif strcmp(object.plate.layout,'1536')
                positionsInt = object.positionXYStr2Int(positionsStr);
            else % unknown
                positionsInt = ones(length(positionsStr),2);                
                positionsInt(:,2) = 1:length(positionsStr);
            end
        end
        
        function [positionsStr] = positionInt2Str(object,positionsInt)
        % convert positions in integer form ([1,1]) to string form ('A1' or 'X01Y01')       
            if strcmp(object.plate.layout,'384')
            	positionsStr = object.positionAlphaInt2Str(positionsInt);
            elseif strcmp(object.plate.layout,'1536')
                positionsStr = object.positionXYInt2Str(positionsInt);
            else % unknown
                positionsStr = object.positionXYInt2Str(positionsInt);                
            end
        end
        
        function [positionsInd] = position2Ind(object,positions,varargin)
        % convert positions in integer or string form to a position index
        % (position can be in specified as: 'A1', 'X01Y01', or [1,1])
        % (multiple inputs can be given as: {'A1','B1'}, {'X01Y01','X02Y01'}, or [1,1; 2,1])
            if length(varargin)>0
                allpos = varargin{1};
            else
                allpos = object.positionStr2Int(object.positions);
            end
            K = 1;
            is_list = false;
            if isnumeric(positions) 
                K = size(positions,1);
                is_list = K>1;                
            elseif iscell(positions)
                K = length(positions);
                is_list = true;                                
            end
            if is_list
                if isnumeric(positions) 
                    positions = num2cell(positions,2);
                end
                positionsInd = zeros(K,1);
                for k = 1:K
                	positionsInd(k) = object.position2Ind(positions{k},allpos);
                end
            else                
                pos = positions;
                if ischar(pos)
                    pos = object.positionStr2Int(pos);
                end
                [ismemb] = ismember(allpos,pos,'rows');
                positionsInd = find(ismemb,1);
                if isempty(positionsInd)
                    positionsInd = nan;
                end
            end
        end                
        
    end
    
    methods % for loading specific datasets
                                           
        function [signal, time, found] = get_data_raw(object,starting_index,consecutive_count)
        % loads the raw data for a single position (specified by index)                                  
            found = object.has.raw;
            if found 
                start = [1,starting_index]; % (x,y) point to start counting start
                count = [object.num.raw,consecutive_count]; % (x,y) number of points to read in
                time = object.get_data(object.location.RawTim,start(1),count(1));
                signal = object.get_data(object.location.RawSig,start,count);
            else
                time = [];
                signal = [];
            end
        end                
        
        function [signal, mass, found] = get_data_spectra(object,starting_index,consecutive_count)
        % loads the spectra data for a single position (specified by index)                                  
            found = object.has.spectra;
            if found 
                start = [1,starting_index]; % (x,y) point to start counting start
                count = [object.num.spectra,consecutive_count]; % (x,y) number of points to read in
                if object.spectra.has.every_mz
                    mass = object.get_data(object.location.SpectraMZ,start,count);
                else
                    mass = object.get_data(object.location.SpectraMZ,start(1),count(1));
                end
                signal = object.get_data(object.location.SpectraSig,start,count);
            else
                mass = [];
                signal = [];
            end
        end    
        
        function [signal, mass, found] = get_data_respectra(object,starting_index,consecutive_count)
        % loads the respectra data for a single position (specified by index)                                  
            found = object.has.respectra;
            if found 
                start = [1,starting_index]; % (x,y) point to start counting start
                count = [object.num.respectra,consecutive_count]; % (x,y) number of points to read in
                mass = object.get_data(object.location.ReSpectraMZ,start(1),count(1));
                signal = object.get_data(object.location.ReSpectraSig,start,count);
            else
                mass = [];
                signal = [];
            end
        end           
        
        function [signal_sum, mass, found] = get_data_respectra_sum(object)
        % loads the respectra data for a single position (specified by index)                                  
            found = object.has.respectra && object.respectra.has.signal_sum;
            if found                 
                start = [1,1]; % (x,y) point to start counting start
                count = [object.num.respectra,1]; % (x,y) number of points to read in                
                mass = object.get_data(object.location.ReSpectraMZ,start(1),count(1));
                signal_sum = object.get_data(object.location.ReSpectraSum,start,count);
            else
                mass = [];
                signal_sum = [];
            end
        end           
        
        
        function [time, found] = get_data_raw_time(object)
        % loads the raw data time vector                               
            found = object.has.raw;
            if found 
                time = object.get_data(object.location.RawTim,1,object.num.raw);
            else
                time = [];
            end
        end   
        
        function [signal, found] = get_data_raw_signal(object,index)
        % loads the raw data signal vector at a position (specified by index)                                
            found = object.has.raw;
            if found 
                signal = object.get_data(object.location.RawSig,[1,index],[object.num.raw,1]);
            else
                signal = [];
            end
        end                        
        
        function [mass, found] = get_data_spectra_mass(object,index)
        % loads the spectra data mass vector at a position (specified by index)                                 
            found = object.has.spectra;
            if found 
                if object.spectra.has.every_mz
                    mass = object.get_data(object.location.SpectraMZ,[1,index],[object.num.spectra,1]);
                else
                    mass = object.get_data(object.location.SpectraMZ,1,object.num.spectra);
                end
            else
                mass = [];
            end
        end   
        
        function [mass, found] = get_data_respectra_mass(object)
        % loads the respectra data mass vector                               
            found = object.has.respectra;
            if found 
                mass = object.get_data(object.location.ReSpectraMZ,1,object.num.respectra);
            else
                mass = [];
            end
        end       
        
        function [signal, found] = get_data_spectra_signal(object,index)
        % loads the spectra data signal vector at a position (specified by index)                                
            found = object.has.spectra;
            if found 
                signal = object.get_data(object.location.SpectraSig,[1,index],[object.num.spectra,1]);
            else
                signal = [];
            end
        end            
                        
        function [signal, found] = get_data_respectra_signal(object,index)
        % loads the respectra data signal vector at a position (specified by index)                                
            found = object.has.respectra;
            if found 
                signal = object.get_data(object.location.ReSpectraSig,[1,index],[object.num.respectra,1]);
            else
                signal = [];
            end
        end                                                    
        
    end
        
    methods % for visualizing data
        
        
        function [fig] = plateview(object,name,variable,varargin)
        % function that plots the vector variable (of the same length and order as positions) 
        % as the color and height of spots arranged on the virtual MALDI plate        
            % options 
            positions = object.positions;
            xyFlipped = false;
            logscale = false;
            linewidth = 0.25;
            marker = 'o';
            markersize = 55;   
            markeredgecolor = 'none';
            zLimits = [];
            Nvarargin = length(varargin);
            custom_cmap = false;
            new_figure = true;
            for n = 1:Nvarargin
                if Nvarargin>=n+1
                    if ischar(varargin{n}) || isscalar(varargin{n})
                        key = lower(varargin{n});
                        value = varargin{n+1};
                        switch key
                            case 'xyflipped'
                                xyFlipped = value;
                            case 'logscale'
                                logscale = value;
                            case 'linewidth'
                                linewidth = value;
                            case 'marker'
                                marker = value;
                            case 'markersize'
                                markersize = value;
                            case 'markeredgecolor'
                                markeredgecolor = value;
                            case 'zlimits'
                                zLimits = value;
                            case 'colormap'
                                cmap = value;
                                custom_cmap = true;
                            case 'positions'
                                positions = value;
                            case 'newfigure'
                                new_figure = value;                                 
                        end
                    end
                end
            end
            % get plate info
            plate_type = object.plate.layout;
            plate_dim = object.plate.dimensions;
            Npos = length(positions);
            % format data as grid
            [XGrid, YGrid, ZGrid] = object.plate_vector_to_grid(variable,'xyflipped',xyFlipped,'positions',positions);
            % plot the resulting plate (z dimension and color correspond to the variable values)
            if new_figure
                fig = figure('color','w','name',name);
            end
            vectorize = @(matrix) reshape(matrix,[1,numel(matrix)]);
            X = vectorize(XGrid);
            Y = vectorize(YGrid);
            Z = vectorize(ZGrid);
            C = vectorize(ZGrid);
            %whos X Y Z C
            scatter3(X,Y,Z,markersize,C, 'filled',marker,'markeredgecolor',markeredgecolor,'linewidth',linewidth)
            if custom_cmap
                cmap = colormap(gca,cmap);                
            else
                cmap = colormap(gca,'parula'); % 'summer'
                cmap = colormap(gca,flipud(cmap)); % reverse the coloring
            end            
            view(0,-90)
            xlimits = [1,plate_dim(2)];
            ylimits = [1,plate_dim(1)];
            xticknames = object.positionInt2Str([xlimits(1) ylimits(1); xlimits(2) ylimits(1)]);
            yticknames = object.positionInt2Str([xlimits(1) ylimits(1); xlimits(1) ylimits(2)]);
            set(gca,...
                'XLim',[xlimits(1)-1,xlimits(2)+1],... 
                'YLim',[ylimits(1)-1,ylimits(2)+1],...                
                'XTick',xlimits,... 
                'YTick',ylimits,...
                'XTickLabel',xticknames,... 
                'YTickLabel',yticknames);                      
            zlabel(name)
            if isempty(zLimits)
                if Npos>1 && length(unique(variable(~isnan(variable))))>1
                    zLimits = [min(variable),max(variable)];
                    zlim(zLimits)  
                    caxis(zLimits)      
                end
            else
                if zLimits(1)==-inf
                    zLimits(1) = min(variable);
                end
                if zLimits(2)==inf
                    zLimits(2) = max(variable);
                end
                zlim(zLimits)                            
                caxis(zLimits)
            end            
            cb = colorbar('location','southoutside');
            cb.Label.String = name;
            set(gca,'fontsize',12)
            if logscale
                set(gca,'zscale','log');
                if isprop(gca,'colorscale') % colorscale property is only supported in newer MATLAB versions
                    set(gca,'colorscale','log');
                end
            end
            offwhite = [1,1,1]*0.975;
            set(gca,'Color',offwhite);
            box on; grid off;
            title([upper(plate_type) ' MALDI Plate'],'fontweight','normal');   
        end
        
        function [XGrid, YGrid, ZGrid] = plate_vector_to_grid(object,variable,varargin)
        % convert a vector to a matrix (if optional argument 'positions' is not provided, the vector will assume be aligned with and of the same length as object.positions)
            % options
            positions = object.positions;
            xyFlipped = false;
            Nvarargin = length(varargin);
            for n = 1:Nvarargin
                if Nvarargin>=n+1
                    if ischar(varargin{n}) || isscalar(varargin{n})
                        key = lower(varargin{n});
                        value = varargin{n+1};
                        switch key
                            case 'xyflipped'
                                xyFlipped = value;
                            case 'positions'
                                positions = value; 
                        end
                    end
                end
            end
            % format positions
            Npos = length(positions);
            Nvar = length(variable);
            if Npos~=Nvar
                error('aMALDIms:lengthMismatch','Error: Input vector must have the same number of elements as positions')    
            end
            positions = object.positionStr2Int(positions); % get position numbers
            if xyFlipped
                positions = fliplr(positions);
            end
            % remove duplicates (first occurence kept)
            ind_dups = object.duplicates;
            positions(ind_dups,:)=[];
            variable(ind_dups)=[];
            % get plate dimensions 
            plate_dim = object.plate.dimensions;
            % generate 2D grid
            [XGrid, YGrid] = meshgrid(1:plate_dim(2),1:plate_dim(1));
            XGrid = XGrid';
            YGrid = YGrid';
            % get the variable value at each grid point
            ZGrid = nan(size(XGrid));
            for n = 1:Npos
                pos = positions(n,:);
                ZGrid(pos(1),pos(2)) = variable(n);
            end                
        end        
        
    end
    
    methods % for analyzing data
                
        function [intensities,window,windowIntensities,windowMasses] = mass_shift_for_intensity(object,mass_to_query,varargin)
        % function to shift for the intensities in a resampled spectra matrix (at all positions) that occur over a small, contiguous range of masses
            % options 
            mass_vector = []; % vector of masses (if left empty, will be loaded from file but to save on time it can be preloaded and provided as an input to this function)
            window_halfwidth = 0.000; % 0.005; % width of the m/z window [Da/C]
            window_function = @(x) max(x,[],1); % function to perform on the m/z window
            Nvarargin = length(varargin);
            for n = 1:Nvarargin
                if Nvarargin>=n+1
                    if ischar(varargin{n}) || isscalar(varargin{n})
                        key = lower(varargin{n});
                        value = varargin{n+1};
                        switch key
                            case 'windowhalfwidth'
                                window_halfwidth = value;
                            case 'windowfunction'
                                window_function = value;
                            case 'massvector'
                                mass_vector = value;  
                        end
                    end
                end
            end
            intensities = [];
            window = struct();            
            if object.has.respectra 
                spectra_matrix_location = object.location.ReSpectraSig;
                if isempty(mass_vector)
                    mass_vector = object.get_data_respectra_mass();  
                end
                % find index of nearest (center) mass value in mass vector as well as window upper and lower bound indices            
                [~, window.index.low] = min(abs(mass_vector-(mass_to_query-window_halfwidth)));
                [~, window.index.center] = min(abs(mass_vector-(mass_to_query)));    
                [~, window.index.high] = min(abs(mass_vector-(mass_to_query+window_halfwidth)));
                window.points = window.index.high - window.index.low + 1;            
                % get the masses at each of the above found indices 
                window.mass.center = mass_vector(window.index.center);
                window.mass.low = mass_vector(window.index.low);
                window.mass.high = mass_vector(window.index.high);    
                window.mass_delta = [window.mass.low, window.mass.high]-window.mass.center;
                windowMasses = mass_vector(window.index.low:window.index.high);
                % retreive the data block containing the desired mass slices
                start = [window.index.low,1]; % (x,y) point to start counting start
                count = [window.points,object.num.positions]; % (x,y) number of points to read in
                windowIntensities = object.get_data(spectra_matrix_location,start,count);
                intensities = window_function(windowIntensities);   
            else
                object.warn('Function requires ReSpectra to be present.');
            end
        end
                        
        
    end
            
    methods % for handling errors and command line notifications

        function notify(object,message)
        % function to display a command line message
            if object.verbose
                disp(['aMALDIms: ' message]);              
            end
        end        
        
        function warn(object,message)
        % function to display warnings
            if object.verbose
                object.notify(['Warning: ' message]); % disp(['Warning: ' message]); % warning(message);                
            end
        end
        
        function handle_error(object,err,location)
        % manually handle or deal with errors encountered during hdf5 file reading
           % handle if location not found
           id = err.identifier;
           message = err.message;
           unable_to_find_location = strcmp(id,'MATLAB:imagesci:h5info:unableToFind') || sum(object.contains(message,{'component not found','traversal operator failed','name doesn''t exist'}));
           if unable_to_find_location
              object.warn(['Location "' location '" not found in HDF5 file.']) 
           else
              rethrow(err)
           end            
        end        
        
        function require(object,fields)
        % throw an error if specified fields are not present
            % if given single input string wrap it in a cell
            if not(iscell(fields))
                fields = {fields};
            end
            % check specified fields
            present = false(size(fields));
            for n = 1:length(fields)
                field = fields{n};
                if isfield(object.has,field)
                    present(n) = object.has.(field);
                end
            end
            % throw error if missing
            if not(all(present))
                fields_missing = fields(not(present));
                missing = [fields_missing{1}];
                for n = 2:length(fields_missing)
                    missing = [missing ', ' fields_missing{n}];
                end
                err.message = ['The following required fields were not found: ' missing];
                err.identifier = 'aMALDIms:requiredFieldsNotFound';
                error(err)                
            end
        end
        
    end
    
    methods (Static) % static 
                
        function name = parse_name(name)
        % parse attribute name (which can contain spaces)
            % replace space with underscore
            space = ' ';
            underscore = '_';
            if aMALDIms.contains(name,space) 
                name = strrep(name,space,underscore);
            end
            % remove other invalid characters
            name(regexp(name,'[ \] \[ , \/ \+ \- \\ ^ \( \) \. \# =]'))=[];
            % clip prefix if not alphebetical
            alphaInd = find(isstrprop(name,'alpha'));
            name = name(alphaInd(1):end);
        end
        
        function value = parse_value(value)
        % parse attribute value (which could be one of many types)
            % convert character cell to character array
            if iscell(value)
                value = char(value);
            end
            % convert logical string to logical
            %{
            if isequal(upper(value),'TRUE')
                value = true;
            elseif isequal(upper(value),'FALSE')
                value = false;                
            end
            %}
        end
        
        function [xy] = positionAlphaStr2Int(str)
        % convert a position from alphabet string ('A1') to [x,y] integers ([1,1]) (accepts multiple positions formatted as an Nx1 cell, such as: {'A1','C10'})
            if iscell(str) || ischar(str)                
                if iscell(str)
                    K = length(str);
                    xy = nan(K,2);
                    for k = 1:K
                        xy(k,:) = aMALDIms.positionAlphaStr2Int(str{k});
                    end
                else                
                    is_let = isstrprop(str,'alpha');
                    is_num = isstrprop(str,'digit');
                    str_let = str(is_let);
                    str_let = upper(str_let); % convert to all capitals
                    len_let = length(str_let);
                    str_num = str(is_num);
                    len_num = length(str_num);
                    if len_let>0 && len_num>0
                        num_let = aMALDIms.letters2number(str_let);
                        num_num = str2double(str_num);
                        xy = [num_let,num_num];
                    else
                        xy = [1,1]*nan;
                    end   
                end
            elseif isnumeric(str)
                xy = str; % assume is already in desired form
            else
                xy = []; % otherwise return empty                
            end
        end
        
        function [xy] = positionXYStr2Int(str)
        % convert a position from XY string ('X01Y01') to [x,y] integers ([1,1]) (accepts multiple positions formatted as an Nx1 cell, such as: {'X01Y01','X32Y47'})
            if iscell(str) || ischar(str)        
                if iscell(str)
                    K = length(str);
                    xy = nan(K,2);
                    for k = 1:K
                        xy(k,:) = aMALDIms.positionXYStr2Int(str{k});
                    end
                else        
                    is_let = isstrprop(str,'alpha');
                    is_num = isstrprop(str,'digit');
                    str_let = str(is_let);
                    str_let = upper(str_let); % convert to all capitals                            
                    if sum(is_let)==2 && aMALDIms.contains(str_let,'X') && aMALDIms.contains(str_let,'Y')
                        ind_num = find(is_num);
                        ind_ind_split = find([diff([ind_num(1)-1,ind_num])]>1);
                        range_num_1 = ind_num(1:(ind_ind_split-1));
                        range_num_2 = ind_num(ind_ind_split:end);   
                        num_1 = str2double(str(range_num_1));
                        num_2 = str2double(str(range_num_2));
                        if strcmp(str_let,'XY') % always put integer for X first
                            xy = [num_1,num_2];
                        else % strcmp(str_let,'YX')
                            xy = [num_2,num_1];
                        end
                    else
                        xy = [1,1]*nan;
                    end  
                end
            elseif isnumeric(str)
                xy = str; % assume is already in desired form
            else
                xy = []; % otherwise return empty                
            end
        end
        
        function [str] = positionAlphaInt2Str(xy)
        % convert a position from [x,y] integers ([1,1]) to alphabet string ('A1') (accepts multiple positions formatted as an Nx2 matrix, such as: [[1,2];[2,2];[1,1];])
            if isnumeric(xy)
                K = size(xy,1);
                if K > 1
                    str = cell(K,1);
                    for k = 1:K
                        str{k} = aMALDIms.positionAlphaInt2Str(xy(k,:));
                    end
                else    
                    num_let = xy(1);
                    str_let = aMALDIms.number2letters(num_let);
                    num_num = xy(2);
                    str_num = sprintf('%01d',num_num);
                    str = [str_let,str_num];
                end  
            elseif ischar(xy)
                str = xy; % assume is already in desired form
            else
                str = ''; % otherwise return empty                                
            end
        end        
        
        function [str] = positionXYInt2Str(xy)
        % convert a position from [x,y] integers ([1,1]) to XY string ('X01Y01') (accepts multiple positions formatted as an Nx2 matrix, such as: [[1,2];[2,2];[1,1];])
            if isnumeric(xy)
                K = size(xy,1);
                if K > 1
                    str = cell(K,1);
                    for k = 1:K
                        str{k} = aMALDIms.positionXYInt2Str(xy(k,:));
                    end
                else     
                    num_1 = xy(1);
                    num_2 = xy(2);
                    str = ['X',sprintf('%02d',num_1),'Y',sprintf('%02d',num_2)];
                end
            elseif ischar(xy)
                str = xy; % assume is already in desired form
            else
                str = ''; % otherwise return empty                                
            end
        end   
        
        function [number] = letters2number(letters)
        % get the number corresponding to the given letters    
            % set base and digit symbols
            symbols = char(double('A'):double('Z'));  
            base = length(symbols);            
            % convert the letters to the number
            N = length(letters);
            number = 0;
            for n = 1:N
                digit = double(letters(n)) - double('A') + 1; % the plus one shifts the digit value back (to start at 0 instead of 1) 
                number = number + base^(N-n)*(digit);
            end                        
            
        end
        
        function [letters] = number2letters(number)
        % get the letters corresponding to the given number    
            % set base and digit symbols
            symbols = char(double('A'):double('Z'));  
            base = length(symbols);
            % determine the expected number of base digits 
            N = max(1,floor(log(number)/log(base))+1); 
            % convert the number to digits in the specified base
            digits = zeros(1,N);
            quotient = number;
            for n=N:-1:1
            digits(n) = rem(quotient,base);      
            quotient = floor(quotient/base);
            end
            % shift the digits so that symbols can be assigned from a digit value of 1 to 26 (as opposed to 0 to 25)
            for n=N:-1:2
            if (digits(n) <= 0)
              digits(n) = base + digits(n);
              digits(n-1) = digits(n-1) - 1;
            end
            end
            % remove zero valued most significant digit (if present)
            ind = find(digits,1,'first'); 
            digits(1:(ind-1)) = [];
            % convert digits to letters
            letters = symbols(digits);            
        end
        
        function [occurrence] = contains(str,pat)
        % determine if a pattern (sub-string) is in string and returns the number of times it occurs (accepts multiple patterns provided in a cell) (to be used in lieu of the "contains" function found in newer MATLAB versions)            
            if not(iscell(pat))
                pat = {pat};
            end
            N = length(pat);
            occurrence = zeros([1,N]);
            for n = 1:N
                pattern = pat{n};
                occurrence(n) = length(strfind(str,pattern));
            end
        end
        
               
    end
    
    
    methods % file conversion (from aMALDIms to aMALDImsRed)
        
        function [compressfiles, duration, success] = compress_multiple_respectra_by_threshold(object,originalfiles,varargin)
        % compress multiple aMALDIms files
            duration = cell(size(originalfiles));
            compressfiles = cell(size(originalfiles));
            success = false(size(originalfiles));
            for n = 1:length(originalfiles)
                amaldims = aMALDIms(originalfiles{n});
                [compressfiles{n},duration{n},success(n)] = object.compress_respectra_by_threshold(amaldims,varargin{:});
                if success(n)
                    object.notify(['Compressed file "' originalfiles{n} '" (' num2str(n) '/' num2str(length(originalfiles)) ')'])
                end
            end  
        end
        
        function [filename,duration,success] = compress_respectra_by_threshold(object,amaldims,varargin)
        % function to compress and export spectra from an aMALDIms (HDF5) file to a reduced size MAT-file
            filename = ''; duration = struct(); success = false;
            if amaldims.loaded 
                if amaldims.has.respectra
                    object.notify('Compression started.')
                    initial_time = tic;
                    
                    % Methods
                    about = struct();
                    if amaldims.has.method
                        about.acquisition = amaldims.method;
                    else
                        about.acquisition = struct();
                    end
                    about.resampling = struct();
                    about.resampling.Info = amaldims.info;                
                    about.resampling.Interpolant = amaldims.respectra.interpolant;
                    about.resampling.Points = amaldims.respectra.points;
                    about.resampling.Count = amaldims.respectra.count;
                    about.acquisition.Parameters = amaldims.respectra.parameters;
                    
                    % Options
                    about.compression.Algorithm = 'SNRThreshold'; % compression method {'SNRThreshold','PeakThreshold'}
                    about.compression.TukeyFactor = 3; % Tukeys's fence hinge step factor, rule of thumb: set to 1.5 (inner fence) for outliers and 3 (outer fence) for extreme outliers
                    about.compression.HeightThreshold = 30; %24; %12; 9; % SNR cutoff in number of background standard deviations
                    about.compression.HeightMode = 'SignalToNoiseRatio'; % {'Intensity','SignalToNoiseRatio'}
                    about.peakfinding.FindPeaks = true;
                    about.peakfinding.WidthReference = 'HalfProm'; % {'HalfProm','HalfHeight'}  
                    about.peakfinding.MinimumSeparation = 0.005; % [Da] 
                    about.peakfinding.HeightReference = 'Prominence'; %  {'Prominence','Height'}
                    about.peakfinding.HeightThreshold = 12;   
                    num_status_update = 50; % 100; % status update rate (the number of spectra to process between giving a status update) (set to "0" for no updates) 
                    about.raw.Include = false; % include raw time data if available
                    if not(isempty(varargin))
                        for n=1:length(varargin)
                            switch lower(varargin{n})
                                case 'algorithm'
                                	about.compression.Algorithm = varargin{n+1};                            
                                case 'threshold'
                                	about.compression.HeightThreshold = varargin{n+1};
                                case 'normalizeintensity'
                                    normalize_intensity = varargin{n+1};
                                    if normalize_intensity
                                        about.compression.HeightMode = 'SignalToNoiseRatio';  
                                    else
                                        about.compression.HeightMode = 'Intensity';                                      
                                    end
                                case 'findpeaks'
                                    about.peakfinding.FindPeaks = varargin{n+1};
                                case 'peakwidthreference'
                                    about.peakfinding.WidthReference = varargin{n+1};                                
                                case 'peakminimumseparation'
                                    about.peakfinding.MinimumSeparation = varargin{n+1};                                
                                case 'peakheightreference'
                                    about.peakfinding.HeightReference = varargin{n+1};                                
                                case 'peakheightthreshold'
                                    about.peakfinding.HeightThreshold = varargin{n+1};    
                                case 'statusupdaterate'
                                    num_status_update = varargin{n+1};     
                                case 'includerawdata'
                                    about.raw.Include = varargin{n+1};                                           
                            end
                        end
                    end
                    if object.verbose
                        object.notify(['Configured to use the "' about.compression.Algorithm '" algorithm.'])                    
                        object.notify(['Threshold set to ' num2str(about.compression.HeightThreshold) ' background standard deviations.'])
                    end                                                
                    
                    % Compression Pass 1
                    if object.verbose
                        object.notify('Compression Pass 1: collecting statistics to determine features to keep.') %  based on a SNR threshold
                    end
                    summary = table(); % compression info
                    summary.Position = amaldims.positions;
                    grid = amaldims.positionStr2Int(amaldims.positions);
                    summary.GridXPosition = grid(:,1);
                    summary.GridYPosition = grid(:,2);
                    clear grid
                    if amaldims.has.scan
                        summary.ElapsedTime = amaldims.scan.elapsedMinutes;
                        summary.TIC = amaldims.scan.TIC;
                    end
                    blank = zeros(amaldims.num.positions,1);
                    summary.Maximum = blank;
                    summary.Mininum = blank;
                    summary.Mean = blank;
                    summary.StdDev = blank;
                    summary.Mode = blank;
                    summary.Median = blank;
                    summary.Quartile1 = blank;
                    summary.Quartile3 = blank;
                    summary.Sum = blank;
                    summary.Area = blank;
                    summary.TukeyLower = blank;
                    summary.TukeyUpper = blank;
                    summary.BackgroundMean = blank;
                    summary.BackgroundStdDev = blank;
                    summary.BackgroundLinearSlope = blank;
                    summary.BackgroundLinearOffset = blank;                
                    duration.Pass1 = blank;
                    duration.Pass2 = blank;
                    duration.PeakFinding = blank;
                    if about.peakfinding.FindPeaks
                        peaks = cell(amaldims.num.positions,1);
                        summary.PeakCount = blank;
                    end
                    clear blank                
                    updateMinutes = @(duration) mean(duration(duration>0))/60;
                    indent = '  ';
                    start_index = 1;
                    spectra_count = 1;
                    [~, mass_to_charge, ~] = amaldims.get_data_respectra(start_index,spectra_count);     
                    keep = false(1,amaldims.num.respectra);
                    for n = 1:amaldims.num.positions
                        % start timer
                        timer = tic;
                        % load resampled spectra data
                        intensity = amaldims.get_data_respectra(n,spectra_count); 
                        % get statistics
                        mask = not(isnan(intensity));
                        summary.Maximum(n) = max(intensity(mask));
                        summary.Mininum(n) = min(intensity(mask));
                        summary.Mean(n) = mean(intensity(mask));   
                        summary.StdDev(n) = std(intensity(mask));   
                        summary.Mode(n) = mode(intensity(mask));
                        summary.Median(n) = median(intensity(mask));        
                        summary.Quartile1(n) = prctile(intensity(mask),25); 
                        summary.Quartile3(n) = prctile(intensity(mask),75); 
                        summary.Sum(n) = sum(intensity(mask));  
                        summary.Area(n) = trapz(mass_to_charge(mask),intensity(mask)); 
                        clear mask
                        % perform background estimation (using Tukey's fences)  
                        summary.TukeyLower(n) = summary.Quartile1(n)-about.compression.TukeyFactor*(summary.Quartile3(n)-summary.Quartile1(n));    
                        summary.TukeyUpper(n) = summary.Quartile3(n)+about.compression.TukeyFactor*(summary.Quartile3(n)-summary.Quartile1(n));
                        fenced = intensity < summary.TukeyUpper(n); % & intensity > fence.bounds(1) % spectra are onesided so no need for first limit
                        summary.BackgroundMean(n) = mean(intensity(fenced));
                        summary.BackgroundStdDev(n) = std(intensity(fenced));
                        linear_fit_coeff = polyfit(mass_to_charge(fenced),intensity(fenced),1);
                        summary.BackgroundLinearSlope(n) = linear_fit_coeff(1);
                        summary.BackgroundLinearOffset(n) = linear_fit_coeff(2);
                        clear fenced linear_fit_coeff
                        % perform peak finding
                        if about.peakfinding.FindPeaks
                            timerpeak = tic;
                            found = table();                        
                            [found.height,found.mass_to_charge,found.width,found.prominence] = findpeaks((intensity - summary.BackgroundMean(n))/summary.BackgroundStdDev(n), mass_to_charge,...
                                ['MinPeak' about.peakfinding.HeightReference],about.peakfinding.HeightThreshold,...
                                'MinPeakDistance',about.peakfinding.MinimumSeparation,...
                                'WidthReference',about.peakfinding.WidthReference); % get height, locations, widths, and prominences of peaks            
                            if strcmp(about.compression.HeightMode,'Intensity')
                                found.height = summary.BackgroundStdDev(n)*found.height+summary.BackgroundMean(n);
                                found.prominence = summary.BackgroundStdDev(n)*found.prominence+summary.BackgroundMean(n);
                            end
                            peaks{n} = found;
                            summary.PeakCount(n) = length(found.height);   
                            %[~,indices] = intersect(mass_to_charge,found.mass_to_charge,'stable');
                            clear found;  
                            about.peakfinding.HeightMode = about.compression.HeightMode; % these are equivalent so one could just look at "about.compression.HeightMode" for both
                            duration.PeakFinding(n) = toc(timerpeak);
                        end
                        % update indices to keep
                        switch about.compression.Algorithm
                            case 'PeakThreshold' 
                                % keep indices of peaks above a certain threshold
                                [~,mask,~,~] = findpeaks((intensity - summary.BackgroundMean(n))/summary.BackgroundStdDev(n),...
                                    ['MinPeak' about.peakfinding.HeightReference],about.compression.HeightThreshold,...
                                    'MinPeakDistance',about.peakfinding.MinimumSeparation,...
                                    'WidthReference',about.peakfinding.WidthReference); % get height, locations, widths, and prominences of peaks            
                            otherwise % case 'SNRThreshold'
                                % keep indices of points above a certain threshold
                                about.compression.Type = 'SNRThreshold'; % indicate using default
                                mask = (intensity' - summary.BackgroundMean(n))/summary.BackgroundStdDev(n) > about.compression.HeightThreshold;
                        end                                    
                        keep(mask) = true; 
                        clear mask;
                        % end timer
                        duration.Pass1(n) = toc(timer);                    
                        if not(mod(n,num_status_update)) && object.verbose
                            object.notify([indent 'Compression Pass 1 Progress: ' sprintf('%0.0f',n/amaldims.num.positions*100) '%, ' sprintf('%0.2f',(amaldims.num.positions-n)*updateMinutes(duration.Pass1)) ' minutes remaining.'])
                        end
                    end
                    compressed = struct();
                    compressed.indices = find(keep);                
                    about.compression.CompressionRatio = length(compressed.indices)/length(intensity);                
                    clear keep                
                    if object.verbose
                        object.notify([indent 'Time elapsed: ' sprintf('%0.2f',sum(duration.Pass1)/60) ' minutes.'])
                        object.notify(['Compression Ratio: ' sprintf('%0.2f',about.compression.CompressionRatio*100) '% (' sprintf('%0.2e',length(intensity)) ' to ' sprintf('%0.2e',length(compressed.indices)) ' points)'])
                    end
    
                    % Compression Pass 2
                    if object.verbose        
                        object.notify('Compression Pass 2: constructing the compressed feature matrix.')
                    end
                    compressed.mass_to_charge = mass_to_charge(compressed.indices)';
                    compressed.signal = zeros(amaldims.num.positions,length(compressed.indices));
                    for n = 1:amaldims.num.positions
                        % start timer
                        timer = tic;
                        % load resampled spectra data
                        intensity = amaldims.get_data_respectra(n,spectra_count); 
                        % collect values in feature matrix
                        compressed.signal(n,:) = intensity(compressed.indices);
                        if strcmp(about.compression.HeightMode,'SignalToNoiseRatio')
                            compressed.signal(n,:) = (compressed.signal(n,:)-summary.BackgroundMean(n))/summary.BackgroundStdDev(n);
                        end        
                        % end timer
                        duration.Pass2(n) = toc(timer);
                        if not(mod(n,num_status_update)) && object.verbose
                            object.notify([indent 'Compression Pass 2 Progress: ' sprintf('%0.0f',n/amaldims.num.positions*100) '%, ' sprintf('%0.2f',(amaldims.num.positions-n)*updateMinutes(duration.Pass2)) ' minutes remaining.'])
                        end    
                    end
                    if object.verbose        
                        object.notify([indent 'Time elapsed: ' sprintf('%0.2f',sum(duration.Pass2)/60) ' minutes.'])
                    end
    
                    % Load Raw Data   
                    duration.Raw = 0;                                
                    if about.raw.Include
                        if object.verbose        
                            object.notify('Raw: copying raw time traces.')
                        end
                        % start timer
                        timer = tic;
                        % load raw data
                        if amaldims.has.raw
                            about.raw.Count = amaldims.raw.count;
                            about.raw.Points = amaldims.raw.points;
                            about.raw.Duration = amaldims.raw.duration; % [s]
                            about.raw.SampleRate = amaldims.raw.sample_rate; % [Hz] 
                            raw = struct();
                            [raw.signal, raw.time] = amaldims.get_data_raw(1,amaldims.num.positions);
                            raw.signal = raw.signal';
                            raw.time = raw.time';
                        else
                            about.raw.Include = false;
                            object.notify([indent 'Raw data not found.' ]);
                        end                    
                        % end timer
                        duration.Raw = toc(timer);
                        if not(mod(n,num_status_update)) && object.verbose
                            object.notify([indent 'Raw. Progress update: ' sprintf('%0.0f',n/amaldims.num.positions*100) '%, ' sprintf('%0.2f',(amaldims.num.positions-n)*updateMinutes(duration.Pass2)) ' minutes remaining.'])
                        end    
                        if object.verbose        
                            object.notify([indent 'Time elapsed: ' sprintf('%0.2f',updateMinutes(duration.Raw)) ' minutes.'])
                        end                    
                        
                    end                
    
                    % export to file
                    duration.Overall = toc(initial_time);                
                    about.compression.TimeStamp = datestr(now,'mmm DD, YYYY - HH:SS PM');
                    about.compression.DurationMinutes.Overall = duration.Overall/60;
                    about.compression.DurationMinutes.Compression = (sum(duration.Pass1)+sum(duration.Pass2)-sum(duration.PeakFinding))/60;
                    about.compression.DurationMinutes.PeakFinding = sum(duration.PeakFinding)/60;
                    about.compression.DurationMinutes.Raw = duration.Raw/60;
                    %about.signal = 'compressed intensity matrix (positions by feature intensities)';
                    %about.mass_to_charge = 'compressed mass-to-charge ratio vector';
                    %about.indices = 'the indices of features in the original mass-to-charge ratio vector';
                    %about.summary = 'a table summarizing the properties of each spectrum'; 
                    %about.methods = 'a structure containing the methods and instrument settings used to record the original spectra/resampled-spectra and the options used for compression';
                    filename = strrep(amaldims.filename,'.hdf5','_compressed.mat');
                    save(filename,'compressed','summary','about','-v7.3')
                    if about.peakfinding.FindPeaks
                        save(filename,'peaks','-append')
                    end
                    if about.raw.Include
                        save(filename,'raw','-append')
                    end
                    success = true;
                    object.notify(['Compression complete. (time elapsed: ' num2str(duration.Overall) ' minutes)'])
                else
                    warning('Compression failed. Respectra not found.')
                end 
           else
                    warning('Compression failed. Data not present.')                
           end           
        end
        
        
    end
    

end


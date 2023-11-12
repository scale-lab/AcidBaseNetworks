
classdef aMALDImsRed < handle
    % @title: AMALDIMSRED.m
    % @description: Class for accessing the data in an reduced (compressed) aMALDIms file (mat file)
    % @author: Chris Arcadia (christopher_arcadia@brown.edu)
    % @created: 2020/01/22
    
    properties % class variables        
        verbose = true; % turn on or off verbose warnings and messages      
        
        filename = ''; % filename of the aMALDIms file
        loaded = false; % boolean indicating if file was successfully loaded
        
        about = struct(); % info about the measurement and file
        summary = table(); % measurement summary
        peaks = cell(0); % peak lists for each spectra
        compressed = struct(); % compressed spectra data (contains signal matrix and mz values and mz indices from original spectra)
        
        mode = ''; % MS signal height mode
        positions = {}; % MALDI spot positions
        duplicates = []; % position index (sample number) of repeated measurements (if present)        
        plate = struct(); % MALDI plate layout info {'1536', '384', 'unknown'}
        
        has = struct(); % structure indicating what datasets and groups have been found and loaded
        num = struct(); % stucture providing the number of points in (length of) the datasets             
    end
    
    properties (Constant) % class constants        
                
        versions = {1}; % compatible versions of aMALDImsRed files, those supported by this class
                
    end
    
    methods % main
                
        function object = aMALDImsRed(filename,layout) 
        % construct an instance of this class   
            % set filename
            object.filename = filename;
            object.plate.layout = layout;            
            if exist(filename, 'file') == 2
                object.load();
            else
                warning('File not found.')
            end
        end
        
        function load(object)
        % load the contents of a file (excluding datasets)
                        
            % load file data
            temp = load(object.filename); % load compressed mass spectra data
            fields = {'about','summary','peaks','compressed'}; 
            for n = 1:length(fields)
                object.has.(fields{n}) = false;
                if isfield(temp,fields{n}) 
                    object.has.(fields{n}) = true;
                    object.(fields{n}) = temp.(fields{n});
                end
            end
            if object.has.compressed
                object.num.points = length(object.compressed.mass_to_charge);
            end
            object.mode = load_signal_mode(object);
            
            % load position info
            [object.positions, object.has.positions, object.duplicates] = object.load_positions(); % get positions while checking that the various versions of the position vector are equivalent
            object.has.duplicates = not(isempty(object.duplicates));
            object.num.positions = length(object.positions);    
            object.plate = object.determine_plate();                                   
            
            % ensure that required fields are present 
            object.require({'about','summary','compressed','positions'}); % a valid reduced aMALDIms file must have these fields      
            object.loaded = true;
        end
                
        function print(object) 
        % display structure of the file     
            fields = {'filename','about','compressed','plate','has','num'};
            for n = 1:length(fields)
                disp(['' fields{n} ':'])
                disp(object.(fields{n}))                    
            end            
        end        
        
    end
    
    methods % for loading data
                                                              
        
        function [signal_mode] = load_signal_mode(object) 
            switch object.about.compression.HeightMode
                case 'SignalToNoiseRatio'
                    signal_mode = 'SNR'; %'Signal-to-Noise Ratio';
                case 'Intensity'
                    signal_mode = 'Intensity';
                otherwise
                    signal_mode = 'Signal';    
            end
        end
        
        function [positions, found, indDuplicates, equivalent] = load_positions(object) 
        % load position vector (while checking if the various positions vectors are equivalent)   
        
            % initialize outputs
            positions = {};   
            indDuplicates = []; % index of duplicate positions
            
            % assign positions
            if object.has.summary
               positions = object.summary.Position;
            end
            
            % check that at least one of the position vectors was present in the file
            num_found = length(positions);
            found = num_found>0;
            equivalent = true; % initially assume all vectors are equivalent
            
            % check that each position is unique
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
        
        function [plate] = determine_plate(object)  
        % determine plate information
            plate = object.plate;
            if not(isempty(object.positions))
                % set plate dimensions (dimensions: [<rows>,<columns>], standard: <position format>)
                if ischar(plate.layout) % specify layout
                    switch plate.layout
                        case '24'
                            plate.dimensions = [4,6];  
                            plate.standard = 'Alpha';
                        case '96'
                            plate.dimensions = [8,12];       
                            plate.standard = 'Alpha';                        
                        case '384'
                            plate.dimensions = [16,24];  
                            plate.standard = 'Alpha';
                        case '1536'
                            plate.dimensions = [32,48];
                            plate.standard = 'XY';  
                        case '6144'
                            plate.dimensions = [64,96];
                            plate.standard = 'XY';          
                        otherwise
                            plate.dimensions = [1,1]*ceil(sqrt(length(object.positions)));
                            %plate.dimensions = [1,length(object.positions)];
                            plate.standard = 'XY';  
                    end
                else % specified dimension
                    plate.dimensions = plate.layout;
                    plate.layout = [num2str(plate.dimensions(1)) 'x' num2str(plate.dimensions(2))];
                    plate.standard = 'XY';                     
                end
                
                % identify position format
                letter_pattern = @(str) num2str(isstrprop(str,'alpha'));
                pos = object.positions{1};
                pos_pattern = letter_pattern(pos);
                pos_count = length(pos);
                switch pos_pattern
                    case letter_pattern('X00Y00')
                        plate.format = 'XY';
                    case letter_pattern('A00')
                        plate.format = 'Alpha';
                    otherwise
                        plate.format = 'unknown';
                end
            end
        end     
        
        function [positionsInt] = positionStr2Int(object,positionsStr)
        % convert positions in string form ('A1' or 'X01Y01') to integer form ([1,1])
            switch object.plate.format
                case 'Alpha'
                	positionsInt = object.positionAlphaStr2Int(positionsStr);
                case 'XY'
                    positionsInt = object.positionXYStr2Int(positionsStr);
                otherwise
                    positionsInt = ones(length(positionsStr),2);                
                    positionsInt(:,2) = 1:length(positionsStr);
            end
        end
        
        function [positionsStr] = positionInt2Str(object,positionsInt)
        % convert positions in integer form ([1,1]) to string form ('A1' or 'X01Y01')    
            switch object.plate.format
                case 'Alpha'
                    positionsStr = object.positionAlphaInt2Str(positionsInt);                    
                case 'XY'
                    positionsStr = object.positionXYInt2Str(positionsInt);
                otherwise
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
        
        
        function [values] = match_data_to_positions(object,data_values,data_positions)
        % match data values (a 1-D matrix or cell) for each given positions to their measurement positions (assumes 1-to-1 mapping)
            if isnumeric(data_values)
                values = nan(size(object.positions));
                for n=1:object.num.positions    
                    index = find(contains(data_positions,object.positions{n}));
                    if not(isempty(index))
                        values(n) = data_values(index);
                    end
                end
            elseif iscell(data_values)
                values = cell(size(object.positions));
                for n=1:object.num.positions    
                    index = find(contains(data_positions,object.positions{n}));
                    if not(isempty(index))
                        values{n} = data_values{index};
                    end
                end                
            else
                values = [];
                object.throw_error('unsupportedType');
            end
            if object.has.duplicates || length(unique(data_positions))<length(data_positions)
                warning('Duplicate positions present.')
            end
        end
        
        function throw_error(object, type)
        % method for throwing an errors by types 
            switch type
                case 'inputLengthMismatch'
                    message = 'Input vector must have the same number of elements as positions.';
                case 'unsupportedType'
                    message = 'Provided data type is unsupported.';                    
                case 'positionMissing'
                    message = 'Some or all positions provided were not found in the spectra file.';  
                case 'invalidInput'
                    message = 'Given input was invalid.';                                                                                               
                otherwise
                    message = type;
                    type = 'SpecificError';
            end
            error(['aMALDImsRed:' type],['Error: ' message])
        end        
        
%         function [indicesDesired] = check_positions(object,positionsDesired)
%         % check that positions in aMALDIms object are matched by the positions in the key and if not throw an error
%             indicesDesired = object.position2Ind(positionsDesired);
%             if isequal(object.positions,positionsDesired)
%                 disp(['Positions in key match those in the aMALDIms spectra file.'])
%             else
%                 if all(ismember(positionsDesired,object.positions))
%                     if isequal(sort(object.positions),sort(positionsDesired))
%                         warning(['Warning : positionOutOfOrder' newline 'The order of positions in key do not match their counterparts in aMALDIms spectra file. Try resorting key positions.)'])                                                
%                     else
%                         warning(['Warning : positionSubSet' newline 'Positions in key are a subset of those present in the aMALDIms spectra file.'])                                                                                                                             
%                     end
%                 else
%                     object.throw_error('positionMissing'); 
%                 end
%             end
%             clear positionsInKey
%         end        
        
        
    end
    
    methods % for analyzing data
        
        function [index, distance] = mz_shift(object,mass_to_charge,varargin)
            % get the indices of the nearest mass to charge ratio (mz) to the desired one or if a tolerance is specified the indices within the desired window
            if nargin<3 % find the nearest mz value
                [distance,index] = min(abs(object.compressed.mass_to_charge-mass_to_charge));
            else % find the mz values within the specified tolerance
                tolerance = varargin{1};
                distance = abs(object.compressed.mass_to_charge-mass_to_charge);
                index = find(distance<=tolerance);
                distance = distance(index);
            end
        end
        
    end

    
                                                                    
    methods % for visualizing data   

        function [Ncurrent] = get_color_count(object)
        % determine the number of colors to generate
        
            % Default number of colors to use
            fig = get(groot,'CurrentFigure');                        
            if isempty(fig)
               Ncurrent = size(get(groot,'DefaultFigureColormap'),1);
            else
               Ncurrent = size(fig.Colormap,1);
            end   
        end        
        
        function [cmap] = color(object,N)
        % default sequential color scheme (for intensity values)
        
            % RGB grid points to interpolate between
            %{
            swatch = [ 247   251   255
                       222   235   247
                       198   219   239
                       158   202   225
                       107   174   214
                        66   146   198
                        33   113   181
                         8    81   156
                         8    48   107 ]; % Blues  
            %}                         
            swatch = [ 255   245   240
                       254   224   210
                       252   187   161
                       252   146   114
                       251   106    74
                       239    59    44
                       203    24    29
                       165    15    21
                       103     0    13 ]; % Reds
            Nswatch = size(swatch,1);
            
            % Input handling
            if nargin < 2
                N = object.get_color_count();
            elseif isempty(N)
                N = object.get_color_count();
            elseif isnan(N)
                N = Nswatch;
            end            
                        
            % Color map generation
            cmap = nan(0,3);            
            if N~=0
                cmap = interp1(1:Nswatch,swatch/255,linspace(1,Nswatch,N),'makima'); % Y = interp1(Xo,Yo,X);
                %cmap = interp1(1:Nswatch,swatch/255,linspace(1,Nswatch,N),'spline'); % Y = interp1(Xo,Yo,X);                
                %cmap(cmap>1) = 1; 
                %cmap(cmap<0) = 0;
            end
            
            % cmap = flipud(bone); % winter % bone % parula % originally using native color schemes
        end        
        
        function [cmap] = colors(object,N)
        % default qualitative color scheme (for distinct traces)
        
            % RGB grid points to repeat
            swatch = [ 228    26    28
                        55   126   184
                        77   175    74
                       152    78   163
                       255   127     0
                       255   200    55 % 255   255    51
                       166    86    40
                       247   129   191
                       153   153   153 ] * 0.9; % Set 1 ()
            Nswatch = size(swatch,1);
            
            % Input handling
            if nargin < 2
                N = object.get_color_count();
            elseif isempty(N)
                N = object.get_color_count();
            elseif isnan(N)
                N = Nswatch;
            end            
                       
            % Color map generation
            cmap = nan(0,3);            
            if N~=0
                cmap = repmat(swatch/255,[ceil(N/Nswatch),1]);
                cmap = cmap(1:N,:);
            end
        end        
        
        
        function [fig] = plateview(object,label,variable,varargin)
        % function that plots the vector variable (of the same length and order as positions) 
        % as the color and height of spots arranged on the virtual MALDI plate        
            % options 
            positions = object.positions;
            xyFlipped = false;
            logscale = false;
            fontsize = 12;
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
                            case 'fontsize'
                                fontsize = value;
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
                fig = figure('color','w','name',['Plate View - ' label]);
            else
                fig = gcf;
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
                cmap = colormap(gca,object.color); % 'parula', 'summer', 'bone'
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
            zlabel(label)
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
            cb.Label.String = label;
            set(gca,'fontsize',fontsize)
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
                object.throw_error('InputLengthMismatch');
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
        
        function simple_plateview(object,vector,label,varargin)
        % simple plate view for subplots (such as in summary data figures)
            % options 
            logscale = false;
            Nvarargin = length(varargin);
            custom_cmap = false;
            new_figure = false;            
            linewidth = 0.5;
            marker = 's';
            markersize = 50;
            markeredgecolor = 'none';                       
            for n = 1:Nvarargin
                if Nvarargin>=n+1
                    if ischar(varargin{n}) || isscalar(varargin{n})
                        key = lower(varargin{n});
                        value = varargin{n+1};
                        switch key
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
                            case 'colormap'
                                cmap = value;
                                custom_cmap = true;
                            case 'newfigure'
                                new_figure = value;                                 
                        end
                    end
                end
            end        
            if new_figure
                fig = figure('color','w','name',['Simple Plate View - ' label]);
            else
                fig = gcf;
            end            
            %figure('color','w','name','Plate View for Summary Data')
            scatter3(object.summary.GridXPosition,object.summary.GridYPosition,vector,markersize,vector,...
                'filled',marker,'markeredgecolor',markeredgecolor,'linewidth',linewidth)
            view(0,-90)
            if custom_cmap
                cmap = colormap(gca,cmap);                
            else
                cmap = colormap(gca,object.color); % 'parula', 'summer', 'bone'
            end            
            cb = colorbar('location','southoutside');
            %cb.Label.String = label;
            title(label,'fontweight','normal')
            grid off; box on;
            xlimits = [min(object.summary.GridXPosition),max(object.summary.GridXPosition)];
            xlim(xlimits+[-1,1])
            ylimits = [min(object.summary.GridYPosition),max(object.summary.GridYPosition)];
            ylim(ylimits+[-1,1])
            if logscale
                set(gca,'zscale','log');
                if isprop(gca,'colorscale') % colorscale property is only supported in newer MATLAB versions
                    set(gca,'colorscale','log');
                end
            end
            
        end  
        % example: object.plateviewforsummary('TIC')
        
        function summary_plateview(object)
        % plot summary data in a plateview
            raw_height_mode = 'Intensity [arb.]';   
            figure('color','w','name','Summary Plateview','units','normalized','position',[0,0,1,1])
            n = 1; nx = 2; ny = 3; %legpos = 'eastoutside'; marker = '.'; markersize = 10; linewidth = 1;
            subplot(nx,ny,[n]); n = n+1;
                object.simple_plateview(object.summary.TIC,'TIC [arb.]')
            subplot(nx,ny,[n]); n = n+1;
                object.simple_plateview(object.summary.Mean,['Mean ' raw_height_mode]); %'ElapsedTime')                        
            subplot(nx,ny,[n]); n = n+1;
                object.simple_plateview(object.summary.PeakCount,'Peak Count')
            subplot(nx,ny,[n]); n = n+1;
                object.simple_plateview(object.summary.BackgroundMean,'Mean Background [arb.]')
            subplot(nx,ny,[n]); n = n+1;
                object.simple_plateview(object.summary.BackgroundStdDev,'Std. Dev. Background [arb.]')
            subplot(nx,ny,[n]); n = n+1;
                object.simple_plateview(object.summary.TukeyUpper,['Tukey Upper ' raw_height_mode])
        end
        
        function summary_plot(object)
        % plot all summary data 
            raw_height_mode = 'Intensity [arb.]';
            figure('color','w','name','Summary Plot','units','normalized','position',[0,0,1,1])
            n = 1; nx = 4; ny = 4; legpos = 'eastoutside'; marker = '.'; markersize = 10; linewidth = 1; cmap = object.colors(nan); ax = [];
            ax(1) = subplot(nx,ny,[n,n+1]); n = n+2;
                hold on
                plot(object.summary.GridXPosition,marker,'markersize',markersize,'color',cmap(1,:))
                plot(object.summary.GridYPosition,marker,'markersize',markersize,'color',cmap(2,:)) % [cmap(3) cmap(2) cmap(1)]
                hold off
                ylabel('Grid Position')
                legend('X (column)','Y (row)','location',legpos)
                title('Sample Number','fontweight','normal')
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
            ax(2) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot([NaN;diff(object.summary.ElapsedTime)]*60,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel(['Time Delay Between' newline 'Measurements [sec.]'])
                %legend('X','Y','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))              
            ax(3) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot(object.summary.TIC,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel('TIC [arb.]')
                %legend('X','Y','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
                %set(gca,'yscale','log')    
            ax(4) = subplot(nx,ny,[n,n+1]); n = n+2;
                hold on
                plot(object.summary.Mean,marker,'markersize',markersize,'color',cmap(1,:))
                plot(object.summary.Mode,marker,'markersize',markersize,'color',cmap(2,:))
                plot(object.summary.Median,marker,'markersize',markersize,'color',cmap(3,:))
                plot(object.summary.Quartile1,marker,'markersize',markersize,'color',cmap(4,:))
                plot(object.summary.Quartile3,marker,'markersize',markersize,'color',cmap(5,:))  
                plot(object.summary.TukeyUpper,marker,'markersize',markersize,'color',cmap(6,:))      
                hold off
                ylabel(raw_height_mode)
                legend('Mean','Mode','Median','Quartile1','Quartile3','TukeyUpper','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
            ax(5) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot(object.summary.Maximum,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel(['Maximum' newline raw_height_mode])
                %legend('Maximum','Mininum','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
                %set(gca,'yscale','log')   
            ax(6) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot(object.summary.Mininum,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel(['Minimum' newline raw_height_mode])
                %legend('Maximum','Mininum','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
                %set(gca,'yscale','log')  
            ax(7) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot(object.summary.StdDev,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel(['Standard Deviation' newline raw_height_mode])
                %legend('Maximum','Mininum','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
                %set(gca,'yscale','log')       
            ax(8) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot(object.summary.Area,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel(['Area Under the' newline 'Spectrum [arb.*Da]'])
                %legend('Maximum','Mininum','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
                %set(gca,'yscale','log')       
            ax(9) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot(object.summary.Sum,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel(['Sum ' raw_height_mode])
                %legend('Maximum','Mininum','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
                %set(gca,'yscale','log')       
            ax(10) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot(object.summary.PeakCount,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel(['Peak Count'])
                %legend('Maximum','Mininum','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
                %set(gca,'yscale','log')    
            ax(11) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot(object.summary.BackgroundMean,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel(['Mean' newline 'Background [arb.]'])
                %legend('Maximum','Mininum','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
                %set(gca,'yscale','log')    
            ax(12) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot(object.summary.BackgroundStdDev,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel(['Standard Deviation' newline 'Background [arb.]'])
                %legend('Maximum','Mininum','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
                %set(gca,'yscale','log')    
            ax(13) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot(object.summary.BackgroundLinearSlope,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel(['Line Slope' newline 'Background [arb./Da]'])
                %legend('Maximum','Mininum','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
                %set(gca,'yscale','log')    
            ax(14) = subplot(nx,ny,[n]); n = n+1;
                hold on
                plot(object.summary.BackgroundLinearOffset,marker,'markersize',markersize,'color',cmap(1,:))
                hold off
                ylabel(['Line Offset' newline 'Background [arb.]'])
                %legend('Maximum','Mininum','location',legpos)
                box on; grid off; axis tight; ylim(ylim+0.1*[-1,1]*diff(ylim)); xlim(xlim+0.025*[-1,1]*diff(xlim))
                %set(gca,'yscale','log')  
            linkaxes(ax,'x')
        end
        
        
        
        function [fig] = heatmap(object,varargin)
        % function that plots a spectra heat map
            % options 
            logscale = true;
            zLimits = [];
            Nvarargin = length(varargin);
            custom_cmap = false;
            new_figure = true;
            %positions = object.positions;                        
            %Npos = length(positions);  
            Npos = length(object.positions);              
            label = 'Spectrum Number';            
            variable = (1:Npos)';
            figpos = [0 0 1 1];
            for n = 1:Nvarargin
                if Nvarargin>=n+1
                    if ischar(varargin{n}) || isscalar(varargin{n})
                        key = lower(varargin{n});
                        value = varargin{n+1};
                        switch key
                            case 'logscale'
                                logscale = value;
                            case 'colormap'
                                cmap = value;
                                custom_cmap = true;
                            case 'newfigure'
                                new_figure = value; 
                            case 'label'
                                label = value;  
                            case 'variable'
                                variable = value;
                            case 'figureposition'
                                figpos = value;
                            case 'zlimits'
                                zLimits = value;
                        end
                    end
                end
            end
            % plot the resulting plate
            if new_figure
               fig = figure('color','w','name',['Spectrum Heat Map - ' label],'units','normalized','position',figpos);
            else
               fig = gcf;
            end       
            
            % format data
            mz = object.compressed.mass_to_charge;
            sig = object.compressed.signal;
            
            % clip snr values below zero
            sig(sig<0) = nan;

            % remove spectra of matrix alone
            mask = isnan(variable);
            variable(mask) = [];
            sig(:,mask) = [];
            % OR
            % include matrix alone samples at end of time series
            % variable(isnan(variable)) = max(variable)+mean(abs(diff(sort(variable(~isnan(variable))))))*[1:sum(isnan(variable))];

            % sort spectra by variable
            [variable,mask] = sort(variable);
            sig = sig(mask,:);
            
            % get limits
            %sigmax = max(sig(:));
            %sigmin = min(sig(:));            

            % plot heatmap
            if not(custom_cmap)            	
                cmap = object.color; % 'parula', 'summer', 'bone'
            end             
            %Ncolor = size(cmap,1);
            if logscale
                heatfunc = @(x) log10(x);                                
                heatlabel = ['log10(' object.mode ')'];
            else
                heatfunc = @(x) x;
                heatlabel = [object.mode];
            end    
%             sig = round(object.linear_map(heatfunc(sig),heatfunc([snrmin,snrmax]),[1,Ncolor]))';
            %image(sig,'CDataMapping','scaled'); 
            image(double(mz), double(variable), heatfunc(sig),'CDataMapping','scaled')            
            %mz_resolution = 0.01; % [Da] 
            %fig = msheatmap(mz,variable,heatfunc(sig),'Resolution',mz_resolution); 
            %set(fig,'color','w','name','Spectra Heat Map','units','normalized','position',[0 0 1 1]); % [0.2 0.2 0.5 0.625])
            %cb = get(gca,'Colorbar');
            ylabel(label)
            xlabel('m/z [Da]')
            colormap(gca,cmap);             
            cb = colorbar;            
            cb.Label.String = heatlabel;   
            if isempty(zLimits)
                zLimits = [0,ceil(heatfunc(max(sig(:))))]; % [0,ceil(heatfunc(sigmax))];
            end
            caxis(zLimits)            
            box on; grid off;
            title(['Spectra Heat Map'],'fontweight','normal');   
        end
        
        
        function simple_boxplot(object,X,Y,varargin) %width,color,linewidth)
        % function to make a boxplot            
            % options 
            Nvarargin = length(varargin);    
            new_figure = false;            
            linewidth = 0.75;
            cmap = object.colors();
            color = cmap(1,:);
            width = min(abs(diff(sort(X))))/2; %mean(abs(diff(X)))/2;
            figpos = [0.3307    0.1870    0.3646    0.6130];
            for n = 1:Nvarargin
                if Nvarargin>=n+1
                    if ischar(varargin{n}) || isscalar(varargin{n})
                        key = lower(varargin{n});
                        value = varargin{n+1};
                        switch key
                            case 'linewidth'
                                linewidth = value;
                            case 'color'
                                color = value;
                            case 'newfigure'
                                new_figure = value;  
                            case 'figureposition'
                                figpos = value;   
                            case 'width'
                                width = value;
                        end
                    end
                end
            end
            
            % make plots
            if new_figure
                fig = figure('color','w','name',['Simple Box Plots'],'units','normalized','position',figpos);
            else
                fig = gcf;
            end
            
            % plot the box for each value of X
            N = length(X);
            if isempty(width)
                width = abs(X)*1/10;
            end
            for n=1:N
                % stats
                Q1 = quantile(Y(:,n),0.25); % 1st quartile
                Q2 = quantile(Y(:,n),0.50); % median
                Q3 = quantile(Y(:,n),0.75); % 3rd quartile
                Max = max(Y(:,n));
                Min = min(Y(:,n));
                %Mean = mean(Y(:,n));
                % box
                x = X(n) - width/2;
                y = Q1;
                h = Q3-Q1;
                rectangle('Position',[x y width h],'LineWidth',linewidth,'EdgeColor',color,'FaceColor','none') % x & y define the lower left corner of the rectangle
                % center
                plot(X(n)+width/2*[-1,1],Q2*[1,1],'-','LineWidth',linewidth,'color',color)
                %plot(X(n),Mean,'o','LineWidth',boxlinewidth,'color',boxcolor)
                % whiskers
                plot(X(n)*[1,1],[Min,Q1],'-','LineWidth',linewidth,'color',color)
                plot(X(n)*[1,1],[Q3,Max],'-','LineWidth',linewidth,'color',color)
                % whisker ends
                plot(X(n)+width/2*[-1,1],Min*[1,1],'-','LineWidth',linewidth,'color',color)
                plot(X(n)+width/2*[-1,1],Max*[1,1],'-','LineWidth',linewidth,'color',color)
            end
        end        
        
        function [fig] = view_mass_window(object,mass_to_charge,tolerance,varargin)
        % function that plots the spectra signals within a given mass window (mass_to_charge +/ tolerance) as a platemap, boxplot, and series plot
            % options 
            Nvarargin = length(varargin);                        
            Npos = length(object.positions);   
            label = 'Spectrum Number';            
            variable = (1:Npos)';            
            cmap = object.colors();
            markersize = 6; 
            pointersize = 5;
            linewidth = 0.75; 
            new_figure = true;
            figpos = [0.3307    0.1870    0.3646    0.6130];
            numformat = '%0.5f';
            selection = 'max'; % {'max','near'}
            use_relative_mz = true;   
            show_selected = false;
            logscale = true; % false;            
            for n = 1:Nvarargin
                if Nvarargin>=n+1
                    if ischar(varargin{n}) || isscalar(varargin{n})
                        key = lower(varargin{n});
                        value = varargin{n+1};
                        switch key
                            case 'linewidth'
                                linewidth = value;
                            case 'markersize'
                                markersize = value;
                            case 'colormap'
                                cmap = value;
                            case 'newfigure'
                                new_figure = value;  
                            case 'figureposition'
                                figpos = value;                                
                            case 'label'
                                label = value;                                   
                            case 'variable'
                                variable = value;  
                            case 'pointselection'
                                selection = value;  
                            case 'pointersize'
                                pointersize = value;                                
                            case 'showselection'
                                show_selected = value;
                            case 'relative'
                                use_relative_mz = value;  
                            case 'logscale'
                                logscale = value;
                        end
                    end
                end
            end
            
            % format data
            mz_target = mass_to_charge;
            mz_tol = tolerance;
            %snr_threshold = 100; %30;

            % get data
            [ind,dist] = object.mz_shift(mz_target,mz_tol);
            mz = object.compressed.mass_to_charge(ind);
            sig = object.compressed.signal(:,ind);
            switch selection
                case 'max'
                    selectLab = 'Maximum';
                    [selectSig,selectInd] = max(sig,[],2);
                case 'near'
                    selectLab = 'Nearest';
                    [~,nearInd] = min(dist);
                    selectSig = sig(:,nearInd);
                    selectInd = nearInd*ones(size(selectSig));
                otherwise
                    object.throw_error('invalidInput')
            end
            mz_max_avg = mean(mz(selectInd));

            % make plots
            if new_figure
                fig = figure('color','w','name',['Mass Search - ' label],'units','normalized','position',figpos);
            else
                fig = gcf;
            end
            
            % plate view
            subplot(4,4,1:12)
            object.simple_plateview(selectSig,[selectLab ' In Range ' object.mode]);
            if logscale
                set(gca,'zscale','log');
                if isprop(gca,'colorscale') % colorscale property is only supported in newer MATLAB versions
                    set(gca,'colorscale','log');
                end
            end
            % object.plateview([],maxsnr,'newfigure',false,'logscale',true,'fontsize',10)
            % title([selectLab ' In Range ' object.mode])
            % object.plateview([],maxsnr>snr_threshold,'newfigure',false,'logscale',false,'fontsize',10); 
            % title(['Maximum in Range SNR > ' num2str(snr_threshold)])

            % shifted data
            subplot(4,4,13:14)
            hold on
            % normal plot
            % plt = plot(mz-use_relative_mz*mz_target,sig,'.','markersize',markersize);
            % set(plt, {'color'}, num2cell(object.colors(size(sig,1)),2));
            % selected points
            if show_selected
                plot(mz(selectInd)-use_relative_mz*mz_target,selectSig,'.k','markersize',markersize)
            end
            % box plot
            object.simple_boxplot(mz - use_relative_mz*mz_target,sig,'color',cmap(1,:),'linewidth',linewidth)            
            % set(gca,'yscale','log')
            axis tight
%             plot([1,1]*mz_target-use_relative_mz*mz_target,ylim,'-','color','k','linewidth',linewidth)
            pltline = plot(mz_target-use_relative_mz*mz_target,max(ylim),'v','color','k','markersize',pointersize,'linewidth',linewidth,'markerfacecolor','k');
            uistack(pltline,'bottom')            
            % plot([1,1]*(mz_target+mz_tol-use_relative_mz*mz_target),ylim,'--k','linewidth',linewidth)
            % plot([1,1]*(mz_target-mz_tol -use_relative_mz*mz_target),ylim,'--k','linewidth',linewidth)
            hold off
            xlim((mz_target-use_relative_mz*mz_target)+mz_tol*[-1,1])
            title([sprintf(numformat,mz_target) ' +/- ' sprintf(numformat,mz_tol) ' Da'],'fontweight','normal')
            ylabel(object.mode)
            if logscale
                set(gca,'yscale','log');
                yticks(10.^(floor(log10(min(sig(:)))):1:ceil(log10(max(sig(:))))))
            end
            
            xlabel('m/z - Target [Da]')
            % view([90, -90])
            box on; grid off;

            % target mass across samples
            mz_drift = abs(mz_target-mz_max_avg);
            subplot(4,4,15:16)
            plot(variable,selectSig,'.','color',cmap(1,:),'markersize',markersize)
            axis tight
            title(['\mu = ' sprintf(numformat,mz_max_avg) ' Da (\Delta = ' sprintf(numformat,mz_drift)  ' Da' ')'],'fontweight','normal')
            ylabel([selectLab ' in' newline 'Range ' object.mode])
            set(gca,'YAxisLocation','right')
            if logscale
                set(gca,'yscale','log');
                yticks(10.^(floor(log10(min(selectSig(:)))):1:ceil(log10(max(selectSig(:))))))
            end
            xlabel(label)
            box on; grid off;

            % print message
            if object.verbose
                object.notify('Mass Search')
                object.notify([' Target mz: ' sprintf(numformat,mz_target) ' (+/-' sprintf(numformat,mz_tol) ')' ' Da'])
                object.notify([' Average mz of Points with the Maximum SNR: ' sprintf(numformat,mz_max_avg) ' Da'])
                object.notify([' Distance between the Target and Found mz: ' sprintf(numformat,mz_drift) ' Da'])
            end
        end
        
        
        
        
        
        
        function [fig] = view_multiple_mass_windows(object,mass_to_charge,tolerance,varargin)
        % function that plots the spectra signals within a multiple given mass window (mass_to_charge +/ tolerance) as a boxplots or scatter plots
            % options 
            Nvarargin = length(varargin);                        
            cmap = object.colors();
            markersize = 10;%6; 
            pointersize = 5;            
            linewidth = 0.75; 
            new_figure = true;
            figpos = [0.3307    0.1870    0.3646    0.6130];
            numformat = '%0.5f';
            use_relative_mz = true;   
            plot_type = 'median'; % 'box'; % {'scatter','mean','median','max','box'}
            label = {};
            name = '';            
            fontsize = 8;
            link_axes = true; % false;
            logscale = true; % false;
            for n = 1:Nvarargin
                if Nvarargin>=n+1
                    if ischar(varargin{n}) || isscalar(varargin{n})
                        key = lower(varargin{n});
                        value = varargin{n+1};
                        switch key
                            case 'linewidth'
                                linewidth = value;
                            case 'markersize'
                                markersize = value;
                            case 'colormap'
                                cmap = value;
                            case 'newfigure'
                                new_figure = value;  
                            case 'figureposition'
                                figpos = value;                                
                            case 'relative'
                                use_relative_mz = value; 
                            case 'plottype'
                                plot_type = value;
                            case 'label'
                                label = value;  
                            case 'name'
                                name = value;                                     
                            case 'fontsize'
                                fontsize = value;
                            case 'pointersize'
                                pointersize = value;
                            case 'link_axes'
                                link_axes = value;
                            case 'logscale'
                                logscale = value;
                        end
                    end
                end
            end
            
            % make plots
            K = length(mass_to_charge);
            if new_figure
                if isempty(name)
                   name = ['Mass Searches - ' upper(plot_type(1)) plot_type(2:end)];
                else
                   name = ['Mass Searches - ' upper(plot_type(1)) plot_type(2:end) ' - ' name];
                end
                fig = figure('color','w','name',name,'units','normalized','position',figpos);
            else
                fig = gcf;
            end
            nsq = ceil(sqrt(K));
            ax = [];
            for k=1:K
                % format data
                mz_target = mass_to_charge(k);
                mz_tol = tolerance;
                %snr_threshold = 100; %30;
                
                % get data
                ind = object.mz_shift(mz_target,mz_tol);
                mz = object.compressed.mass_to_charge(ind);
                sig = object.compressed.signal(:,ind);
                
                % plot the shifted data
                ax(k) = subplot(nsq,nsq,k);
                hold on
                switch plot_type 
                    case 'scatter' % scatter plot                   
                        plt = plot(mz-use_relative_mz*mz_target,sig,'.','markersize',markersize);
                        set(plt, {'color'}, num2cell(object.colors(size(sig,1)),2));
                    case 'mean' % averages plot                   
                        plot(mz-use_relative_mz*mz_target,mean(sig,1),'.','markersize',markersize);
                    case 'median' % medians plot                   
                        plot(mz-use_relative_mz*mz_target,median(sig,1),'.','markersize',markersize);        
                    case 'max' % maxima plot                   
                        plot(mz-use_relative_mz*mz_target,max(sig,[],1),'.','markersize',markersize);                                                
                    case 'box' % box plot
                        object.simple_boxplot(mz-use_relative_mz*mz_target,sig,'color',cmap(1,:),'linewidth',linewidth)            
                end
                if logscale
                    set(gca,'yscale','log')
                    yticks(10.^(0:9))
                    %yticks(10.^(floor(log10(min(sig(:)))):1:ceil(log10(max(sig(:))))))
                end
                axis tight
                %plot([1,1]*mz_target-use_relative_mz*mz_target,ylim,'-','color','k','linewidth',linewidth)
%                 pltline = plot(mz_target-use_relative_mz*mz_target,max(ylim),'v','color','k','markersize',pointersize,'linewidth',linewidth,'markerfacecolor','k'); % markersize/2, % markersize
%                 uistack(pltline,'bottom')
                % plot([1,1]*(mz_target+mz_tol-use_relative_mz*mz_target),ylim,'--k','linewidth',linewidth)
                % plot([1,1]*(mz_target-mz_tol -use_relative_mz*mz_target),ylim,'--k','linewidth',linewidth)
                xlim((mz_target-use_relative_mz*mz_target)+mz_tol*[-1,1])
                %title([sprintf(numformat,mz_target) ' +/- ' sprintf(numformat,mz_tol) ' Da'],'fontweight','normal','fontsize',fontsize)
                if isempty(label)
                    title([sprintf(numformat,mz_target) ' Da'],'fontweight','normal','fontsize',fontsize)
                else
                    title(label{k},'fontweight','normal','fontsize',fontsize)
                end
                ylabel(object.mode)
                xlabel('m/z - Target [Da]')
                % view([90, -90])
                hold off; box on; grid off; 
                
            end
            if link_axes
                linkaxes(ax,'y')
            end
            
        end        
        
        
        
        
       
        
        
        
    end
                
    methods % for handling errors and command line notifications

        function notify(object,message)
        % function to display a command line message
            if object.verbose
                disp(['aMALDImsRed: ' message]);              
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
                err.identifier = 'aMALDImsRed:requiredFieldsNotFound';
                error(err)                
            end
        end
        
    end        
    
    methods % data export
        
        function capture_figures(object,output_directory,figure_handles,figure_format)
            % capture all provided figures
            %figure_handles = findobj('Type', 'figure');
            for n = 1:length(figure_handles)
                fig = figure_handles(n);
                filename = fullfile(output_directory,[sprintf('%03d',n) '_' matlab.lang.makeValidName(strrep(fig.Name,' ','_'))]);
                fig.PaperPositionMode = 'auto';
                fig_pos = fig.PaperPosition;
                fig.PaperSize = [fig_pos(3) fig_pos(4)];
                % supported figure formats: {'pdf','png','fig'}
                if contains(figure_format,'pdf')        
                    print(fig,[filename '.pdf'],'-dpdf') % as PDF
                end    
                if contains(figure_format,'png')    
                    print(fig,[filename '.png'],'-dpng') % as PNG 
                end
                if contains(figure_format,'fig')
                    savefig(fig,[filename '.fig']) % as FIG
                end
                disp(['captured figure: #' sprintf('%03d',n) ' - '  fig.Name])     
            end
            %clear figures fig filename
        end                        
       
        
    end
    
    
    
    
    methods (Static) % position parsing 
                
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
    
    methods % math functions
        
        function [valmapped, newbounds, oldbounds] = linear_map(val,oldbounds,newbounds)
        % function to linearly map a value from one interval to another interval:
        %   val in oldbounds = [A, B] --to--> valmapped in newbounds = [a, b]
            oldMax = max(oldbounds(:));
            oldMin = min(oldbounds(:));
            oldbounds = [oldMin, oldMax];
            newMax = max(newbounds(:));
            newMin = min(newbounds(:));
            newbounds = [newMin, newMax];
            valmapped = (val - oldMin)*(newMax-newMin)/(oldMax-oldMin) + newMin; 
            % the linear mapping is given by : valmapped = (val - A)*(b-a)/(B-A) + a
        end
        
    end
    
    methods (Static) % importers
                
        function reagent = importReagentReview(filename)
            % load reagent review csv file
                        
            % import csv
            importedData = readtable(filename);  

            % collect positions
            position = unique(importedData.position);  
            reagent = struct();
            reagent.Position = aMALDIms.positionXYInt2Str(fliplr(aMALDIms.positionAlphaStr2Int(position))); % convert from alpha to XY notation                        
            reagent.Count.Position = length(reagent.Position);            
            grid = aMALDIms.positionXYStr2Int(reagent.Position); 
            
            % make position grid
            reagent.Grid = table();
            reagent.Grid.X = grid(:,1);
            reagent.Grid.Y = grid(:,2);
                       
            % collect chemicals
            reagent.Chemical = table();
            reagent.Chemical.ID = unique(importedData.id);
            reagent.Count.Chemical = length(reagent.Chemical.ID);
            reagent.Chemical.Mass = nan(size(reagent.Chemical.ID));
            reagent.Chemical.Name = cell(size(reagent.Chemical.ID));
            reagent.Chemical.Type = cell(size(reagent.Chemical.ID)); 
            reagent.Volume = zeros(reagent.Count.Position,reagent.Count.Chemical);
            for n=1:length(importedData.position)
                chemicalIndex = find(reagent.Chemical.ID == importedData.id(n),1);
                reagent.Chemical.Mass(chemicalIndex) = importedData.mass(n);
                reagent.Chemical.Name{chemicalIndex} = importedData.name{n};
                reagent.Chemical.Type{chemicalIndex} = importedData.type{n};
                positionIndex = find(strcmp(position, importedData.position{n}),1);
                reagent.Volume(positionIndex,chemicalIndex) = importedData.volume(n);
            end
            reagent.Chemical = sortrows(reagent.Chemical,'Type'); % sort by type
            
            % collect chemical types
            reagent.Type.Name = unique(importedData.type);
            reagent.Count.Type = length(reagent.Type.Name);
            reagent.Type.Mask = false(reagent.Count.Type,reagent.Count.Chemical);    
            for n=1:reagent.Count.Chemical
                typeIndex = find(strcmp(reagent.Type.Name, reagent.Chemical.Type{n}),1);
                if not(isempty(typeIndex))
                    reagent.Type.Mask(typeIndex,n) = true;
                end
            end
            
        end
        
        function [plate] = importPlateSheet(filename,varargin) % originally outputted "[platesheet, grid]"
            % load platesheet csv file as a table

            % import csv
            [~,~,importedData] = xlsread(filename);
            num_fields = size(importedData,2);
            num_entries = size(importedData,1);    

            % check for additional inputs
            row_of_variable_names = 6;
            rows_to_read = 7:num_entries;     
            grid_dimensions = [16,24]; % default 384-well plate dimensions
            for n=1:length(varargin)
                if ischar(varargin{n})
                    switch lower(varargin{n})
                        case 'varnamerow'
                            row_of_variable_names = varargin{n+1};
                        case 'readrows'
                            rows_to_read = varargin{n+1};     
                        case 'griddimensions'
                            grid_dimensions = varargin{n+1};
                    end
                end
            end            

            % format data to table
            platesheet = table();
            for n=1:num_fields        
                field = genvarname(importedData{row_of_variable_names,n});
                data = importedData(rows_to_read,n);
                if isequal(field,'CID') || isequal(field,'Concentration') || isequal(field,'Volume')
                    mask = cellfun(@(C) isnumeric(C), data);
                    data(not(mask)) = {nan};   %guessing that "blank cells" means "empty array"
                    data = cell2mat(data);
                end
                platesheet.(field) = data;
            end

            % convert position to grid
            num_chemicals = length(platesheet.Start);
            grid = cell(num_chemicals,1);
            pStart = aMALDIms.positionAlphaStr2Int(platesheet.Start);
            pEnd = aMALDIms.positionAlphaStr2Int(platesheet.End);
            for n = 1:num_chemicals 
                grid{n} = false(grid_dimensions);
                grid{n}(pStart(n,1):pEnd(n,1),pStart(n,2):pEnd(n,2)) = true;
            end
            
            % format output
            plate = struct();
            plate.contents = platesheet;
            plate.layout = grid;
            
        end         
                
        
        function [transfers] = importTransferReport(filename,varargin) 
            % load Echo transfer report csv file as a table

            % check for additional inputs
            row_of_variable_names = 10;
            rows_at_end_to_skip = 4;     
            grid_dimensions = [16,24]; % default 384-well plate dimensions
            stampformat = 'yyyy-mm-dd HH:MM:SS.FFF'; % eg. '2019-11-25 10:39:35.765'
            for n=1:length(varargin)
                if ischar(varargin{n})
                    switch lower(varargin{n})
                        case 'varnamerow'
                            row_of_variable_names = varargin{n+1};
                        case 'tailrowcount'
                            rows_at_end_to_skip = varargin{n+1};     
                        case 'griddimensions'
                            grid_dimensions = varargin{n+1};
                        case 'stampformat'
                            stampformat = varargin{n+1}; 
                    end
                end
            end            

            % format data
            report = readtable(filename,'HeaderLines',row_of_variable_names);
            report([end+1-(1:rows_at_end_to_skip)],:) = [];
            transfers = table();
            transfers.TimeStamp = report.DateTimePoint;
            transfers.SourceWell = report.SourceWell;
            transfers.DestinationWell = report.DestinationWell;            
            transfers.TransferVolume = report.TransferVolume;
            transfers.TransferredVolume = report.ActualVolume;
            transfers.SurveyFluidVolume = report.SurveyFluidVolume;
            transfers.ExpectedFluidVolume = report.CurrentFluidVolume;
            transfers.FluidComposition = report.FluidComposition;
            
            % generate elapsed time vector (in seconds)            
            transfers.TimeElapsedInSeconds = 0*transfers.TransferredVolume;
            lastvec = datevec(transfers.TimeStamp{1},stampformat);
            for n=2:length(transfers.TimeElapsedInSeconds)
                currvec = datevec(transfers.TimeStamp{n},stampformat);
                transfers.TimeElapsedInSeconds(n) = transfers.TimeElapsedInSeconds(n-1) + etime(currvec,lastvec);
                lastvec = currvec;
            end
                        
        end                        
               
    end

end


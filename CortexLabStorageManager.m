% Class that provides default Suite2P storage schema
%
% This is an implementation of BaseStorageManager class that provides the default Suite2P file
% storage schema. In addition to fields that are used by BaseStorageManager, this classes uses
% the following fields for db entries (can be provided via options as well):
% rootStorage       - string. Path to folder where data is located. May contain a root path, which
%                     is will be appended by '/mouseName/date/experiment' or it may point to
%                     a folder, which contains Tiff files. rootStorage is used as is if it points
%                     to folder with Tiff files. Must be provided for the first db entry!
%                     If following entries doesn't specify rootStorage, then the value is inherited.
%                     Thus, it is possible to assign a different storage folder to each db entry.
% regFilePath       - string. (A former RegFileRoot from Suite2P). Specifies the location for
%                     binary file, preferable on an SSD drive. Must be provided for the first
%                     db entry!
% resultsSavePath     - string. Path to the folder where result files will be stored. Based on
%                     implementation, files might be stored not in resultsSavePath directly,
%                     but in sub folders. If not provided, rootStorage is used.
%                     CortexLabStorageManager creates a typical file structure:
%                     rootStorage/mouseName/date/experiment
% RegFileBinLocation  - string. Location of binary file on slow storage (i.e. on a remote storage).
%                       If empty, then the registration file is not preserved between Suite2P runs.
%                       The value is inherited through db entries.
% RegFileTiffLocation - string. If not empty, saves registration results as tiff files in
%                       that folder. Saving tiff files is slow.
%                       The value is inherited through db entries.
% temp_tiff         - string. Full path for a temporary tiff file. If specified, each remote tiff
%                     file is copied to temp_tiff before being worked on.
%                     The value is inherited through db entries.
%
classdef CortexLabStorageManager < BaseStorageManager
    properties (Access=public)
    end

    properties (Access=private)
    end

    methods(Static=true)
        % Convert db in original Suite2P format to new format understandable by BaseStorageManager.
        % This function can be used to automatically convert db structure.
        function newEntries = convertDbToStorageManager(db)
            for i = 1:length(db)
                oldEntry = db(i);

                newEntries(i).mouseName = oldEntry.mouse_name;
                newEntries(i).experiments = oldEntry.expts;
                if isfield(oldEntry, 'RootDir')
                    newEntries(i).rootStorage = oldEntry.RootDir;
                    oldEntry = rmfield(oldEntry, 'RootDir');
                end
                oldEntry = rmfield(oldEntry, 'mouse_name');
                oldEntry = rmfield(oldEntry, 'expts');

                leftFields = fieldnames(oldEntry);
                for j = 1:length(numel(leftFields))
                    curField = leftFields{j};
                    newEntries(i).(curField) = oldEntry.(curField);
                end
            end
        end
    end

    methods
        function this = CortexLabStorageManager(options)
            this = this@BaseStorageManager(options);
        end

        function isDone = isRegistrationDone(this)
            isDone = false;

            options = this.loadRegistrationOptions();
            if isempty(options)
                return;
            end

            for i = 1:length(options)
                if exist(options{i}.RegFile, 'file') == 0
                    % file doesn't exist, we have nothing more to do here.
                    return;
                end
            end

            isDone = true;
        end

        function saveRegistrationOptions(this, options) %#ok<INUSD>
            if this.currentEntry == 0
                warning('Suite2P:notInited', 'Storage manager is not initialized, will not save registration file');
                return;
            end

            filePath = this.registrationBinaryPath(this.db(this.currentEntry));
            fileFolder = fileparts(filePath);
            if ~exist(fileFolder, 'dir')
                mkdir(fileFolder);
            end
            save(filePath, 'options');
        end

        function opt = loadRegistrationOptions(this)
            if this.currentEntry == 0
                error('Suite2P:notInited', 'Storage manager is not initialized, can not load registration file');
            end
            opt = [];
            filePath = this.registrationBinaryPath(this.db(this.currentEntry));
            if exist(filePath, 'file') == 0
                return;
            end

            load(filePath, 'options');
            opt = options; %#ok<*PROP>
        end

        function filePath = registrationTiffPath(this, planeInd, frameInd, fileType)
            if nargin < 4
                isRedChannel = false;
            else
                if strcmpi(fileType, 'red')
                    isRedChannel = true;
                else
                    error('Unknown value of parameter ''type'' (given value is %s). See description of function registrationTiffPath.', fileType);
                end
            end
            filePath = '';
            if this.currentEntry == 0 || isempty(this.db)
                return;
            end
            dbEntry = this.db(this.currentEntry);

            [~, experimentDir] = fileparts(dbEntry.allTiffFiles{this.currentExperiment, 1});
            folder = fullfile(this.getOr('RegFileTiffLocation'), dbEntry.mouseName{1}, ...
                dbEntry.date{1}, experimentDir, sprintf('Plane%d', planeInd));
            if ~exist(folder, 'dir')
                mkdir(folder);
            end

            namePart = sprintf('%s_%s_%s_2P_plane%d_%d', dbEntry.date{1}, experimentDir, ...
                dbEntry.mouseName{1}, planeInd, frameInd);
            if isRedChannel
                namePart = [namePart '_RED'];
            end
            fileName = [namePart '.tif'];
            filePath = fullfile(folder, fileName);
        end

        function filePath = registrationFileOnFastStorage(this, planeInd, fileType)
            isRedChannel = false;
            isInterpolation = false;

            if nargin < 3
                isRedChannel = false;
                isInterpolation = false;
            else
                if strcmpi(fileType, 'red')
                    isRedChannel = true;
                elseif strcmpi(fileType, 'interp')
                    isInterpolation = true;
                else
                    error('Unknown value of parameter ''type'' (given value is %s). See description of function registrationFileOnFastStorage.', fileType);
                end
            end
            filePath = '';
            folder = this.registrationResultsFastFolder();
            if isempty(folder) && ~isInterpolation
                % return unless we have interpolation mode on. A different folder is
                % used for interpolation.
                return;
            end
            dbEntry = this.db(this.currentEntry);
            experiments_str = sprintf('%d_', dbEntry.experiments{:});
            experiments_str(end) = []; % remove last '_'
            % experiments_str at this point is CharSubDirs from original Suite2P

            if isRedChannel
                fileName = sprintf('%s_%s_%s_plane%d_RED.bin', dbEntry.mouseName{1}, ...
                    dbEntry.date{1}, experiments_str, planeInd);
            else
                fileName = sprintf('%s_%s_%s_plane%d.bin', dbEntry.mouseName{1}, ...
                    dbEntry.date{1}, experiments_str, planeInd);
            end
            interpolationFolder = this.registrationResultsSlowFolder();
            if isInterpolation
                if isempty(interpolationFolder)
                    filePath = '';
                    return;
                end
                folder = fullfile(interpolationFolder, dbEntry.mouseName{1}, ...
                    dbEntry.date{1}, experiments_str, 'interpolated');
                if ~exist(folder, 'dir')
                    mkdir(folder);
                end
            end

            if ~exist(folder, 'dir')
                mkdir(folder);
            end
            filePath = fullfile(folder, fileName);
        end

        function filePath = registrationFileOnSlowStorage(this, planeInd, fileType)
            isRedChannel = false;
            isInterpolation = false;

            if nargin < 3
                isRedChannel = false;
                isInterpolation = false;
            else
                if strcmpi(fileType, 'red')
                    isRedChannel = true;
                elseif strcmpi(fileType, 'interp')
                    isInterpolation = true;
                else
                    error('Unknown value of parameter ''type'' (given value is %s). See description of function registrationFileOnSlowStorage.', fileType);
                end
            end
            filePath = '';
            folder = this.registrationResultsSlowFolder();
            if isempty(folder)
                return;
            end
            dbEntry = this.db(this.currentEntry);
            experiments_str = sprintf('%d_', dbEntry.experiments{:});
            experiments_str(end) = []; % remove last '_'
            % experiments_str at this point is CharSubDirs from original Suite2P

            if isRedChannel
                fileName = sprintf('%s_%s_%s_plane%d_RED.bin', dbEntry.mouseName{1}, ...
                    dbEntry.date{1}, experiments_str, planeInd);
            else
                fileName = sprintf('%s_%s_%s_plane%d.bin', dbEntry.mouseName{1}, ...
                    dbEntry.date{1}, experiments_str, planeInd);
            end
            interpolationFolder = folder;
            if isInterpolation
                if isempty(interpolationFolder)
                    filePath = '';
                    return;
                end
                folder = fullfile(interpolationFolder, dbEntry.mouseName{1}, ...
                    dbEntry.date{1}, experiments_str, 'interpolated');
            end

            if ~exist(folder, 'dir')
                mkdir(folder);
            end
            filePath = fullfile(folder, fileName);
        end

        function descr = supportedTypesForPlane(~)
            descr = {};
            descr = cat(1, descr, {'f', 'Results with dF'});
            descr = cat(1, descr, {'f_proc', 'Results with dF'});
            descr = cat(1, descr, {'neuropil', 'Neuropil file'});
            descr = cat(1, descr, {'svd', 'file with computed SVD'});
            descr = cat(1, descr, {'svd_roi', 'File with computed SVD for ROI'});
            if nargout == 0
                fprintf('Supported file type per Tiff plane\n');
                for i = 1:size(descr, 1)
                    fprintf('%-20s - %s\n', descr{i, 1:2})
                end
            end
        end

        function isSupported = isTypeForPlaneSupported(this, fileType)
            descr = this.supportedTypesForPlane();
            isSupported = any(ismember(descr(:, 1), fileType));
        end

        function filePath = getFileForPlane(this, fileType, planeInd)
            if ~this.isTypeForPlaneSupported(fileType)
                error('Given file type is not supported. See function supportedTypesForPlane for a list of supported types.');
            end
            filePath = '';
            if isempty(this.db) || this.currentEntry == 0
                return;
            end
            dbEntry = this.db(this.currentEntry);

            folder = this.savePathFromRoot();
            switch lower(fileType)
                case 'f'
                    fileName = sprintf('F_%s_%s_plane%d.mat', dbEntry.mouseName{1}, dbEntry.date{1}, planeInd);

                case 'f_proc'
                    fileName = sprintf('F_%s_%s_plane%d_proc.mat', dbEntry.mouseName{1}, dbEntry.date{1}, planeInd);

                case 'svd'
                    fileName = sprintf('SVD_%s_%s_plane%d.mat', dbEntry.mouseName{1}, dbEntry.date{1}, planeInd);

                case 'svd_roi'
                    fileName = sprintf('SVDroi_%s_%s_plane%d.mat', dbEntry.mouseName{1}, dbEntry.date{1}, planeInd);
            end

            filePath = fullfile(folder, fileName);
        end

        % Return path where result files will be stored.
        % Part of the interface
        function name = savePathFromRoot(this, entryInd)
            name = '';
            if isempty(this.db)
                error('Suite2P:notInited', 'Storage manager is not initialized. Please add a db entry first.');
            end

            if nargin < 2
                entryInd = this.currentEntry;
                if this.currentEntry == 0
                    warning('Suite2P:notInited', 'Storage manager is not initialized.');
                    return;
                end
            end

            if entryInd < 1 || entryInd > length(this.db)
                error('Suite2P:arg', 'Value for argument entryInd is out of range. Possible range is [1; %d], value given is %d\n', length(this.db), entryInd);
            end
            dbEntry = this.db(entryInd);
            name = sprintf('%d_', dbEntry.experiments{:});
            name(end) = []; % remove last '_'
            % name at this point is CharSubDirs from original Suite2P

            if ~isfield(dbEntry, 'resultsSavePath') || isempty(dbEntry.resultsSavePath)
                folder = dbEntry.rootStorage;
            else
                folder = dbEntry.resultsSavePath;
            end

            name = fullfile(folder, dbEntry.mouseName{1}, ...
                dbEntry.date{1}, name);
            if ~exist(name, 'dir')
                mkdir(name);
            end
        end

    end

    methods (Access=protected)
        function checkRequiredFields(this, options, entryInd)
            % let's check that all necessary db fields are specified.
            required = {'regFilePath', 'rootStorage', 'experiments'};
            for i = 1:length(required)
                fieldName = required{i};
                [v1, assigned] = this.getOrInternal(this.db(entryInd), fieldName, []);
                [v2, assigned2] = this.getOrInternal(options, fieldName, []);
                if assigned
                    v = v1;
                else
                    v = v2;
                end
                notAssigned = ~(assigned || assigned2);
                if notAssigned || isempty(v)
                    error('Suite2P:incompleteDb', 'Can not add db entry %d (%s) as required field ''%s'' is missing. Please add it in order to continue.', entryInd, this.db(entryInd).mouseName{1}, fieldName);
                end
            end
        end

        % Initialize current entry with given options. Options are not necessary db related.
        % For more details see initEntry in BaseStorageManager.
        function [inited, info] = initEntry(this, options, entryInd)
            dbEntry = this.db(entryInd);
            %tifPaths = {};
            % cell array that contains list of cell files. Has
            % the same length as tifPaths. Each
            % Resembles fsroot from original Suite2P
            allTiffFiles = {};
            entryNumFiles = 0;

            dbEntry = this.initInheritedValues(dbEntry, {'rootStorage', 'RegFileBinLocation', 'RegFileTiffLocation', ...
                'temp_tiff'});
            this.checkRequiredFields(options, entryInd);

            rootStorage = dbEntry.rootStorage;
            tiffFiles = dir(fullfile(rootStorage, '*.tif*'));
            directPath = ~isempty(tiffFiles);

            nSessions = length(dbEntry.mouseName);
            for i = 1:nSessions
                curExperiments = dbEntry.experiments{i};
                if isempty(curExperiments) && directPath
                    curExperiments = [0]; % the value is not important
                end

                %curSubDirs
                for iExp = 1:length(curExperiments)
                    experiment_string = num2str(curExperiments(iExp));
                    sessionRelativePath = fullfile(dbEntry.mouseName{i}, dbEntry.date{i}, experiment_string);
                    if directPath
                        tiffFolder = rootStorage;
                    else
                        tiffFolder = fullfile(rootStorage, sessionRelativePath);
                    end

                    % Find both .tif and .tiff files
                    files1 = dir(fullfile(tiffFolder, '*.tif'));
                    files2 = dir(fullfile(tiffFolder, '*.tiff'));
                    tiffFiles = cat(1, files1, files2);
                    allFileSizes = [tiffFiles(:).bytes];

                    [tiffFiles, sortIndex] = sort_nat({tiffFiles(:).name});
                    allFileSizes = allFileSizes(sortIndex);
                    tiffFiles = cellfun(@(x) fullfile(tiffFolder, x), tiffFiles, 'uniformoutput', false);
                    allTiffFiles = cat(1, allTiffFiles, {tiffFolder, tiffFiles, allFileSizes});
                    entryNumFiles = entryNumFiles + length(tiffFiles);
                end
            end
            this.numFiles(entryInd) = entryNumFiles;

            options.nplanes = this.getOr('nplanes', 1);
            options.nchannels = this.getOr('nchannels', 1);
            options.readTiffHeader = this.getOr('readTiffHeader', true);

            if options.readTiffHeader
                try
                    [~, header] = loadFramesBuff(allTiffFiles{1, 2}(1), 1, 1, 1);
                    hh = header{1};
                    verStr = 'SI.VERSION_MAJOR = ''2016b''';
                    if contains(hh, verStr) % For scanImage 2016b
                        str = hh(strfind(hh, 'channelSave = '):end);
                        ind = strfind(str, 'SI');
                        ch = str2double(str(15:ind(1)-1));
                        options.nchannels = length(ch);

                        fastZEnable = sscanf(hh(strfind(hh, 'hFastZ.enable = '):end), 'hFastZ.enable = %s');
                        fastZEnable = strcmpi(fastZEnable, 'true');
                        fastZDiscardFlybackFrames = sscanf(hh(strfind(hh, 'hFastZ.discardFlybackFrames = '):end), 'hFastZ.discardFlybackFrames = %s');
                        fastZDiscardFlybackFrames = strcmpi(fastZDiscardFlybackFrames, 'true');
                        stackNumSlices = sscanf(hh(strfind(hh, 'hStackManager.numSlices = '):end), 'hStackManager.numSlices = %d');

                        if fastZEnable
                            options.nplanes = stackNumSlices + fastZDiscardFlybackFrames;
                        end

                        str = hh(strfind(hh, 'scanZoomFactor = '):end);
                        ind = strfind(str, 'SI');
                        options.zoomMicro = str2double(str(18 : ind(1)-1));
                        options.imageRate = sscanf(hh(strfind(hh, 'scanFrameRate = '):end), 'scanFrameRate = %f');

                        % if isfield(db, 'expred') && ~isempty(db.expred) && ...
                        %         (~isfield(db, 'nchannels_red') || isempty(db.nchannels_red))
                        %     [~, header] = loadFramesBuff(options.fsred(1).name, 1, 1, 1);
                        %     hh=header{1};
                        %     str = hh(strfind(hh, 'channelSave = '):end);
                        %     ind = strfind(str, 'SI');
                        %     ch = str2num(str(15 : ind(1)-1));
                        %     options.nchannels_red = length(ch);
                        % end
                    end

                    % Old scanimage
                    str = hh(strfind(hh, 'channelsSave = '):end);
                    ind = strfind(str, 'scanimage');
                    ch = str2double(str(16 : ind(1)-1));
                    options.nchannels = length(ch);

                    fastZEnable = sscanf(hh(strfind(hh, 'fastZEnable = '):end), 'fastZEnable = %d');
                    fastZDiscardFlybackFrames = sscanf(hh(strfind(hh, 'fastZDiscardFlybackFrames = '):end), 'fastZDiscardFlybackFrames = %d');
                    if isempty(fastZDiscardFlybackFrames)
                        fastZDiscardFlybackFrames = 0;
                    end
                    stackNumSlices = sscanf(hh(strfind(hh, 'stackNumSlices = '):end), 'stackNumSlices = %d');

                    options.nplanes = 1;
                    if fastZEnable
                        options.nplanes = stackNumSlices+fastZDiscardFlybackFrames;
                    end

                    str = hh(strfind(hh, 'scanZoomFactor = '):end);
                    ind = strfind(str, 'scanimage');
                    options.zoomMicro = str2double(str(18 : ind(1)-1));

                    options.imageRate = sscanf(hh(strfind(hh, 'scanFrameRate = '):end), 'scanFrameRate = %f');
                catch
                    fprintf('There was an exception during loading of tiff file header. Probably your files do not have headers.\nFile was: %s\n', allTiffFiles{1, 2}(1));
                end
            end % end read tiff header

            if isfield(options, 'zoomMicro')
                options.zoom = getOr(options, 'zoom', options.zoomMicro);
            else
                options.zoom = getOr(options, 'zoom', 1);
            end

            options.planesToProcess = 1:options.nplanes;

            options.resultsSavePath = this.savePathFromRoot(entryInd);

            info = options;
            this.options = options;

            this.db(entryInd).allTiffFiles = allTiffFiles;
            inited = ~isempty(allTiffFiles);
        end
    end

    methods (Access=private)
        function folder = registrationResultsFastFolder(this)
            %folder = this.getOr('RegFileRoot', '');
            folder = this.getOr('regFilePath', '');
        end

        function folder = registrationResultsSlowFolder(this)
            folder = this.getOr('RegFileBinLocation', '');
        end

        function filePath = registrationBinaryPath(this, options)
            filePath = '';
            fileName = this.registrationBinaryName(options);
            folder = this.savePathFromRoot();
            if isempty(folder)
                return;
            end
            filePath = fullfile(folder, fileName);
            filePath = strrep(filePath, '\', '/');
        end

        function name = registrationBinaryName(~, options)
            name = sprintf('regops_%s_%s.mat', options.mouseName{1}, options.date{1});
        end
    end
end
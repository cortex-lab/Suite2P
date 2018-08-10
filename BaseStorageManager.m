% An interface for providing access to raw files
%
% This interfaces hides details of storage schema. It is
% provided as an abstract class. Implementation classes
% are free to define different file schema.
%
% The interface is heavily influenced by original Suite2P
% code, which supported a single file storage and didn't
% abstract things. BaseStorageManager does not abstract all the things
% that are possible to abstract. For example, it uses Suite2P format
% for database entries. This is done in order to minimize the difference
% between BaseStorageManager and original Suite2P.
%
% BaseStorageManager works with file sets, which can be organized
% in partitions. In original Suite2P terms, partition is experiment.
% Location of files and additional information is provided by db
% entries, which are Matlab structures. BaseStorageManager uses
% the following fields in each db entry:
% *) mouseName      - cell array of strings, name/id of the animal.
% *) date           - cell array of strings that represent dates, 'yyyy-mm-dd'.
% *) experiments    - cell array of integer arrays. Must be provided!
% *) expred         - array, points to partitions that contain red channel.
% *) nchannels_red  - how many channels did the red block have in total (assumes red is last).
% *) nchannels      - integer, number of channels. Default is 1.
%
% Please see the documentation of CortexLabStorageManager class for details of how a dbEntry
% can be defined and added.
%
% General usage of BaseStorageManager is:
% 1. Create an instance of BaseStorageManager implementation class.
% 2. Describe things that are needed to be analyzed. This is
%    done by adding 'database entries' to class. Each entry
%    is a Matlab's structure. Implementation classes can use different
%    structures.
% 3. Query BaseStorageManager for data. Queries a performed by calling
%    method getFile().
% 4. When needed, use methods that construct file paths. Example of such
%    method is registrationTiffPath().
%
classdef (Abstract) BaseStorageManager < handle
    properties (Access=public)
        options = struct()
    end

    properties (Access=protected)
        % DB structure. It is used to store information
        % about files to be analyzed.
        db = []
        % Index of a db entry that is currently in use.
        currentEntry = 0
        % Index of current tiff file to read from current db entry
        currentFile = 0
        % Index of current experiment/partition for current db entry
        currentExperiment = 0

        % vector of total number of files per db entry.
        % Index is db entry. Value - number of files.
        numFiles = []
    end

    methods (Access=public)
        % Construct BaseStorageManager and initialized with options.
        % Options is a structure with Suite2P options.
        function this = BaseStorageManager(options)
            this.options = options;
        end

        % Return number of channels for current partition in the current db entry.
        function numChannels = getNumChannels(this)
            numChannels = 0;
            if ~this.isEntryInited
                return;
            end
            numChannels = this.getOr('nchannels', 1);
            redExperiments = this.getOr('expred', []);
            isExperimentRed = ismember(this.currentExperiment, redExperiments);
            if isExperimentRed
                numChannels = this.getOr('nchannels_red');
            end
        end

        % Add new entry to storage manager.
        function addEntry(this, newEntry)
            if ~isfield(newEntry, 'mouseName')
                error('Your DB entry has no field mouseName. Please add it.');
            end
            % make mouseName a cell array
            if ~iscell(newEntry.mouseName)
                newEntry.mouseName = {newEntry.mouseName};
            end

            % make experiments a cell array
            if ~iscell(newEntry.experiments)
                newEntry.experiments = {newEntry.experiments};
            end

            if ~iscell(newEntry.date)
                newEntry.date = {newEntry.date};
            end

            if isempty(this.db)
                this.db = newEntry;
            else
                this.db(end+1) = newEntry;
            end
            [inited, newOptions] = this.initEntry(this.options, length(this.db));
            this.options = newOptions;
            if this.currentEntry == 0
                if inited
                    this.currentFile = 1;
                    this.currentExperiment = 1;
                    this.currentEntry = length(this.db);
                else
                    this.currentFile = 0;
                    this.currentExperiment = 0;
                end
            end
        end

        % Return number of entries in storage manager.
        function num = getNumEntries(this)
            num = length(this.db);
        end

        % Remove all entries from storage manager.
        function removeAllEntries(this)
            this.db = [];
            this.currentEntry = 0;
            this.currentFile = 0;
            this.currentExperiment = 0;
        end

        % Select entry with index entryInd as the current entry.
        % Functions that access data interact with the current entry.
        % This function allows to select current entry.
        function selectEntry(this, entryInd)
            if isempty(this.db)
                error('Suite2P:notInited', 'Storage manager is not initialized, add an entry before calling selectEntry method.');
            end

            if entryInd < 1 || entryInd > length(this.db)
                warning('New entry number is out of range. Should be in [1; %d], %d is given.', ...
                    length(this.db), entryInd);
                return;
            end
            this.currentEntry = entryInd;
        end

        % Reset current entry to state as if it had been just added to storage manager.
        % After calling resetCurrentEntry, hasData() method will return logical 1 (true).
        % This is useful if you want to iterate over the data in the current entry several times.
        function resetCurrentEntry(this)
            if this.currentEntry == 0
                return;
            end
            if this.numFiles(this.currentEntry) > 0
                this.currentFile = 1;
                this.currentExperiment = 1;
            end
        end

        %% The following methods are data-related

        % Determine if data is available to read.
        % Returns logical 1 (true) if there is data available
        % to read from the current entry. If partition argument
        % is specified, then returns logical 1 if that partition
        % has data to read.
        function has = hasData(this, partition)
            has = false;
            if ~this.isEntryInited
                return;
            end
            if this.numFiles(this.currentEntry) == 0
                return;
            end

            if nargin < 2
                partition = [];
            end
            if ~isempty(partition)
                if partition < 0 || partition > size(this.db(this.currentEntry).allTiffFiles, 1)
                    % partition value is invalid
                    return
                end
                if partition ~= this.currentExperiment
                    this.currentExperiment = partition;
                    this.currentFile = 1;
                end

                if isempty(this.db(this.currentEntry).allTiffFiles{partition, 2})
                    return;
                end
            end

            has = true;
        end

        % Return the very first file name for current trial.
        % This method does not change the internal state of storage manager,
        % i.e. multiple calls to it will return the same information.
        %  USAGE
        %   [name, info] = getFirstFile(sm, partition)
        %   sm          - instance of storage manager class
        %   partition   - positive integer. If provided,
        %                 then returns the first file from that partition.
        %                 Default is 1.
        %   name        - full file path.
        %   info        - file information, structure with fields:
        %                 bytes     - file size in bytes.
        %                 partition - partition index.
        %
        function [name, info] = getFirstFile(this, partition)
            name = '';
            if ~this.isEntryInited
                return;
            end

            if nargin < 2
                partition = 1;
            end
            if partition < 0 || partition > size(this.db(this.currentEntry).allTiffFiles, 1)
                partition = 1;
            end

            name = this.db(this.currentEntry).allTiffFiles{partition, 2}{1};
            if nargout > 1
                info.bytes = this.db(this.currentEntry).allTiffFiles{partition, 3}(1);
                info.partition = partition;
            end
        end

        % Get file from current entry of storage manager.
        % Returns file information from the current entry of storage manager.
        % Subsequent calls to the getFile function continue getting files from
        % the endpoint of the previous call.
        % See method getFirstFile for the description of name and info variables.
        function [name, info] = getFile(this)
            name = '';
            info = [];
            if ~this.isEntryInited
                return;
            end

            name = this.db(this.currentEntry).allTiffFiles{this.currentExperiment, 2}{this.currentFile};
            if nargout > 1
                info.bytes = this.db(this.currentEntry).allTiffFiles{this.currentExperiment, 3}(this.currentFile);
                info.partition = this.currentExperiment;
            end
            this.increaseFileCounter();
        end

        % Get next file relative to the current file (last read file).
        % This function does not change the internal state of storage manager,
        % i.e. multiple/subsequent calls to getNextFile return the same information.
        % See method getFirstFile for the description of name and info variables.
        function [name, info] = getNextFile(this)
            name = '';
            info = [];
            if ~this.isEntryInited
                return;
            end

            numFiles = length(this.db(this.currentEntry).allTiffFiles{this.currentExperiment, 2}); %#ok<PROP>
            if this.currentFile >= numFiles %#ok<PROP>
                return;
            end

            fileInd = this.currentFile + 1;
            name = this.db(this.currentEntry).allTiffFiles{this.currentExperiment, 2}{fileInd};
            if nargout > 1
                info.bytes = this.db(this.currentEntry).allTiffFiles{this.currentExperiment, 3}(fileInd);
                info.partition = this.currentExperiment;
            end
        end

        % Return number of files in current db entry or partition.
        % If partition argument is omitted, then returns the total number of files
        % for the current entry. If partition is provided, then returns number
        % of files for that partition.
        function num = getNumFiles(this, partition)
            num = 0;
            if this.currentEntry == 0
                return;
            end
            if nargin == 2
                num = length(this.db(this.currentEntry).allTiffFiles{partition, 2});
            else
                num = this.numFiles(this.currentEntry);
                if isempty(num)
                    num = 0;
                end
            end
        end

        % Return number of partitions in current entry.
        function num = getNumPartitions(this)
            num = 0;
            if this.currentEntry == 0 || isempty(this.db)
                return;
            end
            num = size(this.db(this.currentEntry).allTiffFiles, 1);
        end

        % Prepare partition for reading.
        % Subsequent calls to the getFile function will return file
        % from that partition.
        function seekToPartition(this, partition)
            if this.currentEntry == 0 || isempty(this.db)
                return;
            end
            if partition < 0 || partition > size(this.db(this.currentEntry).allTiffFiles, 1)
                warning('No valid partition %u, will fall back to the first one.', partition);
                partition = 1;
            end
            this.currentExperiment = partition;
            this.currentFile = 1;
        end

        % Return attribute or default value of current db entry or options.
        % The function first checks if attribute with name field exists
        % for the current db entry. If it does, then getOr() returns it's value.
        % If it does not, then getOr() checks if field exists in storage
        % manager's options.
        %
        %  USAGE
        %   v = getOr(sm, field, value)
        %   sm          Instance of storage manager
        %   field       string, attribute name
        %   value       If attribute is not present, then value
        %               is returned.
        %
        function v = getOr(this, field, defaultValue)
            if nargin < 3
                defaultValue = [];
            end
            if this.currentEntry > 0
                s = this.db(this.currentEntry);
                [v, assigned] = this.getOrInternal(s, field, defaultValue);
                if assigned
                    return;
                end
            end
            v = this.getOrInternal(this.options, field, defaultValue);
        end
    end

    methods (Abstract)
        % Determine if registration is done for the current entry.
        tf = isRegistrationDone(this)

        % Save registration options to disc.
        saveRegistrationOptions(this, options)

        % Load registration options from disc.
        options = loadRegistrationOptions(this)

        %%

        % Return relative path from root where result files will be stored.
        % If entryInd is not provided, then the current db entry is used.
        % Root is a user-provided path where results should be stored.
        folder = savePathFromRoot(this, entryInd)

        % Return full file path to registration file for given plane and frame.
        % Optional argument fileType can be set to 'red' in order to get
        % file path for red channel.
        filePath = registrationTiffPath(this, planeInd, frameInd, fileType)

        % Return full file path to registration binary file on fast storage
        % for a given plane.
        % Fast storage is typically a local SSD. SSD have limited capacity and
        % not suited for long-term storage of large binary files.
        %  USAGE
        %   filePath = registrationFileOnFastStorage(sm, planeInd, fileType)
        %   sm          instance of storage manager.
        %   planeInd    index of a plane for which file path is returned.
        %   fileType    optional type of returned file, string.
        %               Possible values are:
        %               1. 'red' indicates that filePath contains path
        %               to red channel.
        %               2. 'interp' indicates that filePath contains path
        %               to files with results of interpolation across frames.
        filePath = registrationFileOnFastStorage(this, planeInd, fileType)

        % Return full file path to registration binary file on slow storage
        % for a given plane.
        % Slow storage is a place where files get copied from fast storage.
        %  USAGE
        %   filePath = registrationFileOnSlowStorage(sm, planeInd, fileType)
        %   sm          instance of storage manager.
        %   planeInd    index of a plane for which file path is returned.
        %   fileType    optional type of returned file, string.
        %               Possible values are:
        %               1. 'red' indicates that filePath contains path
        %               to red channel.
        %               2. 'interp' indicates that filePath contains path
        %               to files with results of interpolation across frames.
        filePath = registrationFileOnSlowStorage(this, planeInd, fileType)

        % Return full file path for file of type fileType for plane with index planeInd.
        % Different files can be produced per Tiff plane.
        %  USAGE
        %   filePath = getFileForPlane(this, fileType, planeInd)
        %   sm          instance of storage manager.
        %   fileType    One of the supported file types. For a complete list
        %               see method supportedTypesForPlane.
        %   planeInd    Index of a plane.
        %   filePath    Full file path.
        %
        filePath = getFileForPlane(this, fileType, planeInd)

        % Return description of supported file types per plane.
        % Returns a Nx2 cell array. First column contains variables
        % that can be passed to the function registrationFileOnSlowStorage as
        % values for fileType argument. Second column contains human description.
        % Example
        % 1x2 cell array
        %   {'svd'} {'File with computed SVD.'}
        %
        % If output variable is omitted, then the function prints out the information.
        % Minimum supported types are: 'f', 'f_proc', 'neuropil', 'svd', 'svd_roi'.
        descr = supportedTypesForPlane(this)

        % Determines if provided fielType is supported as a valid file type for a Tiff plane.
        tf = isTypeForPlaneSupported(this, fileType)
    end

    methods (Access=protected, Abstract=true)
        % Initialize new db entry when it is added to storage manager.
        % Arguments:
        % options       storage manager options.
        % entryInd      index of an entry that should be initialized.
        % Output:
        % inited        boolean flag indicating if entry has been initialized successfully.
        % info          New version of options. Initialization function may change values or add
        %               new fields to structure options. In order to keep track of it, a modified
        %               version is returned via info.
        %
        [inited, info] = initEntry(this, options, entryInd)
    end

    methods (Access=protected)
        % Increase internal file pointer.
        function increaseFileCounter(this)
            if ~this.isEntryInited
                return;
            end
            this.currentFile = this.currentFile + 1;
            if this.currentFile > length(this.db(this.currentEntry).allTiffFiles{this.currentExperiment, 2})
                if this.currentExperiment == size(this.db(this.currentEntry).allTiffFiles, 1)
                    this.currentExperiment = 0;
                    this.currentFile = 0;
                    return;
                end

                this.currentExperiment = this.currentExperiment + 1;
                this.currentFile = 1;
            end
        end

        % Determine if current db entry is initialized.
        % Returns logical 1 (true) if entry is initialized.
        function res = isEntryInited(this)
            res = ~(this.currentEntry == 0 || this.currentFile == 0 || isempty(this.db));
        end

        function [v, assigned] = getOrInternal(~, s, field, default)
            if nargin < 3
                default = [];
            end
            assigned = false;

            fieldExists = isfield(s, field);
            if any(fieldExists)
                assigned = true;
                if iscellstr(field) || isstring(field)
                    v = s.(field{find(fieldExists, 1)});
                else
                    v = s.(field);
                end
            else
                v = default;
            end
        end

        % Initialize values that could be inherited from a previous db entry.
        % fieldNames is a cellArray of field names, which have their values
        % inherited from previous db entry.
        function dbEntry = initInheritedValues(this, dbEntry, fieldNames)
            if ~iscell(fieldNames)
                fieldNames = {fieldNames};
            end
            for i = 1:length(fieldNames)
                fieldName = fieldNames{i};
                if ~isfield(dbEntry, fieldName)
                    if length(this.db) > 1 && this.currentEntry > 1
                        if isfield(this.db(this.currentEntry-1), fieldName)
                            % inherit from previous record
                            this.db(this.currentEntry).(fieldName) = this.db(this.currentEntry-1).(fieldName);
                            dbEntry.(fieldName) = this.db(this.currentEntry).(fieldName);
                        else
                            error('Value for field %s if not specified, animal %s', fieldName, dbEntry.mouseName{1});
                        end
                    end
                end
            end
        end

    end

end
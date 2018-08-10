function options = build_ops3(db, options)
    options.nplanes = getOr(options, 'nplanes', 1);
    options.nchannels = getOr(options, 'nchannels', 1);
    options.readTiffHeader = getOr(options,'readTiffHeader',1);
    if isfield(db, 'expred')
        if ~isempty(db.expred)
            options.nchannels_red = getOr(options, 'nchannels_red', 2);
        end
    end

    % options = db;
    if ~iscell(db.mouse_name)
        % this is the usual case where we have a simple single session recording
        options = addfields(options, db);

        for k = 1:length(db.expts)
            options.SubDirs{k}    = num2str(db.expts(k));
        end
        if isempty(db.expts)
            options.SubDirs{1} = [];
        end

        if ~isfield(options, 'RootDir')
            options.RootDir = fullfile(options.RootStorage, options.mouse_name, options.date);
        end

        % build file list
        options.fsroot = [];
        for j = 1:length(options.SubDirs)
            ffile = dir(fullfile(options.RootDir, options.SubDirs{j}, '*.tif'));
            fname = struct2cell(ffile);
            fname = fname(1,:)';
            [~,index] = sort_nat(fname);
            options.fsroot{j} = ffile(index);
            options.fsroot{j} = cat(1, options.fsroot{j}, ...
                dir(fullfile(options.RootDir, options.SubDirs{j}, '*.tiff')));
            for k = 1:length(options.fsroot{j})
                options.fsroot{j}(k).name = fullfile(options.RootDir, options.SubDirs{j}, options.fsroot{j}(k).name);
            end
        end

        if isfield(db, 'expred') && ~isempty(db.expred) && ...
                (~isfield(db, 'nchannels_red') || isempty(db.nchannels_red))
            ffile = dir(fullfile(options.RootDir, num2str(db.expred), '*.tif'));
            fname = struct2cell(ffile);
            fname = fname(1,:)';
            [~,index] = sort_nat(fname);
            options.fsred = ffile(index);
            for k = 1:length(options.fsroot{j})
                options.fsred(k).name = fullfile(options.RootDir, num2str(db.expred), options.fsred(k).name);
            end
        end
    else
        % here we might have multiple sessions, which we want to be analyzed
        % together (exactly the same FOV)
        nSessions = length(db.mouse_name);
        % a backwards compatible version of db
        dbCompat = db;
        dbCompat.mouse_name = db.mouse_name{1};
        dbCompat.date = db.date{1};
        dbCompat.expts = cell2mat(db.expts(:)');
        options = addfields(options, dbCompat);
        options.db_orig = db;

        options.fsroot = cell(0);
        options.SubDirs = cell(0);
        for iSession = 1:nSessions
            options.RootDir = fullfile(options.RootStorage, db.mouse_name{iSession}, db.date{iSession});
            for iExp = 1:length(db.expts{iSession})
                options.SubDirs{end+1} = num2str(db.expts{iSession}(iExp));
                options.fsroot{end+1} = dir(fullfile(options.RootDir, options.SubDirs{end}, '*.tif'));
                for iFile = 1:length(options.fsroot{end})
                    options.fsroot{end}(iFile).name = fullfile(options.RootDir, options.SubDirs{end}, options.fsroot{end}(iFile).name);
                end
            end
        end
        % this line to be backward compatible (just in case)
        options.RootDir = fullfile(options.RootStorage, options.mouse_name, options.date);
    end

    if options.readTiffHeader
        try
            % MK code for automatically determining number of planes and channels
            [~, header] = loadFramesBuff(options.fsroot{1}(1).name, 1, 1, 1);

            hh=header{1};

            verStr = ['SI.VERSION_MAJOR = ',char(39),'2016b',char(39)];

            if contains(hh, verStr) % For scanimage 2016b, SF
                str = hh(strfind(hh,'channelSave = '):end);
                ind = strfind(str, 'SI');
                ch = str2num(str(15 : ind(1)-1));
                options.nchannels = length(ch);

                fastZEnable = sscanf(hh(strfind(hh,'hFastZ.enable = '):end), 'hFastZ.enable = %s');
                fastZEnable = strcmp(fastZEnable,'true');
                fastZDiscardFlybackFrames = sscanf(hh(strfind(hh, 'hFastZ.discardFlybackFrames = '):end), 'hFastZ.discardFlybackFrames = %s');
                fastZDiscardFlybackFrames = strcmp(fastZDiscardFlybackFrames,'true');
                stackNumSlices = sscanf(hh(strfind(hh, 'hStackManager.numSlices = '):end), 'hStackManager.numSlices = %d');

                options.nplanes = 1;

                if fastZEnable
                    options.nplanes = stackNumSlices+fastZDiscardFlybackFrames;
                end

                str = hh(strfind(hh, 'scanZoomFactor = '):end);
                ind = strfind(str, 'SI');
                options.zoomMicro = str2double(str(18 : ind(1)-1));

                options.imageRate = sscanf(hh(strfind(hh, 'scanFrameRate = '):end), 'scanFrameRate = %f');

                if isfield(db, 'expred') && ~isempty(db.expred) && ...
                        (~isfield(db, 'nchannels_red') || isempty(db.nchannels_red))
                    [~, header] = loadFramesBuff(options.fsred(1).name, 1, 1, 1);
                    hh=header{1};
                    str = hh(strfind(hh, 'channelSave = '):end);
                    ind = strfind(str, 'SI');
                    ch = str2num(str(15 : ind(1)-1));
                    options.nchannels_red = length(ch);
                end

            end

            % Old scanimage

            str = hh(strfind(hh, 'channelsSave = '):end);
            ind = strfind(str, 'scanimage');
            ch = str2num(str(16 : ind(1)-1));
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

            % get number of channels of red experiment
            if isfield(db, 'expred') && ~isempty(db.expred) && ...
                    (~isfield(db, 'nchannels_red') || isempty(db.nchannels_red))
                [~, header] = loadFramesBuff(options.fsred(1).name, 1, 1, 1);
                hh=header{1};
                str = hh(strfind(hh, 'channelsSave = '):end);
                ind = strfind(str, 'scanimage');
                ch = str2num(str(16 : ind(1)-1));
                options.nchannels_red = length(ch);
            end
        catch
        end
    end
    if isfield(options, 'zoomMicro')
        options.zoom = getOr(options, 'zoom', options.zoomMicro);
    else
        options.zoom = getOr(options, 'zoom', 1);
    end


    if ~(isfield(options, 'planesToProcess') && ~isempty(options.planesToProcess))
        options.planesToProcess = 1:options.nplanes;
    else
        % planesToProcess is not working right now
        options.planesToProcess = 1:options.nplanes;
    end

    CharSubDirs = '';
    for i = 1:length(options.SubDirs)
        CharSubDirs = [CharSubDirs options.SubDirs{i} '_'];
    end
    CharSubDirs = CharSubDirs(1:end-1);
    options.CharSubDirs = CharSubDirs;

    options.ResultsSavePath = sprintf('%s//%s//%s//%s//', options.ResultsSavePath, options.mouse_name, options.date, ...
        CharSubDirs);

    if isempty(db.expts)
       options.expts = [];
    end

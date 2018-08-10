%% Unit-tests for CortextLabStorageManager
classdef CortexLabStorageManagerTest < matlab.unittest.TestCase

    properties (Access=protected)

        % Folder where data for this test is stored
        % If empty, then some sample files will be downloaded
        % from Google Drive. You can download the files yourself
        % and then provide a path to the files. DataFolder
        % should point to root folder, i.e. without
        % animal/date/experiment information.
        % Alternatively to setting the path in this file, you may set environment
        % variable SUITE2P_TEST_DATA.
        % If borh SUITE2P_TEST_DATA and DataFolder have values, then DataFolder has priority and
        % will be used.
        DataFolder

        Options % initial static options
    end

    properties (Access=private)
        downloadFiles = true
    end

    methods
        function this = CortexLabStorageManagerTest()
            testDataDir = getenv('SUITE2P_TEST_DATA');
            if isempty(this.DataFolder) && isempty(testDataDir)
                this.DataFolder = fullfile(tempdir, 'Suite2P-tests-CortexLabStorageManagerTest');
                if ~exist(this.DataFolder, 'dir')
                    mkdir(this.DataFolder);
                end
            else
                if ~isempty(testDataDir)
                    candidate = testDataDir;
                end
                if ~isempty(this.DataFolder)
                    candidate = this.DataFolder;
                end
                this.DataFolder = candidate;

                if exist(this.DataFolder, 'dir') ~= 0
                    this.downloadFiles = false;
                end
            end
        end
    end

    methods (TestClassSetup)
        function prepareData(this)
            if ~this.downloadFiles
                return;
            end

            % Make DB structure
            outputFolder = fullfile(this.DataFolder, 'M0', '2017-10-13', '4');
            if ~exist(outputFolder, 'dir')
                mkdir(outputFolder);
            end

            % download couple of files from GDrive
            fprintf('Downloading files. This takes time...\n');
            files = {'file_00002_00001.tif', '0B649boZqpYG1Mk9UM1FCZ0FfaFE';...
                'file_00002_00002.tif', '0B649boZqpYG1OEZnV21ncDVNcVk';...
                'file_00002_00003.tif', '0B649boZqpYG1Z3FkdmZYaEZlTDg'; ...
                'file_00002_00004.tif', '0B649boZqpYG1MzI5Y1g5UXZneGs'; ...
                'file_00002_00005.tif', '0B649boZqpYG1MW9ubURmTk9SQXM'; ...
                };
            options = weboptions('ContentType', 'raw');

            for i = 1:size(files, 1)
                fileName = files{i, 1};
                fileUrl = sprintf('https://drive.google.com/uc?export=download&id=%s', files{i, 2});

                fprintf('File %d/%d...', i, size(files, 1));

                request = matlab.net.http.RequestMessage();
                % First request will be redirected to information page about virus scanning
                % We can get a confirmation code from an associated cookie file
                [~, infos] = sendRequest(matlab.net.URI(fileUrl), request);
                confirmCode = '';
                for j = 1:length(infos('drive.google.com').cookies)
                    if ~isempty(strfind(infos('drive.google.com').cookies(j).Name, 'download'))
                        confirmCode = infos('drive.google.com').cookies(j).Value;
                        break;
                    end
                end
                newUrl = strcat(fileUrl, sprintf('&confirm=%s', confirmCode));
                % We now need to send another request to get the file.
                % However, Matlab doesn't download the whole Tiff file, but
                % only one frame.
                [~, ~, history] = sendRequest(matlab.net.URI(newUrl), request);

                % Thus we must use log file information to find out a
                % direct link and downalod it as raw file
                ind = arrayfun(@(x) ~isempty(strfind(x.URI.Host, 'googleusercontent')), history);
                ind = find(ind, 1);

                imgData = webread(history(ind).URI.EncodedURI, options);
                fid = fopen(fullfile(outputFolder, fileName), 'wb');
                fwrite(fid, imgData);
                fclose(fid);

                fprintf('done\n');
            end
        end
    end

    methods (TestClassTeardown)
        function cleanDataFolder(testCase)
        end
    end

    methods(TestMethodSetup)
        function initializeOptions(testCase)
            testCase.Options = CortexLabStorageManagerTest.getDefaultOptions();
        end
    end

    methods(TestMethodTeardown)
        function removeOptions(testCase)
            testCase.Options = [];
        end
    end

    methods (Test)
        function testConstruction(testCase)
            sm = CortexLabStorageManager(testCase.Options);
        end

        function testDefaultDataset(testCase)
            sm = CortexLabStorageManager(testCase.Options);

            dbEntry.mouseName = 'M0';
            dbEntry.date = '2017-10-13';
            dbEntry.experiments = [4];
            dbEntry.diameter = 16;
            dbEntry.nplanes = 1;
            dbEntry.resultsSavePath = strcat(testCase.DataFolder, '-results');
            dbEntry.rootStorage = testCase.DataFolder;
            dbEntry.regFilePath = dbEntry.resultsSavePath;
            dbEntry.RegFileBinLocation = strcat(testCase.DataFolder, '-out-registration-slow');
            dbEntry.RegFileTiffLocation = '';
            dbEntry.temp_tiff = tempname;

            sm.addEntry(dbEntry);

            testCase.verifyEqual(sm.getNumEntries(), 1);

            sm.selectEntry(1);
            run_pipeline(sm);
            add_deconvolution(sm);

            % now let's check that we have files in proper places
            binFile = 'M0_2017-10-13_4_plane1.bin';
            binExist = exist(fullfile(dbEntry.RegFileBinLocation, binFile), 'file') ~= 0;
            testCase.verifyTrue(binExist);

            % check bin is still on 'fast' storage
            binExist = exist(fullfile(dbEntry.regFilePath, binFile), 'file') ~= 0;
            testCase.verifyTrue(binExist);

            % check output directory structure
            trueOutputFolder = fullfile(dbEntry.resultsSavePath, 'M0', '2017-10-13', '4');
            testCase.verifyTrue(exist(fullfile(dbEntry.resultsSavePath, 'M0'), 'dir') ~= 0);
            testCase.verifyTrue(exist(fullfile(dbEntry.resultsSavePath, 'M0', '2017-10-13'), 'dir') ~= 0);
            testCase.verifyTrue(exist(trueOutputFolder, 'dir') ~= 0);

            % check results files are in place
            testCase.verifyTrue(exist(fullfile(trueOutputFolder, 'F_M0_2017-10-13_plane1.mat'), 'file') ~= 0);
            testCase.verifyTrue(exist(fullfile(trueOutputFolder, 'regops_M0_2017-10-13.mat'), 'file') ~= 0);

            % delete files
            rmdir(dbEntry.resultsSavePath, 's');
            rmdir(dbEntry.RegFileBinLocation, 's');
            delete(dbEntry.temp_tiff);
        end

        function emptyRequiredFields(testCase)
            sm = CortexLabStorageManager(testCase.Options);

            dbEntry.mouseName = 'M0';
            dbEntry.date = '2017-10-13';
            dbEntry.experiments = 4;
            dbEntry.regFilePath = [];
            dbEntry.rootStorage = testCase.DataFolder;

            testCase.verifyError(@()sm.addEntry(dbEntry), 'Suite2P:incompleteDb');
        end

        function missingFields(testCase)
            sm = CortexLabStorageManager(testCase.Options);

            dbEntry.mouseName = 'M0';
            dbEntry.date = '2017-10-13';
            dbEntry.experiments = 4;
            dbEntry.rootStorage = testCase.DataFolder;

            testCase.verifyError(@()sm.addEntry(dbEntry), 'Suite2P:incompleteDb');
        end

        function multipleMouseNames(testCase)
            sm = CortexLabStorageManager(testCase.Options);

            dbEntry.mouseName    = {'MK020', 'M150416_MK020'};
            dbEntry.date          = {'2015-07-30', '2015-07-30'};
            dbEntry.experiments         = {[2010 2107], [1 2 3]};
            dbEntry.diameter      = 12;
            % empty means the folder with data
            dbEntry.resultsSavePath = [];
            dbEntry.rootStorage = fullfile(testCase.DataFolder, 'multiMouse');
            dbEntry.regFilePath = fullfile(dbEntry.rootStorage, 'regFilePath');

            %% prepare files
            mk20_2010 = fullfile(testCase.DataFolder, 'multiMouse', 'MK020', '2015-07-30', '2010');
            mk20_2107 = fullfile(testCase.DataFolder, 'multiMouse', 'MK020', '2015-07-30', '2107');
            mkdir(mk20_2010);
            mkdir(mk20_2107);

            m_1 = fullfile(testCase.DataFolder, 'multiMouse', 'M150416_MK020', '2015-07-30', '1');
            m_2 = fullfile(testCase.DataFolder, 'multiMouse', 'M150416_MK020', '2015-07-30', '2');
            m_3 = fullfile(testCase.DataFolder, 'multiMouse', 'M150416_MK020', '2015-07-30', '3');
            mkdir(m_1);
            mkdir(m_2);
            mkdir(m_3);

            sourceFolder = fullfile(testCase.DataFolder, 'M0', '2017-10-13', '4');
            copyfile(fullfile(sourceFolder, '*.tif'), mk20_2010);

            copyfile(fullfile(sourceFolder, '*1.tif'), mk20_2107);
            copyfile(fullfile(sourceFolder, '*4.tif'), mk20_2107);

            copyfile(fullfile(sourceFolder, '*2.tif'), m_1);
            copyfile(fullfile(sourceFolder, '*3.tif'), m_1);
            copyfile(fullfile(sourceFolder, '*4.tif'), m_1);

            copyfile(fullfile(sourceFolder, '*4.tif'), m_2);

            copyfile(fullfile(sourceFolder, '*4.tif'), m_3);
            copyfile(fullfile(sourceFolder, '*5.tif'), m_3);
            copyfile(fullfile(sourceFolder, '*1.tif'), m_3);
            oneFile = dir(fullfile(m_3, '*5.tif'));
            movefile(fullfile(m_3, oneFile(1).name), fullfile(m_3, 'one_tiff.tiff'));

            %% run the pipeline
            sm.addEntry(dbEntry);
            testCase.verifyEqual(sm.getNumEntries(), 1);

            sm.selectEntry(1);
            run_pipeline(sm);
            add_deconvolution(sm);

            % now let's check that we have files in proper places
            binFile = 'MK020_2015-07-30_2010_2107_1_2_3_plane1.bin';
            binExist = exist(fullfile(dbEntry.regFilePath, binFile), 'file') ~= 0;
            testCase.verifyTrue(binExist);

            resFolder = fullfile(dbEntry.rootStorage, dbEntry.mouseName{1}, dbEntry.date{1}, '2010_2107_1_2_3');
            resultsFolderExist = exist(resFolder, 'dir') ~= 0;
            testCase.verifyTrue(resultsFolderExist);

            % check results files are in place
            testCase.verifyTrue(exist(fullfile(resFolder, 'F_MK020_2015-07-30_plane1.mat'), 'file') ~= 0);
            testCase.verifyTrue(exist(fullfile(resFolder, 'regops_MK020_2015-07-30.mat'), 'file') ~= 0);

            % delete files
            rmdir(dbEntry.rootStorage, 's');
        end

        function datasetWithoutFolderStructure(testCase)
            sm = CortexLabStorageManager(testCase.Options);

            dbEntry.mouseName = 'not_important';
            dbEntry.date = '2016';
            dbEntry.experiments = [];
            dbEntry.diameter = 12;
            dbEntry.rootStorage = fullfile(testCase.DataFolder, 'M0', '2017-10-13', '4');
            dbEntry.resultsSavePath = fullfile(testCase.DataFolder, 'results');
            dbEntry.regFilePath = fullfile(dbEntry.resultsSavePath, 'regFilePath');

            sm.addEntry(dbEntry);

            sm.selectEntry(1);
            run_pipeline(sm);
            add_deconvolution(sm);

            testCase.verifyTrue(exist(dbEntry.regFilePath, 'dir') ~= 0);
            binFile = dir(fullfile(dbEntry.regFilePath, 'not_important_2016__plane1.bin'));
            testCase.verifyTrue(~isempty(binFile));

            testCase.verifyTrue(exist(fullfile(dbEntry.resultsSavePath, ...
                dbEntry.mouseName, dbEntry.date, 'F_not_important_2016_plane1.mat'),...
                'file') ~= 0);

            testCase.verifyTrue(exist(fullfile(dbEntry.resultsSavePath, ...
                dbEntry.mouseName, dbEntry.date, 'regops_not_important_2016.mat'),...
                'file') ~= 0);

            rmdir(dbEntry.resultsSavePath, 's');
        end
    end

    methods (Static=true)
        function options = getDefaultOptions()
            options.useGPU = gpuDeviceCount > 0;
            options.fig = false;
            options.readTiffHeader = false;
            options.DeleteBin = false;

            % registration options
            options.doRegistration = true;
            options.showTargetRegistration = true;
            options.PhaseCorrelation = true;
            options.SubPixel = Inf;
            options.NimgFirstRegistration = 500;
            options.dobidi = true;

            % cell detection options
            options.ShowCellMap = true;
            options.sig = 0.5;
            options.nSVDforROI = 1000;
            options.NavgFramesSVD = 5000;
            options.signalExtration = 'surround';

            % neuropil options
            options.innerNeuropil = 1;
            options.outerNeuropil = Inf;
            if isinf(options.outerNeuropil)
                options.minNeuropilPixels = 400;
                options.ratioNeuropil = 5;
            end

            % spike deconvolution and neuropil subtraction options
            options.imageRate = 7;
            options.sensorTau = 2;
            options.maxNeurop = 1;

            % red channel options
            options.redthres = 1.5;
            options.redmax = 1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Alberto Battistel
%   7 October 2015
%   CollectData
%   It takes the cells saved from picoscope and convert them into a single
%   vector suitable for next analysis. It does not change the data type.
%   Data are stored in 2 int16 vectors, this spare half of the storage space.
%   Lauch int2float to convert to single float precision.	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%% Parameters
% name of the file to transform
file2import = 'xxx';
% name of the file to save
File2save = [file2import,'_collected'];
% name cell to save
cellA  = 'A';
cellB  = 'B';   
% conversion 1 to V
% names of the cells with the data
nameCellA = 'bufferChAmV';
nameCellB = 'bufferChBmV';
% version
FunctionVer = 'v5';

%% import and tests
% collect info to run in for
cells2save = {cellA, cellB};
NameCell = {nameCellA, nameCellB};

% matfile
fprintf('File to import %s\n', file2import)
Mat2imp = matfile(file2import);
InfoMat = whos(Mat2imp);

% check names of the cells in Mat2imp
celltest = [false, false];
for i1 = 1:length(NameCell)
    for i2 = 1:length(InfoMat)
        if strcmp(InfoMat(i2,1).name, NameCell{i1})
            fprintf('Cell %d: %s ok\n', i1, NameCell{i1})
            celltest(i1,1) = true;
            break
        end
    end
end

if any(~celltest)
    disp('wrong cells name in the file!')
    disp('Variables inside:')
    varNames = fieldnames(InfoMat);
    for i = 1:length(varNames)
        fprintf('\t%s\n', varNames{i,1})
    end
    return
end
    
% import
Cell2import = cell(size(NameCell));
for i = 1:length(NameCell)
    Cell2import{i,1} = Mat2imp.(NameCell{i});
end

% check lengths of the data
Ldata = 0;
for i = 1:length(Cell2import{1,1})
    l = length(Cell2import{1,1}{i,1});
    if l ~= 0
        Ldata = Ldata+l;
    else
        break
    end
end
% number of full cells
nCell = i-1;
fprintf('Total length of the data: %d\n', Ldata);
fprintf('Number of full cells: %d\n', nCell);

%% allocate memory
% the starting data are in Cell2import{i,1}
% Cell2save is the final data container
% fBuf is a buffer 
Cell2save = cell(1,1);

%% create file2save
fprintf('File to save %s\n', File2save)
assignin('caller', cells2save{1}, []);
save(File2save, cells2save{1}, '-v7.3');
Mat2save = matfile(File2save, 'writable', true);

%% import and save cell of data
% the starting data are in Cell2import{i,1}
% Cell2save is the final data container

for iCell2import = 1:length(NameCell)
    % allocate matfile
    Mat2save.(cells2save{iCell2import}) = zeros(Ldata,1,'int16');
    % fill
    l1 = 1;
    for i = 1:nCell
        lCell = length(Cell2import{iCell2import,1}{i,1});
        l2 = lCell + l1 - 1;
        % fill buffer in single precision
        Mat2save.(cells2save{iCell2import})(l1:l2,1) = Cell2import{iCell2import,1}{i,1};
        l1 = l2+1;
    end
       
    % save data
    if iCell2import == 1
        % general information
        Info.importedFile = file2import;
        Info.cellImported = nCell;
        Info.importedCell = NameCell;
        Info.savedCell = cells2save;
        Info.CollectData_version = FunctionVer;
        Info.savedDataType = 'int16';
        Mat2save.Info = Info;
        % length of data
        Mat2save.Ldata = Ldata;
    end
    fprintf('Dai Cazzo %d\n', iCell2import)
end

%% clean up
clear

%% done
disp('Saved data are in int16 without conversion')
% to spare space

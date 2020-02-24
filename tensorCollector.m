% ATTENTION: This script must be executed before running 'angleVariation.m'
%
% This script loads the Segmentation Data folder containing the data for a
% particular movie and then iterates the 'singleFrame_tensorCalculator'
% function over each movie frame while collecting the relevant data. 
%
% The appropriate path to (and including) the 'SegmentationData' folder
% must be specified below.

%% Preparation:
% Clear workspace and start timer
clear all;
tic;

% Specify version name for resulting inertiaTensorCollection file
movieName = strcat('c3');
version = strcat(movieName, '_v003');

% Specify movie frame (time) step increment
timeInc = 1;

%% Specify path to (and including) the 'Segmentation Data' folder for desired movie:
% SegmentationData = '/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Movies/wtYiActin_movie2/movie2/20171230_Gap43mCh_moeGFP_B1_1/SegmentationData';
% SegmentationData = '/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Movies/wtLong_movie2/20141123_117_95_A2_2/SegmentationData';
% SegmentationData = '/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Movies/chlorpromazine_movie1/071716chlorpromazineinjection_ecadgap43_gbe_I1_E1_2/SegmentationData';
% SegmentationData = '/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Movies/chlorpromazine_movie2/2016ecadgap43_10mmchlor_1/SegmentationData';
SegmentationData = '/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Movies/chlorpromazine_movie3/20170417ecadgap43_10mMchlorpromazine_1/SegmentationData';

% Extract relevant data from Segmentation Data folder
dirContent = dir(SegmentationData);
isub = [dirContent(:).isdir]; %# returns logical vector
folderNames = {dirContent(isub).name}';
folderNames(ismember(folderNames,{'.','..'})) = [];
numFrames = length(folderNames);

% Specify starting parameters
startFrame = 1;
temp_percent = 0;

% Create container for collection of inertia tensors for each cell in every
% movie frame
inertiaTensorCollection = cell(numFrames, 1);

% Specify file name for saving results
% filename = strcat('InertiaTensorCollection_frames_', int2str(startFrame), '_to_', int2str(numFrames), '_', version, '.mat');
filename = strcat('InertiaTensorCollection_', version, '.mat');

%% Collect moment of inertia tensors for each cell in every movie frame:
for frameNum = startFrame:timeInc:numFrames

    %%% Extract vertex locations:
    if ispc
        file = load(strcat(SegmentationData, '\', string(folderNames(frameNum)), '\dstruct_nodes.mat'));
        file_CellMat = load(strcat(SegmentationData, '\', string(folderNames(frameNum)), '\dstruct_nodeCellMat.mat'));
    else
        file = load(strcat(SegmentationData, '/', string(folderNames(frameNum)), '/dstruct_nodes.mat'));
        file_CellMat = load(strcat(SegmentationData, '/', string(folderNames(frameNum)), '/dstruct_nodeCellMat.mat'));
    end
    
    %%% This option (NN) is always equal to 1 for wild-type movies, but can
    %%% be up to 5 for chlorpromazine movies.
    NN = 1;
%     if strcmp('c3', movieName) || strcmp('c2', movieName)  || strcmp('c1', movieName) 
%         NN = 5;
%     else
%         NN = 1;
%     end
    
    row_verts_list = file.dstruct_nodes(NN).positions(:, 1); % Row
    col_verts_list = file.dstruct_nodes(NN).positions(:, 2); % Column
    nodeCellMat = file_CellMat.dstruct_nodeCellMat(NN).matrix;

    x_verts = {};
    y_verts = {};

    for cellNum = 1:size(nodeCellMat,2)
%         x_verts{end+1} = x_verts_list(find(nodeCellMat(:,cellNum) == 1));
%         y_verts{end+1} = y_verts_list(find(nodeCellMat(:,cellNum) == 1));
        x_verts{end+1} = col_verts_list(find(nodeCellMat(:,cellNum) == 1));
        y_verts{end+1} = row_verts_list(find(nodeCellMat(:,cellNum) == 1));
    end

    % Get number of cells in this movie frame
    numCells = size(nodeCellMat,2);
    
    % Calculate the inertia tensor for each cell in this movie frame
    inertiaTensorCollection{frameNum} = singleFrame_tensorCalculator( x_verts, y_verts );

    % Report progress
    if ((frameNum)/numFrames)*100 >= temp_percent
        disp(['Current Frame Number: ', int2str(frameNum)]);
        disp(['Computation of Moment of Inertia Tensor collection is ', int2str(((frameNum)/numFrames)*100), ' % complete ...']);
        temp_percent = temp_percent + 10;
        % Optional incremental saving
%         save(filename, 'cellMax', 'thetas', 'maxOrientationAngles', 'startFrame', 'numFrames', 'thetaInc');
%         save(filename, 'inertiaTensorCollection');
        toc
    end
end

%% Save resulting collection of this movie's moment of inertia tensors:
save(filename, 'inertiaTensorCollection');
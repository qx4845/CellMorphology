clear all;
tic;

movieName = strcat('c3');
version = strcat(movieName, '_v2');

timeInc = 1;

% SegmentationData = '/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Movies/wtYiActin_movie2/movie2/20171230_Gap43mCh_moeGFP_B1_1/SegmentationData';
% SegmentationData = '/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Movies/wtLong_movie2/20141123_117_95_A2_2/SegmentationData';
% SegmentationData = '/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Movies/chlorpromazine_movie1/071716chlorpromazineinjection_ecadgap43_gbe_I1_E1_2/SegmentationData';
% SegmentationData = '/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Movies/chlorpromazine_movie2/2016ecadgap43_10mmchlor_1/SegmentationData';
SegmentationData = '/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Movies/chlorpromazine_movie3/20170417ecadgap43_10mMchlorpromazine_1/SegmentationData';

dirContent = dir(SegmentationData);
isub = [dirContent(:).isdir]; %# returns logical vector
folderNames = {dirContent(isub).name}';
folderNames(ismember(folderNames,{'.','..'})) = [];
numFrames = length(folderNames);
startFrame = 1;

temp_percent = 0;

inertiaTensorCollection = cell(numFrames, 1);

% filename = strcat('InertiaTensorCollection_frames_', int2str(startFrame), '_to_', int2str(numFrames), '_', version, '.mat');
filename = strcat('InertiaTensorCollection_', version, '.mat');

for frameNum = startFrame:timeInc:numFrames

    %%% Extract vertex locations:
    if ispc
        file = load(strcat(SegmentationData, '\', string(folderNames(frameNum)), '\dstruct_nodes.mat'));
        file_CellMat = load(strcat(SegmentationData, '\', string(folderNames(frameNum)), '\dstruct_nodeCellMat.mat'));
    else
        file = load(strcat(SegmentationData, '/', string(folderNames(frameNum)), '/dstruct_nodes.mat'));
        file_CellMat = load(strcat(SegmentationData, '/', string(folderNames(frameNum)), '/dstruct_nodeCellMat.mat'));
    end
    x_verts_list = file.dstruct_nodes(1).positions(:, 1);
    y_verts_list = file.dstruct_nodes(1).positions(:, 2);
    nodeCellMat = file_CellMat.dstruct_nodeCellMat.matrix;

    x_verts = {};
    y_verts = {};

    for cellNum = 1:size(nodeCellMat,2)
        x_verts{end+1} = x_verts_list(find(nodeCellMat(:,cellNum) == 1));
        y_verts{end+1} = y_verts_list(find(nodeCellMat(:,cellNum) == 1));
    end

    numCells = size(nodeCellMat,2);
    inertiaTensorCollection{frameNum} = singleFrame_tensorCalculator( x_verts, y_verts );

    %%% Report progress:
    if ((frameNum)/numFrames)*100 >= temp_percent
        disp(['Current Frame Number: ', int2str(frameNum)]);
        disp(['Computation of Moment of Inertia Tensor collection is ', int2str(((frameNum)/numFrames)*100), ' % complete ...']);
        temp_percent = temp_percent + 10;
%         save(filename, 'cellMax', 'thetas', 'maxOrientationAngles', 'startFrame', 'numFrames', 'thetaInc');
%         save(filename, 'inertiaTensorCollection');
        toc
    end
end
save(filename, 'inertiaTensorCollection');

% firstMaxFrame = find(maxOrientationAngles(1,:) == max(maxOrientationAngles(1,:)), 1, 'first');
% meanPostMaxFrame = mean(maxOrientationAngles(1, firstMaxFrame:end));
% 
% cardinalAxesOffset = [meanPostMaxFrame - 90, meanPostMaxFrame];
% 
% save(filename, 'cardinalAxesOffset', 'meanPostMaxFrame', 'firstMaxFrame', 'cellMax', 'thetas', 'maxOrientationAngles', 'startFrame', 'numFrames', 'thetaInc', 'cutoff_percentage');
% 
% figure();
% fig_01 = subplot(1,1,1);
% domain = [1:size(maxOrientationAngles, 2)];
% plot(domain, maxOrientationAngles(1,:), 'r', domain, maxOrientationAngles(2,:), 'b');
% xlabel('time (frames)');
% ylabel('Angle (degrees)');
% title(strcat(movieName, ' Max Angle of Rotation'));

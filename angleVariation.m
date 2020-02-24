% ATTENTION: This script requires that 'tensorCollector.m' be run with its 
% output successfully saved before executing this script. 
% 
% This script uses the output file resulting from 'tensorCollector.m': 
% specifically the string corresponding to the version suffix of the 
% 'tensorCollector.m' output file to load this file and incorporate the 
% contained inertia tensors. 
%
% First, this script calculates the Angular Momentum (L) for every cell in 
% every frame at each angle in a prescribed range of angles.  Resulting 
% angular momenta are stored in the 'theta_vs_time' variable and plotted in 
% Subplot #1 as a heatmap where the color corresponds to the average 
% angular momentum for a particular time and rotation angle. 
%
% Next, the eigenvalues and eigenvectors are also calculated from the
% loaded inertia tensors and are used to deteremine the angles which the
% major and minor axes of each cell form with the positive x-axis (a.k.a.
% the extreme angles of elogation, or extreme directions). These results
% are plotted in Subplot #2.
%
% Then, the side-lengths, in pixels, of the equivalent rectangle that 
% corresponds to each cell are computed according to the formulas derived
% here: https://www.dropbox.com/s/w6l65mm77n8dau8/AMequivalentRectangleMeasuremts.pdf?dl=0
%
% Finally, the extreme directions are averaged over the entire movie length
% to determine the cardinal axes of the sample. The cardinal axes of the
% microscope (say, 0 and 90 degrees) are intended to align with those of
% the sample, but this is often not the case. To this end, this script also
% outputs the cardinal axes of the movie. For example, if the cardinal axes
% of the sample are rotated 5 degrees counterclockwise from those of the
% microscope viewing window, then the sample axes would be 5 and 95
% degrees, and the 'cardinalAxesOffset' variable below would equal [5; 95].
%
% All measurements are averaged in time over a frame-window of width 'ww' 
% defined below. Relevant results are saved and plotted, with several
% alternative plot configurations left commented out for future variations.

%% Preparation:
% Clear workspace and start timer
clear all;
tic;

% Specify version for resulting file name suffix
movieName = strcat('c3');
version = strcat(movieName, '_4');

% Specify version of 'tensorCollection' to be loaded
%%% ('tensorCollector.m' must be executed before running this script.
%%% This will collect the inertia tensors required by this script.)
tensorCollection_version = strcat(movieName, '_v003');

% Specify angular resolution for resulting plot in degrees
thetaInc = 5;

% Specify limits of desired angle range in degrees
startAngle = 0;
stopAngle = 180;

% Create vector of angles from 'startAngle' to 'stopAngle' degrees in 
% increments of 'thetaInc' degrees, and then convert the entire vector to 
% units of RADIANS:
thetas = ((startAngle:thetaInc:stopAngle)/360)*2*pi;

% Load the inertia tensor collection version that is specified above
inertiaTensorCollection = load(strcat('InertiaTensorCollection_', tensorCollection_version, '.mat'));
% inertiaTensorCollection = load(strcat('/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Cell_Measurement/InertiaTensor/InertiaTensorCollection_', tensorCollection_version, '.mat'));

% Get numbers of frames and cells from imported data
numFrames = size(inertiaTensorCollection.inertiaTensorCollection, 1);
numCells = size(inertiaTensorCollection.inertiaTensorCollection{end,1}, 1) - 1;

% Create containers for storing new data to be generated
theta_vs_time = NaN(length(thetas), numFrames, 1);
angularMomenta = NaN(numCells, numFrames, 3);
maxOrientationAngles = NaN(2, numFrames);
equivalentRectangleMeasurements = NaN(numCells, numFrames, 3);
extremeDirections = NaN(numCells, numFrames, 2);
celAvgExtremeDirections = NaN(2, numFrames);
cellCounts = NaN(1, numFrames);

% Collect eigenvalues and eigenvectors for each intertia tensor in the
% collection
[ eigenVals, eigenVecs ] = eigenCollector( tensorCollection_version );

% Specify window width (ww) for rolling average in units of frames
ww = 10;

%% Compute and collect angles of extreme elongation and associated statistics:
for frameNum = 1:numFrames
    cellCounts(1, frameNum) = length(inertiaTensorCollection.inertiaTensorCollection{frameNum, 1});
    
    % For case before window width has been reached
    if frameNum <= ww
        % At this frame, go through all angles and find associated angular momenta
        for theta = 1:length(thetas)
            % Create a unit vector pointing forming an angle theta with the
            % positive x-axis
            unitVec = [cos(thetas(theta)); sin(thetas(theta))];
            
            % Cycle thorugh all cells in this frame
            for cellNum = 1:length(inertiaTensorCollection.inertiaTensorCollection{frameNum, 1})
                
                % Get relevant tensor from tensor collection
                tensor = inertiaTensorCollection.inertiaTensorCollection{frameNum,1}(cellNum);
                
                % Compute Angular Momentum for this cell's tensor at this
                % angle
                angularMomentum = tensor{1} * unitVec;
                % if unitVec is composed of angular frequency components, then
                % angularMomentum is the angular momentum "L" in particular direction.
                
                % Organize these angular momentua into collection
                angularMomenta(cellNum, frameNum, 1) = angularMomentum(1);
                angularMomenta(cellNum, frameNum, 2) = angularMomentum(2);
                angularMomenta(cellNum, frameNum, 3) = sqrt(angularMomentum(1)^2 + angularMomentum(2)^2);
            end
            
            % Record cell-averaged angular momenta for each angle and fram
            theta_vs_time(theta, frameNum, 1) = nanmean(angularMomenta(:, frameNum, 3));
        end

        % Once the angular momenta for this frame are calculated, proceed
        % to cycle through each cell in this frame and use the eigenvalues
        % and eigenvectors to determine the extreme elongation directions
        % of each cell in this frame. 
        for cellNum = 1:length(inertiaTensorCollection.inertiaTensorCollection{frameNum, 1})
            
            % Calculate equivalent rectangle measurements in pixels
            % according to the prescription outlined here:
            % https://www.dropbox.com/s/w6l65mm77n8dau8/AMequivalentRectangleMeasuremts.pdf?dl=0
            equivalentRectangleMeasurements(cellNum, frameNum, 3) = nthroot(mean(abs(eigenVals(cellNum, 1:frameNum, 2)), 2) / mean(abs(eigenVals(cellNum, 1:frameNum, 1)), 2), 2); % 'e': elongation factor (a/b)
            equivalentRectangleMeasurements(cellNum, frameNum, 1) = nthroot((18 * mean(abs(eigenVals(cellNum, 1:frameNum, 2)), 2) / equivalentRectangleMeasurements(cellNum, frameNum, 3)) , 4); % 'a': length along long axis
            equivalentRectangleMeasurements(cellNum, frameNum, 2) = nthroot((18 * mean(abs(eigenVals(cellNum, 1:frameNum, 1)), 2) / equivalentRectangleMeasurements(cellNum, frameNum, 3)) , 4); % 'b': length along short axis
%             equivalentRectangleMeasurements(cellNum, frameNum, 2) = nthroot((18 * mean(abs(eigenVals(cellNum, 1:frameNum, 1)), 2) / (equivalentRectangleMeasurements(cellNum, frameNum, 3) * equivalentRectangleMeasurements(cellNum, frameNum, 1)^2)) , 2); % 'b': length along short axis

            % Create temporary variable to store extreme direction results
            % for this early (pre-window-width) averaging
            eD = NaN(2, frameNum);
            
            % Cycle through frames that have already been completed (for
            % this cell only) and compute extreme directions for displaying
            % this frame's average extreme elongations. 
            for wFrame = 1:frameNum
                eVec = eigenVecs{cellNum, wFrame};
                a = [eVec(1, 2), eVec(2, 2), 0];
                b = [1, 0, 0];
                eD(1, wFrame) = atan2(norm(cross(a,b)), dot(a,b)) * (360/(2*pi));
                c = [eVec(1, 1), eVec(2, 1), 0];
                d = [1, 0, 0];
                eD(2, wFrame) = atan2(norm(cross(c,d)), dot(c,d)) * (360/(2*pi));
            end
            
            % Average each cell's extreme directions over that cell's
            % history since the start of the movie.
            extremeDirections(cellNum, frameNum, 1) = nanmean(eD(1,:), 2);
            extremeDirections(cellNum, frameNum, 2) = nanmean(eD(2,:), 2);
        end
        
    % For general case averaging over window width
    else
        
        % At this frame, go through all angles and find associated angular momenta
        for theta = 1:length(thetas)
            
            % Create a unit vector pointing forming an angle theta with the
            % positive x-axis
            unitVec = [cos(thetas(theta)); sin(thetas(theta))];
            
            % Cycle thorugh all cells in this frame
            for cellNum = 1:length(inertiaTensorCollection.inertiaTensorCollection{frameNum, 1})
                
                % Get relevant tensor from tensor collection
                tensor = inertiaTensorCollection.inertiaTensorCollection{frameNum,1}(cellNum);
                
                % Compute the angular momentum for this cell's inertia 
                % tensor at this angle
                angularMomentum = tensor{1} * unitVec;
                % if unitVec is composed of angular frequency components, then
                % angularMomentum is the angular momentum "L" in particular direction.
                
                % Organize these angular momentua into a collection
                angularMomenta(cellNum, frameNum, 1) = angularMomentum(1);
                angularMomenta(cellNum, frameNum, 2) = angularMomentum(2);
                angularMomenta(cellNum, frameNum, 3) = sqrt(angularMomentum(1)^2 + angularMomentum(2)^2);
            end

            % Record cell-averaged angular momenta for each angle and frame
            theta_vs_time(theta, frameNum, 1) = nanmean(angularMomenta(:, frameNum, 3));
        end

        % Once the angular momenta for this frame are calculated, proceed
        % to cycle through each cell in this frame and use the eigenvalues
        % and eigenvectors to determine the extreme elongation directions
        % of each cell in this frame. 
        for cellNum = 1:min(cellCounts(1, frameNum - ww: frameNum))

            % Calculate equivalent rectangle measurements in pixels
            % according to the prescription outlined here:
            % https://www.dropbox.com/s/w6l65mm77n8dau8/AMequivalentRectangleMeasuremts.pdf?dl=0
            equivalentRectangleMeasurements(cellNum, frameNum, 3) = nthroot(mean(abs(eigenVals(cellNum, frameNum - ww: frameNum, 2)), 2) / mean(abs(eigenVals(cellNum, frameNum - ww: frameNum, 1)), 2), 2); % 'e': elongation factor (a/b)
            equivalentRectangleMeasurements(cellNum, frameNum, 1) = nthroot((18 * mean(abs(eigenVals(cellNum, frameNum - ww: frameNum, 2)), 2) / equivalentRectangleMeasurements(cellNum, frameNum, 3)) , 4); % 'a': length along long axis
            equivalentRectangleMeasurements(cellNum, frameNum, 2) = nthroot((18 * mean(abs(eigenVals(cellNum, frameNum - ww: frameNum, 1)), 2) / equivalentRectangleMeasurements(cellNum, frameNum, 3)) , 4); % 'b': length along short axis
%             equivalentRectangleMeasurements(cellNum, frameNum, 2) = nthroot((18 * mean(abs(eigenVals(cellNum, frameNum - ww: frameNum, 1)), 2) / (equivalentRectangleMeasurements(cellNum, frameNum, 3) * equivalentRectangleMeasurements(cellNum, frameNum, 1)^2)) , 2); % 'b': length along short axis

            % Create temporary variable to store extreme direction / elongation results
            % for averaging this cell over the previous window-width ('ww')
            % number of frames, and also specify this domain of frames.
            eD2 = NaN(2, ww + 1);
            eD2domain = frameNum - ww: frameNum;
            
            % Cycle through previous 'ww'/ window-width # of frames (for
            % this cell only) and compute extreme directions for displaying
            % this frame's window-averaged extreme elongations. 
            for wFrame = 1:length(eD2domain)
                eVec = eigenVecs{cellNum, eD2domain(wFrame)};
                a = [eVec(1, 2), eVec(2, 2), 0];
                b = [1, 0, 0];
                eD2(1, wFrame) = atan2(norm(cross(a,b)), dot(a,b)) * (360/(2*pi)); 
                c = [eVec(1, 1), eVec(2, 1), 0];
                d = [1, 0, 0];
                eD2(2, wFrame) = atan2(norm(cross(c,d)), dot(c,d)) * (360/(2*pi));
            end
            
            % Average each cell's extreme directions over that cell's
            % history over the previous window ('ww' #  of frames).
            extremeDirections(cellNum, frameNum, 1) = nanmean(eD2(1,:), 2);
            extremeDirections(cellNum, frameNum, 2) = nanmean(eD2(2,:), 2);
        end
    end
end

%% Find cardinal axes of movie:
cardinalAxesOffset = NaN(2, 1);
cardinalAxesOffset(1, 1) = nanmean(nanmean(extremeDirections(:, :, 1), 1));

if cardinalAxesOffset(1, 1) >= 90
    cardinalAxesOffset(2, 1) = cardinalAxesOffset(1, 1) - 90;
else
    cardinalAxesOffset(2, 1) = cardinalAxesOffset(1, 1) + 90;
end

%% Master Save:
filename = strcat('angle_vs_time_equivRectangle_', version, '.mat');
save(filename, 'theta_vs_time', 'angularMomenta', 'movieName', 'version', 'thetas', 'numFrames', 'equivalentRectangleMeasurements', 'extremeDirections', 'celAvgExtremeDirections', 'cardinalAxesOffset');

%% Plots (with alternate configuration options):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PLOTS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
maxPlotTime = numFrames;
timeInc = 1;

suptitle({strcat('Movie: ', movieName, ' - Cardinal Axes: ', ' ', num2str(round(cardinalAxesOffset(1, 1), 1)), ' x ', ' ', num2str(round(cardinalAxesOffset(2, 1), 1)), ' degrees'), ''});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subplot #1 (with options):

% fig10 = subplot(1,1,1);
% imagesc([theta_vs_time(:, 1:timeInc:maxPlotTime, 1); theta_vs_time(:, 1:timeInc:maxPlotTime, 2)]);
% colorbar;
% title('Prevalence of Rotation Angles');
% ylabel('Angle of Rotation (degrees)');
% set(gca,'YTick', 1:2*length(thetas) ); %This are going to be the only values affected.
% set(gca,'YTickLabel', [(thetas*360)/(2*pi),90+(thetas*360)/(2*pi)] ); %This is what it's going to appear in those places.
% xlabel('time (frames)');

fig10 = subplot(4,1,1);
imagesc(theta_vs_time(:, 1:timeInc:maxPlotTime, 1));
colorbar;
title(strcat('Cell-Averaged Angular Momenta, movie: ', movieName));
ylabel('Angle (degrees)');
% set(gca,'YTick', 1:length(thetas) ); %This are going to be the only values affected.
% set(gca,'YTickLabel', [(thetas*360)/(2*pi),90+(thetas*360)/(2*pi)] ); %This is what it's going to appear in those places.
set(gca,'YTick', 1:2:length(thetas) ); %This are going to be the only values affected.
A = (thetas*360)/(2*pi);
set(gca,'YTickLabel', A(1:2:length(thetas)) ); %This is what it's going to appear in those places.xlabel('time (frames)');
set(gca,'YDir','normal');
% hold on;
% domain = [1:size(maxOrientationAngles, 2)];
% plot(domain, maxOrientationAngles(1,:), 'r', domain, maxOrientationAngles(2,:), 'b');
% legend('Max Orientation Angle', 'Min Orientation Angle');
axis normal;

%%% custom colormap plot:
mymap = ones(180*2+1, 3);
mymap(1:end, 1) = zeros(size(mymap, 1), 1) + 0.5;
mymap(1:end, 2) = abs(linspace(1, -1, 180*2+1));
Atemp = linspace(0,1,181);
mymap(1:end, 3) = [Atemp, linspace(Atemp(end-1), 0, 180)];
mymap(1,1:3) = [0,0,0]; %%% (with zeros being black)
mymap(2,1:3) = [0,0,0]; %%% (with near-zeros being black)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subplot #2:
fig11 = subplot(4,1,2);
% domain = 1:size(extreme_locations, 2);
domain = 1:numFrames;
% plot(domain, nanmean(extremeDirections(:, :, 1), 1), 'r', domain, nanmean(extremeDirections(:, :, 2), 1), 'b');
plot(domain, nanmean(extremeDirections(:, :, 1), 1), 'r');
A = (thetas*360)/(2*pi);
ylim([0 180]);
yticks(A(1:2:length(thetas)));
% % set(gca,'YTickLabel', A(1:2:length(thetas)), 'YAxisLocation', 'right');
ylabel('Angle (degrees)');
% box off
% axes('xlim', [0 numFrames], 'ylim', [0 180], 'color', 'none',...
% 'YTick',[], 'YAxisLocation', 'right')
% xlabel('time (frames)');
title(strcat('Maximum Angles of Rotation, movie: ', movieName));
% legend('max', 'min');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subplot #3:
fig12 = subplot(4,1,3);
% plot(domain, extreme_measurements(1,:), 'r', domain, extreme_measurements(2,:), 'b');
plot(domain, nanmean(equivalentRectangleMeasurements(:, :, 1), 1), 'r', domain, nanmean(equivalentRectangleMeasurements(:, :, 2), 1), 'b');
% xlabel('time (frames)');
% ylabel('Angular Momentum (px^{4}/s)');
% ylabel('Angular Momentum');
% ylabel('L = I*\omega', 'FontWeight', 'bold');
% ylabel('L = I\cdot\omega (px^{4}/s)');
ylabel('Length (pixels)');
title(strcat('Equiv. Rectangle Measurements Extreme Angles, movie: ', movieName));
legend('max', 'min');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subplot #4:
fig13 = subplot(4,1,4);
% plot(domain, extreme_measurements(1,:) ./ extreme_measurements(2,:), 'k');
plot(domain, nanmean(equivalentRectangleMeasurements(:, :, 3), 1), 'k');
xlabel('time (frames)');
ylabel('Elongation');
if strcmp(movieName, 'r4')
    ylim([0 4]);
end
title(strcat('Elongation at Extreme Angles, movie: ', movieName));

%% As yet unsuccessful method to link all above subplots with a common time axis:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linkaxes([fig10 fig11 fig12 fig13],'xy');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Alternative to Subplot #2 (fig11) that instead plots RELATIVE orientation
% %  angles, and also ignores cells that persist for < 5 frames:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% temp_percent_2 = 1;
% maxPlotTime = numFrames;
% relativeOrientation = NaN(numCells, maxPlotTime);
% for row = 1:numCells %%% row -> cell index
%     for col = 2:maxPlotTime            %%% col -> time
%         if ~isnan(ellipseProperties(row,col,4))
%             relativeOrientation(row, col) = (ellipseProperties(row,col,4) - ellipseProperties(row,find(~isnan(ellipseProperties(row,:,4)) == 1, 1, 'first'),4));
%         end
%     end
%     if ~isnan(ellipseProperties(row,col,4))
%         relativeOrientation(row, 1) = ellipseProperties(row,find(~isnan(ellipseProperties(row,:,4)) == 1, 1, 'first'),4);
%     end
%     
%     %%% Output progress:
%     if ((row)/numCells)*100 > temp_percent_2
%         disp(['Orientation plot computation is ', int2str(((row)/numCells)*100), ' % complete ...']);
%         temp_percent_2 = temp_percent_2 + 1;
%         toc
%     end
% end
% 
% %%% Go through and ignore cells that persist for <5 frames
% for row = 1:size(relativeOrientation, 1)
%     reals = ~isnan(relativeOrientation(row, :));
% %     consecGroups = diff(reals);
% %     CM = diff(reals) == 1; %%% (Consecutive Markerrs): number of Markers for Consecutive non-NaN groups
% %     for group = 1:length(CM(CM==1))
% %         find()
% %     end
%     counter = 0;
% %     ind = 1;
%     for ind = 1:length(reals) - 1
% %     while ind < length(reals) - 1
%         
%         if reals(ind) == 1
%             counter = counter + 1;
%             if reals(ind + 1) == 0 && counter <= 5
%                 relativeOrientation(row, ind-counter+1:ind) = NaN(1, counter);
%                 counter = 0;
%             elseif  reals(ind + 1) == 0 && counter > 5
%                 counter = 0;
%             end
%         end
% %         ind = ind + 1;
%     end
%     
% end
% min_preShift_relativeOrientation = min(min(relativeOrientation));
% 
% relativeOrientation = relativeOrientation + 181;
% % shifted_orientation = ellipseProperties(:,:,4) + 91;
% relativeOrientation(isnan(relativeOrientation)) = 0;
% % shifted_orientation(isnan(shifted_orientation)) = 0;
% 
% fig11 = subplot(2,1,2);
% imagesc(relativeOrientation(1:numCells, 2:end))
% ylabel('Cell Index');
% xlabel('time (frames)');
% title('Relative (Current - Initial) Oritentation Angles');
% % cbh=colorbar('h');
% % set(cbh,'YTick',['NaN', -90:30:90]);
% colorbar;
% colormap(gca,mymap);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
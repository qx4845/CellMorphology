clear all;
tic;

movieName = strcat('r3');
version = strcat(movieName, '_v24');

tensorCollection_version = strcat(movieName, '_v12');

thetaInc = 5;
% thetas = ((0:thetaInc:90)/360)*2*pi;
% full_thetas = ((0:thetaInc:180)/360)*2*pi;
thetas = ((0:thetaInc:180)/360)*2*pi;

inertiaTensorCollection = load(strcat('/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Cell_Measurement/InertiaTensor/InertiaTensorCollection_', tensorCollection_version, '.mat'));
% inertiaTensorCollection = load(strcat('/Users/stephenwedekind/Dropbox/School/RESEARCH/Loerke_group/Cell_Measurement/InertiaTensor/InertiaTensorCollection_', tensorCollection_version, '_intermediate2.mat'));

numFrames = size(inertiaTensorCollection.inertiaTensorCollection, 1);
% numFrames = 250;
numCells = size(inertiaTensorCollection.inertiaTensorCollection{end,1}, 1) - 1;
% theta_vs_time = NaN(length(thetas), numFrames, 2);
theta_vs_time = NaN(length(thetas), numFrames, 1);
angularMomenta = NaN(numCells, numFrames, 3);
maxOrientationAngles = NaN(2, numFrames);
% max_combo_range = NaN(length(thetas), numFrames);
% min_combo_range = NaN(length(thetas), numFrames);
equivalentRectangleMeasurements = NaN(numCells, numFrames, 3);
extremeDirections = NaN(numCells, numFrames, 2);
celAvgExtremeDirections = NaN(2, numFrames);
cellCounts = NaN(1, numFrames);

[ eigenVals, eigenVecs ] = eigenCollector( tensorCollection_version );

ww = 10; % Window Width in units of frames
for frameNum = 1:numFrames
    cellCounts(1, frameNum) = length(inertiaTensorCollection.inertiaTensorCollection{frameNum, 1});
    if frameNum <= ww
        for theta = 1:length(thetas)
            unitVec = [cos(thetas(theta)); sin(thetas(theta))];
            for cellNum = 1:length(inertiaTensorCollection.inertiaTensorCollection{frameNum, 1})
                tensor = inertiaTensorCollection.inertiaTensorCollection{frameNum,1}(cellNum);
                angularMomentum = tensor{1} * unitVec;
                % if unitVec is composed of angular frequency components, then
                % angularMomentum is the angular momentum "L" in particular direction.
                angularMomenta(cellNum, frameNum, 1) = angularMomentum(1);
                angularMomenta(cellNum, frameNum, 2) = angularMomentum(2);
                angularMomenta(cellNum, frameNum, 3) = sqrt(angularMomentum(1)^2 + angularMomentum(2)^2);
            end
    %         theta_vs_time(theta, frameNum, 1) = nanmean(inertiaMeasurements(:, frameNum, 1));
    %         theta_vs_time(theta, frameNum, 2) = nanmean(inertiaMeasurements(:, frameNum, 2));
            theta_vs_time(theta, frameNum, 1) = nanmean(angularMomenta(:, frameNum, 3));
        end

        for cellNum = 1:length(inertiaTensorCollection.inertiaTensorCollection{frameNum, 1}) -1
    %         [ eigenVals, eigenVecs ] = eigenCollector( tensorCollection_version );
            equivalentRectangleMeasurements(cellNum, frameNum, 3) = nthroot(abs(nanmean(eigenVals(cellNum, 1:frameNum, 2), 2) / nanmean(eigenVals(cellNum, 1:frameNum, 1), 2)), 2); % 'e': elongation factor (a/b)
            equivalentRectangleMeasurements(cellNum, frameNum, 1) = nthroot(abs( 18 * nanmean(eigenVals(cellNum, 1:frameNum, 2), 2) / nanmean(equivalentRectangleMeasurements(cellNum, 1:frameNum, 3), 2) ), 4); % 'a': length along long axis
            equivalentRectangleMeasurements(cellNum, frameNum, 2) = nthroot(abs( 18 * nanmean(eigenVals(cellNum, 1:frameNum, 1), 2) / nanmean(equivalentRectangleMeasurements(cellNum, 1:frameNum, 3), 2) ), 4); % 'b': length along short axis
%             equivalentRectangleMeasurements(cellNum, frameNum, 2) = nthroot(( 18 * nanmean(eigenVals(cellNum, 1:frameNum, 1), 2) / (nanmean(equivalentRectangleMeasurements(cellNum, 1:frameNum, 3), 2) * nanmean(equivalentRectangleMeasurements(cellNum, 1:frameNum, 1), 2)^2) ), 4); % 'b': length along short axis

%             eVec = eigenVecs{cellNum, frameNum};
            eD = NaN(2, frameNum);
            for wFrame = 1:frameNum
                eVec = eigenVecs{cellNum, wFrame};
%                 eD(1, wFrame) = atan2(eVec(1, 2), eVec(2, 2)) * (360/(2*pi));
                a = real([eVec(1, 2), eVec(2, 2), 0]);
                b = [1, 0, 0];
                eD(1, wFrame) = atan2(norm(cross(a,b)), dot(a,b)) * (360/(2*pi));
                c = real([eVec(1, 1), eVec(2, 1), 0]);
                d = [1, 0, 0];
                eD(2, wFrame) = atan2(norm(cross(c,d)), dot(c,d)) * (360/(2*pi));
%                 eD(1, wFrame) = atan2(eVec(2, 2), eVec(1, 2)) * (360/(2*pi)) + 90;
%                 eD(1, wFrame) = atan2(eVec(1, 1), eVec(2, 1)) * (360/(2*pi)) + 90;
%                 extremeDirections(cellNum, frameNum, 1) = atan2(eVec(1, 2), eVec(2, 2)) * (360/(2*pi)) + 90;
%         %         extremeDirections(cellNum, frameNum, 2) = atan2(eVec(1, 1), eVec(2, 1)) * (360/(2*pi)) + 90;
            end
            extremeDirections(cellNum, frameNum, 1) = nanmean(eD(1,:), 2);
            extremeDirections(cellNum, frameNum, 2) = nanmean(eD(2,:), 2);
%             extremeDirections(cellNum, frameNum, 1) = atan2(eVec(1, 1), eVec(2, 1)) * (360/(2*pi)) + 90;
% %             extremeDirections(cellNum, frameNum, 1) = atan2(eVec(1, 2), eVec(2, 2)) * (360/(2*pi)) + 90;
        end
    else
        for theta = 1:length(thetas)
            unitVec = [cos(thetas(theta)); sin(thetas(theta))];
            for cellNum = 1:length(inertiaTensorCollection.inertiaTensorCollection{frameNum, 1}) - 1
                tensor = inertiaTensorCollection.inertiaTensorCollection{frameNum,1}(cellNum);
                angularMomentum = tensor{1} * unitVec;
                % if unitVec is composed of angular frequency components, then
                % angularMomentum is the angular momentum "L" in particular direction.
                angularMomenta(cellNum, frameNum, 1) = angularMomentum(1);
                angularMomenta(cellNum, frameNum, 2) = angularMomentum(2);
                angularMomenta(cellNum, frameNum, 3) = sqrt(angularMomentum(1)^2 + angularMomentum(2)^2);
            end
    %         theta_vs_time(theta, frameNum, 1) = nanmean(inertiaMeasurements(:, frameNum, 1));
    %         theta_vs_time(theta, frameNum, 2) = nanmean(inertiaMeasurements(:, frameNum, 2));
            theta_vs_time(theta, frameNum, 1) = nanmean(angularMomenta(:, frameNum, 3));
        end

        for cellNum = 1:min(cellCounts(1, frameNum - ww: frameNum)) - 1
    %         [ eigenVals, eigenVecs ] = eigenCollector( tensorCollection_version );
            equivalentRectangleMeasurements(cellNum, frameNum, 3) = nthroot(abs(nanmean(eigenVals(cellNum, frameNum - ww: frameNum, 2), 2) / nanmean(eigenVals(cellNum, frameNum - ww: frameNum, 1), 2)), 2); % 'e': elongation factor (a/b)
            equivalentRectangleMeasurements(cellNum, frameNum, 1) = nthroot(abs(18 * nanmean(eigenVals(cellNum, frameNum - ww: frameNum, 2), 2) / nanmean(equivalentRectangleMeasurements(cellNum, frameNum - ww: frameNum, 3), 2) ), 4); % 'a': length along long axis
            equivalentRectangleMeasurements(cellNum, frameNum, 2) = nthroot(abs(18 * nanmean(eigenVals(cellNum, frameNum - ww: frameNum, 1), 2) / nanmean(equivalentRectangleMeasurements(cellNum, frameNum - ww: frameNum, 3), 2) ), 4); % 'b': length along short axis
%             equivalentRectangleMeasurements(cellNum, frameNum, 2) = nthroot((18 * nanmean(eigenVals(cellNum, frameNum - ww: frameNum, 1), 2) / (nanmean(equivalentRectangleMeasurements(cellNum, frameNum - ww: frameNum, 3), 2) * nanmean(equivalentRectangleMeasurements(cellNum, frameNum - ww: frameNum, 1), 2)^2) ), 4); % 'b': length along short axis

%             eVec = eigenVecs{cellNum, frameNum};
            eD2 = NaN(2, ww + 1);
            eD2domain = frameNum - ww: frameNum;
            for wFrame = 1:length(eD2domain)
                eVec = eigenVecs{cellNum, eD2domain(wFrame)};
%                 a = abs([eVec(1, 2), eVec(2, 2), 0]);
                a = real([eVec(1, 2), eVec(2, 2), 0]);
                b = [1, 0, 0];
                eD2(1, wFrame) = atan2(norm(cross(a,b)), dot(a,b)) * (360/(2*pi)); 
%                 c = abs([eVec(1, 1), eVec(2, 1), 0]);
                c = real([eVec(1, 1), eVec(2, 1), 0]);
                d = [1, 0, 0];
                eD2(2, wFrame) = atan2(norm(cross(c,d)), dot(c,d)) * (360/(2*pi));
%                 eD2(1, wFrame) = atan2(eVec(1, 1), eVec(2, 1)) * (360/(2*pi)) + 90;
%                 extremeDirections(cellNum, frameNum, 1) = atan2(eVec(1, 2), eVec(2, 2)) * (360/(2*pi)) + 90;
%         %         extremeDirections(cellNum, frameNum, 2) = atan2(eVec(1, 1), eVec(2, 1)) * (360/(2*pi)) + 90;
            end
            extremeDirections(cellNum, frameNum, 1) = nanmean(eD2(1,:), 2);
            extremeDirections(cellNum, frameNum, 2) = nanmean(eD2(2,:), 2);
%             extremeDirections(cellNum, frameNum, 1) = atan2(eVec(1, 1), eVec(2, 1)) * (360/(2*pi)) + 90;
% %             extremeDirections(cellNum, frameNum, 1) = atan2(eVec(1, 2), eVec(2, 2)) * (360/(2*pi)) + 90;
        end
    end
end
% firstMaxFrame = find(maxOrientationAngles(1,:) == max(maxOrientationAngles(1,:)), 1, 'first');
% firstMinFrame = find(maxOrientationAngles(2,:) == min(maxOrientationAngles(2,:)), 1, 'first');
% meanPostMaxFrame = mean(maxOrientationAngles(1, firstMaxFrame:end));
% meanPostMinFrame = mean(maxOrientationAngles(2, firstMinFrame:end));
% 
% cardinalAxesOffset = [meanPostMaxFrame, meanPostMinFrame];

% filename = strcat('angle_vs_time_', version, '_frames_1_to_250.mat');
% filename = strcat('angle_vs_time_', version, '.mat');
% save(filename, 'theta_vs_time', 'angularMomenta', 'movieName', 'version', 'thetas', 'numFrames', 'max_combo_range', 'min_combo_range', 'maxOrientationAngles', 'cardinalAxesOffset');

%% Sine-Fitting raw data (obsolete - SW 10/29/2019)
% extreme_AngularMomenta = NaN(2, numFrames);
% % extreme_Elongations = NaN(2, numFrames);
% extreme_Elongation = NaN(1, numFrames);
% 
% ww = 10; % Window Width in units of frames
% fitFunc = NaN(length(thetas), numFrames);
% extreme_locations = NaN(2, numFrames); %% Line Max. Height (row 1) and Line of Min. Width (row 2) measurements
% extreme_measurements = NaN(2, numFrames); %% Max. Height (row 1) and Min. Width (row 2) measurements
% fit_offsets = NaN(2, numFrames); %% vertical (y) offset (row 1) and phase (phi) offset (row 2)
% % sinfit = NaN(length(thetas));
% 
% for frameNum = 1:numFrames
%     if frameNum <= ww
% 
%         bg_guess = mean(mean(theta_vs_time(:, 1:frameNum, 1), 2));
%         
%             dMax = max(mean(theta_vs_time(:, 1:frameNum, 1), 2));  % Max value in the data
%             dMin = min(mean(theta_vs_time(:, 1:frameNum, 1), 2));  % Min value in the data
%         
%         A_guess = ( dMax - dMin ) / 2;
%         omega_guess = 1;
% 
%             min_peak_location = min(thetas(fitFunc(:, frameNum) == dMax), thetas(fitFunc(:, frameNum) == dMin));
%             max_peak_location = max(thetas(fitFunc(:, frameNum) == dMax), thetas(fitFunc(:, frameNum) == dMin));
%             peak_separation = max_peak_location - min_peak_location;
%             midpoint = min_peak_location + peak_separation / 2;
%         
%         phi_guess = abs(pi - midpoint);
%         
% %         guess = [bg_guess, A_guess, omega_guess, phi_guess];
%         guess = [bg_guess A_guess phi_guess omega_guess];
%         fix = [0 0 0 0];
% %         lb = 
% %         constrain = [lb ub];
%         
%         [sine_fit, fitFunc(:, frameNum)] = fitSine(thetas', mean(theta_vs_time(:, 1:frameNum, 1), 2), guess, fix);
%         
%         fit_offsets(1, frameNum) = sine_fit(1);
%         fit_offsets(2, frameNum) = sine_fit(3) * 360 / (2*pi); %%% Phase offset in Degrees
% %         fit_offsets(2, frameNum) = sine_fit(3); %%% Phase offset in Radians
% 
%             fMax = max(fitFunc(:, frameNum));
%             fMin = min(fitFunc(:, frameNum));
%         
% %         extreme_locations(1, frameNum) = thetas(fitFunc(:, frameNum) == fMax) * 360 / (2*pi);
% %         extreme_locations(2, frameNum) = thetas(fitFunc(:, frameNum) == fMin) * 360 / (2*pi);
%         extreme_locations(1, frameNum) = sine_fit(3) * 360 / (2*pi);
%         extreme_locations(2, frameNum) = (sine_fit(3) * 360 / (2*pi)) + 90;
%         
% % %         extreme_measurements(1, frameNum) = fMax;
% % %         extreme_measurements(2, frameNum) = fMin;
%         extreme_measurements(1, frameNum) = dMax;
%         extreme_measurements(2, frameNum) = dMin;
% % 
% %         ellipseThetaInc = 5;
% %         [~, ~, ~, ~, stdDevs] = ellipseStandard(dMax, dMin, ellipseThetaInc);
% %         extreme_measurements(1, frameNum) = stdDevs(1);
% %         extreme_measurements(2, frameNum) = stdDevs(2);
%         
%         
%     else
%         
%         bg_guess = mean(mean(theta_vs_time(:, frameNum - ww: frameNum, 1), 2));
%         
%             dMax = max(mean(theta_vs_time(:, frameNum - ww: frameNum, 1), 2));  % Max value in the data
%             dMin = min(mean(theta_vs_time(:, frameNum - ww: frameNum, 1), 2));  % Min value in the data
%         
%         A_guess = ( dMax - dMin ) / 2;
%         omega_guess = 1;
% 
%             min_peak_location = min(thetas(fitFunc(:, frameNum) == dMax), thetas(fitFunc(:, frameNum) == dMin));
%             max_peak_location = max(thetas(fitFunc(:, frameNum) == dMax), thetas(fitFunc(:, frameNum) == dMin));
%             peak_separation = max_peak_location - min_peak_location;
%             midpoint = min_peak_location + peak_separation / 2;
%         
%         phi_guess = abs(pi - midpoint);
%         
% %         guess = [bg_guess, A_guess, omega_guess, phi_guess];
%         guess = [bg_guess A_guess phi_guess omega_guess];
%         fix = [0 0 0 0];
%         constrain = [];
%         
%         [sine_fit, fitFunc(:, frameNum)] = fitSine(thetas', mean(theta_vs_time(:, frameNum - ww: frameNum, 1), 2), guess, fix, constrain);
%         
%         fit_offsets(1, frameNum) = sine_fit(1);  %%% Background fit paramater
%         fit_offsets(2, frameNum) = sine_fit(3) * 360 / (2*pi); %%% Phase offset in Degrees
% %         fit_offsets(2, frameNum) = sine_fit(3); %%% Phase offset in Radians
% 
%             fMax = max(fitFunc(:, frameNum)); %% Max value of the fitted function
%             fMin = min(fitFunc(:, frameNum)); %% Min value of the fitted function
%         
%         extreme_locations(1, frameNum) = sine_fit(3) * 360 / (2*pi);
%         extreme_locations(2, frameNum) = (sine_fit(3) * 360 / (2*pi)) + 90;
% %         extreme_locations(1, frameNum) = thetas(fitFunc(:, frameNum) == fMax) * 360 / (2*pi);
% %         extreme_locations(2, frameNum) = thetas(fitFunc(:, frameNum) == fMin) * 360 / (2*pi);
% %         extreme_locations(1, frameNum) = thetas(fitFunc(:, frameNum) == sine_fit(1) + sine_fit(2)) * 360 / (2*pi);
% %         extreme_locations(2, frameNum) = thetas(fitFunc(:, frameNum) == sine_fit(1) - sine_fit(2)) * 360 / (2*pi);
%         
% % % %         extreme_measurements(1, frameNum) = fMax;
% % % %         extreme_measurements(2, frameNum) = fMin;
% % %         extreme_measurements(1, frameNum) = sine_fit(1) + sine_fit(2);
% % %         extreme_measurements(2, frameNum) = sine_fit(1) - sine_fit(2);
%         extreme_measurements(1, frameNum) = dMax;
%         extreme_measurements(2, frameNum) = dMin;
% %         ellipseThetaInc = 5;
% %         
% % %         [~, ~, ~, ~, extreme_measurements(1, frameNum)] = circleStandard(dMax, circleThetaInc);
% % %         [~, ~, ~, ~, extreme_measurements(2, frameNum)] = circleStandard(dMin, circleThetaInc);
% %         
% % 
% % %% Need to Load Actual PixelMap Data and use the maximum measurements in each direction as dMax and dMin, with the larger being max and smaller the min.
% %        [~, ~, ~, ~, stdDevs] = ellipseStandard(dMax, dMin, ellipseThetaInc);
% %        extreme_measurements(1, frameNum) = stdDevs(1);
% %        extreme_measurements(2, frameNum) = stdDevs(2);
%         
%     end
% end

%% Find Cardinal Axes of Movie:
cardinalAxesOffset = NaN(2, 1);
% cardinalAxesOffset(1, 1) = mean(fit_offsets(2, :));
cardinalAxesOffset(1, 1) = nanmean(nanmean(extremeDirections(:, :, 1), 1));

if cardinalAxesOffset(1, 1) >= 90
    cardinalAxesOffset(2, 1) = cardinalAxesOffset(1, 1) - 90;
else
    cardinalAxesOffset(2, 1) = cardinalAxesOffset(1, 1) + 90;
end

%% Master Save:
filename = strcat('angle_vs_time_equivRectangle_', version, '.mat');
% save(filename, 'theta_vs_time', 'angularMomenta', 'movieName', 'version', 'thetas', 'numFrames', 'fitFunc', 'extreme_locations', 'extreme_measurements', 'fit_offsets', 'cardinalAxesOffset');
save(filename, 'theta_vs_time', 'angularMomenta', 'movieName', 'version', 'thetas', 'numFrames', 'equivalentRectangleMeasurements', 'extremeDirections', 'celAvgExtremeDirections', 'cardinalAxesOffset');
%%% If file is too big to save, use following 2 lines in addition to above2
% filename1 = strcat('angle_vs_time_', version, '_newSineFit2.mat');
% save(filename1, 'theta_vs_time', 'angularMomenta', 'movieName', 'version', 'thetas', 'numFrames', 'fitFunc', 'extreme_locations', 'extreme_measurements', 'fit_offsets', 'cardinalAxesOffset', '-v7.3');
% % filename2 = strcat('angle_vs_time_', version, '_Lonly_newSineFit1.mat');
% % save(filename2, 'angularMomenta', '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Plots
% figure(); plot(thetas', theta_vs_time(:, frameNum, 1)); hold on; plot(fourier_fit);
% hold on; plot(thetas', fitFunc(:, frameNum));
% hold on; scatter(thetas(fitFunc(:, frameNum) == fMax), fMax); hold on; scatter(thetas(fitFunc(:, frameNum) == fMin), fMin)
% hold on; scatter(fit_offsets(2, frameNum), fit_offsets(1, frameNum))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PLOTS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
maxPlotTime = numFrames;
timeInc = 1;

suptitle({strcat('Movie: ', movieName, ' - Cardinal Axes: ', ' ', num2str(round(cardinalAxesOffset(1, 1), 1)), ' x ', ' ', num2str(round(cardinalAxesOffset(2, 1), 1)), ' degrees'), ''});

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



% figure();
% fig_01 = subplot(1,1,1);
fig11 = subplot(4,1,2);
% domain = 1:size(extreme_locations, 2);
domain = 1:numFrames;

% % plot(domain, extreme_locations(1,:), 'r', domain, extreme_locations(2,:), 'b');
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
fig13 = subplot(4,1,4);
% plot(domain, extreme_measurements(1,:) ./ extreme_measurements(2,:), 'k');
plot(domain, nanmean(equivalentRectangleMeasurements(:, :, 3), 1), 'k');
xlabel('time (frames)');
ylabel('Elongation');
if strcmp(movieName, 'r4')
    ylim([0 4]);
end
title(strcat('Elongation at Extreme Angles, movie: ', movieName));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linkaxes([fig10 fig11 fig12 fig13],'xy');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % numCells = size(ellipseProperties, 1);        %%% Plot all the cells
% numCells = 231;                        %%% 1_wtYiActin_movie2 ONLY         %%% Number of cells to plot
% 
% fig11 = subplot(2,1,2);
% imagesc(ellipseProperties(1:numCells, :, 4));
% ylabel('Cell Index');
% xlabel('time (frames)');
% title('Oritentation Angles');
% % cbh=colorbar('h');
% % set(cbh,'YTick',['NaN', -90:30:90]);
% colorbar;
% colormap(gca,mymap);


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
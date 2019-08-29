function [ inertiaTensors ] = singleFrame_tensorCalculator( x_verts, y_verts )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     numCells = size(nodeCellMat,2);
    inertiaTensors = cell(length(x_verts), 1);
    centroids = zeros(length(x_verts), 2);
    for cellNum = 1:length(x_verts)
        % Find this cell's centroid:
        centroids(cellNum, :) = [mean(x_verts{cellNum}), mean(y_verts{cellNum})];
        % Construct the Moment of Inertia Tensor for this polygon-cell by
        % adding the contributin of each triangle:
        %%% First, find centroids of each triangle contructed by connecting
        %%% each vertex with the centroid. An n-gon has n triangles.
        %%%%%% Then, calculate m_alpha: the effective  "mass" or area of 
        %%%%%% each triangle. 
        tri_centroids = zeros(length(x_verts{cellNum}), 2);
        tri_areas = zeros(length(x_verts{cellNum}), 1);
        inertiaTensor = zeros(2,2);
        for vertNum = 1:length(x_verts{cellNum})
%             base = sqrt((x_verts{cellNum}(vertNum) - centroids(cellNum, 1))^2 + (y_verts{cellNum}(vertNum) - centroids(cellNum, 2))^2);
            A = [x_verts{cellNum}(vertNum) - centroids(cellNum, 1), y_verts{cellNum}(vertNum) - centroids(cellNum, 2)];
            base = sqrt(A(1)^2 + A(2)^2);
            
            if vertNum < length(x_verts{cellNum})
                tri_centroids(vertNum, 1) = mean([centroids(cellNum, 1), x_verts{cellNum}(vertNum), x_verts{cellNum}(vertNum + 1)]);
                tri_centroids(vertNum, 2) = mean([centroids(cellNum, 2), y_verts{cellNum}(vertNum), y_verts{cellNum}(vertNum + 1)]);

                B = [x_verts{cellNum}(vertNum + 1) - centroids(cellNum, 1), y_verts{cellNum}(vertNum + 1) - centroids(cellNum, 2)];
                projection = (dot(A,B)/norm(B)^2)*B;
                rejection = A - projection;
                height = sqrt(rejection(1)^2 + rejection(2)^2);
                tri_areas(vertNum) = base * height / 2;
                
                R = [tri_centroids(vertNum, 1) - centroids(cellNum, 1), tri_centroids(vertNum, 2) - centroids(cellNum, 2)];
                
                inertiaTensor = inertiaTensor + [
                    tri_areas(vertNum) * R(2)^2, -tri_areas(vertNum) * R(1) * R(2);
                    -tri_areas(vertNum) * R(2) * R(1), tri_areas(vertNum) * R(1)^2 ];
                
            elseif vertNum == length(x_verts{cellNum})
                tri_centroids(vertNum, 1) = mean([centroids(cellNum, 1), x_verts{cellNum}(vertNum), x_verts{cellNum}(1)]);
                tri_centroids(vertNum, 2) = mean([centroids(cellNum, 2), y_verts{cellNum}(vertNum), y_verts{cellNum}(1)]);

                B = [x_verts{cellNum}(1) - centroids(cellNum, 1), y_verts{cellNum}(1) - centroids(cellNum, 2)];
                projection = (dot(A,B)/norm(B)^2)*B;
                rejection = A - projection;
                height = sqrt(rejection(1)^2 + rejection(2)^2);
                tri_areas(vertNum) = base * height / 2;
                
                R = [tri_centroids(vertNum, 1) - centroids(cellNum, 1), tri_centroids(vertNum, 2) - centroids(cellNum, 2)];
                
                inertiaTensor = inertiaTensor + [
                    tri_areas(vertNum) * R(2)^2, -tri_areas(vertNum) * R(1) * R(2);
                    -tri_areas(vertNum) * R(2) * R(1), tri_areas(vertNum) * R(1)^2 ];
            else
                disp(['vertex error @ cellNum = ', num2str(cellNum), ", vertNum = ", num2str(vertNum)])
            end
        end
        inertiaTensors{cellNum} = inertiaTensor;
    end
end


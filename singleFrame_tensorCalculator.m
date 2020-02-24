function [ inertiaTensors ] = singleFrame_tensorCalculator( x_verts, y_verts )
% This function computes the inertia tensor for all the cells in a single
% movie frame. 
%
%   First, for each cell, the centroid is calculated and N triangles are
%   specified, where N is the number of vertices for the cell. The triangles
%   are created by connecting the cell-centroid to each vertex. 
%
%   Next, the triangle-centroids are triangle surface areas are calculated. 
%   In this formalism, the surface area stands in for the mass of that 
%   triangle, effectively concentrated at the appropriate triangle-centroid. 
%
%   Finally, the vector R connecting cell-centroid to triangle-centroid
%   is used to calculate the inertial contribution of each triangle to the
%   inertia tensor of the entire cell. 

    inertiaTensors = cell(length(x_verts), 1);
    centroids = zeros(length(x_verts), 2);
    for cellNum = 1:length(x_verts)
        
        % Find this cell's centroid
        centroids(cellNum, :) = [mean(x_verts{cellNum}), mean(y_verts{cellNum})];
        
        % Create containers for this cell's triangle-centroids, triangle 
        % surface areas, and inertia tensor. 
        tri_centroids = zeros(length(x_verts{cellNum}), 2);
        tri_areas = zeros(length(x_verts{cellNum}), 1);
        inertiaTensor = zeros(2,2);
        
        % Cycle through triangles by cycling through vertices
        for vertNum = 1:length(x_verts{cellNum})
            A = [x_verts{cellNum}(vertNum) - centroids(cellNum, 1), y_verts{cellNum}(vertNum) - centroids(cellNum, 2)];
            base = sqrt(A(1)^2 + A(2)^2);
            
            if vertNum < length(x_verts{cellNum})
                % Find triangle centroids
                tri_centroids(vertNum, 1) = mean([centroids(cellNum, 1), x_verts{cellNum}(vertNum), x_verts{cellNum}(vertNum + 1)]);
                tri_centroids(vertNum, 2) = mean([centroids(cellNum, 2), y_verts{cellNum}(vertNum), y_verts{cellNum}(vertNum + 1)]);

                % Find triangle surface areas (~ mass)
                B = [x_verts{cellNum}(vertNum + 1) - centroids(cellNum, 1), y_verts{cellNum}(vertNum + 1) - centroids(cellNum, 2)];
                projection = (dot(A,B)/norm(B)^2)*B;
                rejection = A - projection;
                height = sqrt(rejection(1)^2 + rejection(2)^2);
                tri_areas(vertNum) = base * height / 2;
                
                % Construct vector pointing from cell-centroid to
                % triangle-centroid
                R = [tri_centroids(vertNum, 1) - centroids(cellNum, 1), tri_centroids(vertNum, 2) - centroids(cellNum, 2)];
                
                % Construct inertia tensor
                inertiaTensor = inertiaTensor + [
                    tri_areas(vertNum) * R(2)^2, -tri_areas(vertNum) * R(1) * R(2);
                    -tri_areas(vertNum) * R(2) * R(1), tri_areas(vertNum) * R(1)^2 ];
                
            % Alternatively, the final vertex forms a triangle with the
            % first vertex, hence the slightly different notation below:
            elseif vertNum == length(x_verts{cellNum})
                
                % Find triangle centroids
                tri_centroids(vertNum, 1) = mean([centroids(cellNum, 1), x_verts{cellNum}(vertNum), x_verts{cellNum}(1)]);
                tri_centroids(vertNum, 2) = mean([centroids(cellNum, 2), y_verts{cellNum}(vertNum), y_verts{cellNum}(1)]);

                % Find triangle surface areas (~ mass)
                B = [x_verts{cellNum}(1) - centroids(cellNum, 1), y_verts{cellNum}(1) - centroids(cellNum, 2)];
                projection = (dot(A,B)/norm(B)^2)*B;
                rejection = A - projection;
                height = sqrt(rejection(1)^2 + rejection(2)^2);
                tri_areas(vertNum) = base * height / 2;
                
                % Construct vector pointing from cell-centroid to
                % triangle-centroid
                R = [tri_centroids(vertNum, 1) - centroids(cellNum, 1), tri_centroids(vertNum, 2) - centroids(cellNum, 2)];
                
                % Construct inertia tensor
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


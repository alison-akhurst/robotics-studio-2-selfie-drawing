clear; close all; clc; rng('default');

% Setup webcam and capture image
webcamlist;
cam = webcam('Integrated Camera');
preview(cam);
cam.AvailableResolutions;
cam.Resolution = '960x540';
pause(2);
I = snapshot(cam);
imshow(I);

% Face detection
faceDetector = vision.CascadeObjectDetector('FrontalFaceCART'); 
faceDetector.MinSize = [210 210];
faceDetector.MaxSize = [600 600];
bboxes = step(faceDetector, I);
Iout = insertObjectAnnotation(I, 'rectangle', bboxes, 'Face');
figure, imshow(Iout);

% Crop the face in the image
Iface = imcrop(I, bboxes(1, :));
figure, imshow(Iface);

% Edge detection
[BW, ~] = edge(rgb2gray(Iface), 'Canny', [0.0813 0.1281]);
figure, imshow(BW);

% Skeletonization
BW2 = bwmorph(BW, 'skel', Inf);
figure, imshow(BW2);

% Remove branches
BW3 = bwmorph(BW2, 'spur', 3);
figure, imshow(BW3);

% Mask branch points
branchPoints = bwmorph(BW3, 'branch', 1); 
branchPoints = imdilate(branchPoints, strel('disk', 1));
BW3 = BW3 & ~branchPoints;
figure, imshow(BW3);

% Remove small regions
BWseg = bwareaopen(BW3, 10);
figure, imshow(BWseg);

% Detect start and end points and prepare for G-code conversion
[B, L] = bwboundaries(BWseg, 'noholes');
figure, imshow(label2rgb(L, @jet, [.5 .5 .5]));
hold on;
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'LineWidth', 3);
end

% Convert to 3D coordinates for G-code
origin = [0, 0];  % Adjust the origin to (0,0)
delta = 0.001; 
B2 = cell(size(B));
for i = 1:length(B)
    b = B{i};
    bx = -b(:,2) * delta + origin(1);
    by = b(:,1) * delta + origin(2);
    bz = zeros(length(bx), 1);
    B2{i} = [bx, by, bz];
    plot3(bx, by, bz, 'LineWidth', 2); hold on;
end
grid on;
hold off;

% Apply 180-degree rotation
theta = pi; % 180 degrees in radians
rotationMatrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];

for i = 1:length(B2)
    B2{i}(:,1:2) = (rotationMatrix * B2{i}(:,1:2)')';
end

% Calculate translation to set bottom left to origin
all_points = vertcat(B2{:});
min_x = min(all_points(:, 1));
min_y = min(all_points(:, 2));

for i = 1:length(B2)
    B2{i}(:,1) = B2{i}(:,1) - min_x;
    B2{i}(:,2) = B2{i}(:,2) - min_y;
end

% Scale to 80% of A4 size
scaleFactor = 0.8;
a4_width_cm = 21.0 * scaleFactor;
a4_height_cm = 29.7 * scaleFactor;

figure;
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, a4_width_cm, a4_height_cm]);
set(gca, 'Units', 'centimeters', 'Position', [1, 1, a4_width_cm - 2, a4_height_cm - 2]);

% Initialize cell array to store boundary points
boundaryPoints = {};

% Plot the boundaries to fit within the scaled A4 dimensions and store points
for i = 1:length(B2)
    boundary = B2{i};
    scaledBoundary = [boundary(:,1) * 100 * scaleFactor, boundary(:,2) * 100 * scaleFactor]; % Convert to centimeters and scale
    plot(scaledBoundary(:,1), scaledBoundary(:,2), 'LineWidth', 2); hold on;
    boundaryPoints{i} = scaledBoundary; % Store points in centimeters
end
axis([0 a4_width_cm 0 a4_height_cm]);
axis equal;
grid on;

% Save the scaled figure as an image
saveas(gcf, 'output_image.png');

% Save the boundary points to a .mat file
save('boundary_points.mat', 'boundaryPoints');

% Generate G-code from the processed boundaries
generateGCode(boundaryPoints, scaleFactor);

% Function to generate G-code from boundary data
function generateGCode(boundaryPoints, scaleFactor)
    fileName = 'output.gcode';
    fid = fopen(fileName, 'w');
    
    fprintf(fid, 'G21 ; Set units to mm\nG90 ; Use absolute positioning\n');
    liftHeight = 5;
    engravingDepth = -0.5;
    feedRate = 1000;
    
    fprintf(fid, 'G0 Z%f\n', liftHeight);
    
    % Scale factor to convert from cm to mm
    cmToMm = 10;
    
    for i = 1:length(boundaryPoints)
        boundary = boundaryPoints{i};
        fprintf(fid, 'G0 X%f Y%f\n', boundary(1,1) * cmToMm, boundary(1,2) * cmToMm);
        fprintf(fid, 'G1 Z%f F%d\n', engravingDepth, feedRate / 2);

        for j = 1:size(boundary, 1)
            fprintf(fid, 'G1 X%f Y%f F%d\n', boundary(j,1) * cmToMm, boundary(j,2) * cmToMm, feedRate);
        end
        
        fprintf(fid, 'G0 Z%f\n', liftHeight);
    end

    fprintf(fid, 'M02 ; End of program\n');
    fclose(fid);
    disp(['G-code saved to ', fileName]);
end

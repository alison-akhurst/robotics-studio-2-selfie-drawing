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
Iout = insertObjectAnnotation(I,'rectangle',bboxes,'Face');
figure, imshow(Iout);

% Crop the face in the image
Iface = imcrop(I,bboxes(1,:));
figure, imshow(Iface);

% Edge detection
[BW,~] = edge(rgb2gray(Iface),'Canny',[0.0813 0.1281]);
figure, imshow(BW);

% Skeletonization
BW2 = bwmorph(BW,'skel',Inf);
figure, imshow(BW2);

% Remove branches
BW3 = bwmorph(BW2,'spur',3);
figure, imshow(BW3);

% Mask branch points
branchPoints = bwmorph(BW3,'branch',1); 
branchPoints = imdilate(branchPoints,strel('disk',1));
BW3 = BW3 & ~branchPoints;
figure, imshow(BW3);

% Remove small regions
BWseg = bwareaopen(BW3,10);
figure, imshow(BWseg);

% Detect start and end points and prepare for G-code conversion
[B,L] = bwboundaries(BWseg,'noholes');
figure, imshow(label2rgb(L, @jet, [.5 .5 .5]));
hold on;
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'LineWidth', 3);
end

% Convert to 3D coordinates for G-code
origin = [0.48, -0.229/2]; 
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

% Generate G-code from the processed boundaries
generateGCode(B2);

% Function to generate G-code from boundary data
function generateGCode(B2)
    fileName = 'output.gcode';
    fid = fopen(fileName, 'w');
    
    fprintf(fid, 'G21 ; Set units to mm\nG90 ; Use absolute positioning\n');
    liftHeight = 5;
    engravingDepth = -0.5;
    feedRate = 1000;
    
    fprintf(fid, 'G0 Z%f\n', liftHeight);
    
    for i = 1:length(B2)
        boundary = B2{i};
        fprintf(fid, 'G0 X%f Y%f\n', boundary(1,1), boundary(1,2));
        fprintf(fid, 'G1 Z%f F%d\n', engravingDepth, feedRate / 2);

        for j = 1:size(boundary, 1)
            fprintf(fid, 'G1 X%f Y%f F%d\n', boundary(j,1), boundary(j,2), feedRate);
        end
        
        fprintf(fid, 'G0 Z%f\n', liftHeight);
    end

    fprintf(fid, 'M02 ; End of program\n');
    fclose(fid);
    disp(['G-code saved to ', fileName]);
end

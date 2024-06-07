clc;
close all;
clear all;

% Setup webcam and capture image
webcamlist;
cam = webcam('Integrated Camera');
preview(cam);
cam.AvailableResolutions;
cam.Resolution = '960x540';
pause(5);
img = snapshot(cam);
imshow(img);

%Computer Vision System toolbox

faceDetector = vision.CascadeObjectDetector('FrontalFaceCART'); 
faceDetector.MinSize = [210 210];
faceDetector.MaxSize = [315 315];
bboxes = step(faceDetector, img);
Iout = insertObjectAnnotation(img,'rectangle',bboxes,'Face');
figure, imshow(Iout);

Iface = imcrop(img,bboxes(1,:));
figure, imshow(Iface); 

% Convert to grayscale and smooth the image
grayImg = rgb2gray(Iface);
smoothedImg = imgaussfilt(grayImg, 0.5); % Apply Gaussian filter for smoothing
imshow(smoothedImg);

% Edge detection on the smoothed grayscale image
edgeImg = edge(smoothedImg, 'canny', [0.0813 0.1281]); % Adjust the Canny threshold as needed

%skeletonisation

BW2 = bwmorph(edgeImg,'skel',Inf);
figure, imshow(BW2);

BW3 = bwmorph(BW2,'spur',3);
figure, imshow(BW3);

branchPoints = bwmorph(BW3,'branch',1); 
branchPoints = imdilate(branchPoints,strel('disk',1));
BW3 = BW3 & ~branchPoints;
figure, imshow(BW3)

edgeImgDilated = bwareaopen(BW3,10);
figure, imshow(edgeImgDilated)

% Optional: Enhance edges with morphological operation
% se = strel('line', 1, 50);
% edgeImgDilated = imdilate(edgeImg, se);

% Get connected components from the binary edge image
 cc = bwconncomp(edgeImgDilated);

% Initialize figure for plotting
% figure;
% imshow(edgeImg);
% hold on;

% Initialize a string to hold the G-code
% Setup and initial image processing code here...
% [Omitted for brevity - include the webcam capture, edge detection, and processing logic before this]

% Open a G-code file for writing
gcodeFileName = 'output.gcode';
fid = fopen(gcodeFileName, 'w');

% Initial setup commands (adjust as needed)
fprintf(fid, 'G90 ; Set to Absolute Positioning\n');
fprintf(fid, 'G21 ; Set to Millimeter Units\n');
fprintf(fid, 'G92 X0 Y0 Z0 ; Set Current Position to Origin\n');

% Set the feed rate (adjust as needed)
feedRate = 100; % in mm/min
fprintf(fid, 'F%d ; Set Feed Rate\n', feedRate);

% Assuming a starting position might be necessary before the loop
% For simplicity, this example starts from the origin

% Process connected components for G-code generation
for idx = 1:cc.NumObjects
    [y, x] = ind2sub(size(edgeImg), cc.PixelIdxList{idx});

    if length(x) > 1
        try
            ft = fittype('smoothingspline');
            curve = fit(x, y, ft, 'SmoothingParam', 0.2);
            xx = linspace(min(x), max(x), 100); % Fewer points for G-code simplicity
            yy = feval(curve, xx);
            
            % Rapid move to the starting point of the curve
            fprintf(fid, 'G0 X%.2f Y%.2f ; Rapid Move to Start of Curve\n', xx(1), yy(1));
            
            % Move through the curve points
            for i = 1:length(xx)
                fprintf(fid, 'G1 X%.2f Y%.2f ; Move through Curve\n', xx(i), yy(i));
            end
        catch ME
            disp(['Error processing component ', num2str(idx), ': ', ME.message]);
        end
    elseif length(x) == 1
        % Handle a single point by moving directly to it
        fprintf(fid, 'G0 X%.2f Y%.2f ; Move to Single Point\n', x, y);
    end
end

% Optionally return to origin or lift the tool at the end
fprintf(fid, 'G0 Z1 ; Lift tool above the workpiece\n');
fprintf(fid, 'G0 X0 Y0 ; Return to origin\n');
fprintf(fid, 'M2 ; End of Program\n');

% Close the file
fclose(fid);

disp(['G-code generated and saved as ', gcodeFileName]);


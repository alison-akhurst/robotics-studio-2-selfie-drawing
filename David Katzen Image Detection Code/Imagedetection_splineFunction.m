
clc;
close all;
clear all;

% Setup webcam and capture image
webcamlist;
cam = webcam('Integrated Camera');
preview(cam);
cam.AvailableResolutions;
cam.Resolution = '960x540';
pause(3);
img = snapshot(cam);
imshow(img);

%Computer Vision System toolbox

faceDetector = vision.CascadeObjectDetector('FrontalFaceCART'); 
faceDetector.MinSize = [210 210];
faceDetector.MaxSize = [500 500];

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

branchPoints = bwmorph(BW3,'branchpoints',1); 
branchPoints = imdilate(branchPoints,strel('disk',2));
BW3 = BW3 & ~branchPoints;
figure, imshow(BW3);

edgeImgDilated = bwareaopen(BW3,10);
figure, imshow(edgeImgDilated);

% Get connected components of the binary edge image
cc = bwconncomp(edgeImgDilated);

% Assuming previous steps for setting up webcam and image processing remain unchanged

% Assuming previous steps for setting up webcam and image processing remain unchanged

% Loop through each connected component to fit smooth curves or handle small segments appropriately
figure;
imshow(edgeImgDilated);
hold on;

for idx = 1:cc.NumObjects
    % Get pixel indices and convert them to subscripts
    [y, x] = ind2sub(size(edgeImgDilated), cc.PixelIdxList{idx});

    if length(x) >= 2 % Ensure there are enough points for fitting
        % Fit a smooth curve to the points
        ft = fittype('smoothingspline');
        try
            curve = fit(x, y, ft, 'SmoothingParam', 1); % Adjust SmoothingParam as needed

            % Generate points from the fitted curve for plotting
            xx = linspace(min(x), max(x), 1000); % Generate a range of x values for smooth plotting
            yy = feval(curve, xx); % Evaluate the fitted curve at these x values

            plot(xx, yy, 'r', 'LineWidth', 1); % Draw the smooth curve
        catch
            % Handle cases with insufficient data points (less than 2) or other fitting errors
            fprintf('Skipping a component with insufficient points or fitting error at index %d.\n', idx);
        end
    else
        % For a single point, you could plot it or decide to skip
        % plot(x, y, 'r.', 'MarkerSize', 15); % Example: Plotting the point
        fprintf('Skipping a single-point component at index %d.\n', idx);
    end
end

hold off;

% Since the next steps involve generating G-code based on the processed image, and the handling of connected components
% and the drawing of curves have been adjusted, ensure the method of interpreting these curves into G-code
% accounts for the nature of the data being processed (e.g., smooth curves or straight lines between points).

% The G-code generation part would then proceed here, with considerations for how to translate the image analysis
% results into movements and operations for a CNC machine, 3D printer, or similar device.

% Note: This section is context-dependent and requires additional information about the target machine and the
% intended use of the image data for precise instructions.

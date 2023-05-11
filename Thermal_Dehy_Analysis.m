clearvars
close all

% Identify the folder with data
path = uigetdir(path);
% Find all .dat files in folder
pathcontent = dir([path '/*.dat']);

% Number of .dat files in folder
curvelength = size(pathcontent,1);

% Create matrices based on number of .dat files
finalcurve = zeros(curvelength,1);
dehydration_Lesion = zeros(curvelength,1);
dehydration_Control = zeros(curvelength,1);
dQ_Lesion = strings([curvelength,1]);
dQ_Control = strings([curvelength,1]);
Diff_dQ = strings([curvelength,1]);

% Setting up calling string
name = char({pathcontent(1).name});
ext = name(1:end-5);

% Find Absolute Minimum
disp('Calculating Absolute Minimum...');
img_size = load([path '/' ext num2str(1) '.dat']);
stack_raw = zeros(size(img_size,1),size(img_size,2),curvelength);
for loopfolder = 1:curvelength
    % load individual .dat file
    stack_raw(:,:,loopfolder) = load([path '/' ext num2str(loopfolder) '.dat']);
end
absmin = min(stack_raw,[],'all');
disp(['Absolute Minimum: ' num2str(absmin)]);

% Make Image Visible, Draw ROI on 1st Image
img_raw = load([path '/' ext num2str(1) '.dat']);
img_org = (img_raw-absmin)*100; % dehydration intensity was gather with *0.01
img8 = uint8(img_org);
first_img = img8;

%Draw ROI in imfreehand and get ROI info
fontSize = 16;
imshow(img8, []);
axis on;
title('Original Image', 'FontSize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

%%%%% DRAW TOOTH %%%%%
dehydration_tooth = zeros(size(img8,1),size(img8,2));
% Ask user to draw freehand mask.
message_tooth = sprintf('Left click and hold to draw out TOOTH.\nSimply lift the mouse button to finish.');
uiwait(msgbox(message_tooth));
hFH_tooth = imfreehand(); % Actual line of code to do the drawing.
% Create a binary image ("mask") from the ROI object.
binaryImage_tooth = hFH_tooth.createMask();
xy_tooth = hFH_tooth.getPosition;

% Calculate the area, in pixels, that they drew.
numberOfPixels1_tooth = sum(binaryImage_tooth(:));
% Another way to calculate it that takes fractional pixels into account.
numberOfPixels2_tooth = bwarea(binaryImage_tooth);

% Get coordinates of the boundary of the freehand drawn region.
structBoundaries_tooth = bwboundaries(binaryImage_tooth);
xy_tooth = structBoundaries_tooth{1}; % Get n by 2 array of x,y coordinates.
x_tooth = xy_tooth(:, 2); % Columns.
y_tooth = xy_tooth(:, 1); % Rows.

close % Closes drawing window.

% Calculate Temperature Values for Each Pixel
stack = zeros(size(img8,1),size(img8,2),curvelength);

for loopfolder = 1:curvelength

disp(['Analyzing frame ' num2str(loopfolder) ' ...']);
% load individual .dat file
img_raw = load([path '/' ext num2str(loopfolder) '.dat']);
temp_data = (img_raw*100)*160/16384; %%% Revert 0.01 by x100, 160K Range (248K - 408K) on 14-Bit Intensity %%%

% Mask the images outside the mask, and display it.
% Will keep only the part of the image that's inside the mask, zero outside mask.
blackMaskedImage_tooth = temp_data;
blackMaskedImage_tooth(~binaryImage_tooth) = 0;

% Insert image into stack
stack(:,:,loopfolder) = blackMaskedImage_tooth;

end

% Find max of each pixel
tempmax = max(stack,[],3);
% Difference between Max throughout dehydration and curve for each pixel
tempdiff = tempmax - stack;
% Integral of difference for each pixel
sumtempdiff = sum(tempdiff,3);

temp_norm = sumtempdiff/max(max(sumtempdiff));%%%Normalize to max intensity
temp_norm16 = temp_norm*2^16-1; %%%16bit (0-65535)
temp_Plot = uint16(temp_norm16);

% Plot Map
fontSize = 16;
imshow(temp_Plot, []);
% colorbar;
% colormap(jet(65536));
axis on;
title('Dehydration Map', 'FontSize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

%%%%% DRAW LESION %%%%%
% Ask user to draw freehand mask.
message_Lesion = sprintf('Left click and hold to draw for LESION.\nSimply lift the mouse button to finish.');
uiwait(msgbox(message_Lesion));
hFH_Lesion = imfreehand(); % Actual line of code to do the drawing.
% Create a binary image ("mask") from the ROI object.
binaryImage_Lesion = hFH_Lesion.createMask();
xy_Lesion = hFH_Lesion.getPosition;

% Calculate the area, in pixels, that they drew.
numberOfPixels1_Lesion = sum(binaryImage_Lesion(:));
% Another way to calculate it that takes fractional pixels into account.
numberOfPixels2_Lesion = bwarea(binaryImage_Lesion);

% Get coordinates of the boundary of the freehand drawn region.
structBoundaries_Lesion = bwboundaries(binaryImage_Lesion);
xy_Lesion = structBoundaries_Lesion{1}; % Get n by 2 array of x,y coordinates.
x_Lesion = xy_Lesion(:, 2); % Columns.
y_Lesion = xy_Lesion(:, 1); % Rows.

%%%%% DRAW CONTROL %%%%%
% Ask user to draw freehand mask.
message_Control = sprintf('Left click and hold to draw for CONTROL.\nSimply lift the mouse button to finish.');
uiwait(msgbox(message_Control));
hFH_Control = imfreehand(); % Actual line of code to do the drawing.
% Create a binary image ("mask") from the ROI object.
binaryImage_Control = hFH_Control.createMask();
xy_Control = hFH_Control.getPosition;

% Calculate the area, in pixels, that they drew.
numberOfPixels1_Control = sum(binaryImage_Control(:));
% Another way to calculate it that takes fractional pixels into account.
numberOfPixels2_Control = bwarea(binaryImage_Control);

% Get coordinates of the boundary of the freehand drawn region.
structBoundaries_Control = bwboundaries(binaryImage_Control);
xy_Control = structBoundaries_Control{1}; % Get n by 2 array of x,y coordinates.
x_Control = xy_Control(:, 2); % Columns.
y_Control = xy_Control(:, 1); % Rows.

close % Closes drawing window.

% Calculate Average Temperature Values for Each Image
for loopfolder = 1:curvelength

disp(['Analyzing frame ' num2str(loopfolder) ' ...']);
% load individual .dat file
img_raw = load([path '/' ext num2str(loopfolder) '.dat']);
temp_data = (img_raw*100)*160/16384; %%% Revert 0.01 by x100, 160K Range (248K - 408K) on 14-Bit Intensity %%%

% Mask the images outside the mask, and display it.
% Will keep only the part of the image that's inside the mask, zero outside mask.
blackMaskedImage_Lesion = temp_data;
blackMaskedImage_Lesion(~binaryImage_Lesion) = 0;
blackMaskedImage_Control = temp_data;
blackMaskedImage_Control(~binaryImage_Control) = 0;

% Calculate the temperature means
dehydration_Lesion(loopfolder,1) = mean(blackMaskedImage_Lesion(binaryImage_Lesion));
dehydration_Control(loopfolder,1) = mean(blackMaskedImage_Control(binaryImage_Control));
end

dehydration_Lesion_min = dehydration_Lesion - min([min(dehydration_Lesion),min(dehydration_Control)]);
dehydration_Control_min = dehydration_Control - min([min(dehydration_Lesion),min(dehydration_Control)]);

%%%%% Calculate Delta Q %%%%%
dQ_L = round(max(dehydration_Lesion)*curvelength-sum(dehydration_Lesion),0);
dQ_C = round(max(dehydration_Control)*curvelength-sum(dehydration_Control),0);
DdQ = dQ_L-dQ_C;
dQ_Lesion(1,1) = num2str(dQ_L);
dQ_Control(1,1) = num2str(dQ_C);
Diff_dQ(1,1) = num2str(DdQ);

%%%%% Generating Figure %%%%%
% Plot 1st Image
subplot(2, 3, 1);
imshow(first_img, []);
axis on;
drawnow;
xlabel(['Pixel (' num2str(size(first_img,2)) ')']); 
ylabel(['Pixel (' num2str(size(first_img,1)) ')']);
fontSize = 16;
title([strrep(ext,'_',' ') 'Frame ' num2str(1)]);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

% Plot Last Image
img_raw_Last = load([path '/' ext num2str(curvelength) '.dat']);
img_org_Last = (img_raw_Last-absmin)*100;
Last_img = uint8(img_org_Last);

subplot(2, 3, 2);
imshow(Last_img, []);
axis on;
drawnow;
xlabel(['Pixel (' num2str(size(first_img,2)) ')']); 
ylabel(['Pixel (' num2str(size(first_img,1)) ')']);
title([strrep(ext,'_',' ') 'Frame ' num2str(curvelength)]);

% Plot Dehydration Map
map = subplot(2, 3, 3);
imshow(temp_Plot, []);
colorbar;
colormap(map,jet(65536));
axis on;
drawnow;
xlabel(['Pixel (' num2str(size(first_img,2)) ')']); 
ylabel(['Pixel (' num2str(size(first_img,1)) ')']);
title([strrep(ext,'_',' ') 'Dehydration Map']);

% Plot Lesion and Control ROIs
subplot(2, 3, 4);
imshow(Last_img, []);
axis on;
drawnow;
xlabel(['Pixel (' num2str(size(first_img,2)) ')']); 
ylabel(['Pixel (' num2str(size(first_img,1)) ')']);
title([strrep(ext,'_',' ') 'Frame ' num2str(curvelength) ' ROIs']);

subplot(2, 3, 4); % Plot Lesion over original image.
hold on; 
plot(x_Lesion, y_Lesion, 'r-', 'LineWidth', 2);
drawnow; % Force it to draw immediately.

subplot(2, 3, 4); % Plot Control over original image.
hold on; 
plot(x_Control, y_Control, 'b-', 'LineWidth', 2);
drawnow; % Force it to draw immediately.

% Plot Dehydration Curves
xaxis = 0:curvelength-1;
subplot(2, 3, 5); % Plot Control over original image.
hold on; 
title([strrep(ext,'_',' ') 'Dehydration Curve']);
plot(xaxis, dehydration_Lesion_min, 'r-', 'LineWidth', 2);
hold on;
plot(xaxis, dehydration_Control_min, 'b-', 'Linewidth', 2);
legend({'Lesion','Control'},'Location','southwest')
xlabel('Frames'); 
ylabel('Change in Temperature (K)');
%%% Text data onto Plot
str = {['dQ Lesion: ' num2str(str2double(dQ_Lesion(1,1)))],['dQ Control: ' num2str(str2double(dQ_Control(1,1)))],['Diff. dQ: ' num2str(str2double(Diff_dQ(1,1)))]};
xtext = [curvelength/2, curvelength/2, curvelength/2+100];
yL = dehydration_Lesion_min(round(curvelength/2,0),1);
yC = dehydration_Control_min(round(curvelength/2,0),1);
yDiff = (max([max(dehydration_Control_min),max(dehydration_Lesion_min)])-min([min(dehydration_Lesion_min),min(dehydration_Control_min)]))/2;
ytext = [yL,yC,yDiff];
text(xtext,ytext,str);

% Plot Dehy Map + ROIs
map = subplot(2, 3, 6);
imshow(temp_Plot, []);
colorbar;
colormap(map,jet(65536));
axis on;
drawnow;
xlabel(['Pixel (' num2str(size(first_img,2)) ')']); 
ylabel(['Pixel (' num2str(size(first_img,1)) ')']);
title([strrep(ext,'_',' ') 'Dehydration Map + ROIs']);

subplot(2, 3, 6); % Plot Lesion over original image.
hold on; 
plot(x_Lesion, y_Lesion, 'r-', 'LineWidth', 2);
drawnow; % Force it to draw immediately.

subplot(2, 3, 6); % Plot Control over original image.
hold on; 
plot(x_Control, y_Control, 'b-', 'LineWidth', 2);
drawnow; % Force it to draw immediately.

%%%%% Output Data and Save Data & Figure %%%%%
% Write data to be exported to Excel
data = horzcat(dehydration_Lesion,dehydration_Control);
Lesion = dehydration_Lesion;
Control = dehydration_Control;
finaldata = table(Lesion,Control,dQ_Lesion,dQ_Control,Diff_dQ);
% savepath = uigetdir;
savepath = '/Users/nai-yuannicholaschang/Desktop/Aim2 Dehy Results';
filename = [savepath '/' ext 'Results.xlsx'];
writetable(finaldata,filename);

% Save Figure
figname = [savepath '/' ext 'Results.png'];
saveas(gcf,figname);

disp(['****** Analysis of Sample ' num2str(ext(1:end-1)) ' Complete ******']);
disp(' ');
disp(['Dehydration Value (Lesion): ' num2str(str2double(dQ_Lesion(1,1))) ' a.u.']);
disp(' ');
disp(['Dehydration Value (Control): ' num2str(str2double(dQ_Control(1,1))) ' a.u.']);
disp(' ');
disp(['Dehydration Value Difference: ' num2str(str2double(Diff_dQ(1,1))) ' a.u.']);
disp(' ');
disp('************************');
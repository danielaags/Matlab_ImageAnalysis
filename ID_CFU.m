clear

%Info
day = '190921';
plate = '1';
%pixel_size=1/DistancePix;
%1inch x 96 pixels; 1inch = 2.54cm
pixel_size=2.54/96; %cm

%Read file
file = 'agar1mMp';
%16-bit file, RGB range 0-2500
I = imread(file, 'tif');
%8-bit, RGB range 0-255
I = uint8(I/257);

%%
%%Remove uninterested region from the image
%get the image size to remove the rim of the petri dish
imageSize = size(I);
%center and radius of circle ([c_row, c_col, r]). The images are a bit
%off. Center doesn't seem to be at 1024, 1024
ci = [1044, 1024, 850];     
[xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
mask = uint8((xx.^2 + yy.^2)<ci(3)^2);
%mask = uint16((xx.^2 + yy.^2)<ci(3)^2);
croppedImage = uint8(ones(size(I)));
%croppedImage = uint16(ones(size(I)));
croppedImage(:,:,1) = I(:,:,1).*mask;
croppedImage(:,:,2) = I(:,:,2).*mask;
croppedImage(:,:,3) = I(:,:,3).*mask;

%imshow(croppedImage)

%Find circles without segmentation
%[centers, radii] = imfindcircles(croppedImage,[20, 60], 'Sensitivity',0.92, 'EdgeThreshold', 0.04);


%%
%%Segmentation
grayImage = rgb2gray(croppedImage);
%Trasformation to get into 0-255 RGB range and plot a hist to see the
%distribution
figure
imhist(grayImage)
%imhist(uint8(grayImage/257))
print(strcat(file,'hist'),'-dpng');
close;

%For brighter colonies than the background
grayImage_1 = grayImage > 30& grayImage < 135; %imshow(grayImage_1)
%grayImage_1 = grayImage > 30*257 & grayImage < 135*257; %imshow(grayImage_1)
figure
imshow(grayImage_1)
print(strcat(file,'seg'),'-dpng');
close;
%Two iterations for better detection
[centers_1, radii_1] = imfindcircles(grayImage_1,[20 60], 'ObjectPolarity', 'dark');
[centers_2, radii_2] = imfindcircles(grayImage_1,[70 110], 'ObjectPolarity', 'dark');
%[centers_2, radii_2] = imfindcircles(grayImage_1,[90 150], 'Sensitivity', 0.98)

% More than one iteration of imfindcircles
centers = [centers_1 ; centers_2]; 
radii = [radii_1 ; radii_2];

%%
%%Plot the colonies found
%Labels
colony = cell(1,length(radii)); 
label = cell(1,length(radii));
for i = 1:length(radii)
    colony{i} = int2str(i);
    label{i}= strcat(day,'-',plate,'-',int2str(i));
end

%Plot the colony number
figure
imshow(croppedImage);
viscircles(centers,radii,'EdgeColor','b');
text(centers(:,1), centers(:,2), colony)
print(strcat(file,'IDs'),'-dpng');
close;

%%
%Extract information. 
%test
c1 = I(round(centers(1,2)-(radii(1)+50):centers(1,2)+(radii(1)+50)), round(centers(1,1)-(radii(1)+50):centers(1,1)+(radii(1)+50)));

x = floor(centers(1,1)-radii(1):centers(1,1)+radii(1));
y = repmat(floor(centers(1,2)), 1, length(x));

figure
subplot(1,2,1)
imshow(c1)

subplot(1,2,2)
improfile(I,x,y); grid on 
close;

c9 = I(round(centers(9,2)-(radii(9)+50):centers(9,2)+(radii(9)+50)), round(centers(9,1)-(radii(9)+50):centers(9,1)+(radii(9)+50)));

x = floor(centers(9,1)-radii(9):centers(9,1)+radii(9));
y = repmat(floor(centers(9,2)), 1, length(x));

figure
subplot(1,2,1)
imshow(c9)

subplot(1,2,2)
improfile(I,x,y); grid on 
close;


%%
%Save data
data = table(label', centers(:,1), centers(:,2), radii, radii*pixel_size, 'VariableNames', {'Label', 'x', 'y', 'r_px', 'r_cm'});
writetable(data,strcat(file,'.csv'),'Delimiter',',');

%%
%Garbage%

%For colonies similar to the background
% grayImage_2 = grayImage > 75; %imshow(grayImage_2)
% [centers_2, radii_2] = imfindcircles(grayImage_2,[20, 60], 'ObjectPolarity', 'bright');

% c = figure;
% subplot(1,2,1)
% imshow(croppedImage);
% 
% subplot(1,2,2)
% imshow(gray_image)
% viscircles(centers,radii,'EdgeColor','b');

% Raw seems too noisy
% red = I(round(centers(1,1)-radii(1):centers(1,1)+radii(1)),round(centers(1,2)),1);
% green = I(round(centers(1,1)-radii(1):centers(1,1)+radii(1)),round(centers(1,2)),2);
% blue = I(round(centers(1,1)-radii(1):centers(1,1)+radii(1)),round(centers(1,2)),3);
% 
% figure
% plot(red, 'red');
% hold on
% plot(green, 'green');
% plot(blue, 'blue');
% hold off
% 
% define a disk filter for image convolution. It integrates the fluoresce 
% over a radius of 80%
% diskfilter = fspecial('disk',0.8*20); 
% 
% convolve the primary image with the diskfilter. This creates a new images
% that locally integrates fluorescence.
% Iconvolved = imfilter(I,diskfilter,'replicate');
% Extract data
% red = Iconvolved(round(centers(1,1)-radii(1):centers(1,1)+radii(1)),round(centers(1,2)),1);
% green = Iconvolved(round(centers(1,1)-radii(1):centers(1,1)+radii(1)),round(centers(1,2)),2);
% blue = Iconvolved(round(centers(1,1)-radii(1):centers(1,1)+radii(1)),round(centers(1,2)),3);
% 
% figure
% plot(red, 'red');
% hold on
% plot(green, 'green');
% plot(blue, 'blue');
% hold off

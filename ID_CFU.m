%Read file
file = 'day5-40µl-nr_5';
I = imread(strcat(file,'.tif'));

%%%Remove uninterested region from the image%%%
%get the image size to remove the rim of the petri dish
imageSize = size(I);
%center and radius of circle ([c_row, c_col, r]). The images are a bit
%off. Center doesn't seem to be at 1024, 1024
ci = [1044, 1024, 850];     
[xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
mask = uint16((xx.^2 + yy.^2)<ci(3)^2);
croppedImage = uint16(ones(size(I)));
croppedImage(:,:,1) = I(:,:,1).*mask;
croppedImage(:,:,2) = I(:,:,2).*mask;
croppedImage(:,:,3) = I(:,:,3).*mask;

imshow(croppedImage)

%Find circles without segmentation
%[centers, radii] = imfindcircles(croppedImage,[20, 60], 'Sensitivity',0.92, 'EdgeThreshold', 0.04);

%Segmentation%%%
grayImage = rgb2gray(croppedImage);
%bg = gaussf(grayImage,1,'best');
imshow(grayImage);
%imshow(bg);
%Trasformation to get into 0-255 RGB range
imhist(uint8(gray_image/257))

%For brighter colonies than the background
grayImage_fd = grayImage > 30*257 & grayImage < 120*257; %imshow(grayImage_f)
[centers_fd, radii_fd] = imfindcircles(grayImage_f,[20, 60], 'ObjectPolarity', 'dark');

%For colonies similar to the background
grayImage_fb = grayImage > 75; %imshow(grayImage_f)
[centers_fb, radii_fb] = imfindcircles(grayImage_f,[20, 60], 'ObjectPolarity', 'bright');

centers = [centers_fd ; centers_fb];
data = table(centers);
writetable(data,strcat(file,'_centers.txt'),'Delimiter','\t');

radii = [radii_fd ; radii_fb];
data = table(radii);
writetable(data,strcat(file,'_radii.txt'),'Delimiter','\t');


% %Plot the colonies found
% c = figure;
% subplot(1,2,1)
% imshow(croppedImage);
% 
% subplot(1,2,2)
% imshow(gray_image)
% viscircles(centers,radii,'EdgeColor','b');

%Plot the colony number
imshow(croppedImage);
viscircles(centers,radii,'EdgeColor','b');

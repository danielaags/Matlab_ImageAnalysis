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

i=1;
%Red channel
rgbI = imadjust(croppedImage(:,:,i)); %imshow(RGBImage);

%Correct non-uniform ilumination
se = strel('disk',175);
background = imopen(rgbI,se);
%imshow(background)
rgbIbackground = rgbI - background; %imshow(rgbIbackground);

%For brighter colonies than the background
%imhist(RGBImage_)
rgbI_Filter = rgbIbackground > 50;%imshow(rgbI_Filter);
%to fill up holes
rgbI_fill = imfill(rgbI_Filter,'holes');%imshow(rgbI_fill)
%remove connected objest
rgbI_nobord = imclearborder(rgbI_fill,4);%imshow(rgbI_nobord)
%smooth object
seD = strel('diamond',1);
rgbI_final = imerode(rgbI_nobord,seD); rgbI_final = imerode(rgbI_final,seD);imshow(rgbI_final)

%Find colonies using boundary
[B,L,n,A] = bwboundaries(rgbI_final,'noholes');
imshow(L)
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
end

%get stats
stats=  regionprops(L, 'Centroid', 'Area', 'Perimeter');
Centroid = cat(1, stats.Centroid);
Perimeter = cat(1,stats.Perimeter);
Area = cat(1,stats.Area);
CircleMetric = (Perimeter.^2)./(4*pi*Area);  %circularity metric
Radii = sqrt(Area/pi);

% %Threshold for metric
% isCircle =   (CircleMetric < 1.1);
% %assign shape to each object
% whichShape = cell(n,1);  
% whichShape(isCircle) = {'Circle'};

imshow(rgbI_final);
viscircles(Centroid,Radii,'EdgeColor','r');


% now label with results
% RGB = label2rgb(L);
% imshow(RGB); hold on;
% for k=1:n
%    display metric values and which shape next to object
%    text( Centroid(k,1), Centroid(k,2), 'a');
%    text( Centroid(k,1)-20, Centroid(k,2)+20, whichShape{k});
% end

%%
%%Plot the colonies found
%Labels
colony = cell(1,length(radii)); 
label = cell(1,length(radii));
for i = 1:length(radii)
    colony{i} = int2str(i);
    label{i}= strcat(day,'-',plate,'-',int2str(i));
end

%Plot 
figure
subplot(2,2,1)
imshow(rgbI_Filter)

subplot(2,2,2)
imshow(rgbI_fill)

subplot(2,2,3)
imshow(rgbI_nobord)

subplot(2,2,4)
imshow(rgbI_final)
print(strcat(file,'flow'),'-dpng');
close;

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

%Other ways to highlight colonies
% imshow(labeloverlay(grayImage,grayImage_1))
% 
% BWoutline = bwperim(grayImage_1);
% Segout = grayImage; 
% Segout(BWoutline) = 100; 
% imshow(Segout)


% figure
imhist(grayImage)
% %imhist(uint8(grayImage/257))
% print(strcat(file,'hist'),'-dpng');
% close;

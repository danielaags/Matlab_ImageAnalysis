function[] = ID_CFU_rgb(day0, plateN, filename)
%This script identifies bacteria colonies in a petri dish. It generates a
%mask first to avoid the outer part of the plate. Later, segmentation and
%label identification allows to target even those colonies that are not
%completely circular. It retrieves a .jpg file with the identifies colonies
%and a .csv file with a label and basic information.

%clear

%Info
day = day0;
plate = plateN;
%pixel_size=1/DistancePix;
%1inch x 96 pixels; 1inch = 2.54cm
pixel_size=2.54/96; %cm

%Read file
file = filename;
%16-bit file, RGB range 0-2500
I = imread(file, 'tif');
%8-bit, RGB range 0-255
I = uint8(I/257);

%%
%%Remove uninterested region from the image
%get the image size to remove the rim of the petri dish
imageSize = size(I);
%center and radius of circle ([c_col, c_row, r]). The images are a bit
%off. Center doesn't seem to be at 1024, 1024
ci = [1044, 1024, 750];     
[xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
mask = uint8((xx.^2 + yy.^2)<ci(3)^2);
%mask = uint16((xx.^2 + yy.^2)<ci(3)^2);
croppedImage = uint8(ones(size(I)));
%croppedImage = uint16(ones(size(I)));
croppedImage(:,:,1) = I(:,:,1).*mask;
croppedImage(:,:,2) = I(:,:,2).*mask;
croppedImage(:,:,3) = I(:,:,3).*mask;


%%
%%Segmentation
c1 = cell(1,3);
c2 = cell(1,3);
r = cell(1,3);

for i = 1:3
    %Different channels
    %i=1;
    rgbI = imadjust(croppedImage(:,:,i)); %imshow(rgbI);
    
    %Correct non-uniform ilumination. Disk size should be smaller than colonies
    %size.
    %%#important parameter to take into consideration#%%
    se = strel('disk',200);
    background = imopen(rgbI,se);
    rgbIbackground = rgbI - background; %imshow(rgbIbackground);
    
    %For brighter colonies than the background
    %imhist(RGBImage_)
    %%#important parameter to take into consideration#%%
    rgbI_filter = rgbIbackground > 65;%imshow(rgbI_Filter);
    %to fill up holes
    rgbI_fill = imfill(rgbI_filter,'holes');%imshow(rgbI_fill)
    %remove connected objects
    rgbI_nobord = imclearborder(rgbI_fill,4);%imshow(rgbI_nobord)
    %smooth object
    seD = strel('diamond',1);
    rgbI_final = imerode(rgbI_nobord,seD); rgbI_final = imerode(rgbI_final,seD);
    
    %Find colonies using boundary
    [B,L,n,A] = bwboundaries(rgbI_final,'noholes');
%     imshow(L)
%     hold on
%     for k = 1:length(B)
%         boundary = B{k};
%         plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
%     end
    
    %get stats
    %stats=  regionprops(L, 'Centroid', 'Area', 'Perimeter');
    stats=  regionprops(L, 'Centroid', 'Area');
    Centroid = cat(1, stats.Centroid);
    %Perimeter = cat(1,stats.Perimeter);
    Area = cat(1,stats.Area);
    %radii of shapes
    Radii = sqrt(Area/pi);
    
    %Filter by Radii bigger than n
    filter = Radii > 6 & Radii < 100;
    % imshow(rgbI_final);
    % viscircles(Centroid(filter,:),Radii(filter),'EdgeColor','r');
    
%    c1{i} = Centroid(filter,1);
%    c2{i} = Centroid(filter,2);
   c1{i} = floor(Centroid(filter,1));
   c2{i} = floor(Centroid(filter,2));
   r{i} = Radii(filter);
   %p{i} = Perimeter(filter);
end

%Save last boundaries, blue channel
figure
imshow(rgbI_final)
print(strcat(file,'bw_boundaries'),'-dpng');
close;

%Further filter will stablish that only centroids found more than one will
%be taken into account in the final output.

%Concatenate the values from cell to array
center1 = cat(1,c1{1},c1{2},c1{3});
center2 = cat(1,c2{1},c2{2},c2{3});
radii = cat(1,r{1},r{2},r{3});

%Find unique elements in x
tol = 2/max(abs(center1));
C = uniquetol(center1, tol);
idx_C = zeros(length(C),1);

for n = 1:length(C)
    idx_C(n) = find(center1 == C(n),1);
end

%Find unique elements in y
tol = 2/max(abs(center2(idx_C)));
D = uniquetol(center2(idx_C), tol);
idx_D = zeros(length(D),1);

for n = 1:length(D)
    idx_D(n) = find(center2 == D(n),1);
end

%Filter out ID colonies close/outside the plate
%Get indexes from mask
x = center1(idx_D);
y = center2(idx_D);
r = radii(idx_D);

% centers = [x, y];
% radii = r;

%Find the elements with nonzero elements. The pixels 
[X_mask, Y_mask] = find(mask);
idx_E = zeros(length(x),1);

%calculate indexes to plot a circle and compare with the indexes in mask
th = 0:pi/5:2*pi;
for i = 1:length(x)
%     imshow(croppedImage);
    %Fit a circle
    xunit = floor(r(i) * cos(th) + x(i));
    yunit = floor(r(i) * sin(th) + y(i));
    %Find within the boundaries. Check ci variable
    mean_x = mean(xunit>400 & xunit<1700);
    mean_y = mean(yunit>400 & yunit<1700);
%     hold on
%     h = plot(xunit, yunit);
%     hold off
    if mean_x == 1 && mean_y == 1
       idx_E(i) = 1;
    end
end

idx = find(idx_E);
centers = [x(idx), y(idx)];
radii = r(idx);

%%
%%Plot the colonies found
%Labels
colony = cell(1,length(radii)); 
label = cell(1,length(radii));
for i = 1:length(radii)
    colony{i} = int2str(i);
    label{i}= strcat(day,'-',plate,'-',int2str(i));
end

figure
imshow(croppedImage);
viscircles(centers,radii,'EdgeColor','b');
text(centers(:,1), centers(:,2), colony)
print(strcat(file,'IDs'),'-dpng');
close;


%%
%Save data
data = table(label', centers(:,1), centers(:,2), round(radii), fprintf(radii*pixel_size), 'VariableNames', {'Label', 'x', 'y', 'r_px', 'r_cm'});
writetable(data,strcat(file,'.csv'),'Delimiter',',');


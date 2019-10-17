function[] = ID_CFU_bright_dark(day0, plateN, filename)
%This script identifies bacteria colonies in a petri dish. It generates a
%mask first to avoid the outer part of the plate. Later, segmentation and
%label identification allows to target even those colonies that are not
%completely circular. It retrieves a .jpg file with the identifies colonies
%and a .csv file with a label and basic information.

% clear
 
% %to test without function
% day0 = '190921';
% plateN = '1';
% filename = 'd20-20µl-nr5';

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
%%Segmentation just in grayscale
c1 = cell(2,1);
c2 = cell(2,1);
r = cell(2,1);

rgbI = imadjust(rgb2gray(croppedImage)); %imshow(rgbI);

%Correct non-uniform ilumination. Disk size should be smaller than colonies
%size.#important parameter to take into consideration#
se = strel('disk',100);
background = imopen(rgbI,se);

for i = 1:2
    if i == 1
        rgbIbackground = rgbI - background; %imshow(rgbIbackground);
        %For brighter colonies than the background
        %%#important parameter to take into consideration#%%
        %Filter out most of the plate
        rgbI_filter = rgbIbackground > 45; %imshow(rgbI_filter);
        %to fill up holes
        rgbI_fill = imfill(rgbI_filter,'holes');%imshow(rgbI_fill)
        %remove connected objects
        rgbI_nobord = imclearborder(rgbI_fill,8);%imshow(rgbI_nobord)
        %smooth object
        seD = strel('diamond',1);
        rgbI_final = imerode(rgbI_nobord,seD); %rgbI_final = imerode(rgbI_final,seD);

        %Find colonies using boundary
        [B,L,n,A] = bwboundaries(rgbI_final,'noholes');
        %figure
        %imshow(L)
        %hold on
        %for k = 1:length(B)
        %    boundary = B{k};
        %    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
        %end

        %get stats
        %stats=  regionprops(L, 'Centroid', 'Area', 'Perimeter');
        %Perimeter = cat(1,stats.Perimeter);
        stats=  regionprops(L, 'Centroid', 'Area');
        Centroid = cat(1, stats.Centroid);
        Area = cat(1,stats.Area);
        %radii of shapes
        Radii = sqrt(Area/pi);

        %Filter by Radii bigger than n
        filter = Radii > 10 & Radii < 150;

        if isempty(filter) == 0
            c1{i} = floor(Centroid(filter,1));
            c2{i} = floor(Centroid(filter,2));
            r{i} = Radii(filter);
        end
    end
    if i ==2
        %For darker colonies than the background
        rgbIbackground = rgbI + background; %imshow(rgbIbackground);
        %%#important parameter to take into consideration#%%
        %Filter out most of the plate
        rgbI_filter = rgbIbackground < 75; %imshow(rgbI_filter);
        rgbI_filter = rgbI_filter.*double(mask);
        %to fill up holes
        rgbI_fill = imfill(rgbI_filter,'holes');%imshow(rgbI_fill)
        %remove connected objects
        rgbI_nobord = imclearborder(rgbI_fill,8);%imshow(rgbI_nobord)
        %smooth object
        seD = strel('diamond',1);
        rgbI_final = imerode(rgbI_nobord,seD); %rgbI_final = imerode(rgbI_final,seD);

        %Find colonies using boundary
        [B,L,n,A] = bwboundaries(rgbI_final,'noholes');
%         figure 
%         imshow(L)
%         hold on
%         for k = 1:length(B)
%             boundary = B{k};
%             plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
%         end
        
        %get stats
        %stats=  regionprops(L, 'Centroid', 'Area', 'Perimeter');
        %Perimeter = cat(1,stats.Perimeter);
        stats=  regionprops(L, 'Centroid', 'Area');
        Centroid = cat(1, stats.Centroid);
        Area = cat(1,stats.Area);
        %radii of shapes
        Radii = sqrt(Area/pi);

        %Filter by Radii bigger than n
        filter = Radii > 10 & Radii < 150;

        if isempty(filter) == 0
            c1{i} = floor(Centroid(filter,1));
            c2{i} = floor(Centroid(filter,2));
            r{i} = Radii(filter);
        end

    end    

end

x = cat(1,c1{1},c1{2});
y = cat(1,c2{1},c2{2});
r = cat(1,r{1},r{2});

%Find the elements with nonzero elements.
idx = zeros(length(x),1);

%Calculate indexes to plot a circle and compare with the indexes in mask
th = 0:pi/5:2*pi;
for i = 1:length(x)
%     imshow(croppedImage);
    %Fit a circle
    xunit = floor(r(i) * cos(th) + x(i));
    yunit = floor(r(i) * sin(th) + y(i));
    %Find within the boundaries. Check ci variable
    mean_x = mean(xunit>300 & xunit<1700);
    mean_y = mean(yunit>300 & yunit<1700);
%     hold on
%     h = plot(xunit, yunit);
%     hold off
    if mean_x == 1 && mean_y == 1
       idx(i) = 1;
    end
end

idx_in = find(idx);
centers = [x(idx_in), y(idx_in)];
radii = r(idx_in);

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
text(centers(:,1), centers(:,2), colony);
print(strcat(file,'IDs'),'-dpng');
close;


%Save data
data = table(label', centers(:,1), centers(:,2), round(radii), round(radii*pixel_size,2), 'VariableNames', {'Label', 'x', 'y', 'r_px', 'r_cm'});
writetable(data,strcat(file,'.csv'),'Delimiter',',');
 
end


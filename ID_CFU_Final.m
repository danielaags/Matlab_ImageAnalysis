function[] = ID_CFU_Final(day0, plateN, filename)
%This script identifies bacteria colonies in a petri dish. It generates a
%mask first to avoid the outer part of the plate. Later, segmentation and
%label identification allows to target even those colonies that are not
%completely circular. It retrieves a .jpg file with the identifies colonies
%and a .csv file with a label and basic information.

% clear
 
% %to test without function
% day0 = '190921';
% plateN = '1';
% filename = '_d30-10µl-nr3';

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
%[y0,x0]
ci = [1040, 1015, 850];     
[xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
mask = uint8((xx.^2 + yy.^2)<ci(3)^2);
croppedImage = uint8(true(size(I)));
croppedImage(:,:,1) = I(:,:,1).*mask;
croppedImage(:,:,2) = I(:,:,2).*mask;
croppedImage(:,:,3) = I(:,:,3).*mask;
% figure
% imshow(croppedImage)


%%
%%Segmentation just in grayscale
c1 = cell(2,1);
c2 = cell(2,1);
r = cell(2,1);
e = cell(2,1);


%Correct non-uniform ilumination. Adjust specific range
rgbI = imadjust(croppedImage(:,:,1), [0 0.60],[]);%imshow(rgbI);

for i = 1:2
    if i == 1
        %Filter bright colonies
        rgbBW = rgbI >=190;%imshow(rgbBW)
        %remove connected objects
        rgbI_nobord = imclearborder(rgbBW,8);%imshow(rgbI_nobord)
        %to fill up holes
        rgbI_final = imfill(rgbI_nobord,'holes');%imshow(rgbI_final)       
        %smooth object
%         seD = strel('diamond',1);
%         rgbI_final = imerode(rgbI_nobord,seD); %rgbI_final = imerode(rgbI_final,seD);

        %Find colonies using boundary
        [B,L1,n,A] = bwboundaries(rgbI_final,'noholes');      
%         figure
%         imshow(L)
%         hold on
%         for k = 1:length(B)
%            boundary = B{k};
%            plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
%         end

        %get stats
        stats=  regionprops(L1, 'Centroid', 'Area', 'Perimeter', 'Eccentricity');
        %stats=  regionprops(L, 'Centroid', 'Area');
        Eccentricity = cat(1,stats.Eccentricity);
        Perimeter = cat(1,stats.Perimeter);
        Centroid = cat(1, stats.Centroid);
        Area = cat(1,stats.Area);
        Radii = sqrt(Area/pi);

        %Filter by Area bigger than n
        filter = Area > 200 & Area < 70000 & Eccentricity < 0.7; 

        if isempty(filter) == 0
            c1{i} = floor(Centroid(filter,1));
            c2{i} = floor(Centroid(filter,2));
            r{i} = Radii(filter);
            e{i} = Eccentricity(filter);
        end
    else
        rgbBW = rgbI < 50;%imshow(rgbBW)
        %remove connected objects
        rgbI_nobord = imclearborder(rgbBW,8);%imshow(rgbI_nobord)
        %smooth object
        seD = strel('diamond',1);
        rgbI_final = imerode(rgbI_nobord,seD); %rgbI_final = imerode(rgbI_final,seD);

        %Find colonies using boundary
        [B,L2,n,A] = bwboundaries(rgbI_final,'noholes');
%         figure 
%         imshow(L2)
%         hold on
%         for k = 1:length(B)
%             boundary = B{k};
%             plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
%         end
        
        %get stats
        stats=  regionprops(L2, 'Centroid', 'Area', 'Perimeter', 'Eccentricity');
        %stats=  regionprops(L, 'Centroid', 'Area');
        Eccentricity = cat(1,stats.Eccentricity);
        Perimeter = cat(1,stats.Perimeter);
        Centroid = cat(1, stats.Centroid);
        Area = cat(1,stats.Area);
        Radii = sqrt(Area/pi);

        %Filter by Area bigger than n
        filter = Area > 200 & Area < 70000 & Eccentricity < 0.7; 
        if isempty(filter) == 0
            c1{i} = floor(Centroid(filter,1));
            c2{i} = floor(Centroid(filter,2));
            r{i} = Radii(filter);
            e{i} = Eccentricity(filter);

        end

    end 
    
end

 %BW image
    L = (L1 | L2);  
    figure
    imshow(L)
    print(strcat('BW-', file),'-dpng');
    close;

x = cat(1,c1{1},c1{2});
y = cat(1,c2{1},c2{2});
r = cat(1,r{1},r{2});
e = cat(1,e{1},e{2});

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
    mean_xy = mean(yunit>225 & yunit<1875 & xunit>225 & xunit<1875);
    if mean_xy == 1
       idx(i) = 1;
    end
end

idx_in = find(idx);
centers = [x(idx_in), y(idx_in)];
radii = r(idx_in);
eccentricity = e(idx_in);


%%
%%Plot the colonies found
%Labels
colony = cell(1,length(radii)); 
label = cell(1,length(radii));
for i = 1:length(radii)
    colony{i} = int2str(i);
    label{i}= strcat(day,'-',plate,'-',int2str(i));
end
% 
% figure
% imshow(I);
% viscircles(centers,radii,'EdgeColor','b');
% text(centers(:,1), centers(:,2), colony);
% print(strcat(file,'IDs'),'-dpng');
% close;
% 
% 
% %Save data
% data = table(label', centers(:,1), centers(:,2), round(radii), round(radii*pixel_size,2), eccentricity, 'VariableNames', {'Label', 'x', 'y', 'r_px', 'r_cm', 'eccentricity'});
% writetable(data,strcat(file,'.csv'),'Delimiter',',');
 
end

%[a,b,c] = impixel()

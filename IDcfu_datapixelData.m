%function[] = ID_CFU_Final(day0, plateN, filename)
%This script identifies bacteria colonies in a petri dish. It generates a
%mask first to avoid the outer part of the plate. Later, segmentation and
%label identification allows to target even those colonies that are not
%completely circular. It retrieves a .jpg file with the identifies colonies
%and a .csv file with a label and basic information.

clear
 
%to test without function
% day0 = '190921';
% plateN = '1';
%filename = 'i0-d30-20µl-nr2';

%Info
day = '190921'; 
%Replace plateN it with n counter
%plate = plateN;
%pixel_size=1/DistancePix;
%1inch x 96 pixels; 1inch = 2.54cm
pixel_size=2.54/96; %cm

%Read files. Batch of pictures from a plate taken along several days
images = dir('*.tif');

for n = 1:length(images)
    file= images(n).name;
    filename = split(file, '.');
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
    ci = [1040, 1015, 845];     
    [xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
    mask = uint8((xx.^2 + yy.^2)<ci(3)^2);
    croppedImage = uint8(true(size(I)));
    croppedImage(:,:,1) = I(:,:,1).*mask;
    croppedImage(:,:,2) = I(:,:,2).*mask;
    croppedImage(:,:,3) = I(:,:,3).*mask;
%     figure
%     imshow(croppedImage)
    
    %Only take information from the red channel
    
    %Correct non-uniform ilumination. Adjust specific range
    rgbI = imadjust(croppedImage(:,:,1), [0 0.60],[]);%imshow(rgbI);
    
    %top-hat filtering 
%     se = strel('disk',100);
%     rgbI = imtophat(croppedImage(:,:,1), se);imshow(rgbI);
%     

    %Filter bright colonies
    rgbBW = rgbI >=190;%imshow(rgbBW)
    %remove connected objects
    rgbI_nobord = imclearborder(rgbBW,8);%imshow(rgbI_nobord)
    %to fill up holes
    rgbI_final = imfill(rgbI_nobord,'holes');%imshow(rgbI_final)       

    %Find colonies using boundary. Save the data collected.
    %B, returns an array of pixel locations
    %L, label matrix of objects
    %n, number of objects (labels)
    %A, adjacent matrix
    [B1,L1,n1,A1] = bwboundaries(rgbI_final,'noholes');      
        figure
        imshow(L)
        hold on
        for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
        end

    rgbBW = rgbI < 50;%imshow(rgbBW)
    %remove connected objects
    rgbI_nobord = imclearborder(rgbBW,8);%imshow(rgbI_nobord)
    %smooth object
    seD = strel('diamond',1);
    rgbI_final = imerode(rgbI_nobord,seD);

    %Find colonies using boundary
    [B2,L2,n2,A2] = bwboundaries(rgbI_final,'noholes');

    if isempty(n2)
        B = B1;
        L = L1;
        A = A1;
    else
        %BW image
        B = [B1; B2];
        L = (L1 | L2);
        A = [A1; A2];
    end 

    figure
    imshow(L)
    print(strcat(filename{1},'-BW'),'-dpng');
    close;

    %get morphological stats
    stats=  regionprops(L, 'Centroid', 'Area', 'EquivDiameter', 'Perimeter', 'Circularity', 'Eccentricity');
    Eccentricity = cat(1,stats.Eccentricity);
    Area = cat(1,stats.Area);
    Diameter = cat(1,stats.EquivDiameter);
    %Radii = sqrt(Area/pi); %calculated from area
    
    %%
    %Get pixel information per channel
    for k = 1:3
        rgbI = I(:,:,k);
        %Save according to RGB channel
        if k == 1
            %Better to have a different name for the different structures
            statsPixel_red = regionprops(L,rgbI,{'Centroid','PixelValues','MeanIntensity'});
            %save(strcat(filename{1},'-pixel_red.mat'), 'statsPixel_red');
        elseif k == 2
            %Better to have a different name for the different structures
            statsPixel_green = regionprops(L,rgbI,{'Centroid','PixelValues', 'MeanIntensity'});
            %save(strcat(filename{1},'-pixel_green.mat'), 'statsPixel_green');
        else
            %Better to have a different name for the different structures
            statsPixel_blue = regionprops(L,rgbI,{'Centroid','PixelValues', 'MeanIntensity'});
            %save(strcat(filename{1},'-pixel_blue.mat'), 'statsPixel_blue');
        end

    end

    %%
    %Filter by area bigger than n and eccentricity, the closes to 1 the more
    %line segment, the closes to 0 the more circular shape
    filter1 = find(Area > 200 & Area < 70000 & Eccentricity < 0.7); 

    %Find the elements close the the petri dish walls
    filter2 = zeros(length(filter1),1);

    %Calculate indexes to plot a circle and compare with the indexes in mask
    th = 0:pi/5:2*pi;
    for i = 1:length(filter1)
        %imshow(croppedImage);
        %Fit a circle
        xunit = floor(Diameter(filter1(i)) * cos(th) + Centroid(filter1(i),1));
        yunit = floor(Diameter(filter1(i)) * sin(th) + Centroid(filter1(i),2));
        %Find within the boundaries. Check ci variable
        mean_xy = mean(yunit>225 & yunit<1875 & xunit>225 & xunit<1875);
        if mean_xy == 1
           filter2(i) = filter1(i);
        end
    end

    filter2 = nonzeros(filter2);
    
    diameter = floor(Diameter(filter2));

    %%
    %%Plot the colonies found and save data
    %Labels
    colony = cell(1,length(diameter)); 
    label = cell(length(diameter),1);
    ID = cell(length(diameter),1);
    %Give labels
    for i = 1:length(diameter)
        colony{i} = int2str(i);
        label{i}= strcat(day,'-',int2str(n),'-',int2str(i));
        ID{i}=strcat(int2str(n),'-',int2str(i));
    end
    % 
    centroid= floor(cat(1,stats(filter2).Centroid));
    figure
    imshow(croppedImage);
    viscircles(centroid,diameter/2,'EdgeColor','b');
    text(centroid(:,1), centroid(:,2), colony);
    print(strcat(filename{1},'IDs'),'-dpng');
    close;
%%
%Check how to save the data, all the variables bwboundaries, regiongprops
%morphological and pixel values
%Save data
    statsData = struct('label', label, 'ID', ID, 'centroid', centroid, 'perimeter', perimeter, 'area', area, 'radii', radii, 'eccentricity', eccentricity, 'circularity', circularity);
    save(strcat(filename{1},'-data.mat'), 'statsData');


end 
%end

%[a,b,c] = impixel()

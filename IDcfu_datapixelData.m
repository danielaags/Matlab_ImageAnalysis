%function[] = ID_CFU_Final(day0, plateN, filename)
%This script identifies bacteria colonies in a petri dish. It generates a
%mask first to avoid the outer part of the plate. Later, segmentation and
%label identification allows to target even those colonies that are not
%completely circular. It retrieves a .jpg file with the identifies colonies
%and a .csv file with a label and basic information.

clear
 
%to test without function
day0 = '190921';
plateN = '1';
%filename = 'i0-d30-20µl-nr2';

%Info
day = day0; 
plate = plateN;
%pixel_size=1/DistancePix;
%1inch x 96 pixels; 1inch = 2.54cm
pixel_size=2.54/96; %cm

%Read file
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
            [B,L1,n1,A] = bwboundaries(rgbI_final,'noholes');      
    %         figure
    %         imshow(L)
    %         hold on
    %         for k = 1:length(B)
    %            boundary = B{k};
    %            plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
    %         end
        else
            rgbBW = rgbI < 50;%imshow(rgbBW)
            %remove connected objects
            rgbI_nobord = imclearborder(rgbBW,8);%imshow(rgbI_nobord)
            %smooth object
            seD = strel('diamond',1);
            rgbI_final = imerode(rgbI_nobord,seD); %rgbI_final = imerode(rgbI_final,seD);

            %Find colonies using boundary
            [B,L2,n2,A] = bwboundaries(rgbI_final,'noholes');
        end     
    end

    if isempty(n2)
        L = L1;
    else
        %BW image
        L = (L1 | L2);
    end 

    % figure
    % imshow(L)
    % imwrite(L, strcat('BW-', file),'-png');
    % close;

    %get stats
    stats=  regionprops(L, 'Centroid', 'Area', 'Perimeter', 'Circularity', 'Eccentricity');
    Eccentricity = cat(1,stats.Eccentricity);
    Circularity = cat(1,stats.Circularity);
    Centroid = round(cat(1, stats.Centroid)); %x,y
    Perimeter = cat(1,stats.Perimeter);
    Area = cat(1,stats.Area);
    Radii = sqrt(Area/pi); %calculated from area

    %Filter by area bigger than n and eccentricity, the closes to 1 the more
    %line segment, the closes to 0 the more circular shape
    filter1 = find(Area > 200 & Area < 70000 & Eccentricity < 0.7); 

    %Find the elements close the the petri dish walls
    filter2 = zeros(length(filter1),1);

    %Calculate indexes to plot a circle and compare with the indexes in mask
    th = 0:pi/5:2*pi;
    for i = 1:length(filter1)
    %     imshow(croppedImage);
        %Fit a circle
        xunit = floor(Radii(filter1(i)) * cos(th) + Centroid(filter1(i),1));
        yunit = floor(Radii(filter1(i)) * sin(th) + Centroid(filter1(i),2));
        %Find within the boundaries. Check ci variable
        mean_xy = mean(yunit>225 & yunit<1875 & xunit>225 & xunit<1875);
        if mean_xy == 1
           filter2(i) = filter1(i);
        end
    end

    filter2 = nonzeros(filter2);

    eccentricity = Eccentricity(filter2);
    circularity = Circularity(filter2);
    perimeter = Perimeter(filter2);
    area = Area(filter2);
    radii = Radii(filter2);

    %%
    %Compare two data points and get the change in radii
    if n>1
        %Save the i-1 and the i to compared
        centroidsO = centroid;
        centroid = Centroid(filter2,:);

        %Check if there are colonies in i-1 not seen in i
        isDay = ones(length(centroidsO),1);

        for i = 1:length(centroidsO)
            x = find(centroidsO(i,1)-3 <= centroid(:,1) & centroidsO(i,1)+3 >= centroid(:,1),1);
            if isempty(x)
                y = find(centroidsO(i,2)-3 <= centroid(:,2) & centroidsO(i,2)+3 >= centroid(:,2),1);
                if isempty(y)
                    %Not found in i
                    isDay(i) = false;
                else
                    %Found in i
                    isDay(i) = y;
                end
            else
                isDay(i) = x;
            end     
        end    
    else
        %It will save the first time point
        centroid = round(Centroid(filter2,:)); %x,y    
    end

    %%
    %%Plot the colonies found and save data
    %Labels
    colony = cell(1,length(radii)); 
    label = cell(length(radii),1);
    ID = cell(length(radii),1);
    %Give labels
    for i = 1:length(radii)
        colony{i} = int2str(i);
        label{i}= strcat(day,'-',plate,'-',int2str(i));
        ID{i}=strcat(plate,'-',int2str(i));
    end
    % 
    figure
    imshow(I);
    viscircles(centroid,radii,'EdgeColor','b');
    text(centroid(:,1), centroid(:,2), colony);
    print(strcat(filename{1},'IDs'),'-dpng');
    close;

    %Save data
    statsData = struct('label', label, 'ID', ID, 'centroid', centroid, 'perimeter', perimeter, 'area', area, 'radii', radii, 'eccentricity', eccentricity, 'circularity', circularity);
    save(strcat(filename{1},'-data.mat'), 'statsData');

    %%
    %Get pixel information per channel
    for k = 1:3
        rgbI = I(:,:,k);
        %Create BW canvas
        imageSize = size(rgbI);
        BW = false(imageSize);

        %The number of identified colonies don't stay the same. I will better
        %save the whole struct with centroids included to filter later.
        imshow(rgbI);
        hold on
        for i = 1:length(radii)
            h = drawcircle(gca,'Center',centroid(i,:),'Radius',(radii(i)));
            mask = h.createMask();
            BW = BW | mask;
        end
        hold off
        close

        % figure
        % imshow(BW)
        %Save according to RGB channel
%         if k == 1
%             %Better to have a different name for the different structures
%             statsPixel_red = regionprops(BW,rgbI,{'Centroid','PixelValues', 'MaxIntensity', 'MinIntensity', 'MeanIntensity'});
%             save(strcat(filename{1},'-pixel_red.mat'), 'statsPixel_red');
%         elseif k == 2
%             %Better to have a different name for the different structures
%             statsPixel_green = regionprops(BW,rgbI,{'Centroid','PixelValues', 'MaxIntensity', 'MinIntensity', 'MeanIntensity'});
%             save(strcat(filename{1},'-pixel_green.mat'), 'statsPixel_green');
%         else
%             %Better to have a different name for the different structures
%             statsPixel_blue = regionprops(BW,rgbI,{'Centroid','PixelValues', 'MaxIntensity', 'MinIntensity', 'MeanIntensity'});
%             save(strcat(filename{1},'-pixel_blue.mat'), 'statsPixel_blue');
%         end

    end

end 
%end

%[a,b,c] = impixel()

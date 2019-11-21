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
day = day0; 
plate = plateN;
%pixel_size=1/DistancePix;
%1inch x 96 pixels; 1inch = 2.54cm
pixel_size=2.54/96; %cm

%Read files. Batch of pictures from a plate taken along several days
images = dir('*.tif');

for n = 0:length(images)-1

    file= images(n+1).name;
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
    
    %Correct non-uniform ilumination. Adjust specific range
    rgbI = imadjust(croppedImage(:,:,1), [0 0.60],[]);%imshow(rgbI);

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

    rgbBW = rgbI < 50;%imshow(rgbBW)
    %remove connected objects
    rgbI_nobord = imclearborder(rgbBW,8);%imshow(rgbI_nobord)
    %smooth object
    seD = strel('diamond',1);
    rgbI_final = imerode(rgbI_nobord,seD); %rgbI_final = imerode(rgbI_final,seD);

    %Find colonies using boundary
    [B,L2,n2,A] = bwboundaries(rgbI_final,'noholes');

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
        %imshow(croppedImage);
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

    eccentricity{n+1} = Eccentricity(filter2);
    circularity{n+1} = Circularity(filter2);
    perimeter{n+1} = Perimeter(filter2);
    area{n+1} = Area(filter2);
    radii{n+1} = Radii(filter2);
    centroid{n+1} = round(Centroid(filter2,:));
    %test = round(Centroid(filter2,:));
    %%
    %Compare two data points and get the change in radii
    if n>0
        %Check if there are colonies in i-1 not seen in i
        %match will store the colonies found in t-1 and t. It takes the
        %indexes from t
        %new will save the colony from t-1 that wasn't found in t
        %update_radii = zeros(length(centroidOld),1);
        for i = 1:length(centroid{n})
            x = find(centroid{n}(i,1) - 10 <= centroid{n+1}(:,1) & centroid{n}(i,1) + 10 >= centroid{n+1}(:,1));
            y = find(centroid{n}(i,2) - 10 <= centroid{n+1}(:,2) & centroid{n}(i,2) + 10 >= centroid{n+1}(:,2));
            if isempty(x)
                if isempty(y)                 
                    %index from centroid, no colony found in previous
                    %search
                    new_centroid{i} = centroid{n}(i,:);
                    %new_area{i}=area{n}(i);
                    new_radii{i} = radii{n}(i);
                else
                    if length(y) == 1
                        if centroid{n}(i,2)-1 <= centroid{n+1}(y,2) && centroid{n}(i,2)+1 >= centroid{n+1}(y,2)
                            new_centroid{i} = centroid{n+1}(y,:);
                            new_area{i}=area{n+1}(y);
                            new_radii{i} = radii{n+1}(y);
                            match(i) = y;
                        else
                            new_centroid{i} = centroid{n}(i,:);
                            new_area{i}=area{n}(i);
                            new_radii{i} = radii{n}(i);
                        end
                    end
                end
            else
                if isempty(y)
                    new_centroid{i} = centroid{n+1}(x,:);
                    new_area{i}=area{n+1}(x);
                    new_radii{i} = radii{n+1}(x);
                    match(i) = x;
                else
                    if length(x) == length(y)
                        %found in i
                        new_centroid{i} = centroid{n+1}(x,:);
                        new_area{i}=area{n+1}(x);
                        new_radii{i} = radii{n+1}(x);
                        match(i) = x;
                    else
                        in = intersect(x,y);
                        new_centroid{i} = centroid{n+1}(in,:);
                        new_area{i}=area{n+1}(in);
                        new_radii{i} = radii{n+1}(in);
                        match(i) = in;
                    end
                end

            end 
        end 
        
        new = setdiff(1:length(radii{n+1}), match);
        
%         centroid_new{n} = [cat(1,new_centroid{:});centroid{n+1}(new,:)];
%         radii_new{n} = [new_radii';radii{n+1}(new)];
        
        centroid{n+1} = [cat(1,new_centroid{:});centroid{n+1}(new,:)];
        radii{n+1} = [cat(1,new_radii{:});radii{n+1}(new)];
    end  

    %%
    %%Plot the colonies found and save data
    %Labels
    colony = cell(1,length(radii{n+1})); 
    label = cell(length(radii{n+1}),1);
    ID = cell(length(radii{n+1}),1);
    %Give labels
    for i = 1:length(radii{n+1})
        colony{i} = int2str(i);
        label{i}= strcat(day,'-',plate,'-',int2str(i));
        ID{i}=strcat(plate,'-',int2str(i));
    end
    % 
    figure
    imshow(croppedImage);
    viscircles(centroid{n+1},radii{n+1},'EdgeColor','b');
    text(centroid{n+1}(:,1), centroid{n+1}(:,2), colony);
    print(strcat(filename{1},'IDs'),'-dpng');
    close;
% 
%     %Save data
    statsData = struct('label', label, 'ID', ID, 'centroid', centroid, 'perimeter', perimeter, 'area', area, 'radii', radii, 'eccentricity', eccentricity, 'circularity', circularity);
    save(strcat(filename{1},'-data.mat'), 'statsData');

%     %%
%     %Get pixel information per channel
    for k = 1:3
        rgbI = I(:,:,k);
        %Create BW canvas
        imageSize = size(rgbI);
        BW = false(imageSize);

        %The number of identified colonies don't stay the same. I will better
        %save the whole struct with centroids included to filter later.
        imshow(rgbI);
        hold on
        for i = 1:length(radii{n+1})
            h = drawcircle(gca,'Center',centroid{n+1}(i,:),'Radius',(radii{n+1}(i)));
            mask = h.createMask();
            BW = BW | mask;
        end
        hold off
        close

        figure
        imshow(BW)
        Save according to RGB channel
        if k == 1
            %Better to have a different name for the different structures
            statsPixel_red = regionprops(BW,rgbI,{'Centroid','PixelValues', 'MaxIntensity', 'MinIntensity', 'MeanIntensity'});
            save(strcat(filename{1},'-pixel_red.mat'), 'statsPixel_red');
        elseif k == 2
            %Better to have a different name for the different structures
            statsPixel_green = regionprops(BW,rgbI,{'Centroid','PixelValues', 'MaxIntensity', 'MinIntensity', 'MeanIntensity'});
            save(strcat(filename{1},'-pixel_green.mat'), 'statsPixel_green');
        else
            %Better to have a different name for the different structures
            statsPixel_blue = regionprops(BW,rgbI,{'Centroid','PixelValues', 'MaxIntensity', 'MinIntensity', 'MeanIntensity'});
            save(strcat(filename{1},'-pixel_blue.mat'), 'statsPixel_blue');
        end

    end

end 
%end

%[a,b,c] = impixel()

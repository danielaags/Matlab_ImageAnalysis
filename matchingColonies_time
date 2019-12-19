function[] = matchingColonies_time(folderpath)
%The function is desing to match colonies over time after having extracting
%morphological characteristics from each time point. The tif data has to be
%organised by plate and include the different time points to be used.

%path = 'C:\Users\au648169\Documents\Postdoc_TorringLab\ImagingAnalysis\02-SecondRound\10%NBA(190exp)\PerPlate';
path = folderpath;
folders = dir('*nr*');
for a = 1:length(folders)
    subpath = strcat(path, '\', folders(a).name);
    cd(subpath);

    %clear all
    %close all

    files = dir('*data.mat');
    images = dir('*tif');
    %Better to read from the last image to the first
    lst=sort(1:length(files), 'descend');

    for i = 1:length(files)
        load(files(lst(i)).name);
        centroid{i} = statsData.centroid;
    end 

    %%
    %Find the size of the contents of each array. Maximum value to pre-locate
    %memory
    cellsz = cellfun(@length,centroid,'uni',false);
    m = max(cell2mat(cellsz));
    n = length(centroid);


    for k = 1:n-1
        %Clean the index variable to save the new positions
        index={};
        found = [];
        %j=1;
        %Check the lenght of the current centroids list (k), if it is different
        %from the maximum m then the for loop will take the length of the
        %current list (k) otherwise it will take m.
        if length(centroid{k}) ~= m
            l = length(centroid{k});
        else
            l = m;
        end 
        for i = 1:l
                %Check the y access within a range of +-10
                idy = find(centroid{k}(i,1)-10 <= centroid{k+1}(:,1) & centroid{k}(i,1)+10 >= centroid{k+1}(:,1));
                %If not found check withing a range of +-15
                if isempty(idy)
                    idy = find(centroid{k}(i,1)-15 <= centroid{k+1}(:,1) & centroid{k}(i,1)+15 >= centroid{k+1}(:,1));
                end
            %Check the x access within a range of +-10
                idx = find(centroid{k}(i,2)-10 <= centroid{k+1}(:,2) & centroid{k}(i,2)+10 >= centroid{k+1}(:,2));
                %If not found check withing a range of +-15
                if isempty(idx)
                    idx = find(centroid{k}(i,2)-15 <= centroid{k+1}(:,2) & centroid{k}(i,2)+15 >= centroid{k+1}(:,2));
                end

                %If find() empty, asume is not in the other plate and presetve
                %position. Important when the list k has more elements than k+1
                if isempty(idy) && isempty(idx)
                    %Assign empty value
                    index{i} = NaN;
                else
                    %If the searches are not empty, get the intersection, save
                    %the index only if the intersection length is 1 
                    id = intersect(idy, idx);
                    if length(id) == 1
                        index{i} = id;
                        %Keep track of the indexes that match, later it will be
                        %use to identify the new colonies in k+1
                        found = [found id];
                    else
                        index{i} = NaN;
                    end
                end
        end
        %Turn the cell array with the re-ordered indexes to an array
        index = cell2mat(index);
        %Find the indexes not identified
        new = setdiff(1:length(centroid{k+1}(:,1)), found);
        %Build a new array with the re+-ordered and the new colonies from k+1
        indexN= [index new];
        %Save the indexes to use later
        %ID = indexN;

        %Re-order
        load(files(lst(k+1)).name);
        file = strsplit(files(lst(k+1)).name,'.');
        %Transform into table
        arrayfile = struct2array(statsData);

        n = length(indexN);
        %Create array. Since my data was all converted to strings, an empty
        %array of strings is needed.
        arraynew = strings(n, 22);
        for j = 1:n
            if isnan(indexN(j))
            else
                arraynew(j,:) = arrayfile(indexN(j),:);
            end
        end

        centroid{k+1} = [double(arraynew(:,4)), double(arraynew(:,5))];
        
        %%
        %Print new labels on plate
        colony = string(1:length(centroid{k+1}(:,1)));
        %Find the nan values
        tomap = ~isnan(centroid{k+1}(:,1));
        diameter = double(arraynew(:,7));
        
        %On figure
        I = imread(images(lst(k+1)).name);
        figure
        imshow(I);
        viscircles(centroid{k+1}(tomap,:),diameter(tomap)/2,'EdgeColor','b');
        text(centroid{k+1}(tomap,1), centroid{k+1}(tomap,2), colony(tomap));
        print(strcat(file{1},'-reIDs'),'-dpng');
        close;

        %Save re-ordered data
        dataReordered = struct('label', arraynew(:,1), 'sample', arraynew(:,2),...
            'ID', arraynew(:,3), 'centroids', double(arraynew(:,4:5)),...
            'area', double(arraynew(:,6)), 'diameter', diameter, 'perimeter', double(arraynew(:,8)),...
            'circularity', double(arraynew(:,9)), 'eccentricity', double(arraynew(:,10)),'RGB_mean', double(arraynew(:,11:13)),...
            'RGB_std', double(arraynew(:,14:16)), 'RGBt_mean', double(arraynew(:,17:19)), 'RGBt_std', double(arraynew(:,20:22)));

        save(strcat(file{1},'-Reordered.mat'), 'dataReordered');


    end
    
end    
end


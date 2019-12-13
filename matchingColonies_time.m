clear all
close all

files = dir('*data.mat');
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
%Create a new cell to save the new indexes
ID = cell(1,n-1);

for k = 1:n-1
    %Clean the index variable to save the new positions
    index={};
    found = [];
    %j=1;
    %Check the lenght of the current centroids list (k), if it is different
    %from the maximum m then the for loop will take the lenght od the
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
    new = setdiff(1:length(centroid{k+1}), found);
    %Build a new array with the re+-ordered and the new colonies from k+1
    indexN= [index new];
    %Save the indexes to use later
    ID{k} = indexN;
end

%%
%Re-order
for i = 1:length(ID)
    load(files(lst(i+1)).name);
    file = strsplit(files(lst(i+1)).name,'.');
    %Transform into table
    arrayfile = struct2array(statsData);
    %Lenght new array
    n = length(ID{i});
    %Create array. Since my data was all converted to strings, an empty
    %array of strings is needed.
    arraynew = strings(n, 15);
    for j = 1:n
        if isnan(ID{i}(j))
        else
            arraynew(j,:) = arrayfile(ID{i}(j),:);
        end
    end
    
    dataReordered = struct('label', arraynew(:,1), 'ID', arraynew(:,1), 'centroids', double(arraynew(:,3:4)),...
        'area', double(arraynew(:,5)), 'diameter', double(arraynew(:,6)), 'perimeter', double(arraynew(:,6)),...
        'circularity', double(arraynew(:,7)), 'eccentricity', double(arraynew(:,8)),'R_mean', double(arraynew(:,9)),...
        'R_std', double(arraynew(:,10)), 'G_mean', double(arraynew(:,11)), 'G_std', double(arraynew(:,11)),...
        'B_mean', double(arraynew(:,12)), 'B_std', double(arraynew(:,13)));
    save(strcat(file{1},'-Reordered.mat'), 'dataReordered');
end




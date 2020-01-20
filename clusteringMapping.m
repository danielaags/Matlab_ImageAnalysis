clear all
images = dir('*.tif');
load('pixel_bio-Data.mat', 'data_isolation')

%Extract data to use
%X = cat(1,data_isolation.RGB_mean); %meanArea
X = cat(1,data_isolation.RGBt_std); %stdTrasversal

%R saves the row numbers that were removed
[X, R] = rmmissing(X);

%Consider rows without NaN
nrow = logical(abs(R-1));

%Clustering
C0 = clusterdata(X,'linkage','ward', 8);

%Mapping variables
id = cat(1,data_isolation.ID2);
id = id(nrow);
%Extract the plate number and order
plate = split(id, '-');
plate = double(plate(:,1));
%Extract centroids
centroid = cat(1,data_isolation.centroid);
centroid = centroid(nrow,:);
%Extract diameters
diameter = cat(1,data_isolation.diameter);
diameter = diameter(nrow);

%Colors to use during mapping
%For 8 clusters
col = [61 38 168; 
    71 80 243;
    45 135 247;
    18 176 215;
    56 200 149;
    170 199 57;
    254 195 55;
    248 250 19];

col = col/255;

for j = 1:length(images)
    %j=1;
    file= images(j).name;
    I = imread(file);
    
    %Find index where the plate belongs
    temp_idx = find(plate == j);
    %Find the centroids found in that plate and the raddii for the colonies
    temp_centroid = centroid(temp_idx,:);
    temp_radii = diameter(temp_idx)/2;
    %Find the corresponding cluster the colony belongs to
    %Mean
    temp_cluster_mean = C0(temp_idx);
    %Standard deviation
    %temp_cluster_sd = C1(temp_idx);
    
    %Gets the cluster that correspond to the specific plate
    %clusters = unique(temp_cluster);
    figure
    imshow(I);
    hold on
    %Goes along the different clusters in the plate and plots
%     for i = 1:length(clusters) 
%         idx_c = find(temp_cluster == clusters(i));
%         viscircles(temp_centroid(idx_c,:), temp_radii(idx_c), 'Color', col(clusters(i),:));
%     end  
    for i = 1:length(temp_cluster_mean) 
        y = temp_centroid(i,2);
        x = temp_centroid(i,1);
        r = temp_radii(i);
        %rectangle('Position', [x-r y-r r r], 'EdgeColor', col(temp_cluster_sd(i),:), 'LineWidth', 2);
        rectangle('Position', [x-r y-r r*2 r*2], 'EdgeColor', col(temp_cluster_mean(i),:), 'LineWidth', 2);
        
    end 
    hold off
    
    output = split(file, '.');
    print(output{1},'-dpng');
    close;
end

%Remove missing values and add clustering information
T = struct2table(data_isolation);
dataT = T(nrow, :);
%dataT.sumActivity = sum(table2array(dataT(:,26:38)), 2);
%dataT.cluster = C0;
%writetable(dataT, 'pixel_bio_cluster-Data.xls')

%Save as excel

%Save as .mat
cluster_data_isolation = table2struct(dataT);
%save('pixel_bio_cluster-Data.mat', 'cluster_data_isolation');

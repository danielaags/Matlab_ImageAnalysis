% clear all
% close all
load('pixel_bio-Data.mat', 'data_isolation')

T = struct2table(data_isolation);

%Extract data to use
X = cat(1,data_isolation.RGB_mean); %meanArea
%X = cat(1,data_isolation.RGB_std); %stdArea
%X = cat(1,data_isolation.RGBt_mean); %areaTrasversal
%X = cat(1,data_isolation.RGBt_std); %stdTrasversal

%R saves the row numbers that were removed
[X, R] = rmmissing(X);

%Consider rows without NaN
nrow = logical(abs(R-1));
Tfilter = T(nrow,:);

%Find activity = 1 for each pathogen
for i = 1:13
    idx = i + 26;
    pick = table2array(Tfilter(:,idx)) == 1;
    table = [Tfilter(pick,1), Tfilter(pick,idx)];
    writetable(table, strcat(names{idx}, '-ForActivityCheck.xls'))
end

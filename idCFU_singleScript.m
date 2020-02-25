%Inputs
day0 = '190921';
%singleFolder 1 for analysis of only one folder and 0 for more than one
%folder
singleFolder = 0;
filename_matfile = 'FirstRound-';

%%
%Give the path where your folders are located
path = pwd;

%%
%Single or multiple data points

%Single folder
if singleFolder == 1
    fnames = dir('*tif');
    numfids = length(fnames);

    %Go to each file
    for k = 1:numfids
      %Read the tif file
      file = strsplit(fnames(k).name,'.');
      plate =  num2str(k);
      %Run the desired ID function.
      %IDcfu_datapixelFinal(day, num2str(k), file{1});
      IDcfu_datapixelLabFinal_transversalcolonyPixel(day0, plate, file{1});
    end
    
    dataCollectionRGBLab_matfiles(filename_matfile);

else
    %To check in more than one folder
    folders = dir('*nr*');
    %day of the experiment
    day = day0;
    for i = 1:length(folders)
        subpath = strcat(path, '\', folders(i).name);
        cd(subpath);
        %Batch_idCFU('191123', i)
        %Get all the the tif files in a folder
        fnames = dir('*tif');
        numfids = length(fnames);
        plate =  num2str(i);

        %Go to each one and gets do identification 
        for k = 1:numfids
          %Read the tif file
          file = strsplit(fnames(k).name,'.');
          %Run the desired ID function.
          %IDcfu_datapixelFinal(day, num2str(k), file{1});
          %IDcfu_datapixelFinal_transversalcolonyPixel(day, plate, file{1});
          IDcfu_datapixelLabFinal_transversalcolonyPixel(day, plate, file{1});
        end
    end
    
    cd(path);
    matchingColoniesRGBLab_timegrowth(path);
    
    %%
    %Collect all the .mat data
    cd(path);
    %dataCollectionRGBLab_matfiles(filename_matfile);
    dataCollectionRGBLab_growth_matfiles(filename_matfile);
end

%%
%The script requires three different functions: i)ID_ function, ii)
%matching colonies time, iii)data_collection_matfiles. Depending
%on the identificyation function use, the two other functions need to be 
%tailored to the specific output. Since the number of columns might change
%depending on first function used.

%%
%Inputs
day0 = '191101';
%singleFolder 1 for analysis of only one folder and 0 for more than one
%folder
singleFolder = 0;
filename_matfile = 'test';

%%
%Give the path where your folders are located
folderpath = pwd;

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
      IDcfu_Final(day, plate, file{1});
    end   
else
    %To check in more than one folder
    folders = dir('*nr*');
    %day of the experiment
    day = day0;
    for i = 1:length(folders)
        subpath = strcat(folderpath, '\', folders(i).name);
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
          IDcfu_Final(day, plate, file{1});
        end
    end
    
    cd(folderpath);
    matchColonies_Final(folderpath);

end

%Collect all the .mat data
cd(folderpath);
dataCollection_Final(filename_matfile);

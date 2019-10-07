file = 'agar1mMp';
I = imread(strcat(file,'.tif'));

%Get the image size to remove the rim of the petri dish
imageSize = size(I);
%center and radius of circle ([c_row, c_col, r]). The images are a bit
%off. Center doesn't seem to be at 1024, 1024
ci = [1044, 1024, 850];     
[xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
mask = uint16((xx.^2 + yy.^2)<ci(3)^2);
croppedImage = uint16(zeros(size(I)));
croppedImage(:,:,1) = I(:,:,1).*mask;
croppedImage(:,:,2) = I(:,:,2).*mask;
croppedImage(:,:,3) = I(:,:,3).*mask;


%CFUradius=100; % CFU, maximum size to take into account
[centers1, radii1, metric1] = imfindcircles(croppedImage,[20, 60],'Sensitivity',0.95);
[centers2, radii2, metric2] = imfindcircles(croppedImage,[250, 300],'Sensitivity',0.97);

%data = table(centers);
%writetable(data,strcat(file,'_centers.txt'),'Delimiter','\t');

%data = table(radii);
%writetable(data,strcat(file,'_radii.txt'),'Delimiter','\t');


imshow(croppedImage);
viscircles(centers,radii,'EdgeColor','b');

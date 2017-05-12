function [BW,area,perimeter,max_l,min_l,orientation,center]=trace_grain_res_manu(file,res,r)

% read image
I=imread(file);
I=imresize(I,r);
% draw outline
BW = roipoly(I);
% plot original image with estimated outline
imshow(I), title('outlined original image');
% properties of image
area=regionprops(BW,'Area');
area.Area1=area.Area*(res/r)^2;
perimeter=regionprops(BW,'Perimeter');
perimeter.Perimeter1=perimeter.Perimeter*(res/r);
max_l=regionprops(BW,'MajorAxisLength');
max_l.MajorAxisLength1=max_l.MajorAxisLength*(res/r);
min_l=regionprops(BW,'MinorAxisLength');
min_l.MinorAxisLength1=min_l.MinorAxisLength*(res/r);
orientation=regionprops(BW,'Orientation');
center=regionprops(BW,'Centroid');
% script to load files from folder to calculate FTs from all possible
% combinations of images taken parallel and perpendicular to the
% crystallographic c-axis

mineral='ap';
R=0.4; % original resolution of images in microns per pixel
res=0.5; % reduce resolution by 0...1; final resultion is R/res
r232_238=0; % ratio of 232Th and 238U, measured in mol
r147_238=0; % ratio of 147Sm and 238U, measured in mol 
rf=pwd; % current folder
root_data=[rf '/Test/']; % folder with images
listing = dir(fullfile(root_data,'*.jpg'));
details=struct2cell(listing);
file_names=details(1,:);
file_name='Test-1'; % string of files that does not change from image to image

Index_par = strncmp(file_names, [file_name '-par'],10); % find images with c-axis parralel
Files_par = file_names(Index_par);
Index_per = strncmp(file_names, [file_name '-per'],10); % find images with c-axis perpendicular
Files_per = file_names(Index_per);
o=0;
for j=1:sum(Index_par)
    for k=1:sum(Index_per)
        o=o+1;
        [Ft(o),mass(o),matrix_3D]=calcFTphoto(mineral,char(strcat(root_data,Files_par(j))),char(strcat(root_data,Files_per(k))),R,res,r232_238,r147_238);
        Ft
        mass
        % show 3D model, you may want to comment this out if not needed
        plot_grain(matrix_3D,0.5)
    end
end

% script to load files from folder to calculate FTs from all possible
% combinations of images taken parallel and perpendicular to the
% crystallographic c-axis
% written by C. Glotzbach, 04/2017

% input parameters, adjust to your sample
mode='normal';  % 'normal' mode for grains or broken grains (example: Test-1), 
                % 'cut' mode for grains grained and polished (example: Test-2)
mineral='ap';
R=0.3937; % original resolution of images in microns per pixel
res=0.3281; % reduce resolution by 0...1; final resultion is R/res
r232_238=0.4444; % ratio of 232Th and 238U, measured in mol; if 0 mean values are used
r147_238=0.0615; % ratio of 147Sm and 238U, measured in mol; if 0 mean values are used
root_data='/Test/'; % folder that contains the images, example folder is /test/
file_name='Test-1'; % string of files that does not change from image to image, example grains: 
                    % 'Test-1' : used with the mode 'normal'
                    % 'Test-2' : used with the mode 'cut'
                    

% do not change from here
rf=pwd; % current folder
root_data=[rf root_data];
listing = dir(fullfile(root_data,'*.jpg'));
details=struct2cell(listing);
file_names=details(1,:);


Index_par = strncmp(file_names, [file_name '-par'],length(file_name)+4); % find images with c-axis parralel
Files_par = file_names(Index_par);
Index_per = strncmp(file_names, [file_name '-per'],length(file_name)+4); % find images with c-axis perpendicular
Files_per = file_names(Index_per);

o=0;
done_par=[];
done_per=[];
for j=1:sum(Index_par)
    for k=1:sum(Index_per)
        o=o+1;
        if any(done_par==j)
            image_par=images.(['image_par' num2str(j)]);
        else
            image_par=[];
        end
        if any(done_per==k)
            image_per=images.(['image_per' num2str(k)]);
        else
            image_per=[];
        end
        switch mode
            case 'normal'
            % calculate Fts of intact or broken grains
            [Ft_ana(o),Ft_num(o),mass_ana(o),mass_num(o),vol_photo,images.(['image_par' num2str(j)]),images.(['image_per' num2str(k)])]=calcFTphoto(mineral,char(strcat(root_data,Files_par(j))),char(strcat(root_data,Files_per(k))),R,res,r232_238,r147_238,image_par,image_per);
            case 'cut'
            % calculate Fts of grains cut by polishing parallel to the c-axis
            [Ft_num(o),mass_num(o),vol_photo]=calcFTphotoCut(mineral,char(strcat(root_data,Files_par(j))),char(strcat(root_data,Files_per(k))),R,res,r232_238,r147_238);
        end
        done_par=[done_par j];
        done_per=[done_per k];
    end
end

% plot grain model
plot_grain(vol_photo,0.5)

% % make video of slices
% min_x=0;
% min_y=0;
% max_x=size(vol_photo,1);
% max_y=size(vol_photo,2);
% k=0;
% hFig=figure('Color',[0 0 0]);
% set(gca,'Color',[0 0 0])
% set(gca,'YColor',[1 1 1])
% set(gca,'XColor',[1 1 1])
% set(gca,'ZColor',[1 1 1])
% set(hFig,'Position',[200,200,600,750])
% hold on
% for i=200:450
%     k=k+1;
%     h=surf([min_x max_x],[min_y max_y],repmat(i, [2 2]),...
%     squeeze(vol_photo(i,:,:))*i,'facecolor','texture')
%     view(45,30);
%     axis([150 500 30 180 200 450])
%     set(h,'FaceAlpha',  'texturemap', 'AlphaDataMapping', 'none', 'AlphaData',~isnan(squeeze(vol_photo(i,:,:))));
%     M(k)=getframe(gcf);
% end
% V=VideoWriter('Model_Grain20.avi');
% V.FrameRate=50;
% open(V)
% writeVideo(V,M)
% close(V)

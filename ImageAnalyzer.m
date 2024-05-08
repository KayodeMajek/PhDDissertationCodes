%Written by Olukayode T. Majekodunmi 
%For a PhD Dissertation titled "Discontinuous Colloidal Clogging in Tapered Microchannels"
%Submitted to the Department of Chemical Engineering, Northeastern University, Boston MA
%Submitted in May, 2024.
%See the accompanying readme.txt file for a functional description

%% Input area
clear; clc; close all;

%% Specify tiling rows and columns
nrows = 7; 
ncols = 21;
%% Image folder directory
f = 'C:\Users\Tilescan\*.tif';
files = dir(f);

%% Preallocation for speed & conversion rate
prelength = 10000; %predefined length of the preallocated arrays
cloglength = zeros(1,prelength); %preallocated array for cloglengths - Major axes of the centroids - excess zeros will deleted or turned into NaN in the final vector for histograming 
clogwidth = zeros(1,prelength);  %preallocated array for clogwidths - Minor axes of the centroids - excess zeros will deleted or turned into NaN in the final vector for histograming
cntrsx = zeros(1,prelength);
cntrsy = zeros(1,prelength);
names = fullfile({files.folder}, {files.name});
numfiles = length(names);
imdata = cell(1,numfiles);
imtotal = cell(1,nrows);
meter2pixel = 2.631726e-3; %%conversion rate from pixels to meter - will change based on camera magnification - always check image metadata (".xml")
img = cell(1,ncols);
%% Image reading
for k = 1:numfiles
    N = names{k};
    im = imread(N);
    %im = flipud(im); % if you need to flip individual images about the horizontal axis 
    im(1839:end,:)=[]; %delete last 210 rows in each image (to remove the auto stitching correction in each image).
    im(:,1839:end)=[]; %delete last 210 columns in each image (to remove the auto stitching correction in each image).
    imdata{k}=im; %collection of all image data
end

%%  Image tiling for individual rows
for row = 1:nrows %tile each row of the tilescan separately
    imtotal{row} = imtile(imdata, 'Frames', (row-1)*ncols+1:row*ncols, 'GridSize', [1 ncols]);
end 

%% Replace the even-numbered rows after reordering because of the microscope tilescan process (left to right - odd row; right to left - even row)
for row = 2:2:nrows 
   dim = flip((row-1)*ncols+1:row*ncols); %index of images on the referenced even-numbered row
    for t = 1:ncols
        img{t} = imdata{dim(t)}; %collects the image data into a new cell for horizontal tiling
    end
    imtotal{row} = horzcat(img{:}); %use horzcat or imtile in one column to horizontally concatenate the referenced images
    %imtotal{row} = imtile(imdata, 'Frames', dim, 'GridSize', [1 ncols]);
end

im = vertcat(imtotal{:}); %final image of clogged channels

%im = fliplr(im); to flip tiled image about the vertical axis
figure 
imshow(im); %tiled image
%imagesc(im); colormap gray; hold on %imagesc for scaled view

%% Clogs identification
 if min(size(im))==3 %depending on the dimension of the original image file (no worries for uint16 RAW images from the microscope). 
    im = im(:,:,1);
 end 

level = graythresh(im);
bwim = imbinarize(im,level); %binarized image
smallestAcceptableArea = 100;  %adjust this to remove single particles (not clogs) from the estimation
bwim = bwareaopen(bwim,smallestAcceptableArea);
shapeprops=regionprops(bwim,'Centroid','MajorAxisLength','MinorAxisLength');
numshapes=length(shapeprops);

figure 
imshow(bwim); hold on %binarized tiled image with clogs identified
%imagesc(bwim); colormap gray; %for scaled view

for i=1:numshapes 
    cloglength(i)=shapeprops(i).MajorAxisLength; %Cloglength
    clogwidth(i)=shapeprops(i).MinorAxisLength; %Clogwidth
    cntrsx(i) =shapeprops(i).Centroid(1); % x-axis of Centroid for image "k"
    cntrsy(i) =shapeprops(i).Centroid(2); % y-axis of Centroid for image "k"
    plot(cntrsx(i),cntrsy(i),'g.','MarkerSize',10); %plot of Centroids on image "k"
end
hold off

%% Remove extra zeros from clog data and also remove identified objects greater than 40 microns
cnvrt = meter2pixel/2048; %data conversion from pixels to centimeter
DataCloglength = ((cnvrt.*cloglength)*1e6)';
DataClogwidth = ((cnvrt.*clogwidth)*1e6)';

Dummy = zeros(length(cloglength),1); %creates a dummy vector to compare with the data vector and delete
DataCloglength(Dummy == DataCloglength) = []; %for loop doesn't work in deleting the zero elements because the vector is shrinking in size upon each iteration but it can change the zeros to NaN.
Dummy = zeros(length(clogwidth),1);
DataClogwidth(Dummy == DataClogwidth) = [];

GreaterThan40 = length(find(DataClogwidth > 40));
 
p = length(DataClogwidth); %To remove values greater than 40um for uniform histogram scaling purposes.
for q = 1:p      
   if DataClogwidth(q) > 40
    DataClogwidth(q) = NaN; %Turn all the values greater than 40um into NaN
   end 
end

TF=isnan(DataClogwidth); %To make the cloglengths corresponding to the "NaN" clogwidths (i.e. > 40um) also "NaN". 
DataCloglength(TF)=NaN;

DataClogwidth(isnan(DataClogwidth))=[]; %delete all elements that are "NaN"
DataCloglength(isnan(DataCloglength))=[];


%% Histogramming Cloglengths and Clogwidths (with NONCLOGS particles)
figure  
edges = 0:50:max(DataCloglength); 
histogram(DataCloglength,edges) 
xlabel('Clog length ({\mu}m)','fontweight','bold'); ylabel('Count','fontweight','bold');
title('With nonclogs')

figure
edges = 0:2:40; %clogs cannot be wider than 4o micrometer
histogram(DataClogwidth,edges)
xlabel('Clog width ({\mu}m)','fontweight','bold'); ylabel('Count','fontweight','bold');
title('With nonclogs')
%% Data collection for centroid locations of clogs: clogs position in the image. 
DataCntrsx = ((cnvrt.*cntrsx)*1e2)'; %convert to centimeter. Length of channel is ~5cm. 
DataCntrsy = ((cnvrt.*cntrsy)*1e2)'; %convert to centimeter. Width of device is ~2cm. 

Dummy = zeros(length(cntrsx),1); %creates a dummy vector to compare with the data vector and delete
DataCntrsx(Dummy == DataCntrsx) = []; %for loop doesn't work in deleting the zero elements because the vector is shrinking in size upon each iteration but it can change the zeros to NaN.
Dummy = zeros(length(cntrsy),1);
DataCntrsy(Dummy == DataCntrsy) = [];

%To make the centroid positions corresponding to the "NaN" clogwidths (i.e. > 40um) also "NaN". 
DataCntrsx(TF)=NaN;
DataCntrsy(TF)=NaN;
DataCntrsx(isnan(DataCntrsx))=[]; %delete all elements that are "NaN"
DataCntrsy(isnan(DataCntrsy))=[];

%% Histogramming the clog locations (with NONCLOGS particles)
figure
edges = 0:0.25:5; %length of each channel is approximately 5cm
histogram(DataCntrsx,edges)
xlabel('Clog centroid x-location (cm)','fontweight','bold'); ylabel('Count','fontweight','bold');
title('With nonclogs')

figure 
edges = 0:0.2:2; %width of device is about 2cm
histogram(DataCntrsy, edges)
xlabel('Clog centroid y-location (cm)','fontweight','bold'); ylabel('Count','fontweight','bold');
title('With nonclogs')
%% Location of clog start and end points (with NONCLOGS particles) - To be used only when needed.
DataCntrsx_start = DataCntrsx - 0.5*(DataCloglength/1e4); %x-location of clog start is the centroid minus half the clog length  %convert for units uniformity
%figure  
%histogram(DataCntrsx_start)
%xlabel('Clog centroid x-location start point (cm)','fontweight','bold'); ylabel('Count','fontweight','bold');           
                                                 
DataCntrsx_end = DataCntrsx + 0.5*(DataCloglength/1e4); %x-location of clog start is the centroid plus half the clog length  %convert for units uniformity
%figure  
%histogram(DataCntrsx_end)
%xlabel('Clog centroid x-location end point (cm)','fontweight','bold'); ylabel('Count','fontweight','bold');   

%% Plotting clog dimensions versus locations (with NONCLOGS particles)
figure
plot(DataCntrsx, DataClogwidth, 'r.')
xlabel('Clog centroid x-location (cm)','fontweight','bold'); ylabel('Clog width ({\mu}m)','fontweight','bold');
axis([0 5 0 40]) 
title('With nonclogs')

figure
plot(DataCloglength, DataClogwidth, 'b.')
xlabel('Clog length ({\mu}m)','fontweight','bold'); ylabel('Clog width ({\mu}m)','fontweight','bold');
axis([0 5000 0 40]) 
title('With nonclogs')

%% Fractional Area of Clogs
AreaClogs = sum(DataClogwidth.*DataCloglength); %Rectangular area occupied by clogs
ChnArea = 165e6; %Total channel area = 165e6 (um)^2; area of one channel = 2*area of right-angle triangle (half channel) = 20 um * 5cm = 1e6 (um)^2.
FractionalAreaClogs = (AreaClogs/ChnArea)*100; %Percentage of total channels area occupied by clogs. 

%% Deleting the black spaces in the y-direction 
cntrsy_p = cntrsy';
Dummy = zeros(length(cntrsy),1);
cntrsy_p(Dummy == cntrsy_p) = [];

cntrsy_p(TF)=NaN;
cntrsy_p(isnan(cntrsy_p))=[]; %delete all elements that are "NaN" (clog widths greater than 40um)
n_bwim = bwim;

rg = 200; %half width of a channel in pixels
mgy1 = round(min(cntrsy_p)); %minimum y-centroid location (first channel midpoint)
mgy2 = round(max(cntrsy_p)) - numel(1:mgy1-rg); %maximum y-centroid location (last channel midpoint)% the substraction is needed to adjust the position of the maximum y-centroid since the deletion of the top blank areas shifts the array upwards.


n_bwim(1:mgy1-rg,:) = [];
n_bwim(mgy2+rg:end,:) = [];

%% Deleting empty spaces in the x-direction (at the outlet only). 

cntrsx_p = cntrsx';
Dummy = zeros(length(cntrsx),1);
cntrsx_p(Dummy == cntrsx_p) = [];

cntrsx_p(TF)=NaN;
cntrsx_p(isnan(cntrsx_p))=[]; %delete all elements that are "NaN" (clog widths greater than 40um)

mgx = round(min(cntrsx_p)); %minimum x-centroid location
n_bwim(:, 1:mgx-rg) = []; %new processed image

%figure
%imshow(n_bwim) %this shows the tiled image of the clogged channels with extra empty spaces at the outlet and sides removed

%% Final image with new x and y-centroids identification and histogramming

[DataCntrsy_n, DataCntrsx_n, DataCloglength_n, DataClogwidth_n] = sppcrop(n_bwim,cnvrt);
DataCntrsx_n(TF)=NaN;
DataCntrsy_n(TF)=NaN;
DataCntrsx_n(isnan(DataCntrsx_n))=[]; %delete all elements that are "NaN" to remove clogs > 40um
DataCntrsy_n(isnan(DataCntrsy_n))=[];

DataCloglength_n(TF)=NaN;
DataClogwidth_n(TF)=NaN;
DataCloglength_n(isnan(DataCloglength_n))=[]; %delete all elements that are "NaN" to remove clogs > 40um
DataClogwidth_n(isnan(DataClogwidth_n))=[];

%% Data Collection (x and y locations have been redetermined. Data contains NONCLOGS but without objects > 40um)
delete('WithNonClogsResults.xlsx') %to enable overwriting of any existing 'Results' file. 

T = table(DataCloglength_n, DataClogwidth_n, DataCntrsx_n, DataCntrsy_n); %Data contains NONCLOGS but without clogs > 40um. 
filename = 'WithNonClogsResults.xlsx';
writetable(T,filename);

%% Removing nonClogs: "clogs" that are smaller than the clogs that do fit the channel width. 
[n_data, n_clogs, excl] = nonclogs(DataCntrsy_n, DataCntrsx_n, DataCloglength_n, DataClogwidth_n);

%Data collection without NONCLOGs is obtained from the "sppcrop" function and available in the "UpdatedResults.xlsx" file

%% Export the final image file
%imwrite(n_bwim,'export.tif')
%%
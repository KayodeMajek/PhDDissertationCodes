%Written by Olukayode T. Majekodunmi 
%For a PhD Dissertation titled "Discontinuous Colloidal Clogging in Tapered Microchannels"
%Submitted to the Department of Chemical Engineering, Northeastern University, Boston MA
%Submitted in May, 2024.
%"sppcrop" is a function called in the main "ImageAnalyzer.m" file. 
%See the accompanying readme.txt file for a functional description

function [DataCntrsy_n, DataCntrsx_n, DataCloglength_n, DataClogwidth_n] = sppcrop(n_bwim,cnvrt)
%SPPCROP takes a tiled image "n_bwim", runs the regionprops function to
%identify all clogs (bright spots or streaks) in the image and convert the
%values of the clog length and width from pixels to micrometer and that of
%the clog locations to centimeter

shapeprops_n=regionprops(n_bwim,'Centroid','MajorAxisLength','MinorAxisLength');
numshapes_n=length(shapeprops_n);

figure
imshow(n_bwim); hold on %processed tile image with clogs identified
%imagesc(bwim); colormap gray; %for scaled view

for i=1:numshapes_n 
    cloglength_n(i)=shapeprops_n(i).MajorAxisLength; %Cloglength
    clogwidth_n(i)=shapeprops_n(i).MinorAxisLength; %Clogwidth
    cntrsx_n(i) =shapeprops_n(i).Centroid(1); % x-axis of Centroid for image "k"
    cntrsy_n(i) =shapeprops_n(i).Centroid(2); % y-axis of Centroid for image "k"
    plot(cntrsx_n(i),cntrsy_n(i),'g.','MarkerSize',10); %plot of Centroids on image "k"
end

title('Blank spaces removed on all sides except the channels inlet area')
hold off

DataCntrsx_n = ((cnvrt.*cntrsx_n)*1e2)'; %convert to centimeter. Length of channel is ~5cm. 
DataCntrsy_n = ((cnvrt.*cntrsy_n)*1e2)'; %convert to centimeter. Width of device is ~2cm. 


DataCloglength_n = ((cnvrt.*cloglength_n)*1e6)';%convert to micrometer. 
DataClogwidth_n = ((cnvrt.*clogwidth_n)*1e6)';


end

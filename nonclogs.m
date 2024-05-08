%Written by Olukayode T. Majekodunmi 
%For a PhD Dissertation titled "Discontinuous Colloidal Clogging in Tapered Microchannels"
%Submitted to the Department of Chemical Engineering, Northeastern University, Boston MA
%Submitted in May, 2024.
%"nonclogs" is a function called in the main "ImageAnalyzer.m" file. 
%See the accompanying readme.txt file for a functional description

function [n_data, n_clogs, excl] = nonclogs(DataCntrsy_n, DataCntrsx_n, DataCloglength_n, DataClogwidth_n)

%NONCLOGS remove clogs with width that are outside a tolerance of two
%particle diameters of the width of the channel at the point of clogging. 
% It also confirms the taper angle from the slope of the plot of clogwidth
% versus cloglocation in x-axis.

    shift = 8.5; % axis shift by two particle diameters
    ClogLength = DataCloglength_n;
    ClogWidth = DataClogwidth_n;
    ClogXlocation = DataCntrsx_n; %this is the clog centroid x-location
    %ClogXlocation = DataCntrsx_n - (0.5*DataCloglength_n)*1e-4; %this is the clog start x-location
    ClogYlocation = DataCntrsy_n;

    initial = length(ClogYlocation); %Initial number of clogs

    [a,b] = max(ClogWidth);
    [c,d]  = min(ClogXlocation);

    slope = (a-ClogWidth(d))/((ClogXlocation(b) - c));
    zx = (linspace(c, max(ClogXlocation), length(ClogXlocation)))';
    line = slope*zx + ClogWidth(d);
    %line = slope*zx + min(ClogWidth);
    

    figure
    plot(ClogXlocation, ClogWidth, 'r.')
    hold on
    plot(zx, line, 'b--', 'LineWidth',1.5)
    plot(zx, line-shift, 'k--', 'LineWidth',1.5) %with tolerance = +/- 2particle diameters
    plot(zx, line+shift, 'k--', 'LineWidth',1.5)
    hold off
    ylim([-10 60])
    xlim([0 5])
    xlabel('x_i (cm)','fontweight','bold'); ylabel('w_i','fontweight','bold');

    p_up = [slope (ClogWidth(d)+shift)];
    y_up = polyval(p_up,ClogXlocation);
    y = y_up - ClogWidth;
        for i = 1:length(y)
            if y(i) < 0
                ClogWidth(i) = NaN;
                ClogLength(i) = NaN;
                ClogXlocation(i) = NaN;
                ClogYlocation(i) = NaN;
            end
        
        end

    p_low = [slope (ClogWidth(d)-shift)];
    y_low = polyval(p_low,ClogXlocation);
    y = y_low-ClogWidth;
        for i = 1:length(y)
            if y(i) > 0
                ClogWidth(i) = NaN;
                ClogLength(i) = NaN;
                ClogXlocation(i) = NaN;
                ClogYlocation(i) = NaN;
            end
        end
        
    ClogWidth(isnan(ClogWidth))=[];
    ClogLength(isnan(ClogLength))=[];
    ClogXlocation(isnan(ClogXlocation))=[];
    ClogYlocation(isnan(ClogYlocation))=[]; 


    %figure
    %plot(ClogXlocation, ClogWidth, 'r.')
    %hold on
    %plot(zx, line, 'b--', 'LineWidth',1.5)
    %plot(zx, line-shift, 'k--', 'LineWidth',1.5) %with tolerance = +/- 2particle diameters
    %plot(zx, line+shift, 'k--', 'LineWidth',1.5)
    %hold off

  
  n_data = [ClogLength, ClogWidth, ClogXlocation, ClogYlocation];

  excl = [zx,line,line-shift,line+shift]; %export data to plot the nonclogs figure
    
 final = length(ClogYlocation); %final number of clogs

n_clogs = [initial final (initial - final)]; %array containing the initial and final number of clogs, and number of nonclogs. 

T = table(ClogLength, ClogWidth, ClogXlocation, ClogYlocation);
delete('WithoutNonClogs.xlsx') %to enable overwriting of any existing 'Results' file. 

filename_n = 'WithoutNonClogs.xlsx';
writetable(T,filename_n);

end
 %% This is used for processing micro-imaging data
% Data from JASCO .jwa (export with coordinate info) >> .csv
% Then, .csv is used as input for this script

% 2023/01/13 v1 By Zijiang Yang
% 2023/01/14 v2 binary plot
% 2023/01/14 v3 binary plot - with detected plastics
% 2023/01/16 v4 final output displays multiple types
% 2023/01/17 v5 Remove DSW method - pruningf is enough (´;ω;`)ヾ(･∀･`)
% 2023/01/17 v5bc using avg(all spectra) to determine background
% 2023/01/19 v6 - change v4 to nearest for spatial interpolation; estimated r = 1/2*rv


%% Prepare workspace
warning off

clc, clear, close all

tic

%% Input data

para.rv = 100; % revolution of measurement


Namesmp = '..\RMS Project\Data\EXP02  ThreeSize PE rv100.csv';
Namestd = '..\RMS Project\Data\Raman Standard Plastics.csv';

T.smp = readtable(Namesmp,'PreserveVariableNames',true); % 
% smp = sample
T.std = readtable(Namestd,'PreserveVariableNames',true); % smp = sample

% ID information
ID.smp     = T.smp.Properties.VariableNames(2:end);         % Spectrum ID
ID.std     = T.std.Properties.VariableNames(2:end); 

%% Parameters
dx = para.rv;    % measurement resolution of x-direction
dy = para.rv;    % measurement resolution of y-direction
para.dd = dx/2;  % revolution for spatial interpolation

para.wd = 25;    % window size; may be 10:5:200 for optimization; window size and step size for baseline Lenz et al., 2015 para.ws = 128 para.ss = 128; para.wd and 5*para.wd are the two window sizes for two potential baselines. The final estimated baseline is based on these two potetial baselines
para.pr = 25;    % moving windowsize for pruning, 25 is an empirical value
para.hq = 0.4;   % Default 0.4/ flowcell with ethanol 0.6/HQI ranges from -1 to 1; may be 0.0:0.05:0.8 
para.dn = 0;     % direvative order: para.dn = 0 no; para.dn = 1 1st order; para.dn = 2 2nd order; 1st order seemed to be the best (Renner et al., 2019)
para.di = 0.00001;   % co-polymer shreshold, default 0.1

% Options -  for validation and checking
para.vd  = 0;    % for model validation: yes = 1, no = 0;
para.ck  = 0;    % the i-th sample spectra for checking baseline correaction and CI % 123
para.ci  = 0;    % calculate CI: yes = 1, no = 0;
para.tp  = 1;    % show seperated identified plots of each plastic type for checking
para.it  = 'nearest'; % interpolation method 'nearest','v4','cubic'

% Fixed parameters
SpectrumType  = "Raman"; % Raman, ATR, or μFTIR in future. This automatically determines para.ct
para.X1  = [1650 1850]; % CI range, based on Almond et al., 2020
para.X2  = [1420 1500]; %CI range, based on Almond et al., 2020
para.nr  = 0;           % for noise removal 1 = yes, 0 = no; * no noise removal seemed better
para.dv  = 0.5;         % default 0.5 alignment interval; accuracy of spectrum 1D-interpolation
para.sp  = 11;          % smoothing range, Nakano et al., 2021 para.sp = 11
para.qt  = 0;           % Quantile for localized baseline; Lenz et al., 2015 para.qt = 0
para.rg  = [400 4000];  % range of measured spectrum, Nakano et al., 2021 para.sp = [400 4000]
para.sm  = 1;           % para.sm = 0 for moving average (Nakano et al., 2021); para.sm = 1 for Savitzky-Golay (Lenz et al., 2015)
para.od  = 2;           % polynomial order for Savitzky-Golay
if SpectrumType == "Raman"
    para.ct  = [3500 4000; 2000 2500; 400 500];          % For Raman [3500 4000; 2000 2500; 400 500]                                    
elseif SpectrumType == "ATR"
    para.ct  = [650 700; 2250 2450; 000 500; 3500 8000]; % For ATR [650 700; 2250 2450; 000 500; 3500 8000], Nakano et al., 2021 
elseif SpectrumType == "μFTIR"
    para.ct  = [400 700; 1000 1350; 2250 2450; 3300 5000]; % Nakano et al., 2021 
end

%% Extract information of  coordinate(xy)

% For coordinates (CD = coordinates related)
CD.xy0 = (T.smp.Properties.VariableNames(2:end))';    % coordinate into column

CD.xy1 = split(CD.xy0,["_" "X" "Y"]);                 % keep numbers of x and y
CD.xy1(:,all(ismissing(CD.xy1)))=[];                  % remove empty colums (from "_", "X", "Y")

xy = cellfun(@str2num,CD.xy1);                        % coordinates of each spectrum
x  = xy(:,1);
y  = xy(:,2);
wd = max(x) - min(x);
ht = max(y) - min(y);
% dx = unique(abs(diff(x)));
% dx = dx(2);
% dy = unique(abs(diff(y)));
% dy = dy(1);

% %% Spectral preprocessing (I) - baseline correction
% % DSW method for baseline correction
% [SNR.smp,XY.pre1smp,~,~]  = Function_DSW_f04(table2array(T.smp),para.wd,0); % bcnr = baseline corrected & noise removed
% [SNR.std,XY.pre1std,~,~]  = Function_DSW_f04(table2array(T.std),para.wd,0); % bc = baseline corrected only
% 
% % Results for spectral preprocessing (I)
% T.pre1smp  = array2table(XY.pre1smp,"VariableNames",['X',ID.smp]);
% T.pre1std  = array2table(XY.pre1std,"VariableNames",['X',ID.std]);
% 
% T.snrsmp   = array2table(SNR.smp(3,:)',"RowNames",ID.smp);  % snr = signal-2-noise ratio
% T.snrstd   = array2table(SNR.std(3,:)',"RowNames",ID.std);
% 
% T.sigmasmp = array2table(SNR.smp(2,:)',"RowNames",ID.smp);   % sigma = std of noise (assume noise is normal)
% T.sigmastd = array2table(SNR.std(2,:)',"RowNames",ID.std);


%% Spectral proprocessing (II) - pruning
% Pruning and differential after pruning
% [XY.pre2smpd0,XY.pre2smpd1,XY.pre2smpd2] = Function_Pruning_f02(table2array(T.pre1smp),para.dv,para.sp,0,para.qt,para.pr,para.pr,para.rg,para.ct,para.sm,para.od);
% [XY.pre2stdd0,XY.pre2stdd1,XY.pre2stdd2] = Function_Pruning_f02(table2array(T.pre1std),para.dv,para.sp,0,para.qt,para.pr,para.pr,para.rg,para.ct,para.sm,para.od);

[XY.pre2smpd0,XY.pre2smpd1,XY.pre2smpd2] = Function_Pruning_f02(table2array(T.smp),para.dv,para.sp,0,para.qt,para.pr,para.pr,para.rg,para.ct,para.sm,para.od);
[XY.pre2stdd0,XY.pre2stdd1,XY.pre2stdd2] = Function_Pruning_f02(table2array(T.std),para.dv,para.sp,0,para.qt,para.pr,para.pr,para.rg,para.ct,para.sm,para.od);
 
% Results for spectral preprocessing (II)
T.pre2smpd0 = array2table(XY.pre2smpd0,"VariableNames",['X',ID.smp]); % pre2 = preprocessing (II)
T.pre2smpd1 = array2table(XY.pre2smpd1,"VariableNames",['X',ID.smp]); % di = the i-th order
T.pre2smpd2 = array2table(XY.pre2smpd2,"VariableNames",['X',ID.smp]);

T.pre2stdd0 = array2table(XY.pre2stdd0,"VariableNames",['X',ID.std]);
T.pre2stdd1 = array2table(XY.pre2stdd1,"VariableNames",['X',ID.std]);
T.pre2stdd2 = array2table(XY.pre2stdd2,"VariableNames",['X',ID.std]);

T.pre2smp = T.pre2smpd0;
T.pre2std = T.pre2stdd0;

%% Identification
% Based on correlation coeffecient (Renner et al., 2019)
[T.HQI,T.Rsmp,T.copolymer] = Function_Identification_f01(T.pre2smp,T.pre2std,para.hq,para.di); % HQI; R-matrix, Identification results

%% Plots for checking
% Input data
if para.ck > 0
    Function_ChekingPlot(T.smp,T.std,T.pre1smp,T.pre1std,T.pre2smp,T.pre2std,T.Rsmp,para.ct,para.rg,para.ck)%para.ck)
    %para.ck)
end

%% Plots - HQI (plot the first 8 types)
clear zc

R.std  =  T.Rsmp{:,:};

R.label = T.copolymer(:,2);         % predicted types of plastics
% R.types = ["PA" "PLC" "PC""PE" "PET" "PP" "PS" "PVC" "unknown"];
R.types = ID.std;
% R.type = unique(T.copolymer(:,2));% unique types

R.index = R.label{:,:} == R.types;         % if plastics = 1, if unkown = 0 (background)
R.Rsmp  = R.index(:,:) .* R.std(:,:);      % assign R to certain plastics

figure
for i = 1:1:8
    subplot(2,4,i)
    zi = R.Rsmp(:,i);
    z = R.std(:,i);

    if sum(zi) ~= 0

        [xc, yc] = meshgrid(min(x):para.dd:max(x),min(y):para.dd:max(y));
        zc(:,:,i) = griddata(x,y,z,xc,yc,para.it);
        [C,h] = contourf(xc,yc, zc(:,:,i),100);
        set(h,'LineColor','none')
        colorbar
    end
        axis equal
        colormap jet    
        xlabel('X(μm)')
        ylabel('Y(μm)')
        title(ID.std(i))    
end

%% Plots - Particle identification (binary plots)
clear zc 
% clc
% close all

% Major and minor axis length are saved
nmax  = length(x); % maximium possible number of particles (assume each point corresponds to 1 particle)
ALmj0 = NaN(nmax,length(ID.std));
ALmn0 = NaN(nmax,length(ID.std));

%
map = [1 1 1
       0.7 0.7 0.7];

for i = 1:1:length(ID.std)

    % z = R.Rsmp(:,i);
    z = R.index(:,i) + 0; % if plastic, 1, not plastic, 0

    if sum(z) ~= 0

        figure

        fig = figure;
        [xc, yc] = meshgrid(min(x):para.dd:max(x),min(y):para.dd:max(y));
        zc(:,:,i) = griddata(x,y,z,xc,yc,para.it);
        [M,c] = contourf(xc,yc, zc(:,:,i),1);
        set(c,'LineColor','none')
        axis equal
        colormap gray
        % xlabel('X(μm)')
        % ylabel('Y(μm)')
        axis off
        axis tight 

        % Convert figure to binary image
        exportgraphics(fig,string(append('plot of ',ID.std(i),'.png')),'Resolution',300); % Export figure, it is related to pixels and scaled space
        close
        I    = imread(string(append('plot of ',ID.std(i),'.png')));
        % imshow(I,'Border','tight')
        BW   = im2bw(I);
        
        [hti,wdi] = size(BW);
        r    = 1/2*(ht/(dy + hti) + wd/(wdi + dy)); % average ratio between pixels and real size. h/hi ≈ w/wi, but there are some small differences due to white space
                      
        % extract information from image
        Area = regionprops(BW,'area');
        Cent = regionprops(BW,'centroid');
        ALmj = regionprops(BW,'MajorAxisLength');
        ALmn = regionprops(BW,'MinorAxisLength');

        Area = cat(1, Area.Area);
        ALmj = cat(1, ALmj.MajorAxisLength) * r; % *100/15 = convertion factor
        ALmn = cat(1, ALmn.MinorAxisLength) * r;
        Cent = cat(1, Cent.Centroid);
        
        index = find(ALmj == max(ALmj));

        Area(index) = [];
        ALmj(index) = [];
        ALmn(index) = [];
        Cent(index,:) = [];

        ALmj0(1:length(ALmj),i) = ALmj;
        ALmn0(1:length(ALmn),i) = ALmn;    

        % Plot with identified particle and data
        imshow(BW,'Border','tight','Colormap',map)
        imshow(BW,map)
        hold on
        text(Cent(:,1),Cent(:,2),strcat('(',strcat(num2str(round(ALmj,0))),', ',num2str(round(ALmn,0)),')'),'Color','red')
        title([ID.std(i),'Major and minor axis length (μm), n = ',num2str(length(Area))])
        title(strcat(ID.std(i),' | Major and minor axis length (μm), n = ',num2str(length(Area))))
        xlabel('X-axis')
        ylabel('Y-axis')
        axis off
        
        % Determine if wanna check the idetified results of each plastic type
        if para.tp == 0
            close
        end

    end
end


% Axis length
idxmj  = isnan(ALmj0);
idxmjS = sum(idxmj,2);
idxmjS(idxmjS==max(idxmjS)) = []; % the length of column that needs to be kept
nKeep = length(idxmjS);
ALmj0(nKeep+1:end,:) = [];
ALmn0(nKeep+1:end,:) = [];

T.LengthMajorAxis = array2table(ALmj0,'VariableNames',ID.std);
T.LengthMinorAxis = array2table(ALmn0,'VariableNames',ID.std);


%% Color image test

% RGB.ColorNames = {'FF7C80';'FFB480';'F8F38D';'42D6A4';'08CAD1';'59ADF6';'9D94FF';'C780E8'};
RGB.colors = 255 - [255, 124, 128; 255, 180, 128; 248, 243, 141; 66, 214, 164;8, 202, 209; 89, 173, 246; 157, 148, 255; 199, 128, 232];


% Seperated color figure
for i = 1:1:length(ID.std)

    z = R.Rsmp(:,i);

    if sum(z) ~= 0

        cfig = figure;

        RGB.rgb = imread(string(append('plot of ',ID.std(i),'.png'))); % Export figure, it is related to pixels and scaled space
        
        if i <= length(RGB.colors)
            RGB.rgb(:,:,1) = RGB.rgb(:,:,1) - RGB.colors(i,1);
            RGB.rgb(:,:,2) = RGB.rgb(:,:,2) - RGB.colors(i,2);
            RGB.rgb(:,:,3) = RGB.rgb(:,:,3) - RGB.colors(i,3);
        end

        % RGB.rgb(:,:,3) = 0;
        imshow(RGB.rgb)
        exportgraphics(cfig,string(append('color plot of ',ID.std(i),'.png')),'Resolution',300); % Export figure, it is related to pixels and scaled space
        close
    end
end

% Legend
figure
for i = 1:1:(height(RGB.colors)+1)
    if i <=height(RGB.colors)
        rectangle('Position',[1 0-i 0.5 0.5],'FaceColor',  -(RGB.colors(i,:) - 255)/255)
        hold on
        text(1.7,0-i + 0.3,ID.std(i));
    else
        rectangle('Position',[1 0-i 0.5 0.5],'FaceColor',  [150 150 150]/255)
        hold on
        text(1.7,0-i + 0.3,'others');
    end
end

axis equal
axis tight
axis off


% Conbined color figure
RGB.pooled = 0;

for i = 1:1:length(ID.std)

    z = R.Rsmp(:,i);

    if sum(z) ~= 0
        RGB.ith = imread(string(append('color plot of ',ID.std(i),'.png')));
        RGB.pooled = imfuse(RGB.ith,RGB.pooled,"blend",Scaling="independent");
    end


end

RGB.layers = round(sum(max(R.Rsmp) ./ max(R.Rsmp),'omitnan')); % number of layers blended - will be used for intensity correction
RGB.pooled = RGB.layers * RGB.pooled;

% figure
% imshow(RGB.pooled)
% title(append('Identified plastics | 1000 μm = ',string(round(width(RGB.pooled)/wd*1000,2)),'pixels'))
% xlabel('X-axis (pixels)')
% ylabel('Y-axis (pixels)')
% axis on

%close all 
figure
imshow(imresize(RGB.pooled,1/(width(RGB.pooled)/wd)))
title('Identified plastics')
xlabel('X-axis (μm)')
ylabel('Y-axis (μm)')
axis on
axis equal
%axis tight

%% Get axis-length
AxisLength = [ALmj ALmn]; % only for one type

toc










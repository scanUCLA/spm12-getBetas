%=====================================================================%
% BETA EXTRACTION
%
%   This script allows the user to select a set of ROIs and a set of
%   contrast images
%
%   Created 1/23/2015 by Jared Torre
%
%=====================================================================%
clear all; home;

javaaddpath('/space/raid8/data/lieber/JTORRE/toolbox/xlwrite/jxl.jar')
javaaddpath('/space/raid8/data/lieber/JTORRE/toolbox/xlwrite/MXL.jar')
addpath(genpath('/space/raid8/data/lieber/JTORRE/toolbox/xlwrite'))

startDIR = pwd;

roifiles = spm_select([1 Inf], {'image'}, 'Select the ROI .hdr or .nii files...', [], pwd, '.(hdr|nii)', 1:99999);
confiles = spm_select([1 Inf], {'image'}, 'Select the contrast/beta .hdr or .nii files...', [], pwd, '.(hdr|nii)', 1:99999);

roifiles = cellstr(roifiles);
confiles = cellstr(confiles);

[roipath roiname roi_e] = cellfun(@(x) fileparts(x), roifiles, 'UniformOutput', false);
[conpath conname con_e] = cellfun(@(x) fileparts(x), confiles, 'UniformOutput', false);

nroi = length(roifiles);
ncon = length(confiles);

% Convert the .hdr ROIs to .nii
roi_hdrs_in_idx = strcmpi(roi_e,'.hdr');
converted_hdr2nii_roifiles = cellfun(@(x,y) fullfile(x,strcat(y,'.nii')), roipath(roi_hdrs_in_idx), roiname(roi_hdrs_in_idx), 'UniformOutput', false);
Nii = cellfun(@(x) load_nii(x), roifiles(roi_hdrs_in_idx));
for i = 1:length(Nii)
    save_nii(Nii(i),converted_hdr2nii_roifiles{i})
end
roifiles(roi_hdrs_in_idx) = converted_hdr2nii_roifiles;

% Create list of new resliced output filenames
rtmp_roifiles = cellfun(@(x,y) fullfile(x,strcat('rtmp_',y,'.nii')),roipath,roiname,'UniformOutput',false);
%cellfun(@(x,y) copyfile(x,y), roifiles, rtmp_roifiles)

% Reslice ROI .nii images into sample contrast space
%==================================================%
% flags = struct('interp', 1, ... % b-spline
%     'mask', 0, ...              % do not mask
%     'mean', 0, ...              % do not write mean image
%     'hold', -1, ...             % i don't think this is used anymore
%     'which', 1, ...             % reslice 2nd-nth only
%     'wrap', [0 0 0]', ...        % the default; don't know what this is
%     'prefix', 'rtmp_' ...
%     );
% P = [confiles(1);rtmp_roifiles];
% spm_reslice(P,flags)

% read in reference image
RefHead = spm_vol(confiles{1});
RefData = spm_read_vols(RefHead);
mat=RefHead.mat;
dim=RefHead.dim;

for r = 1:length(roifiles)
    % read in image to reslice
    SourceHead = spm_vol(roifiles{r});
    SourceData = spm_read_vols(SourceHead);
    
    % do the reslicing
    [x1,x2,x3] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
    d           = [1*[1 1 1]' [1 1 0]'];
    C = spm_bsplinc(SourceHead, d);
    v = zeros(dim);
    M = inv(SourceHead.mat)*mat; % M = inv(mat\SourceHead.mat) in spm_reslice.m
    y1   = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
    y2   = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
    y3   = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
    out    = spm_bsplins(C, y1,y2,y3, d);
    
    %Revised by YAN Chao-Gan 121214. Apply a mask from the source image: don't extend values to outside brain.
    tiny = 5e-2; % From spm_vol_utils.c
    Mask = true(size(y1));
    Mask = Mask & (y1 >= (1-tiny) & y1 <= (SourceHead.dim(1)+tiny));
    Mask = Mask & (y2 >= (1-tiny) & y2 <= (SourceHead.dim(2)+tiny));
    Mask = Mask & (y3 >= (1-tiny) & y3 <= (SourceHead.dim(3)+tiny));
    
    out(~Mask) = 0;
    outmat = mat;
    
    OutHead=SourceHead;
    OutHead.mat      = mat;
    OutHead.dim(1:3) = dim;
%    [p n e] = fileparts(SourceHead.fname);
%    newname = sprintf('%s_%dx%dx%d%s',n,dim,e);
%    OutHead.fname = [p filesep newname];
    OutHead.fname = rtmp_roifiles{r};
    spm_write_vol(OutHead,out);
end

% Set up data matrix
data = cell(length(confiles)+1,length(roifiles)+2);
data(:,1) = ['Full File'; confiles];
data(:,2) = ['Con Name'; conname];
data(1,3:end) = roiname';

% Extract the parameters
for r = 1:nroi
    roi = rtmp_roifiles{r};
    hdr = spm_vol(roi); img = spm_read_vols(hdr);
    roiIDX{r} = find(img);
    
    for c = 1:ncon
        conimg = confiles{c};
        hdr = spm_vol(conimg); img = spm_read_vols(hdr);
        data{c+1,r+2} = nanmean(img(roiIDX{r}));
    end
end

a=datestr(clock,31);   % returns date string of the form 'YYYY-MM-DD HH:MM:SS' e.g., 2006-12-27 15:03:37
time_stamp = [a(6:7) a(9:10) a(3:4) '_' a(12:13) a(15:16)];   % timestamp is a function name, hence the _ in time_stamp
xlwrite(strcat('roidata_',time_stamp,'.xls'),data);

% Cleanup the resliced and converted ROI files
cellfun(@(x) delete(x), [converted_hdr2nii_roifiles;rtmp_roifiles])

%This script is for .png 3D images (e.g. MRI Scans)
%Assumptions: The input image is assumed to have no background noise.
%Output #1: EC is a cell, EC{i} is the *integrated* EC curve for ith direction.
%Output #2: Shapes is the resulting structual array where:
%.Name = Name of the shape and the slice number (e.g. apple_1)
%.EC = A matrix of EC curves for that shape and slice number (#columns = #rotations)

%Load functions to compute Euler Characteristics
addpath /home/lac55/GBM_Scripts/ECFunctions/

%Define where results should be stored
resultdir = '/data/mukherjeelab/GBM/';

%Set the working directory where images are stored
maindir = '/data/mukherjeelab/GBM/TCGAClean/';
cd(maindir)

%Load each of the images in that directory (e.g. apple as example)
TCGA_patients = dir();

%Set up the desired number of rotations
rotstep=72;
theta=-pi:2*pi/rotstep:pi;
d1=cos(theta);d2=sin(theta);
d=[d1;d2];

for k=4:length(TCGA_patients)
   
   %Go into each subdirectory for each patient 
   subdir = strcat(maindir,TCGA_patients(k).name,'/Segmentations'); 
   cd(subdir)
   slices = dir('*.png');
   EC = 0;
   realslice = 0;
   
   %Compute the ECs for each slice of the MRI
   for j=1:length(slices)
       filename=slices(j).name;
       I=imread(filename);
       if sum(sum(I))~=0
           BW=im2bw(I);
           [start_r,start_c] = find(BW,1,'first');
           Z = bwtraceboundary(BW,[start_r start_c],'W',8,Inf,'counterclockwise');
           E=[1:length(Z);[2:length(Z) 1]]';
           complex.V=Z;
           complex.E=E;
           complex.F=[];
           complex.T=[];
           Z=Z-repmat(mean(Z),length(Z),1);
           Z=Z/norm(Z,'fro');
           for i=2:length(d)
               fun=Z*d(:,i);
               C=mao_gEuler(complex,fun,100);
               ec(:,i-1) = C(:,2);
           end
           EC = EC+ec; %Integration Step #1
           realslice = realslice+1;
           clear ec;
       end
   end
       %Record the Results
       MRIs(k-3).Name = TCGA_patients(k).name;
       MRIs(k-3).EC = EC/realslice; %Integration Step #2
end

%Set the working directory where images are stored
cd(resultdir)
 
save('MRI_SECTs.mat','MRIs');

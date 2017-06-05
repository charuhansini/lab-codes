close all;
clear all;
minimumsize=35;
maximumsize=12000;
circularity=[0 0.8];
cutoff=100;
thresh=0.5;
radiusrange=[900 1400];
structuralelement = strel('disk',2);
smooth=0.1;

%RGB threshold values for face/phone region
rmin1=130;
rmax1=201;
gmin1=100;
gmax1=186;
bmin1=55;
bmax1=125;
%RGB threshold values for non- face/phone region
rmin2=180;
rmax2=256;
gmin2=165;
gmax2=256;
bmin2=80;
bmax2=210;

junkratio=1.1;


% r=a(:,:,1);
% g=a(:,:,2);
% b=a(:,:,3);
% figure(2);
% imagesc(r);
% figure(3);
% imagesc(g);
% figure(4);
% imagesc(b);

% figure(1);
% imagesc(a);
% [gx, gy]=gradient(pic);
% figure(5);
% imagesc(gx);
% figure(6);
% imagesc(gy);
% figure(7);
% imagesc(divergence(gx, gy));
% mag=sqrt((gx.*gx)+(gy.*gy));
% figure(8);
% imagesc(mag);
% imwrite(uint16(add),'sum.tif','tiff','WriteMode','append','Compression','none');


%meshgrid

imageset=dir('*.jpg');
% imageset=dir('*.jpg');
for i=25:30
%     i=1:length(imageset)
    a = imread(imageset(i).name);
    add = sum(a,3);
    pic=double(add);
    [sizex,sizey,sizez]=size(a);
    siz=[sizex,sizey];
    %smoothen using gaussian filter
%     a = imgaussfilt(b,smooth);
    
    % find objects of colour with specified rgb values WITH FACE
    
        binaryim1=zeros(siz);
    binaryimr1=zeros(siz);
    binaryimg1=zeros(siz);
    binaryimb1=zeros(siz);
    positionr1=intersect(find(a(:,:,1)>rmin1), find(a(:,:,1)<rmax1));
    positiong1=intersect(find(a(:,:,2)>gmin1), find(a(:,:,2)<gmax1));
    positionb1=intersect(find(a(:,:,3)>bmin1), find(a(:,:,3)<bmax1));
    binaryimr1(positionr1)=1;
    binaryimg1(positiong1)=1;
    binaryimb1(positionb1)=1;
    addition1=binaryimr1+binaryimg1+binaryimb1;
    final1=find(addition1==3);
    binaryim1(final1)=1;
        
    binaryim1=logical(binaryim1);
    filtered1=binaryim1;
        filtered1 = bwpropfilt(filtered1, 'Area', [minimumsize maximumsize]);
 %close holes
    filtered1 = imclose(filtered1,structuralelement);
    
    filtered1 = bwpropfilt(filtered1, 'Eccentricity', circularity);
    labeld1=bwlabel(filtered1);
    
    stats1=regionprops(labeld1,'all');
    notjunk1=zeros(size(labeld1));
    for m=1:max(labeld1(:))
        
        if ((stats1(m).ConvexArea)/(stats1(m).Area))<junkratio
            notjunkpixels1=find(labeld1==m);
            notjunk1(notjunkpixels1)=1;
        end
    end
    
    %end of 1
    
     
    % find objects of colour with specified rgb values WITHOUT FACE
    binaryim2=zeros(siz);
    binaryimr2=zeros(siz);
    binaryimg2=zeros(siz);
    binaryimb2=zeros(siz);
    positionr2=intersect(find(a(:,:,1)>rmin2), find(a(:,:,1)<rmax2));
    positiong2=intersect(find(a(:,:,2)>gmin2), find(a(:,:,2)<gmax2));
    positionb2=intersect(find(a(:,:,3)>bmin2), find(a(:,:,3)<bmax2));
    binaryimr2(positionr2)=1;
    binaryimg2(positiong2)=1;
    binaryimb2(positionb2)=1;
    addition2=binaryimr2+binaryimg2+binaryimb2;
    final2=find(addition2==3);
    binaryim2(final2)=1;
    binaryim=zeros(siz);
    
    binaryim=binaryim1+binaryim2;
    
    binaryim2=logical(binaryim2);
    filtered2=binaryim2;
        filtered2 = bwpropfilt(filtered2, 'Area', [minimumsize maximumsize]);
 %close holes
    filtered2 = imclose(filtered2,structuralelement);
    
    filtered2 = bwpropfilt(filtered2, 'Eccentricity', circularity);
    labeld2=bwlabel(filtered2);
    
    stats2=regionprops(labeld2,'all');
    notjunk2=zeros(size(labeld2));
    for m=1:max(labeld2(:))
        
        if ((stats2(m).ConvexArea)/(stats2(m).Area))<junkratio
            notjunkpixels2=find(labeld2==m);
            notjunk2(notjunkpixels2)=1;
        end
    end
    
    %end of 2
    
    
    
%     %keyboard
% %     figure(2);
% %     imagesc(binaryim);
% %     title('binarised');
%     
%     
%     %      reduce=zeros(size(pic));
%     %     segmentedpixels=find(pic>cutoff);
%     %     reduce(segmentedpixels)=1;
%     %     new=reduce.*pic;
%     %     % to crop plate and threshold based on plate intensity
%     %     [centres, radii] = imfindcircles(pic, radiusrange);
%     %     radius=max(radii);
%     %     coordinate=find(radii==radius);
%     %     cx=centres(coordinate,1);
%     %     cy=centres(coordinate,2);
%     %     [X,Y] = meshgrid(1:2368, 1:4208);
%     %     pixels=find(((X-cx).^2+(Y-cy).^2)<(radius*radius));
%     %     mask=zeros(4208,2368);
%     %     mask(pixels)=1;
%     %     plate=pic.*mask;
%     %     maxi=max(plate(:));
%     %     avg1=mean(add(:));
%     %     thresholdvalue=thresh*maxi;
%     
%     % % threshold using intensity
%     %     binaryim=zeros(size(plate));
%     %     segmentedpixels=find(plate>thresholdvalue);
%     %     binaryim(segmentedpixels)=1;
%     
%     % imagesc(binaryim);
%     % binaryim=bwareaopen(binaryim,minimumsize);
%     
%     % labelledim=bwlabel(binaryim);
%     % imagesc(labelledim);
%     
%     binaryim=logical(binaryim);
%     filtered=binaryim;
%     %circles only
%     %     filtered = bwpropfilt(binaryim, 'EulerNumber', [1 1]);
%     %     keyboard
%     
%     filtered = bwpropfilt(filtered, 'Area', [minimumsize maximumsize]);
%     %     keyboard
%     %     figure(4);
%     %     imagesc(filtered);
%     %     title('binarised+area filtered');
%     
%       %close holes
%     structuralelement = strel('disk',2);
%     filtered = imclose(filtered,structuralelement);
%     
%     filtered = bwpropfilt(filtered, 'Eccentricity', circularity);
%     
%     %     figure(3);
%     %     imagesc(filtered);
%     %     title('binarised+area+eccentricity filtered');
%         
%     labeld=bwlabel(filtered);
%     
%     stats=regionprops(labeld,'all');
%     notjunk=zeros(size(labeld));
%     for m=1:max(labeld(:))
%         
%         if ((stats(m).ConvexArea)/(stats(m).Area))<junkratio
%             notjunkpixels=find(labeld==m);
%             notjunk(notjunkpixels)=1;
%         end
%     end
%     
%     
%     
     
    notjunk=notjunk1+notjunk2;
    filtered=filtered1+filtered2;
    notjunk=logical(notjunk);
    filtered=logical(filtered);
    finalpic=notjunk.*filtered;
    stat=regionprops(bwlabel(finalpic),'all');
    figure(i);
    subplot(2,2,1);
    imagesc(a);
    subplot(2,2,4);
    imagesc(filtered);
    title('filtered'); 
    
   
    subplot(2,2,2);
    imagesc(finalpic);
    title('final');
    subplot(2,2,3);
    imagesc(binaryim);
    title('binarised');
%     saveas(gcf,['\media\charu\alpha\Data\evolution of freeze thaw tolerance\colony counting\Image analysis\analysed\mask_' imageset(i).name]);
%     close
end


% a = imread('image.jpg');
% add = sum(a,3);
% pic=double(add);
% binaryim=zeros(size(pic));
% segmentedpixels=find(pic>thresholdvalue);
% binaryim(segmentedpixels)=1;
% % imagesc(binaryim);
% % binaryim=bwareaopen(binaryim,minimumsize);
%
% % labelledim=bwlabel(binaryim);
% % imagesc(labelledim);
% binaryim=logical(binaryim);
%
% filtered = bwpropfilt(binaryim, 'EulerNumber', [1 1]);
% filtered = bwpropfilt(filtered, 'Area', [minimumsize maximumsize]);

% %watershed copied from documentation

% D = bwdist(~binaryim);
% figure(3);
% imshow(D,[],'InitialMagnification','fit')
% title('Distance transform of ~pic')
%
% D = -D;
% D(~binaryim) = -Inf;
%
% L = watershed(D);
% rgb = label2rgb(L,'jet',[.5 .5 .5]);
% figure(4);
% imshow(rgb,'InitialMagnification','fit')
% title('Watershed transform of D')
% %end of watershed
%
% filtered = bwpropfilt(filtered, 'Eccentricity', circularity);
% stats=regionprops(filtered,'all');
% figure(1);
% subplot(1,2,1);
% imagesc(a);
% subplot(1,2,2);
% imagesc(filtered);
%
% [CENTERS, RADII, METRIC] = imfindcircles(filtered, [4 50]);
%
% imagesc(a)
% hold on;
% plot(CENTERS(:,1),CENTERS(:,2),'o');
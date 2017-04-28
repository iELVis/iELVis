function elec = ntools_elec_savebin(loc,hdr,savename)
% save the depth electrodes into a binary image volume
% usage: ntools_elec_savebin(location_matrix, head_info, output_volume)
% 
% location_matrix: N*3 matrix contains the electrodes locations in ras
% head_info: head information from input volume using load_nifti
% output_volume: output image name

if isempty(loc)
    disp('No location is saved');
    return;
end

fprintf('Saving electrodes in T1 space volume...');
tic;
tm_vox2ras = hdr.vox2ras;
elec_ras = [loc ones(length(loc),1)]';
elec_vox = tm_vox2ras\elec_ras;
elec = elec_vox(1:3,:)';


volume = zeros(hdr.dim(2),hdr.dim(3),hdr.dim(4));
% volume = zeros(hdr.volsize);

elec = floor(elec)+1;
for i=1:length(elec)
    [X, Y, Z] = meshgrid(elec(i,1)-1:elec(i,1)+1,elec(i,2)-1:elec(i,2)+1,elec(i,3)-1:elec(i,3)+1);
    count = 1;
    for k=1:3
        for j=1:3
            for l=1:3
                location1(count,1)=X(l,j,k);
                location1(count,2)=Y(l,j,k);
                location1(count,3)=Z(l,j,k);
                count = count+1;
            end
        end
    end
    for m=1:27
        volume(location1(m,1),location1(m,2),location1(m,3)) = i;
    end
end

hdr.vol = volume;
ntools_elec_save_nifti(hdr,savename);

fprintf('Done. (%f seconds) \n\n', toc);

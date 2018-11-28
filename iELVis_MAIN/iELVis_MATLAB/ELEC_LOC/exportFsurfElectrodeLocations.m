function exportFsurfElectrodeLocations(sub,elecReconPath,elecStems,elecNums,elecType,elecHem,VOX2RAS,ctRAS,leptoRAS,pialRAS);

nElec=length(elecStems);

%%%%%% Output Electrode Names to Text Files %%%%%%%%%
fnameLabels=fullfile(elecReconPath,[sub '.electrodeNames']);
fprintf('Saving electrode labels to: %s\n',fnameLabels);
fidLabels=fopen(fnameLabels,'w');
fprintf(fidLabels,'%s\n',datestr(now));
fprintf(fidLabels,'Name, Depth/Strip/Grid, Hem\n');
for a=1:nElec,
    if elecHem(a)
        hem='L';
    else
        hem='R';
    end
    fprintf(fidLabels,'%s%d %s %s\n',elecStems{a},elecNums(a),elecType{a},upper(hem(1)));
end
fclose(fidLabels);

%%%%%% Output RAS Coordinates to Text Files %%%%%%%%%
% CT RAS COORDINATES
fnameCtRAS = fullfile(elecReconPath,[ sub '.CT']);
fprintf('Saving CT RAS electrode locations to: %s\n',fnameCtRAS);
fidCt=fopen(fnameCtRAS,'w');
fprintf(fidCt,'%s\n',datestr(now));
fprintf(fidCt,'R A S\n');
for a=1:nElec,
    fprintf(fidCt,'%f %f %f\n',ctRAS(a,1),ctRAS(a,2),ctRAS(a,3));
end
fclose(fidCt);

% Lepto RAS COORDINATES
fnameLeptoRAS = fullfile(elecReconPath,[sub '.LEPTO']);
fprintf('Saving Lepto RAS electrode locations to: %s\n',fnameLeptoRAS);
fidLepto=fopen(fnameLeptoRAS,'w');
fprintf(fidLepto,'%s\n',datestr(now));
fprintf(fidLepto,'R A S\n');
for a=1:nElec,
    fprintf(fidLepto,'%f %f %f\n',leptoRAS(a,1),leptoRAS(a,2),leptoRAS(a,3));
end
fclose(fidLepto);

% Pial RAS COORDINATES
fnamePialRAS = fullfile(elecReconPath,[sub '.PIAL']);
fprintf('Saving Pial RAS electrode locations to: %s\n',fnamePialRAS);
fidPial=fopen(fnamePialRAS,'w');
fprintf(fidPial,'%s\n',datestr(now));
fprintf(fidPial,'R A S\n');
for a=1:nElec,
    fprintf(fidPial,'%f %f %f\n',pialRAS(a,1),pialRAS(a,2),pialRAS(a,3));
end
fclose(fidPial);

%%%%%% Output VOX Coordinates to Text Files %%%%%%%%%
% Lepto VOX COORDINATES
RAS2VOX=inv(VOX2RAS);
leptoVOX=(RAS2VOX*[leptoRAS'; ones(1, nElec)])';
fnameLeptoVOX = fullfile(elecReconPath,[sub '.LEPTOVOX']);
fprintf('Saving lepto VOX electrode locations to: %s\n',fnameLeptoVOX);
fidLeptoVox=fopen(fnameLeptoVOX,'w');
fprintf(fidLeptoVox,'%s\n',datestr(now));
fprintf(fidLeptoVox,'X Y Z\n');
for a=1:nElec,
    fprintf(fidLeptoVox,'%f %f %f\n',leptoVOX(a,1),leptoVOX(a,2),leptoVOX(a,3));
end
fclose(fidLeptoVox);

% Pial VOX COORDINATES
pialVOX=(RAS2VOX*[pialRAS'; ones(1, nElec)])';
fnamePialVOX = fullfile(elecReconPath,[sub '.PIALVOX']);
fprintf('Saving pial VOX electrode locations to: %s\n',fnamePialVOX);
fidPialVox=fopen(fnamePialVOX,'w');
fprintf(fidPialVox,'%s\n',datestr(now));
fprintf(fidPialVox,'X Y Z\n');
for a=1:nElec,
    fprintf(fidPialVox,'%f %f %f\n',pialVOX(a,1),pialVOX(a,2),pialVOX(a,3));
end
fclose(fidPialVox);

%% Created text file of Inflated Pial Surface Coordinates (relies on just created text files) 
infRAS=pial2InfBrain(sub,[]);
fnameInfRAS = fullfile(elecReconPath,[sub '.INF']);
fprintf('Saving inflated pial RAS electrode locations to: %s\n',fnameInfRAS);
fidInf=fopen(fnameInfRAS,'w');
fprintf(fidInf,'%s\n',datestr(now));
fprintf(fidInf,'R A S\n');
for a=1:nElec,
    fprintf(fidInf,'%f %f %f\n',infRAS(a,1),infRAS(a,2),infRAS(a,3));
end
fclose(fidInf);

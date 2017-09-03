% This script tests Manu's PTD index function using an SEEG patient

%% Compute PTD index
PTD_idx = getPtdIndex('TWH077');

%% Plot and output PTD index for each contact
figure(1); clf;
plot(1:length(PTD_idx.elec),PTD_idx.PTD,'-o'); hold on;
for a=1:length(PTD_idx.elec),
    fprintf('%d: %s %f\n',a,PTD_idx.elec{a},PTD_idx.PTD(a));
    h=plot(a,PTD_idx.PTD(a),'o');
    clickText(h,PTD_idx.elec{a});
end
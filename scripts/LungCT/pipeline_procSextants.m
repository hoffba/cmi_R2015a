function res = pipeline_procSextants

sv_path = 'Y:\WLabaki_LVRS_2';
opts = load(fullfile(sv_path,'Pipeline_opts.mat'));
D = [dir(fullfile(sv_path,'BLVR0*'));dir(fullfile(sv_path,'WL0*'))];
N = numel(D);

res = [];
for i = 1:N
    fprintf('%d/%d : %s\n',i,N,D(i).name);
    procdir = fullfile(D(i).folder,D(i).name);
    tres = CTlung_Pipeline_sub_sextants(procdir,opts);
    res = addResultToTable(res,tres);
end
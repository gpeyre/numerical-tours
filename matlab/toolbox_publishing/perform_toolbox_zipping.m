function perform_toolbox_zipping(outdir)

if nargin<1
    outdir='../html/';
end


% zip all toolbox files
system('rm toolbox_*.zip');
system(['rm ' outdir 'toolbox_*.zip']);
a = dir('toolbox_*');
for i=1:length(a)
    name = a(i).name;
    disp(['---> Zipping ' name ' ...']);
    if exist(name)==7
        zip([name '.zip'], [name '/']);
    end
    system(['mv ' name '.zip ' outdir ]);
end
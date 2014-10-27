function perform_scilab_conversion(name, outdir, toolbox_dir)

% perform_scilab_conversion - convert a matlab file to scilab
%
%   perform_scilab_conversion(name);
%
%   If name is empty, process all the files.
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin<2
    outdir = '../scilab/';
end
if nargin<3
    toolbox_dir = '../matlab/';
end


if nargin<1 || isempty(name)
    % back processing
    list_ext = {'coding_' 'introduction_' 'image_' 'audio_' 'wavelet_' 'sparsity_' 'cs_' ...
                'denoising_' 'inverse_' 'graphics_' 'multidim_' 'meshproc_' 'meshdeform_' ...
                'meshwav_' 'variational_' 'fastmarching_'};
    a = dir('*_*.m');
    for i=1:length(a)
        name = a(i).name;
        for k=1:length(list_ext)
            if not(isempty( findstr(name, list_ext{k}) ))
                disp(['---> Translating ' name ' ...']);
                perform_scilab_conversion(name(1:end-2));
            end
        end
    end
    return;
end

fid     = fopen([name '.m'], 'rt');
fidout  = fopen([outdir name '.sce'], 'wt');

if fid<0
    error(['Cannot open ' name '.m.']);
end


str1 = {'%' '\' 'getd(''' };
str2 = {'//' '\\'  ['getd(''' toolbox_dir] };

while true
    s = fgets(fid);
    if s<0
        break;
    end
    % remove comments and stufs
    for i=1:length(str1)
        s = strrep(s,str1{i},str2{i});
    end
    makeoutput = 1;
    if length(s)>5 && strcmp(s(1:6), 'getd =')
        makeoutput = 0;
    end
    if length(s)>28 && strcmp(s(1:28), 'perform_toolbox_installation')
        process_header(fidout, s);
        makeoutput = 0;
    end
    % write output
    if makeoutput
        fprintf(fidout, s);
    end
end


fclose(fid);
fclose(fidout);


%%
function process_header(fidout, s)

tbx = {};
if findstr(s, 'signal')
    tbx{end+1} = 'signal';
end
if findstr(s, 'general')
    tbx{end+1} = 'general';
end
if findstr(s, 'graph')
    tbx{end+1} = 'graph';
end
if findstr(s, 'wavelet_meshes')
    tbx{end+1} = 'wavelet_meshes';
end
if findstr(s, 'additional')
    tbx{end+1} = 'additional';
end
for i=1:length(tbx)
    output_line(fidout, ['getd(''toolbox_' tbx{i} '/'');\n']);
end


%%
function output_line(fid,s)

s = strrep(s,'%','%%');
fprintf(fid, s);

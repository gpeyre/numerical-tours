function perform_publishing(name, options)

% perform_publishing - publish a file to HTML format
%
%    perform_publishing(name, options);
%
%   If name is empty, process all the files (and also zip all toolboxes).
%
%   options.rep set output directory (default '../../../numerical-tours-site/matlab/')
%   options.repprivate set output directory for exercices (default '../solutions/')
%   options.stylesheet set XSL stylesheet path.
%   options.format set output format (default 'html').
%
%   To put part of code that are not taken into account during publishing,
%   put them as
%       %CMT
%       ...
%       %CMT
%
%   To put part of code that are release as exercises, put them as
%       %EXO
%       ...
%       %EXO
%
%   Copyright (c) 2008 Gabriel Peyre

path('../toolbox_general/', path);

options.null = 0;

if nargin<1 || isempty(name)
    % batch processing
%    list_ext = {'coding_' 'introduction_' 'image_' 'audio_' 'wavelet_' 'sparsity_' ...
%                'denoising_' 'inverse_' 'graphics_' 'multidim_' 'mesh_' 'variational_' 'fastmarching_'};
    list_ext = {'audio' 'coding' 'cs' 'denoisingadv' 'denoisingsimp' 'denoisingwav' 'fastmarching' ...
                'graphics' 'introduction' 'inverse' ...
                'meshdeform' 'meshproc' 'meshwav' 'multidim' 'numerics' ...
                'optim' 'optimaltransp' 'segmentation' 'shapes' ... 
                'sparsity' 'wavelet'};
    a = dir('*_*.m');
    for i=1:length(a)
        name = a(i).name;
        for k=1:length(list_ext)
            if not(isempty( findstr(name, list_ext{k}) ))
                disp(['---> Publishing ' name ' ...']);
                perform_publishing(name(1:end-2));
            end
        end
    end
    % zip all toolbox files
    perform_toolbox_zipping();
    return;
end

if iscell(name)
    for i=1:length(name)
        disp(['---> Publishing ' name{i} ' ...']);
        perform_publishing(name{i}, options);
    end
    return;
end

if not(isempty(strfind(name, '*')))
    a = dir([name '.m']);
    name = {};
    for i=1:length(a)
        name{end+1} = a(i).name(1:end-2);
    end
    perform_publishing(name, options);
    return;
end

%% set up path
% directory where the wavelet-tour web site is
repweb = getoptions(options, 'rep', '../../../numerical-tours-site/matlab/');
% directory where the specific publishing is made
rep = [repweb name '/'];
if not(exist(rep))
    mkdir(rep);
end
% directory for the exercices to be stored
repprivate = getoptions(options, 'repprivate', '../solutions/');
reppriv = [ repprivate name '/'];
if not(exist(reppriv))
    mkdir(reppriv);
end
% open files
fid     = fopen([name '.m'], 'rt');
name_out = 'index';
fidout  = fopen([repweb name_out '.m'], 'wt');
if fid<0
    error(['Cannot open ' name '.m.']);
end
% copy style files
if 0
    copyfile([repweb 'style.css '],[rep 'style.css']);
else
    str = ['cp ' repweb 'style.css ' rep 'style.css'];
%    system( str );
end

%% First pre-process the file for publishing

exo_num = 0;

% output_line(fidout, '% special publishing token\npublishing_time = 1;\n\n');

while true
    s = fgets(fid);
    if not(isstr(s))
        break;
    end
    
    if exo_mode(s)
        % SPECIAL EXERCICE MODE
        exo_num = exo_num + 1;
        % Create a specific   
        process_exercice(fid, fidout, exo_num, repweb, reppriv, name);
    elseif cmt_mode(s)
        % SPECIAL CMT MODE
        process_comment(fid);
    elseif header_mode(s)
        % SPECIAL header more
        process_header(fidout, s);
    else
        s = strrep(s,'\','\\');
        output_line(fidout, s);
    end
end


fclose(fid);
fclose(fidout);

%% do the publishing
opts.format = getoptions(options, 'format', 'html');
if not(strcmp(opts.format, 'latex'))
    opts.stylesheet = getoptions(options, 'stylesheet', [repweb 'nt.xsl']);
    opts.stylesheet = [pwd '/' opts.stylesheet];
end
%filename = [repweb name_out '.m'];
opts.outputDir = [pwd '/' rep];
% opts.format = 'xml';

path(repweb, path);
file = publish([name_out '.m'],opts);

% web(file); % uncomment for display
delete([repweb 'exo*.m']);
delete([repweb name_out '.m']);
% make the name of the output file index.html
% movefile([opts.outputDir '/' name_out '.html'], [opts.outputDir '/index.html']);

%% do the online publishing 
% disp('Performing online publishing (might take some time) ...');
% perform_online_publishing(name);

%% do the python conversion
if 0
    str = ['python ../../scripts/m2nb_converter.py ./' name '.m ../'];
    system(str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  AUX FUNCTIONS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function output_line(fid,s)

s = strrep(s,'%','%%');
fprintf(fid, s);

%%
function r = cmt_mode(s)
r = (length(s)>4 && strcmp(s(1:4),'%CMT'));

%%
function r = exo_mode(s)
r = (length(s)>4 && strcmp(s(1:4),'%EXO'));

%%
function r = header_mode(s)
r = (length(s)>28 && strcmp(s(1:28),'perform_toolbox_installation'));


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


output_line(fidout, '%% Installing toolboxes and setting up the path.\n\n');

% output_line(fidout, '%%\n%<? include ''../nt.sty''; ?>\n\n');

output_line(fidout, '%%\n');
output_line(fidout, '% You need to download the following files: \n'); 
for i=1:length(tbx)
    output_line(fidout, ['% <../toolbox_' tbx{i} '.zip ' tbx{i} ' toolbox>']);
    if i<length(tbx)-1
        output_line(fidout, ', \n');
    elseif i==length(tbx)-1
        output_line(fidout, ' and \n');
    else
        output_line(fidout, '.\n');
    end
end

output_line(fidout, '\n%%\n');
output_line(fidout, '% You need to unzip these toolboxes in your working directory, so\n');
output_line(fidout, '% that you have \n');
for i=1:length(tbx)
    output_line(fidout, ['% |toolbox_' tbx{i} '|']);
    if i<length(tbx)-1
        output_line(fidout, ', \n');
    elseif i==length(tbx)-1
        output_line(fidout, ' and \n');
    else
        output_line(fidout, '\n');
    end
end
output_line(fidout, '% in your directory.\n\n');

output_line(fidout, '%%\n');
output_line(fidout, '% *For Scilab user:* you must replace the Matlab comment ''%'' by its Scilab\n');
output_line(fidout, '% counterpart ''//''.\n\n');

output_line(fidout, '%%\n');
output_line(fidout, '% *Recommandation:* You should create a text file named for instance |numericaltour.sce| (in Scilab) or |numericaltour.m| (in Matlab) to write all the\n');
output_line(fidout, '% Scilab/Matlab command you want to execute. Then, simply run |exec(''numericaltour.sce'');| (in Scilab) or |numericaltour;| (in Matlab) to run the commands. \n\n');

output_line(fidout, '%%\n');
output_line(fidout, '% Execute this line only if you are using Matlab.\n\n');

output_line(fidout, 'getd = @(p)path(p,path); % scilab users must *not* execute this\n\n');

output_line(fidout, '%%\n');
output_line(fidout, '% Then you can add the toolboxes to the path.\n\n');
for i=1:length(tbx)
    output_line(fidout, ['getd(''toolbox_' tbx{i} '/'');\n']);
end



%%
function process_comment(fid)
while true
    s = fgets(fid);
    if not(isstr(s)) || cmt_mode(s)
        break;
    end
end

%%
function process_exercice(fid, fidout, exo_num, rep, reppriv, curname)
% create a new file
% The m file is created and executed in rep
% and also copied to reppriv.

name = ['exo' num2str(exo_num)];
filename = [name '.m'];
fidexo = fopen([rep filename], 'wt');

output_line(fidout, ['%%\n% _Exercice ' num2str(exo_num) ':_' ]);
output_line(fidout, [' (<../missing-exo/ check the solution>)\n']);

% process files
while true
    s = fgets(fid);
    if not(isstr(s)) || exo_mode(s)
        break;
    end
    s = strrep(s,'#',num2str(exo_num));
    s = strrep(s,'\','\\');
    if length(s)>1 && strcmp(s(1:2), '%%')
        s = s(2:end);
        output_line(fidout, s);
    else
        output_line(fidexo, s);
    end
end

output_line(fidout, ['\n' name ';\n']);

fclose(fidexo);
if 0
    copyfile([rep filename],[reppriv filename]);
else
    system(['cp ' rep filename ' ' reppriv filename]);
end
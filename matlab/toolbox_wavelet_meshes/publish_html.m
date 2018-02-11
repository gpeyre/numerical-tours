function publish_html(filename, outputDir, stylesheet)

% publish_html - publish a file to HTML format
%
%    publish_html(filename, outputDir, stylesheet);
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin<1
    filename = 'content';
end
if nargin<2
    outputDir = 'html';
end
if nargin<3
    stylesheet = [outputDir '/gpeyre.xsl'];
end

opts.outputDir = outputDir;
if not(isempty(stylesheet))
    opts.stylesheet = stylesheet;
end
file = publish(filename,opts);

web(file);

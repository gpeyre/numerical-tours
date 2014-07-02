function [model] = vol3d(varargin)
%H = VOL3D Volume render 3-D data. 
% VOL3D uses the orthogonal plane 2-D texture mapping technique for 
% volume rending 3-D data in OpenGL. Use the 'texture' option to fine 
% tune the texture mapping technique. This function is best used with
% fast OpenGL hardware.
%
% H = vol3d('CData',data) Create volume render object from input 
%                         3-D data. Use interp3 on data to increase volume
%                         rendering resolution. Returns a struct 
%                         encapsulating the pseudo-volume rendering object. 
%
% vol3d(...'Parent',axH) Specify parent axes
%
% vol3d(...,'texture','2D')  Default. Only render texture planes parallel
%                            to nearest orthogonal viewing plane. Requires
%                            doing vol3d(h) to refresh if the view is
%                            rotated (i.e. using cameratoolbar).
%
% vol3d(...,'texture','3D')  Render x,y,z texture planes simultaneously. 
%                            This avoids the need to refresh the view but 
%                            requires faster OpenGL hardware peformance.
%
% vol3d(H)  Refresh view. Updates rendering of texture planes 
%           to reduce visual aliasing when using the 'texture'='2D'
%           option.
%
% NOTES
% Use vol3dtool for editing the colormap and alphamap. 
% Adjusting these maps will allow you to explore your 3-D volume 
% data at various intensity levels. See documentation on 
% alphamap and colormap for more information.
%
% Use interp3 on input date to increase/decrease resolution of data
%
% Examples:
%
% % Visualizing fluid flow
% v = flow(50);
% h = vol3d('cdata',v,'texture','2D');
% view(3); 
% % Update view since 'texture' = '2D'
% vol3d(h);  
% alphamap('rampdown'), alphamap('decrease'), alphamap('decrease')
% 
% % Visualizing MRI data
% load mri.mat
% D = squeeze(D);
% h = vol3d('cdata',D,'texture','2D');
% view(3); 
% % Update view since 'texture' = '2D'
% vol3d(h);  
% axis tight;  daspect([1 1 .4])
% alphamap('rampup');
% alphamap(.06 .* alphamap);
%
% See also vol3dtool, alphamap, colormap, opengl, isosurface

% Copyright Joe Conti, 2004

if isstruct(varargin{1})
    model = varargin{1};
    if length(varargin) > 1
       varargin = {varargin{2:end}};
    end
else
    model = localGetDefaultModel;
end


if length(varargin)>1
  for n = 1:2:length(varargin)
    switch(lower(varargin{n}))
        case 'cdata'
            model.cdata = varargin{n+1};
        case 'parent'
            model.parent = varargin{n+1};
        case 'texture'
            model.texture = varargin{n+1};
    end
    
  end
end

if isempty(model.parent)
    model.parent = gca;
end

% choose default 3-D view
ax = model.parent;
axis(ax,'vis3d');
axis(ax,'tight');

[model] = local_draw(model);


%------------------------------------------%
function [model] = localGetDefaultModel

model.cdata = [];
model.xdata = [];
model.ydata = [];
model.zdata = [];
model.parent = [];
model.handles = [];
model.texture = '2D';

%------------------------------------------%
function [model,ax] = local_draw(model)

cdata = model.cdata; 
siz = size(cdata);

% Define [x,y,z]data
if isempty(model.xdata)
    model.xdata = [0 siz(2)];
end
if isempty(model.ydata)
    model.ydata = [0 siz(1)];
end
if isempty(model.zdata)
    model.zdata = [0 siz(3)];
end

try,
   delete(model.handles);
end

ax = model.parent;
cam_dir = camtarget(ax) - campos(ax);
[m,ind] = max(abs(cam_dir));


h = findobj(ax,'type','surface','tag','vol3d');
for n = 1:length(h)
  try,
     delete(h(n));
  end
end

is3DTexture = strcmpi(model.texture,'3D');
handle_ind = 1;

% Create z-slice
if(ind==3 || is3DTexture )    
  x = [model.xdata(1), model.xdata(2); model.xdata(1), model.xdata(2)];
  y = [model.ydata(1), model.ydata(1); model.ydata(2), model.ydata(2)];
  z = [model.zdata(1), model.zdata(1); model.zdata(1), model.zdata(1)];
  diff = model.zdata(2)-model.zdata(1);
  delta = diff/size(cdata,3);
  for n = 1:size(cdata,3)

   slice = double(squeeze(cdata(:,:,n)));
   h(handle_ind) = surface(x,y,z,'Parent',ax);
   set(h(handle_ind),'cdatamapping','scaled','facecolor','texture','cdata',slice,...
	 'edgealpha',0,'alphadata',double(slice),'facealpha','texturemap','tag','vol3d');
   z = z + delta;
   handle_ind = handle_ind + 1;
  end

end

% Create x-slice
if (ind==1 || is3DTexture ) 
  x = [model.xdata(1), model.xdata(1); model.xdata(1), model.xdata(1)];
  y = [model.ydata(1), model.ydata(1); model.ydata(2), model.ydata(2)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  diff = model.xdata(2)-model.xdata(1);
  delta = diff/size(cdata,2);
  for n = 1:size(cdata,2)

   slice = double(squeeze(cdata(:,n,:)));
   h(handle_ind) = surface(x,y,z,'Parent',ax);
   set(h(handle_ind),'cdatamapping','scaled','facecolor','texture','cdata',slice,...
	 'edgealpha',0,'alphadata',double(slice),'facealpha','texturemap','tag','vol3d');
   x = x + delta;
   handle_ind = handle_ind + 1;
  end
end

  
% Create y-slice
if (ind==2 || is3DTexture)
  x = [model.xdata(1), model.xdata(1); model.xdata(2), model.xdata(2)];
  y = [model.ydata(1), model.ydata(1); model.ydata(1), model.ydata(1)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  diff = model.ydata(2)-model.ydata(1);
  delta = diff/size(cdata,1);
  for n = 1:size(cdata,1)

   slice = double(squeeze(cdata(n,:,:)));
   h(handle_ind) = surface(x,y,z,'Parent',ax);
   set(h(handle_ind),'cdatamapping','scaled','facecolor','texture','cdata',slice,...
	 'edgealpha',0,'alphadata',double(slice),'facealpha','texturemap','tag','vol3d');
   y = y + delta;
   handle_ind = handle_ind + 1;
  end
end

model.handles = h;
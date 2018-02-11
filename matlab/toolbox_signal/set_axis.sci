function set_axis(v)

// draw axis on/off
//  set_axis(0);
//  set_axis(1);

if v==0
a = gca();
a.box = "off";
a.axes_visible = ["off","off","off"];
else
a = gca();
a.box = "on";
a.axes_visible = ["on","on","on"];
end


endfunction
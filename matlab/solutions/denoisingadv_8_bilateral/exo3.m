clf;
vlist = round([.1 .3 .6 .9]*p);
for i=1:4
    v = vlist(i);
    imageplot( F(:,:,v), ['v_i=' num2str((v-1)/(p-1),2)], 2,2,i );
end

n = 256;
freq = [[.02; 0] [0; .03] [.04; .05]];
freq = round(freq*n);
t = 0:n-1;
[Y,X] = meshgrid(t,t);
clf;
for i=1:size(freq,2)
    subplot(1, size(freq,2), i);
    imageplot( cos( 2*pi/n*( freq(1,i)*X+freq(2,i)*Y )) ); 
end

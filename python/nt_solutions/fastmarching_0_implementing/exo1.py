"""
Implement the Dijkstra algorithm by iterating these step while the
stack |I| is non empty.
Display from time to time the front that propagates.
"""
pstart = transpose(array([x0]))
[D,Dsvg,Ssvg] = perform_dijstra_fm(W, pstart, inf,'dijstr', 'sym',n*6)
clf;
for i in arange(0,4):
    subplot(2, 2, i+1)
    d = Dsvg[:,:,i]
    d[d==Inf] = 0
    imageplot(d)
    set_cmap('jet')
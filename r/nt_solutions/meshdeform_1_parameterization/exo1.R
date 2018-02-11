q=floor(p/4)
t = seq(from=0,to=q-1,by=1)/q
t1 = seq(from=0, to=(p-3*q)-1)/(p-3*q)
Z = rbind(c(t,t*0+1,1-t,t1*0),c(t*0, t, t*0+1,1-t1))
ind = c(seq(from=0,to=p-1),c(0))
ax <- list(title="", zeroline=FALSE, showline = FALSE, showticklabels = FALSE, showgrid=FALSE)


plot_ly(x=Z[1,ind+1],y=Z[2,ind+1], type="scatter", mode="lines", name="line", line=list(color="blue"), frame=FALSE)%>%add_trace(x=Z[1,ind+1], y=Z[2,ind+1], type="scatter", mode="markers", name="markers", marker=(list(color="red")))%>%layout(xaxis=ax,yaxis=ax)
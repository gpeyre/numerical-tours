if 0
    clf; imageplot(M(:,:,delta));
    title('Pick starting point');
    start_point = round( ginput(1) );
    start_point = [start_point(2); start_point(1); delta];
end
%EXO
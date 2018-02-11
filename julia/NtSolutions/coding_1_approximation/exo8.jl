fL = zeros(n, n);

for i in 1 : Base.div(n, w)
    for j in 1 : Base.div(n, w)
        seli = (i - 1)*w + 1 : i*w;
        selj = (j - 1)*w + 1 : j*w;
        fL[seli, selj] = plan_dct(f[seli, selj])*f[seli, selj];
    end
end

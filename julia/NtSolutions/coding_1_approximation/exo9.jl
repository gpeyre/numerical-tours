f1 = copy(fL);
for i in 1 : Base.div(n, w)
    for j in 1 : Base.div(n, w)
        seli = (i - 1)*w + 1 : i*w;
        selj = (j - 1)*w + 1 : j*w;
        f1[seli, selj] = plan_idct(f1[seli, selj])*f1[seli, selj];
    end
end

print("Error |f-f1|/|f| =", vecnorm(f - f1)/vecnorm(f))
remap = @(FV,FV1,gamma)exp( ( gamma*FV1 + FV-FV1 ) * (b-a) + a ) - epsilon;
gamma = .1;
FVmapped = remap(FV,FV1,gamma);
imageplot(color_recompose(FVmapped), ['\gamma=' num2str(gamma)]);

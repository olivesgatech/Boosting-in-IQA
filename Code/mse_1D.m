function [mseVal] = mse_1D(img1, img2)
imgRes=img1-img2;
mseMap=(imgRes).^2;
mseVal=sqrt(mean(mseMap));
end
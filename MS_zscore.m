function z = MS_zscore(x, mu, sigma)
%% MS_zscore: same as the build in zscore but lets you specify the mean and std



sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = (x-mu) ./ sigma0;
function obj=obj_rician(x,ydata,A,noise_var,weights)
% Objective function that takes into account the Rician noise
% characteristics
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

predicted_vals=A*x;
temp=(predicted_vals(:).*ydata(:)/noise_var);
obj=sum(weights.*(((predicted_vals(:).^2)/(2*noise_var)) - (log(besseli(0,temp,1))) -temp));% See work by Sijbers et. Al, 1999 (Reference 4 of the article)

return;
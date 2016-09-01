% loess computes a robust locally weighted regression
%
%   vi=loess(x,v,xi,span) computes a robust locally fitted regression. The
%   inputs are:
%       x - The input locations given as an N x D matrix for N input points
%           in D-dimensional space
%       v - The values for the regression at the locations given in x.
%           Vector with dimensions N x 1
%       xi - The locations where to compute the regression values. Has
%           dimensions Ni x D.
%       span - Indicates the span of the regression. It is a scalar value
%              that indicates the fraction of datapoints used for
%              the regression
%
%   vi=loess(x,v,xi,span,niter) to specify the number of robust iterations
%   niter. Default is 0.
%
%   vi=loess(x,v,xi,span,niter,order) to specify the regression order. This
%   should be 1 or 2. Default is 1.
%
%   vi=loess(x,v,xi,span,niter,order,nthreads) to specify the number of
%   threads used in the computation. This value is at least one and at most
%   the number of threads reported by the system.
%
%
%   Example:
%
%   %% Generate some input with outliers and nans
%   n_inputs=100000;
%   n_outliers=100;
%   n_nans=10;
%   x=[1-2*rand(n_inputs,1) 1-2*rand(n_inputs,1)];
%   val=sin(pi*x(:,1).*x(:,2))/0.5+x(:,1).^2-exp(x(:,1));
%   f_outl=ceil(length(x)*rand(n_outliers,1));
%   val(f_outl) = 15-30*rand(n_outliers,1);
%   f_nan=ceil(length(x)*rand(n_nans,1));
%   val(f_nan) = nan;
%
%   %% Generate output locations
%   [XI,YI]=meshgrid(linspace(-1, 1, 100), linspace(-1,1,100));
%   VI=nan(size(XI));
%
%   %% smooth
%   VI(:)=loess(x,val,[XI(:),YI(:)],0.01,2,2,12);
%
%
%   Source:
%
%   Cleveland, W.S. and Grosse, E (1991), "Computational methods for local
%   regression", Statistics and Computing, vol. 1,pp. 47-62
%
%   Cleveland, W.S. (1979), "Robust Locally Weighted Regression and
%   Smoothing Scatterplots", Journal of American Statistical Association,
%   vol. 74, no. 368 , pp. 829-836
%
%
%   Implementation:
%
%   This function is implemented in C++ using the CGAL library for spatial
%   searching (http://www.cgal.org) and the Eigen library for the
%   regressions (http://eigen.tuxfamily.org)
%
%   Author: Bart Vermeulen


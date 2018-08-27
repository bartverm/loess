  %% Generate some input with outliers and nans
  n_inputs=100000;
  n_outliers=100;
  n_nans=10;
  x=[1-2*rand(n_inputs,1) 1-2*rand(n_inputs,1)];
  val=sin(pi*x(:,1).*x(:,2))/0.5+x(:,1).^2-exp(x(:,1));
  f_outl=ceil(length(x)*rand(n_outliers,1));
  val(f_outl) = 15-30*rand(n_outliers,1);
  f_nan=ceil(length(x)*rand(n_nans,1));
  val(f_nan) = nan;

  %% Generate output locations
  [XI,YI]=meshgrid(linspace(-1, 1, 100), linspace(-1,1,100));
  VI=nan(size(XI));

  %% smooth
  VI(:)=loess(x,val,[XI(:),YI(:)],0.01,2,2,12);



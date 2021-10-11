clearvars
close all

%% Generate some input with outliers and nans
  n_inputs=10000;
  n_outliers=80;
%   n_nans=0;
  x=[1-2*rand(n_inputs,1) 1-2*rand(n_inputs,1)];
  val=sin(pi*x(:,1).*x(:,2))/0.5+x(:,1).^2-exp(x(:,1));
  f_outl=ceil(length(x)*rand(n_outliers,1));
  val(f_outl) = -8; %15-30*rand(n_outliers,1);
%   f_nan=ceil(length(x)*rand(n_nans,1));
%   val(f_nan) = nan;

  %% Generate output locations
  [XI,YI]=meshgrid(linspace(-1, 1, 100), linspace(-1,1,100));
  VI=nan(size(XI));

  %% smooth
  VI(:)=loess(x,val,[XI(:),YI(:)],0.01,3,1);
  res=val-loess(x,val,x,0.02,3,1,1);
%% plot
figure
plot3(x(:,1), x(:,2), val, 'k.','MarkerSize',1)
hold on
plot3(x(f_outl,1), x(f_outl,2), val(f_outl), 'rx')


surf(XI,YI,VI)
shading flat
lighting gouraud
camlight('headlight')

%%
figure
filt=true(size(x,1),1);
filt(f_outl)=false;
scatter(x(filt,1),x(filt,2),10,res(filt),'filled')

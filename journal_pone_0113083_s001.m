%% Example Matlab-script illustrating 
%  SHAPE MODE ANALYSIS by PRINCIPAL COMPONENT ANALYSIS
%  for a synthetic data set of sperm flagellar bending waves
%
% BM Friedrich et al., 23.09.2014

%% size of data set =======================================================
n=1024; % number of experimental observations (e.g. frames recorded)
m=41;   % number of features recorded with each observation

%% construct pseudo-data set ==============================================
% (can be ignored on first reading)
ds=1.4; % arc-length step along flagellum [um]
% simple model for principlal shape modes ...
s=1:ds:m*ds; % arc-length position for which tangent angle shall be recorded
L=m*ds; % flagellar length [um]
lambda=50; % flagellar wavelength [um]
v1=cos(2*pi*s/lambda); % prototypical wave modes
v2=sin(2*pi*s/lambda);
% simple model for phase fluctuations
T=0.030; % period of flagellar beat [s]
w0=2*pi/T; % frequency of flagellar beat [1/s]
dt=0.004; % time-step [s]
D0=2; % flagellar phase diffusion coefficient
w1=0.38; % non-isochrony parameter
phi=cumsum(w0*dt+sqrt(2*D0*dt)*randn(n,1)); % flagellar phase angle 
% simple model for amplitude fluctuations ...
A0=1; % mean flagellar bending amplitude
sigma_A=0.01; % variance of amplitude fluctuations
tau_A=0.02; % correlation time of amplitude fluctuations
dA=conv(sqrt(sigma_A*2*dt/tau_A)*randn(n,1),exp(-(0:dt:10*dt)/tau_A),'same'); % amplitude fluctuations
% var(dA)/sigma_A % test
A=A0+dA;
phi=phi-cumsum(2*w1*dA*sqrt(A0)*dt); % non-isochrony
% "ideal" shape scores
B1=A.*cos(phi);
B2=A.*sin(phi);
% assemble data set using phase and amplitude
kappa0=1/L; % mean flagellar curvature [1/um]
psi0=repmat(kappa0*s,n,1); % mean tangent angle
psi1=B1*v1+B2*v2; % oscillatory part of tangent angle 
psi=psi0+psi1; % reconstructed tangent angle data set
psi=psi+0.1*randn(n,m); % add measurement noise

%% plot kymograph of (n x m)-data matrix ==================================
figure(1), clf, hold on
pcolor(s,dt*(1:100)*1e3,psi(1:100,:)) % note: only first 100 observations shown
axis tight, shading flat
colormap jet, caxis(pi*[-1/2 1]), colorbar('YTick',pi/2*(-1:2),'YTickLabel',{'-pi/2','0','pi/2','pi'}) % set color-code
xlabel('s [um]'), ylabel('t [ms]'), title('kymograph tangent angle') % plot annotation

%% linear principal component analysis ====================================

% mean of all n observations (represented as (n x m)-matrix)
psi0=repmat(mean(psi),n,1); 

% mean-corrected (n x m)-observation matrix 
Delta=psi-psi0; 

% feature-feature covariance matrix of dimensions m x m
C=Delta'*Delta;

% standard principal component analysis
[V,D]=eig(C); % V=[V1 V2 ... Vm] ... eigenvectors, D=[d1 d2 ... dm] ... corresponding eigenvalues
[D,ind]=sort(diag(D),'descend'); % sort eigenvalues (and corresponding eigenvectors) by decreasing magnitude
V=V(:,ind);

% principal shape modes
V1=V(:,1);
V2=V(:,2);

% determine shape scores from linear least-square fit
B1=V1\Delta'; B1=B1(:);
B2=V2\Delta'; B2=B2(:);

%% show eigenvalues and shape modes =======================================

figure(1), clf
% show shape modes
subplot(1,2,1), hold on
plot(s,V1,'b.-','LineWidth',2)
plot(s,V2,'r.-','LineWidth',2)
xlabel('arclength s along flagellum [um]')
ylabel('tangent angle [rad]')
title('principal shape modes')
% show eigenvalues
subplot(1,2,2),
semilogy(D/sum(D),'k.-','MarkerSize',20), hold on
plot(1,D(1)/sum(D),'b.','MarkerSize',20), plot(2,D(2)/sum(D),'r.','MarkerSize',20) % mark data points of interest
xlabel('index shape modes')
ylabel('corresponding eigenvalue (normalized)')
title('relative variance explained by shape modes')

% show shape dynamics in reduced shape space ------------------------------

j=randi(n); % pick random observation

figure(2), clf 
% show shape dynamics in reduced shape space
subplot(2,2,[1 3]), hold on
plot(B1,B2,'k.'), plot(B1(j),B2(j),'r.','MarkerSize',20),
axis image
xlabel('B_1'), ylabel('B_2'), title('shape dynamics in reduced shape space')

subplot(2,2,2), hold on
psi_proj=psi0(1,:)'+B1(j)*V1+B2(j)*V2; % projection on reduced shape space
plot(s,psi_proj,'m','LineWidth',3) 
plot(s,psi(j,:),'k.') % original observation
xlabel('s [um]'), ylabel('\psi(s) [rad]'), title('tangent angle profile (original and reconstructed)')
ylim([-2 2])

subplot(2,2,4), hold on
plot(s,B1(j)*V1,'b','LineWidth',2) % contribution first shape mode
plot(s,B2(j)*V2,'r','LineWidth',2) % contribution second shape mode
plot(s,psi0(1,:),'k--','LineWidth',2) % mean shape
xlabel('s [um]'), ylabel('\psi(s) [rad]'), title('individual contributions from shape modes')
h=legend({'V_1','V_2','\psi_0'});
set(h,'Location','EastOutside')
ylim([-2 2])

%% reconstruct limit cycle for shape space dynamics =======================

% compute phase angle from limit cycle representation ---------------------

% The polar angle in the planar shape space defines a so-called proto-phase; 
% its mean time-derivative may vary along the limit cycle.
protophi=unwrap(atan2(B2,B1)); 

% Following [Kralemann, ..., Pikovsky, PRE, 2008], we can define a proper phase
% from the proto-phase that increases at a uniform rate along the limit cycle.
nharm=10;
Sn=nan(2*nharm+1,1);
for k=-nharm:nharm
    Sn(nharm+1+k)=mean(exp(-1i*k*protophi)); % eq. (15)
end;
phi=protophi;
for k=[-nharm:-1,1:nharm]
    phi=phi+Sn(nharm+1+k)*(exp(1i*k*protophi)-1)/(1i*k); % eq. (16)
end;
phi=unwrap(real(phi));

% plot flagellar phase
figure(3), subplot(1,3,2), hold on
plot(dt*(1:n),phi)
xlabel('time [s]')
ylabel('phase \phi [rad]')

% recontruct limit cycle by Fourier average -------------------------------
nharm=1;
phi0=pi*(0:0.01:2);
B1limit=zeros(size(phi0));
B2limit=zeros(size(phi0));
for k=[-nharm:-1,1:nharm]
    B1limit=B1limit+...
        exp(1i*k*phi0)*mean(B1.*exp(-1i*k*phi));
    B2limit=B2limit+...
        exp(1i*k*phi0)*mean(B2.*exp(-1i*k*phi));        
end;
B1limit=real(B1limit); % chop tiny imaginary part
B2limit=real(B2limit);

% plot shape limit cycle with phase information ---------------------------
figure(3), clf 
set(gcf, 'Renderer', 'zbuffer');
subplot(1,2,1), hold on
plot(B1limit,B2limit,'k','LineWidth',10)
scatter(B1,B2,100*ones(n,1),hsv2rgb([mod(phi/(2*pi),1) ones(n,2)]),'.') % color-code shape points according to phase
scatter(B1limit,B2limit,100*ones(length(phi0),1),hsv2rgb([mod(phi0'/(2*pi),1) ones(length(phi0),2)]),'.') % color-coded limit circle
axis image
xlabel('B_1'), ylabel('B_2'), title('shape dynamics in reduced shape space')

subplot(1,2,2), hold on
plot(dt*(1:n),phi)
xlabel('time [s]')
ylabel('\phi [rad]')
title('phase angle')

%% Principal component analysis and singular value decomposition ----------
% ... agrees with singular value decomposition
[Us,S,Vs]=svd(Delta);
[~,ind]=sort(diag(S'*S),'descend');
S=S(:,ind); % S'*S and D are equal
Vs=Vs(:,ind); % V and Vs are equal UP TO SIGN!

% We are free to change the sign of the singular values S(k,k) and the corresponding columns of the matrix Vs. 
% The matrices V (from PCA) and Vs (from SVD) will be found equal if these signs are chosen appropriately.
for k=1:m
    if norm(V(:,k)-Vs(:,k))>norm(V(:,k)+Vs(:,k))
        Vs(:,k)=-Vs(:,k);
        S(k,k)=-S(k,k);
    end;
end;

figure(4), clf
% plot shape mode matrix V as obtained by principal component analysis
subplot(1,2,1), hold on
slist=(1:m)*ds;
smat=repmat(slist,m,1);
pcolor(repmat((1:m),m,1),smat',V), axis tight, shading flat
xlabel('k'), ylabel('s [um]')
colormap jet, caxis([-1 1]), colorbar('YTick',(-1:1),'YTickLabel',{'-1','0','1'}) % set color-code
title('shape modes from PCA')
% plot shape mode matrix V as obtained by singular value decomposition (after adjusting signs)
subplot(1,2,2), hold on
pcolor(repmat((1:m),m,1),smat',Vs), axis tight, shading flat
xlabel('k'), ylabel('s [um]')
colormap jet, caxis([-1 1]), colorbar('YTick',(-1:1),'YTickLabel',{'-1','0','1'}) % set color-code
title('shape modes from SVD')

%% kernel PCA for angular data ============================================

% feature-feature simularity matrix ---------------------------------------
C_kernel=nan(m);
for i=1:m
    for j=1:m
        C_kernel(i,j)=sum(cos(psi1(:,i)-psi1(:,j)));
    end;
end;

% kernel-centering --------------------------------------------------------
C_kernel=C_kernel-sum(C_kernel)'*ones(1,m)/m-ones(m,1)*sum(C_kernel)/m+sum(C_kernel(:))/m^2;

% plot --------------------------------------------------------------------
figure(5), clf
% show feature-feature covariance matrix
subplot(1,2,1), hold on
pcolor(smat,smat',C), shading flat
caxis([-1e3 1e3]), colorbar
xlim([0 60]), ylim([0 60])
xlabel('s [um]'), ylabel('s [um]')
title('feature-feature covariance matrix')
% show feature-feature similarity matrix
subplot(1,2,2), hold on
pcolor(smat,smat',C_kernel), shading flat
caxis([-1e3 1e3]), colorbar
xlim([0 60]), ylim([0 60])
xlabel('s [um]'), ylabel('s [um]')
title('feature-feature similarity matrix')

%% eigenvalue decomposition of similarity matrix ==========================
[V_kernel,D_kernel]=eig(C_kernel);
% sort eigenvalues in desending order
[D_kernel,ind]=sort(diag(D_kernel),'descend');
V_kernel=V_kernel(:,ind);

% shape modes only def'd up to sign: adjust sign to match linear PCA
for i=1:m
    if V(:,i)'*V_kernel(:,i)<0
        V_kernel(:,i)=-V_kernel(:,i);
    end;
end;

%% show plot
figure(6), clf
% first shape mode
subplot(1,2,1), hold on
plot(slist,V(:,1),'b','LineWidth',2)
plot(slist,V_kernel(:,1),'g','color',[0 0.8 0],'LineWidth',2)
title('first shape mode')
legend({'linear PCA','kernel PCA'})
xlabel('s [um]'), ylabel('[rad]')
% second shape mode
subplot(1,2,2), hold on
plot(slist,V(:,2),'r','LineWidth',2)
plot(slist,V_kernel(:,2),'g','color',[0 0.8 0],'LineWidth',2)
title('second shape mode')
legend({'linear PCA','kernel PCA'})
xlabel('s [um]'), ylabel('[rad]')
%% find shape scores ======================================================

% mean flagellar shape using circular mean --------------------------------
psi0_kernel=angle(mean(exp(1i*psi)));

% get shape scores by nonlinear fit ---------------------------------------
B_kernel=nan(n,2); % alloc mem
% loop over frames
for iframe=1:n
    % nonlinear fit to find shape scores that maximize similarity measure
    B_kernel(iframe,:)=fminsearch(@(param) -sum(cos(psi(iframe,:)-psi0_kernel-param(1)*V_kernel(:,1)'-param(2)*V_kernel(:,2)')),[0 0]);
end;

%% show plot of shape scores
figure(7), clf, hold on
tlist=(1:n)*dt;
plot(tlist,B1,'LineWidth',2)
plot(tlist,B_kernel(:,1),'color',[0 0.8 0],'LineWidth',2)
xlabel('time [s]'), xlim([0 0.25])
ylabel('shape scores')
legend({'linear PCA','kernel PCA'})

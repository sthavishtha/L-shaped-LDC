%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to compute flows in a single-sided L-shaped cavity
% Using Lattice Boltzmann Technique (TRT) model
% in a D2Q9 model
% Boundary conditions used:
% Zou-He boundary condition on top moving wall
% Bounce-back on remaining stationary walls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code modified based on Jonas Latt's code cavity2d.m
%available at http://wiki.palabos.org/numerics:matlab_samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% GENERAL FLOW CONSTANTS 
lx = 210; %Grid size of the domain
ly = 210; 

uLid  = 0.2; % horizontal lid velocity 
vLid  = 0;    % vertical lid velocity 
Re    = 500;  % Reynolds number 
nu    = uLid *lx / Re;     % kinematic viscosity 
omega_s = 1. / (3*nu+1./2.); % symmetrical relaxation parameter 
omega_as = 8*(2-omega_s)/(8-omega_s); %anti-symmetrical relaxation parameter
%maxT  = 1e4; % total number of iterations 
%tPlot = 10;    % cycles for graphical output
ar=0.2; %aspect ratio 

% D2Q9 LATTICE CONSTANTS 
t   = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36]; %lattice weights
cx  = [ 0, 1, 0, -1, 0, 1, -1, -1, 1];  %lattice directions
cy  = [ 0, 0, 1, 0, -1, 1, 1, -1, -1]; 
opp = [ 1, 4, 5, 2,  3, 8, 9,  6,  7]; 
lid = [2: (lx-1)]; %lid coordinates

[y,x] = meshgrid(1:ly,1:lx); 
obst = ones(lx,ly); 
obst(lid,2:ly) = 0; 
bbRegion = find(obst); %solid nodes to apply bounce back 
obs=zeros(lx,ly);
obs(1:(1-ar)*lx+1,(1-ar)*ly+1)=1.0;
obs((1-ar)*lx+1,1:(1-ar)*ly+1)=1.0;
region=find(obs); %solid nodes on the corner to apply bounce-back (sim to obstacle)

% INITIAL CONDITION: (rho=0, u=0) ==> fIn(i) = t(i) 
fIn = reshape( t' * ones(1,lx*ly), 9, lx, ly); 

err_max=1; %to keep a track of error for convergence
iter=0; %to keep a track of iterations 

% MAIN LOOP (TIME CYCLES) 
while err_max>1e-12 %Convergence criteria

% MACROSCOPIC VARIABLES 
rho = sum(fIn); 
ux_n = reshape ( (cx * reshape(fIn,9,lx*ly)), 1,lx,ly ) ./rho; 
uy_n = reshape ( (cy * reshape(fIn,9,lx*ly)), 1,lx,ly ) ./rho; 

ux_n(1,[1 lx],:)=0.0; uy_n(1,[1 lx],:)=0.0; %dirichlet boundary conditions
ux_n(1,:,1)=0.0; uy_n(1,:,1)=0.0;
ux_n(1,:,ly)=uLid; uy_n(1,:,ly)=vLid;
ux_n(region)=0.0; uy_n(region)=0.0;

%checking convergence
if(iter~=0 & rem(iter,100)==0)
    err1=reshape(ux_n-ux,[],1); err2=reshape(uy_n-uy,[],1);
    A=ux.^2+uy.^2;
    iter
    err_max=sqrt(err1'*err1+err2'*err2)/sqrt(sum(A(:)))
end
    
ux=ux_n; uy=uy_n;

% MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS 

ux(:,lid,ly) = uLid; %lid x - velocity 
uy(:,lid,ly) = vLid; %lid y - velocity 
rho(:,lid,ly) = 1 ./ (1+uy(:,lid,ly)) .* ( ... 
                sum(fIn([1,2,4],lid,ly)) + 2*sum(fIn([3,6,7],lid,ly)) ); 

% MICROSCOPIC BOUNDARY CONDITIONS: LID (Zou/He BC)
fIn(5,lid,ly) = fIn(3,lid,ly) - 2/3*rho(:,lid,ly).*uy(:,lid,ly); 
fIn(9,lid,ly) = fIn(7,lid,ly) + 1/2*(fIn(4,lid,ly)-fIn(2,lid,ly))+ ... 
                1/2*rho(:,lid,ly).*ux(:,lid,ly) - 1/6*rho(:,lid,ly).*uy(:,lid,ly); 
fIn(8,lid,ly) = fIn(6,lid,ly) + 1/2*(fIn(2,lid,ly)-fIn(4,lid,ly))- ... 
                1/2*rho(:,lid,ly).*ux(:,lid,ly) - 1/6*rho(:,lid,ly).*uy(:,lid,ly); 

% COLLISION STEP 
for i=1:9 
    cu = 3*(cx(i)*ux+cy(i)*uy); 
    fEq(i,:,:) = rho .* t(i) .* ... 
        ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2) ); 
end 

for i=1:9
    fOut(i,:,:) = fIn(i,:,:) - 0.5 * omega_s .* (fIn(i,:,:) + fIn(opp(i),:,:) - fEq(i,:,:) - fEq(opp(i),:,:))...
                - 0.5 * omega_as .* (fIn(i,:,:) - fIn(opp(i),:,:) - fEq(i,:,:) + fEq(opp(i),:,:)); 
end

% MICROSCOPIC BOUNDARY CONDITIONS: NO-SLIP WALLS (bounce-back)
for i=1:9 
    fOut(i,bbRegion) = fIn(opp(i),bbRegion); 
end 

for i=1:9 
    fOut(i,region) = fIn(opp(i),region); 
end 

% STREAMING STEP 
for i=1:9 
    fIn(i,:,: ) = circshift(fOut(i,:,: ), [0,cx(i),cy(i)]); 
end

% VISUALIZATION
% if (mod(cycle,tPlot)==0)
%     u = reshape(sqrt(ux.^2+uy.^2),lx,ly);
%     u(bbRegion) = nan;
%     imagesc(u(:,ly:-1:1)'./uLid);
%     colorbar
%     axis equal off; drawnow
% end

iter=iter+1;

end 
%Code for plotting streamlines in MATLAB
[x,y]=meshgrid(1:lx,1:ly);
streamslice((x-1)/(lx-1),(y-1)/(ly-1),squeeze(ux)',squeeze(uy)',10);
grid on
axis tight

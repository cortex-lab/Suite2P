% align mean planes from multi-plane imaging to z-stack
% allows for rotations and translation
% returns transformation PtoZ, z-stack matrix Aq
% Zpatch are z-stack patches aligned to mean images Bimg using PtoZ
function [PtoZ, A0, Zpatch, Bimg] = affinePlanestoZ(ops)

% discard first plane from multi-plane imaging (flyback)
if ops.nplanes > 1
    wpl = 2:ops.nplanes;
else
    wpl = 1;
end

%% load mean planes from suite2p files
clear Bimg;
for ipl = wpl
    planefile = sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, ipl);
    dat = load(planefile);
    Bimg{ipl} = dat.ops.mimg1(dat.ops.yrange, dat.ops.xrange);
end

%% load z-stack
zfile = sprintf('%s/stack_%s_%s.mat', ops.ZstackSavePath, ops.mouse_name, ops.date);
load(zfile);

Zplanes = ops.Zplanes;

% A is z-stack (discard first two planes -- flyback)
A = Mimg0(:,:,3:Zplanes);

% make A into correct zoom
Aq = [];
if ops.Zzoom ~= ops.zoom
    zoom = ops.zoom / ops.Zzoom;
    [x,y]   = meshgrid([1:size(A,1)]);
    [xq,yq] = meshgrid(linspace(1, size(A,1), size(A,1)*zoom));
    for j = 1:size(A,3)
        Aq(:,:,j) = interp2(x,y,A(:,:,j),xq,yq);
    end
    Ly = dat.ops.Ly;
    Lx = dat.ops.Lx;
    ry = -floor(Ly/2) + floor(size(Aq,1)/2);
    rx = -floor(Lx/2) + floor(size(Aq,2)/2);
    Aq = Aq(ry+[1:Ly], rx+[1:Lx],:);
    
else
    Aq = A;
end
% Aq at same zoom as planes

%% approx angle between z-stack and planes
% default is ratio between number of planes in z-stack vs imaging
% angle in radians
if numel(wpl) > 1
    ang = -1 * atan(ops.Zplanes/ops.nplanes / size(Bimg{wpl(1)},1));
else
    ang = 0;
end

%%
[cx1, ix1, cZ1] = ZtoPlanes(Aq, Bimg, ang, ops.useGPU);
cx1
% use Z-positions of these planes to initialize affine transformation

%% which of these planes are in z-stack?
clf;
hold all;
cmap = colormap('jet');
cmap = flipud(cmap);
cm = cmap(round(linspace(1,64,max(wpl))), :);
for j = wpl
    plot(cZ1(:,j),'color',cm(j,:));
    text(double(ix1(j,3)), double(cx1(j)+.01), num2str(j), 'color', cm(j,:),...
        'fontweight', 'bold','fontsize',14);
end
ylabel('correlation');
xlabel('depth');
drawnow;

if length(wpl) > 1
    prompt = {'First plane within z-stack:','Last plane within z-stack:'};
    dlg_title = 'Input planes';
    num_lines = 1;
    defaultans = {num2str(wpl(1)),num2str(wpl(end))};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
end
try
    ipl  = str2num(answer{1}) : str2num(answer{2});
    % exclude any planes outside range
    ipl(ipl<wpl(1))   = [];
    ipl(ipl>wpl(end)) = [];
    fprintf('computing transformation for planes %d to %d\n', ipl(1), ipl(end));
catch
    disp('your inputs were not increasing numerics, set to default (all planes)');
end

%% find local maxima of patches - divide each plane into 64 patches

x1 = [];
y1 = [];
z1 = [];
for j = 1:length(Bimg)
    PtoZ{j} = [];
end
ny        = 30;
nx        = 30;
A0        = pad3Dzeros(Aq, ny, nx, 0);
for iplane = ipl
    %%
    [Ny, Nx]  = size(Bimg{iplane});
    
    
    
    npatch = 8;
    pixy   = floor(Ny/npatch);
    pixx   = floor(Nx/npatch);
    
    % patch coordinates in plane
    y0     = floor(pixy/2) + 1 + [0:pixy:pixy*(npatch-1)];
    x0     = floor(pixx/2) + 1 + [0:pixx:pixx*(npatch-1)];
    yc      = repmat(y0, npatch, 1);
    yc      = yc(:);
    xc      = repmat(x0, 1, npatch);
    xc      = xc(:);
    
    % make yx patches with centers at (yc, xc) and size (pixy, pixx)
    ipix0 = patchXY(Ny, Nx, yc, xc, pixy, pixx);
    
    % fft of patches
    Bi = Bimg{iplane};
    Bp = reshape(Bi(ipix0(:)), pixy, pixx, npatch^2);
    m2 = fft(fft(Bp,[],1),[],2);
    eps0 = single(1e-20);
    m2 = m2./(eps0 + abs(m2));
    
    % patches of z-stack
    [ny nx Np] = size(A0);
    zspread = 15;
    % recentering
    ry = - floor(Ny/2) + floor(ny/2);
    rx = - floor(Nx/2) + floor(nx/2);
    
    % initialized patch coordinates in z-stack
    yxz    = [(yc+ry)  (xc+rx)  -(yc+ry+ix1(iplane,1))*tan(ang)];
    yxz    = bsxfun(@plus, yxz, ix1(iplane,:));
    
    niter  = 5;
    yxzit  = zeros([size(yxz) niter]);
    epsit  = zeros(niter,1);
    Tit    = zeros(3,3,niter);
    for it = [1:niter]
        yxz0  = yxz;
        epsit0 = epsit;
        ipixz = patchXY(ny, nx, yxz0(:,1), yxz0(:,2), pixy, pixx);
        
        % make patches of z-stack with z depth of 2*zspread + 1
        zz    = yxz0(:,3);
        Ai = reshape(A0, [], Np);
        Ai = reshape(Ai(ipixz(:),:), pixy, pixx, npatch^2, Np);
        Az = NaN*zeros(size(Ai,1),size(Ai,2),size(Ai,3),zspread*2 + 1);
        for j = 1:npatch^2
            zinds      = max(1,floor(zz(j))-zspread) : min(Np,floor(zz(j))+zspread);
            Az(:,:,j,max(0,zspread-floor(zz(j)))+[1:length(zinds)]) = Ai(:,:,j,zinds);
        end
        
        % fft of patches
        m1 = fft(fft(Az,[],1),[],2);
        m1 = m1./(eps + abs(m1));
        
        [iyxz, cc] = xcorr_patches(m1, m2);
        
        % fit lines predicting z-stack position from Y and X
        yxz       = bsxfun(@plus, iyxz, yxz);
        X         = [yc xc ones(size(yc))]';
        Taff      = ((X * X') \ (X * yxz))';
        
        yxz       = round((Taff*X)');
        epsit(it) = median(cc);
        yxzit(:,:,it) = yxz;
        Tit(:,:,it)   = Taff;
    end
    
    [~,ibest] = max(epsit);
    yxz       = yxzit(:,:,ibest);
    Taff      = Tit(:,:,ibest);
    
    PtoZ{iplane} = Taff;
    
end


%%
for j = 1:length(Bimg)
    Zpatch{j} = [];
end

% check out aligned images
for iplane = ipl
    %%
    Taff   = PtoZ{iplane};
    Bi     = Bimg{iplane};
    
    [Ny, Nx]  = size(Bimg{iplane});
    
    npatch = 8;
    pixy   = floor(Ny/npatch);
    pixx   = floor(Nx/npatch);
    
    % patch coordinates in plane
    y0     = floor(pixy/2) + 1 + [0:pixy:pixy*(npatch-1)];
    x0     = floor(pixx/2) + 1 + [0:pixx:pixx*(npatch-1)];
    yc      = repmat(y0, npatch, 1);
    yc      = yc(:);
    xc      = repmat(x0, 1, npatch);
    xc      = xc(:);
    % make yx patches with centers at (yc, xc) and size (pixy, pixx)
    ipixz = patchXY(ny, nx, yc+ry, xc+rx, pixy, pixx);
        
    Zp     = zeros(Ny, Nx);
    for j = 1:npatch^2
        zcoord = round(Taff * [yc(j); xc(j); 1]);
        zpatch = A0(zcoord(1) - floor(pixy/2) + [0:pixy-1],...
            zcoord(2) - floor(pixx/2) + [0:pixx-1],zcoord(3));
        Zp(ipix0(:,j)) = zpatch;
    end
    
    Zpatch{iplane} = Zp;
    
    clf;
    subplot(1,2,1),
    imagesc(Bi);
    title(sprintf('plane %d',iplane));
    colormap('gray')
    subplot(1,2,2),
    imagesc(Zp)
    title('aligned z patches');
    colormap('gray')
    drawnow;
    disp('press enter to go to next plane...');
    pause;    
end




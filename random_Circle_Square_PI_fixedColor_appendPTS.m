
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate mpeg_4 file of randomly generated points in a circle
%AUTHOR: J. DeMarco & chatGPT assistant
%CREATION DATE: 07/27/25
%LAST MODIFIED: 07/27/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Random square sampling vs. circle (fixed colors, cumulative points) + live pi estimate & CI
% - Points are NEVER trimmed or redrawn from scratch.
% - We append new points each frame using animatedline/addpoints.
% - Inside circle: GREEN; Outside circle: BLUE.

rng(1);

%% --- Geometry ---
R = 1;
xlimv = [-R R]; ylimv = [-R R];

%% --- Accuracy / performance controls ---
Npoints     = 1.0e5;    % total samples (increase to reduce error; memory ~ 2*Npoints*8 bytes)
ptsPerFrame = 0.5e3;    % samples processed per frame
fps         = 30;       % video frame rate
updateEvery = 2;        % update the title every k frames
useHalton   = true;     % low-discrepancy sampling if available

% Marker sizes
size_in  = 8;
size_out = 8;

%% --- Output path: next to this script (fallback to pwd) ---
scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptFullPath);
end
outFileName = 'random_points_square_vs_circle_cumulative_pi.mp4';
outFilePath = fullfile(scriptDir, outFileName);

%% --- Generate samples in the square [-R,R]^2 ---
if useHalton && exist('haltonset','file')
    hs = haltonset(2,'Skip',1e3,'Leap',1e2);
    hs = scramble(hs,'RR2');
    U  = net(hs, Npoints);               % [0,1]^2
    xs = (U(:,1)*2 - 1) * R;             % [-R, R]
    ys = (U(:,2)*2 - 1) * R;
else
    warning('Using pseudorandom sampling (haltonset not found or disabled).');
    xs = (2*rand(Npoints,1) - 1) * R;
    ys = (2*rand(Npoints,1) - 1) * R;
end
insideMask = xs.^2 + ys.^2 <= R^2;

%% --- Figure & axes ---
fig = figure('Color','w','Position',[100 100 800 820]);
ax  = axes('Parent',fig);
hold(ax,'on'); axis(ax,'equal'); box(ax,'on');
xlim(ax,xlimv); ylim(ax,ylimv);

% Draw boundaries once
t = linspace(0,2*pi,400);
plot(ax, R*cos(t), R*sin(t), 'k-', 'LineWidth', 1.0);                     % circle
plot(ax, xlimv([1 2 2 1 1]), ylimv([1 1 2 2 1]), 'k-', 'LineWidth', 1.0); % square
xlabel(ax,'x'); ylabel(ax,'y');

% Two-line title (first line updated periodically)
hTitle = title(ax, { ...
    '\pi \approx --, N = --, 95% CI: [--, --]', ...
    'Uniform square sampling: inside circle (green) vs. outside (blue)' ...
}, 'Interpreter','tex', 'FontWeight','bold');

% Animated lines for cumulative points (no lines, marker-only)
green = [0.0 0.6 0.0];
blue  = [0.0 0.4470 0.7410];

hIn  = animatedline(ax, 'Color', green, 'LineStyle','none', 'Marker','.', 'MarkerSize', size_in);
hOut = animatedline(ax, 'Color', blue,  'LineStyle','none', 'Marker','.', 'MarkerSize', size_out);

%% --- Video writer ---
vw = VideoWriter(outFilePath, 'MPEG-4');
vw.FrameRate = fps;
open(vw);

%% --- Monte Carlo accumulators ---
Ntotal  = 0;
Ninside = 0;

nFrames = ceil(Npoints/ptsPerFrame);
for k = 1:nFrames
    i1 = (k-1)*ptsPerFrame + 1;
    i2 = min(k*ptsPerFrame, Npoints);

    xb = xs(i1:i2);
    yb = ys(i1:i2);
    ib = insideMask(i1:i2);

    % Monte Carlo counts
    Nbatch   = numel(xb);
    NinBatch = sum(ib);
    Ntotal   = Ntotal  + Nbatch;
    Ninside  = Ninside + NinBatch;

    % Current estimates
    pi_hat = 4 * (Ninside / Ntotal);
    p_hat  = pi_hat / 4;
    SE     = 4 * sqrt(max(p_hat*(1-p_hat), eps) / Ntotal);
    CI95   = 1.96 * SE;

    % Append new points cumulatively (no trimming)
    if any(ib)
        addpoints(hIn,  xb(ib),  yb(ib));
    end
    if any(~ib)
        addpoints(hOut, xb(~ib), yb(~ib));
    end

    % Update title less frequently
    if mod(k, updateEvery) == 0 || k == nFrames
        set(hTitle,'String', { ...
            sprintf('\\pi \\approx %.8f   N = %d   95%%%% CI: [%.8f, %.8f]', ...
                    pi_hat, Ntotal, pi_hat - CI95, pi_hat + CI95), ...
            'Uniform square sampling: inside circle (green) vs. outside (blue)' ...
        });
    end

    drawnow limitrate;
    writeVideo(vw, getframe(fig));
end

close(vw);
disp(['Saved video to: ' outFilePath]);


%% ---- helpers ----
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function S = fadeSizes(m, sz)
% Linear size fade from oldest -> newest
S = linspace(sz(1), sz(2), m).';
end
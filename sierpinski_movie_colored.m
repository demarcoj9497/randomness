
function sierpinski_movie_colored(outputFile, totalPoints, pointsPerFrame, seed)
% SIERPINSKI_MOVIE_COLORED  Chaos game animation of the Sierpiński triangle with color coding.
% Colors indicate which of the three *corner triangles* (at the first subdivision)
% each point lies in.
%
% Usage:
%   sierpinski_movie_colored('sierpinski_movie_colored.mp4', 200000, 1000, 2025)
%
% Args:
%   outputFile     - MP4 file to write (e.g., 'sierpinski_movie_colored.mp4')
%   totalPoints    - total number of points to plot (excluding burn-in)
%   pointsPerFrame - how many new points per frame
%   seed           - RNG seed for reproducibility
%
% Notes:
% - We classify each point into one of the three large corner triangles:
%       T_A = (A, midAB, midAC)
%       T_B = (B, midAB, midBC)
%       T_C = (C, midAC, midBC)
%   Points in the central (removed) triangle will not be visited by the chaos game.
% - Colors are set to three distinct options; adjust RGB values below as you wish.
%
% Author: ChatGPT (GPT-5 Thinking)

if nargin < 1 || isempty(outputFile);     outputFile = 'sierpinski_movie_colored.mp4'; end
if nargin < 2 || isempty(totalPoints);    totalPoints = 60000; end
if nargin < 3 || isempty(pointsPerFrame); pointsPerFrame = 300; end
if nargin < 4 || isempty(seed);           seed = 2025; end

% ---------- Geometry ----------
A = [0, 0];
B = [1, 0];
h = sqrt(3)/2;
C = [0.5, h];
verts = [A; B; C];

% Midpoints
Mab = (A + B)/2;
Mac = (A + C)/2;
Mbc = (B + C)/2;

% Corner triangles
TA = [A; Mab; Mac];
TB = [B; Mab; Mbc];
TC = [C; Mac; Mbc];

% ---------- Colors (RGB triplets) ----------
colA = [0.85 0.10 0.10];  % red-ish
colB = [0.10 0.45 0.85];  % blue-ish
colC = [0.10 0.65 0.30];  % green-ish

% ---------- RNG & video ----------
rng(seed);
vw = VideoWriter(outputFile, 'MPEG-4');
vw.Quality  = 100;
vw.FrameRate = 30;
open(vw);

burnIn = 1000;
N = totalPoints;
numFrames = ceil(N / pointsPerFrame);

% Figure
hfig = figure('Color','w','Position',[100 100 800 800]);
hax  = axes('Parent', hfig);
axis(hax,'equal'); axis(hax,'off');
xlim(hax,[0 1]); ylim(hax,[0 h]);
hold(hax,'on');

% Fast buffers
bufA = zeros(pointsPerFrame, 2);
bufB = zeros(pointsPerFrame, 2);
bufC = zeros(pointsPerFrame, 2);
nA = 0; nB = 0; nC = 0;

% Start point
p = [0.12, 0.34];
drawn = 0;

for f = 1:numFrames
    nThis = min(pointsPerFrame, N - drawn);
    nA = 0; nB = 0; nC = 0;

    for k = 1:(nThis + burnIn*(f==1))
        v = verts(randi(3), :);
        p = (p + v)/2;
        if f==1 && k <= burnIn
            continue; % skip burn-in
        end
        % Classify p into TA, TB, or TC
        if pointInTri(p, TA(1,:), TA(2,:), TA(3,:))
            nA = nA + 1;
            bufA(nA,:) = p;
        elseif pointInTri(p, TB(1,:), TB(2,:), TB(3,:))
            nB = nB + 1;
            bufB(nB,:) = p;
        elseif pointInTri(p, TC(1,:), TC(2,:), TC(3,:))
            nC = nC + 1;
            bufC(nC,:) = p;
        else
            % Should not happen (central hole not visited),
            % but ignore if encountered due to numerical issues.
        end
    end

    if nA > 0
        scatter(hax, bufA(1:nA,1), bufA(1:nA,2), 1, colA, 'filled');
    end
    if nB > 0
        scatter(hax, bufB(1:nB,1), bufB(1:nB,2), 1, colB, 'filled');
    end
    if nC > 0
        scatter(hax, bufC(1:nC,1), bufC(1:nC,2), 1, colC, 'filled');
    end

    if mod(f,10)==0 || f==1 || f==numFrames
        title(hax, sprintf('Sierpinski Triangle - Colored Chaos Game (frame %d/%d)', f, numFrames), 'FontName','Helvetica','FontSize',14);
    end

    frame = getframe(hfig);
    writeVideo(vw, frame);

    drawn = drawn + nThis;
end

hold(hax,'off');
close(vw);
disp(['Wrote ', outputFile]);

end % function

% ---------- Helper: point-in-triangle via barycentric technique ----------
function inside = pointInTri(P, V1, V2, V3)
% Returns true if P is in triangle V1V2V3 (including edges).
v0 = V3 - V1;
v1 = V2 - V1;
v2 = P  - V1;

% Solve for barycentric coords: v2 = a*v0 + b*v1
dot00 = dot(v0, v0);
dot01 = dot(v0, v1);
dot02 = dot(v0, v2);
dot11 = dot(v1, v1);
dot12 = dot(v1, v2);
invDen = 1 / (dot00*dot11 - dot01*dot01 + eps);

a = (dot11*dot02 - dot01*dot12) * invDen;
b = (dot00*dot12 - dot01*dot02) * invDen;

inside = (a >= -1e-12) && (b >= -1e-12) && (a + b <= 1 + 1e-12);
end

% function stim = mkDotsC(stimSz, stimSz, dotDirection, dotSpeed, stimLoc
%                      stim.density, dotCoherence, dotContrast, dotLuminance,
%                      dotSz, densityStyle, sampleFactor, interpMethod)
%
% MKDOTSC makes a circular aperture drifting dot stimulus.
%
% Required arguments:
% N.B. All relevant parameters have b- counterparts that control the
% background noise patch (regardless of which patch is at the nearer
% disparity).
% stimSz            The dimensions of the entire stimulus, in [Y X T] coordinates;
% maskRad           The dimensions of the stimulus aperture (scalar).
% dotDirection      The direction of movement, in radians, with 0 = rightward.
%                   If dotDirection is just one number, the dots will move
%                   in that direction at all times. If dotDirection is a
%                   vector of the same length as the number of frames in
%                   the stimulus, the direction of movement in each frame
%                   will be specified by the corresponding element of
%                   dotDirection.
% dotSpeed          The speed of movement, in frames/second.
%                   If dotSpeed is just one number, the dots will move
%                   with that speed at all times. If dotSpeed is a
%                   vector of the same length as the number of frames in
%                   the stimulus, the speed of movement in each frame
%                   will be specified by the corresponding element of
%                   dotSpeed.
%
% Optional arguments:
% stimLoc           The location of the stimulus within the surround.
%                   Options are: N, NE, E, SE, S, SW, W, NW, DEFAULT = C
% stimDist          The distance from stimulus center to aperture center. 0
%                   if stimLoc == C.
% density           The density of the dots, which can be from 0 to 1. DEFAULT = .1
% dotCoherence      The coherence of dot motion. DEFAULT = 1.
% dotContrast       The Michelson contrast ratio of the dots with the background.
% dotLuminance      The luminance of the dots, which can be 0 to 1.
%                   DEFAULT = 1 (i.e. 100% Luminance).
% dotSz         The radius of the dots. If dotSz < 0, the dots will
%                   be single pixels. If dotSz > 0, the dots will be
%                   Gaussian blobs with sigma = .5 * dotSz. DEFAULT = -1                    
% dotPlacementStyle The number of dots is calculated by multiplying the
%                   stim.density by the size of the image window. If the dots
%                   are placed randomly, some dots will overlap, so the
%                   actual dot density will be lower than stim.density. If,
%                   however, the dots are placed exactly so that they never
%                   overlap, all the dots will move in register. Both of
%                   these are problems that only occur when you use pixel
%                   dots rather than Gaussian dots. Your choices for
%                   dotPlacementStyle are 'random' (for the first problem)
%                   or 'exact' (for the second problem). DEFAULT = 'random'

function stim = mkDispDotsRotMovie(stim)
ppdeg = 54; %New LCD Monitor
% ppdeg = 21; % Old CRT Monitor

% Parse arguments out of varargin and set default values
if ~isfield(stim, 'loc');        stim.loc = [0 0];               end %center of stimulus (px)
if ~isfield(stim, 'density');    stim.density = 0.002;           end %signal dot density (dots/px)
if ~isfield(stim, 'bdensity');   stim.bdensity = 0.002;          end %noise dot density
if ~isfield(stim, 'disparity');  stim.disparity = [-0.2,0.2];    end %disparity of [signal, noise] (degrees)
if ~isfield(stim, 'coherence');  stim.coherence = 1;             end %signal coherence (0-1)
if ~isfield(stim, 'bcoherence'); stim.bcoherence = 0;            end %noise coherence
if ~isfield(stim, 'contrast');   stim.contrast = -1;             end %signal michelson contrast (-1 if defining luminance)
if ~isfield(stim, 'bcontrast');  stim.bcontrast = -1;            end %noise contrast
if ~isfield(stim, 'dotSz');      stim.dotSz = 2;                 end %dot diameter (px)
if ~isfield(stim, 'luminance');  stim.luminance = [1 0];         end %dot luminance of [signal, signal background] (0-1)
if ~isfield(stim, 'bluminance'); stim.bluminance = [1 0];        end %dot luminance of [noise, noise background]
if ~isfield(stim, 'color');      stim.color = [1 1 0];           end %all dot color (red and green will be set automatically)
if ~isfield(stim, 'lifetime');   stim.lifetime = 0;              end %dot lifetime
if ~isfield(stim, 'background'); stim.background = 0.0;          end %background luminance

if ~isfield(stim, 'size') || ~isfield(stim, 'direction') || ~isfield(stim, 'speed')
    error('stim needs field size, direction, and speed');
end

%set default parameters
if ~isfield(stim, 'maskRad')
    stim.maskRad = floor(min(stim.size(1:2))./4);
    stim.loc = [0 0];
end
if ~isfield(stim, 'loc')
    stim.loc = [0 0];
end
if ~isfield(stim, 'bMaskRad')
    stim.bMaskRad = floor(min(stim.size(1:2))./2);
    stim.bLoc = [0 0];
end
if ~isfield(stim, 'bLoc')
    stim.bLoc = [0 0];
end

% resize stim.size, dotDirection and dotSpeed if necessary
if length(stim.bdirection) == 1
    stim.bdirection = repmat(stim.bdirection, stim.size(3), 1);
end
if length(stim.bspeed) == 1
    stim.bspeed = repmat(stim.bspeed, stim.size(3), 1);
end
color = reshape(stim.color,1,1,length(stim.color)).*stim.luminance(1);

if length(stim.bdirection) ~= stim.size(3)
    error('If stim.direction is a vector, it must have the same number of entries as there are frames in the stimulus.');
end
if length(stim.bspeed) ~= stim.size(3)
    error('If stim.speed is a vector, it must have the same number of entries as there are frames in the stimulus.');
end

% Make sure stim size is odd to center the circular aperture
if ~mod(stim.size(1),2) || ~mod(stim.size(2),2)
    warning('Stimulus size must be odd to center circular aperture.\n Stimulus has been trimmed\n');
    stim.size = stim.size - 1 + mod(stim.size,2);
end
if (stim.bcontrast ~= -1) %calculate values for dots based on background and michelson contrast
    if(~stim.background);   error('Background must be >0 for contrast calculation');    end
    stim.bluminance(1) = stim.background.*(stim.bcontrast(1)+1)./(1-stim.bcontrast(1));
    stim.bluminance(2) = stim.background.*(stim.bcontrast(2)+1)./(1-stim.bcontrast(2));
    stim.bluminance = stim.bluminance./80.1;
    if(isfield(stim,'luminance'));  warning('Overriding luminance input with contrast');    end
elseif (stim.bluminance ~= [1 1])
    stim.bluminance(1) = stim.bluminance(1);
    stim.bluminance(2) = stim.bluminance(2);
end

%%%%%%%%%%%%%%%%% NOW, ON WITH THE CODE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%
bstim = stim;
bstim.density = stim.bdensity;
bstim.coherence = stim.bcoherence;
bstim.luminance = stim.bluminance;
bstim.contrast = stim.bcontrast;
bstim.direction = stim.bdirection;
bstim.maskRad = stim.bMaskRad;
bstim.loc = stim.bLoc;
bstim.speed = stim.bspeed;
noiseStim = mkDotsRotMovie(bstim);
%flipBook(noiseStim.s);
signalStim = mkDotsRotMovie(stim); %create signal patch

gmask = repmat([0 1 0]',stim.size(3),1);
rmask = repmat([1 0 0]',stim.size(3),1);

for i = 1:size(noiseStim.s,3)
    gs(:,:,i) = gmask(i).*noiseStim.s(:,:,i);
    rs(:,:,i) = rmask(i).*noiseStim.s(:,:,i);
    gsig(:,:,i) = gmask(i).*signalStim.s(:,:,i);
    rsig(:,:,i) = rmask(i).*signalStim.s(:,:,i);
end
disp = round(stim.disparity.*ppdeg);
d = zeros(stim.size(1),stim.size(2)+max(abs(disp))+1,stim.size(3)*3);
xcd = ceil(size(d,2)/2)-1;
xc = floor(stim.size(2)/2);
zd = xcd-xc;
offset = round(abs(disp/2));
if disp(1)>=0
   d(:,zd-offset(1)+1:zd-offset(1)+stim.size(2),:) = gsig+d(:,zd-offset(1)+1:zd-offset(1)+stim.size(2),:);
   d(:,zd+offset(1)+1:zd+offset(1)+stim.size(2),:) = rsig+d(:,zd+offset(1)+1:zd+offset(1)+stim.size(2),:);
else
   d(:,zd-offset(1)+1:zd-offset(1)+stim.size(2),:) = rsig+d(:,zd-offset(1)+1:zd-offset(1)+stim.size(2),:);
   d(:,zd+offset(1)+1:zd+offset(1)+stim.size(2),:) = gsig+d(:,zd+offset(1)+1:zd+offset(1)+stim.size(2),:);
end
if disp(2)>=0
    d(:,zd-offset(2)+1:zd-offset(2)+stim.size(2),:) =   gs+d(:,zd-offset(2)+1:zd-offset(2)+stim.size(2),:);
    d(:,zd+offset(2)+1:zd+offset(2)+stim.size(2),:) =   rs+d(:,zd+offset(2)+1:zd+offset(2)+stim.size(2),:);
else
    d(:,zd-offset(2)+1:zd-offset(2)+stim.size(2),:) =   rs+d(:,zd-offset(2)+1:zd-offset(2)+stim.size(2),:);
    d(:,zd+offset(2)+1:zd+offset(2)+stim.size(2),:) =   gs+d(:,zd+offset(2)+1:zd+offset(2)+stim.size(2),:);
end
    
%%Following three lines for fixation point
% fp(:,:,1:2) = ones(3,3,2);
% fp(:,:,3) = zeros(3,3,1);
% d(xcd-1:xcd+1,xcd-1:xcd+1,:) = repmat(fp,1,1,stim.size(3));

d(d>stim.luminance(1)) = stim.luminance(1);
stim.s = d;
% stim=signalStim; 
% bidir = mkBiStim(signalStim.s,noiseStim.s);
% flipBook(stim.s);
end

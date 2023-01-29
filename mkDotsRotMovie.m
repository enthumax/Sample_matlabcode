% function stim = mkDotsC(stimSz, stimSz, dotDirection, dotSpeed, stimLoc
%                      stim.density, dotCoherence, dotContrast, dotLuminance,
%                      dotSz, densityStyle, sampleFactor, interpMethod)
%
% MKDOTSC makes a circular aperture drifting dot stimulus.
%
% Required arguments:
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
% stim.density        The density of the dots, which can be from 0 to 1. DEFAULT = .1
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

function stim = mkDotsRotMovie(stim)


% Parse arguments out of varargin and set default values
if ~isfield(stim, 'loc');        stim.loc = [0 0];               end
if ~isfield(stim, 'density');    stim.density = 0.02;            end
if ~isfield(stim, 'coherence');  stim.coherence = 1;             end
if ~isfield(stim, 'contrast');   stim.contrast = -1;             end
if ~isfield(stim, 'dotSz');      stim.dotSz = 2;                 end
if ~isfield(stim, 'luminance');  stim.luminance = [1 0];         end
if ~isfield(stim, 'color');      stim.color = [1 1 1];           end
if ~isfield(stim, 'lifetime');   stim.lifetime = Inf;            end
if ~isfield(stim, 'background'); stim.background = 0;            end%0.0445;       end


if ~isfield(stim, 'size') || ~isfield(stim, 'direction') || ~isfield(stim, 'speed')
    error('stim needs fields size, direction, and speed');
end
fsize = stim.size(1:2);
stim.size(1:2) = ceil(fsize.*sqrt(2));

if ~isfield(stim, 'maskRad')
    stim.maskRad = floor(min(stim.size(1:2))./2);
    stim.loc = [0 0];
end

% resize stim.size, dotDirection and dotSpeed if necessary
if length(stim.direction) == 1
    stim.direction = repmat(stim.direction, stim.size(3), 1);
end
if length(stim.speed) == 1
    stim.speed = repmat(stim.speed, stim.size(3), 1);
end

if length(stim.direction) ~= stim.size(3)
    error('If stim.direction is a vector, it must have the same number of entries as there are frames in the stimulus.');
end
if length(stim.speed) ~= stim.size(3)
    error('If stim.speed is a vector, it must have the same number of entries as there are frames in the stimulus.');
end

% Make sure stim size is odd to center the circular aperture
if ~mod(stim.size(1),2) || ~mod(stim.size(2),2)
    warning('Stimulus size must be odd to center circular aperture.\n Stimulus has been trimmed\n');
    stim.size = stim.size - 1 + mod(stim.size,2);
end

%%%%%%%%%%%%%%%%% NOW, ON WITH THE CODE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%

if (stim.contrast ~= -1)
    if(~stim.background);   error('Background must be >0 for contrast calculation');    end
    stim.luminance(1) = stim.background.*(stim.contrast(1)+1)./(1-stim.contrast(1));
    stim.luminance(2) = stim.background.*(stim.contrast(2)+1)./(1-stim.contrast(2));
    stim.luminance = stim.luminance./80.1;
    if(isfield(stim,'luminance'));  warning('Overriding luminance input with contrast');    end
elseif (stim.luminance ~= [1 1])
    stim.luminance(1) = stim.luminance(1);
    stim.luminance(2) = stim.luminance(2);
end
color = reshape(stim.color,1,1,length(stim.color)).*stim.luminance(1);

lifetime = stim.lifetime;                     %in frames
pathlength = lifetime.*stim.speed(1);         %in px

% There is a buffer area around the image so that we don't have to worry
% about getting wraparounds exactly right.  This buffer is twice the size
% of a dot diameter.
bufferSize = round(max(stim.speed(:))+stim.dotSz...
    +round(min(stim.size(1:2)/2)*sqrt(2)-min(stim.size(1:2)/2)));
rad = floor(stim.dotSz./2);

% the 'frame' is the field across which the dots drift. The final ouput of
% this function will consist of 'snapshots' of this frame without buffer.
% We store the size of the frame and a coordinate system for it here:
frameSzX = stim.size(2) + 2.*bufferSize;
frameSzY = stim.size(1) + 2.*bufferSize;
frameSz = [frameSzY frameSzX];
[xFrame, yFrame] = meshgrid([1:frameSzX], [1:frameSzY]);
yFrame = flipud(yFrame);


% set up mask
cFrame = ceil(([frameSzY frameSzX]./2));
cStim = cFrame+[-1*stim.loc(1), stim.loc(2)];

% nDots is the number of coherently moving dots in the stimulus.
% nDotsNonCoherent is the number of noncoherently moving dots.

sqOffset = min(stim.size(1:2))/2;

pos = (abs(xFrame-cFrame(2))<=sqOffset) & (abs(yFrame-cFrame(1))<=sqOffset);
% vPos = (abs(xFrame-cFrame(2))<=sqOffset) & (yFrame-cFrame(1)<=sqOffset) & (-yFrame+cFrame(1)<=sqOffset-sin(stim.ds)*pathlength);
% flipBook(vPos);
pos = find(pos == 1);
vPos = pos;
nDots = max([1, round(stim.density.*length(pos))]);

nDotsNonCoherent = round((1-stim.coherence).*stim.density.*length(pos));

% Set the initial dot positions.
% dotPositions is a matrix of positions of the coherently moving dots in
% [y, x] coordinates; each row in dotPositions stores the position of one
% dot.
%% For motion %%
% cd P:\labFolderNew\Anjani\dispstimuli_videos\new_LCD\size3p5deg
% load('randseed_sig_new','dotInds');

% %% For stat %%
dotInds = randperm(length(vPos), nDots);
% save('randseed_sig_new','dotInds')
% % 

%% Common %%
dotInds = vPos(dotInds);
[I,J] = ind2sub(frameSz, dotInds);
dotPositions = [I,J];

nKill = nDots./lifetime;        %dots to kill per frame
nKill = [0:stim.size(3)].*nKill;
nKill = diff(floor(nKill));
dot2kill = 1;
% dotInds2 = sub2ind(frameSz, dotPositions2(:,1), dotPositions2(:,2));


% s will store the output. After looping over each frame, we will trim away
% the buffer from s to obtain the final result.
r = zeros(frameSzY, frameSzX, stim.size(3).*3)+stim.luminance(2);
for t = 1:stim.size(3)
    
    % move the positions of all the dots
    dotVelocity = [0, 1];
    dotVelocity = dotVelocity*stim.speed(t);
    dotPositions = dotPositions + repmat(dotVelocity, size(dotPositions, 1), 1);
    tmpDotPositions = round(dotPositions);
    dotInds = sub2ind(frameSz, tmpDotPositions(:,1), tmpDotPositions(:,2));
    
    % FOR RANDOMIZED NON-COHERENCE IB_EDIT
    if nDotsNonCoherent>0
        dotIndsNonCoherent = randperm(length(vPos), nDotsNonCoherent);
        dotIndsNonCoherent = vPos(dotIndsNonCoherent);
        overlap = ismember(dotIndsNonCoherent,[dotInds]);
        while sum(overlap)>0
            newInds = randperm(length(vPos), sum(overlap));
            dotIndsNonCoherent(overlap) = vPos(newInds);
            overlap = ismember(dotIndsNonCoherent,[dotInds]);
        end

        randDotIndex = randperm(nDots,nDotsNonCoherent);
        dotInds(randDotIndex) = dotIndsNonCoherent(1:nDotsNonCoherent);
        [tmpDotPositions(randDotIndex,1), tmpDotPositions(randDotIndex,2)] = ind2sub(frameSz, dotInds(randDotIndex));
        dotPositions(randDotIndex,:) = tmpDotPositions(randDotIndex,:);
    end
    
    % FOR LIFETIME RESTRICTION
    if nKill(t)>0
        dotIndsRespawn = randperm(length(vPos), nKill(t));
        dotIndsRespawn = vPos(dotIndsRespawn);
        overlap = ismember(dotIndsRespawn,dotInds);
        while sum(overlap)>0
            newInds = randperm(length(vPos), sum(overlap));
            dotIndsRespawn(overlap) = vPos(newInds);
            overlap = ismember(dotIndsRespawn,dotInds);
        end

        kInds = [dot2kill:dot2kill+nKill(t)-1];
        kInds = 1+mod(kInds,nDots);
        dot2kill = dot2kill+nKill(t);
        dotInds(kInds) = dotIndsRespawn(1:nKill(t));
        [tmpDotPositions(kInds,1), tmpDotPositions(kInds,2)] = ind2sub(frameSz, dotInds(kInds));
        dotPositions(kInds,:) = tmpDotPositions(kInds,:);
    end
    
    % wrap around for all dots that have gone past the image borders
    w = tmpDotPositions(:,2)>cStim(2)+sqOffset;
    while sum(w)>0
        tmpDotPositions(w,2) = tmpDotPositions(w,2)-2*sqOffset;
        dotPositions(w,:) = tmpDotPositions(w,:);
        dotInds(w,:) = sub2ind(frameSz, tmpDotPositions(w,1), tmpDotPositions(w,2));
        w = tmpDotPositions(:,2)>cStim(2)+sqOffset;
    end
    
    % rotate dot positions
    rotMat = [cos(stim.direction(t)), -sin(stim.direction(t)); sin(stim.direction(t)), cos(stim.direction(t))];
    tmpDotPositions = round((tmpDotPositions-cStim)*rotMat+cStim);
    
    s.dotPositions(:,:,t) = tmpDotPositions;
    s.dotSz(:,:,t) = repmat(stim.dotSz,1,nDots);
    s.color(:,:,t) = repmat(stim.color,1,nDots);
    
    
    % prepare a matrix of zeros for this frame
    thisFrame = zeros([frameSz,3])+stim.background;
    
    if stim.dotSz>1
        if mod(stim.dotSz,2)
            for j = 1:length(color)
                for i = 1:size(tmpDotPositions,1)
                    thisFrame(tmpDotPositions(i,1)-rad:tmpDotPositions(i,1)+rad, ...
                        tmpDotPositions(i,2)-rad:tmpDotPositions(i,2)+rad,j) = color(j);
                end
            end
        else
            for j = 1:length(color)
                for i = 1:size(tmpDotPositions,1)
                    thisFrame(tmpDotPositions(i,1)-rad+1:tmpDotPositions(i,1)+rad, ...
                        tmpDotPositions(i,2)-rad+1:tmpDotPositions(i,2)+rad,j) = color(j);
                end
            end
        end
    else
        for j = 1:3
            thisFrame(tmpDotPositions,j) = color(j);
        end
    end
    r(:,:,(t-1)*3+1:(t-1)*3+3) = flipud(thisFrame);%.*mask;

end
% Now trim away the buff
s = r(bufferSize+1:end-bufferSize, bufferSize+1:end-bufferSize, :);
if(~strcmp(stim.maskRad,'none'))
    s = patchMask(s,'C',stim.maskRad);
end

s(s>stim.luminance(1)) = stim.luminance(1);
trim = (stim.size(1:2)-fsize)./2;
s = s(trim+1:end-trim,trim+1:end-trim,:);
stim.size(1:2) = fsize;
stim.s = s;
% flipBook(s);
end

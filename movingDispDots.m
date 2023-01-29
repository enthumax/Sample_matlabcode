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

function [stim,rt,resp] = movingDispDots(stim,display,duration)
%ppdeg = 21;


ppdeg=angle2pix(display,1);


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
if ~isfield(stim, 'color');      stim.color = [30 80 0];           end %all dot color (red and green will be set automatically)
if ~isfield(stim, 'lifetime');   stim.lifetime = 0;              end %dot lifetime
if ~isfield(stim, 'background'); stim.background = 0.0;          end %background luminance

if ~isfield(stim, 'size') || ~isfield(stim, 'direction') || ~isfield(stim, 'speed')
    error('stim needs field size, direction, and speed');
end

stim.size(1)=stim.size(1).*ppdeg;
stim.size(2)=stim.size(2).*ppdeg;
stim.maskRad=stim.maskRad.*ppdeg;
stim.bmaskRad=stim.bmaskRad.*ppdeg;
stim.loc=stim.loc.*ppdeg;
stim.speed=stim.speed.*ppdeg;
stim.bspeed=stim.bspeed.*ppdeg;
stim.density=stim.density./ppdeg;
stim.bdensity=stim.bdensity./ppdeg;



% fsize = stim.size(1:2);
% stim.size(1:2) = ceil(fsize.*sqrt(2));

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
if ~isfield(stim, 'bmaskRad')
    stim.bmaskRad = floor(min(stim.size(1:2))./2);
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
%color = reshape(stim.color,1,1,length(stim.color)).*stim.luminance(1);

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


lifetime = stim.lifetime;                     %in frames
% pathlength = lifetime.*stim.speed(1); %in px
% pathlengthn=lifetime.*stim.bspeed(1);

% There is a buffer area around the image so that we don't have to worry
% about getting wraparounds exactly right.  This buffer is twice the size
% of a dot diameter.
bufferSize = round(max(stim.speed(:))+stim.dotSz+round(min(stim.size(1:2)/2)*sqrt(2)-min(stim.size(1:2)/2)));
%bufferSize=0;
%rad = floor(stim.dotSz./2);

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
%cFrame=ceil(display.resolution/2);
%stim.loc=ceil(display.resolution/2)+stim.loc*ppdeg;


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
nDotsN= max([1, round(stim.bdensity.*length(pos))]);
nDotsNonCoherent = round((1-stim.coherence).*stim.density.*length(pos));
nDotsNonCoherentN = round((1-stim.bcoherence).*stim.bdensity.*length(pos));

% Set the initial dot positions.
% dotPositions is a matrix of positions of the coherently moving dots in
% [y, x] coordinates; each row in dotPositions stores the position of one
% dot.
dotInds = randperm(length(vPos), nDots);
dotInds = vPos(dotInds);
[I,J] = ind2sub(frameSz, dotInds);
dotPositions = [I,J];

dotIndsN = randperm(length(vPos), nDotsN);
dotIndsN = vPos(dotIndsN);
[In,Jn] = ind2sub(frameSz, dotIndsN);
dotPositionsN = [In,Jn];

nKill = nDots./lifetime;        %dots to kill per frame
nKill = [0:stim.size(3)].*nKill;
nKill = diff(floor(nKill));
dot2kill = 1;

nKillN = nDotsN./lifetime;        %dots to kill per frame
nKillN = [0:stim.size(3)].*nKillN;
nKillN = diff(floor(nKillN));
dot2killN=1;
% dotInds2 = sub2ind(frameSz, dotPositions2(:,1), dotPositions2(:,2));

 dotVelocity = [0, 1];
  dotVelocityN = [0, 1];
 dotVelocity = dotVelocity*stim.speed(1);
 dotVelocityN = dotVelocityN*stim.bspeed(1); 
 stim.noise_level=0;
 startime=GetSecs;

t=0;
resp=NaN;
rt=NaN;
%for t = 1:stim.size(3)
while GetSecs-startime<duration && t<stim.size(3)
    
%  [keyIsDown, secs, keyCode] = KbCheck; 
%  
%  if keyIsDown
%  resp=find(keyCode);
%   resp=resp(1);
%   rt=secs-startime;
% if resp>0
%     break;
% end
%  else
   t=t+1;
    
    % move the positions of all the dots
   
    dotPositions = dotPositions + repmat(dotVelocity, size(dotPositions, 1), 1);
    tmpDotPositions = round(dotPositions);
    dotInds = sub2ind(frameSz, tmpDotPositions(:,1), tmpDotPositions(:,2));
    
    if nDotsNonCoherent>0
        dotIndsNonCoherent = randperm(length(vPos), nDotsNonCoherent);
        dotIndsNonCoherent = vPos(dotIndsNonCoherent);
        overlap = ismember(dotIndsNonCoherent,dotInds);
        while sum(overlap)>0
            newInds = randperm(length(vPos), sum(overlap));
            dotIndsNonCoherent(overlap) = vPos(newInds);
            overlap = ismember(dotIndsNonCoherent,dotInds);
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
    
     % Put a circular aperture mask for these rectangular stimuli
    if(~strcmp(stim.maskRad,'none'))
     idx=((tmpDotPositions(:,2)-(cStim(2))).^2+(tmpDotPositions(:,1)-(cStim(1))).^2) < stim.maskRad.^2;
     clear tmpDotPositionsC;
     tmpDotPositionsC(:,1)=tmpDotPositions(idx,1);
     tmpDotPositionsC(:,2)=tmpDotPositions(idx,2);
     tmpDotPositions=tmpDotPositionsC;
    end
    sizes=repmat(stim.dotSz*1.6,1,length(tmpDotPositions));
    color=repmat(stim.color'*2,1,length(tmpDotPositions));
    colorG=repmat([0 stim.color(2) 0]',1,length(tmpDotPositions));
    colorR=repmat([stim.color(1) 0 0]',1,length(tmpDotPositions));
     disp = round(stim.disparity.*ppdeg);
    offset = round(abs(disp/2));
    center=[display.resolution(1)/2-cFrame(2),display.resolution(2)/2-cFrame(1)];
    %center=[700 700];
    %center=[display.windowsize(3)/2 display.windowsize(4)/2];
    Screen('BlendFunction',display.windowPtr, 'GL_SRC_ALPHA','GL_ONE_MINUS_SRC_ALPHA');
    stim.dot_type=2;
    
if disp(1)>0
    clear  tmpDotPositionsG
    clear  tmpDotPositionsR
    tmpDotPositionsG(:,2)=tmpDotPositions(:,2)-offset(1);
    tmpDotPositionsR(:,2)=tmpDotPositions(:,2)+offset(1);
    Screen('DrawDots',display.windowPtr,[tmpDotPositionsG(:,2)';tmpDotPositions(:,1)'], sizes, colorG, center,stim.dot_type);
    Screen('DrawDots',display.windowPtr,[tmpDotPositionsR(:,2)';tmpDotPositions(:,1)'], sizes, colorR, center,stim.dot_type);
    display.fixation.color = [255,255,0];
display.fixation.flip=0;
drawFixation(display);
elseif disp(1)<0
    sizes=repmat(stim.dotSz*1.2,1,length(tmpDotPositions));
    clear  tmpDotPositionsG
    clear  tmpDotPositionsR
    tmpDotPositionsG(:,2)=tmpDotPositions(:,2)+offset(1);
    tmpDotPositionsR(:,2)=tmpDotPositions(:,2)-offset(1);
    
    Screen('DrawDots',display.windowPtr,[tmpDotPositionsG(:,2)';tmpDotPositions(:,1)'], sizes, colorG, center,stim.dot_type);
    Screen('DrawDots',display.windowPtr,[tmpDotPositionsR(:,2)';tmpDotPositions(:,1)'], sizes, colorR, center,stim.dot_type);
     display.fixation.color = [255,255,0];
display.fixation.flip=0;
drawFixation(display);
else
   Screen('DrawDots',display.windowPtr,[tmpDotPositions(:,2)';tmpDotPositions(:,1)'], sizes, color, center,stim.dot_type);
    display.fixation.color = [255,255,0];
display.fixation.flip=0;
drawFixation(display);
end
    
    
    % Now plot noise dots
    dotPositionsN = dotPositionsN + repmat(dotVelocityN, size(dotPositionsN, 1), 1);
    tmpDotPositionsN = round(dotPositionsN);
    dotIndsN = sub2ind(frameSz, tmpDotPositionsN(:,1), tmpDotPositionsN(:,2));
    
%     overlapS = ismember(dotIndsN,dotInds);
%        % nS=sum(overlapS);
%         while sum(overlapS)>0
%             newIndsS = randperm(length(vPos), sum(overlapS));
%             dotIndsN(overlapS) = vPos(newIndsS);
%             overlapS = ismember(dotIndsN,dotInds);
%         end
    % FOR RANDOMIZED NON-COHERENCE IB_EDIT
    
    
    if nDotsNonCoherentN>0
        dotIndsNonCoherentN = randperm(length(vPos), nDotsNonCoherentN);
        dotIndsNonCoherentN = vPos(dotIndsNonCoherentN);
        
        overlapS = ismember(dotIndsNonCoherentN,dotInds);
        nS=sum(overlapS);
        while sum(overlapS)>0
            newIndsS = randperm(length(vPos), sum(overlapS));
            dotIndsNonCoherentN(overlapS) = vPos(newIndsS);
            overlapS = ismember(dotIndsNonCoherentN,dotInds);
        end
       if nS<nDotsNonCoherentN 
           
            overlapN = ismember(dotIndsNonCoherentN,dotIndsN);
            while sum(overlapN)>0
            newIndsN = randperm(length(vPos), sum(overlapN));
            dotIndsNonCoherentN(overlapN) = vPos(newIndsN);
            overlapN = ismember(dotIndsNonCoherentN,dotIndsN);
            end
            stim.noise(t)=nDotsNonCoherentN/nDotsN;
            
       elseif nS==nDotsNonCoherentN 
           stim.noise(t)=nDotsNonCoherentN/nDotsN;
       else
           stim.noise(t)=nS/nDotsN;
       end
        

        randDotIndexN = randperm(nDotsN,nDotsNonCoherentN);
        dotIndsN(randDotIndexN) = dotIndsNonCoherentN(1:nDotsNonCoherentN);
        [tmpDotPositionsN(randDotIndexN,1), tmpDotPositionsN(randDotIndexN,2)] = ind2sub(frameSz, dotIndsN(randDotIndexN));
        dotPositionsN(randDotIndexN,:) = tmpDotPositionsN(randDotIndexN,:);
    end
    
   
    
    if nKillN(t)>0
        dotIndsRespawnN = randperm(length(vPos), nKillN(t));
        dotIndsRespawnN = vPos(dotIndsRespawnN);
        overlapN = ismember(dotIndsRespawnN,dotIndsN);
        while sum(overlapN)>0
            newIndsN = randperm(length(vPos), sum(overlapN));
            dotIndsRespawnN(overlapN) = vPos(newIndsN);
            overlapN = ismember(dotIndsRespawnN,dotIndsN);
        end

        kIndsN = dot2killN:dot2killN+nKillN(t)-1;
        kIndsN = 1+mod(kIndsN,nDotsN);
        dot2killN = dot2killN+nKillN(t);
        dotIndsN(kInds) = dotIndsRespawnN(1:nKillN(t));
        [tmpDotPositionsN(kIndsN,1), tmpDotPositionsN(kIndsN,2)] = ind2sub(frameSz, dotIndsN(kIndsN));
        dotPositionsN(kIndsN,:) = tmpDotPositionsN(kIndsN,:);
    end
    
   
    
     wN = tmpDotPositionsN(:,2)>cStim(2)+sqOffset;
    while sum(wN)>0
        tmpDotPositionsN(wN,2) = tmpDotPositionsN(wN,2)-2*sqOffset;
        dotPositionsN(wN,:) = tmpDotPositionsN(wN,:);
        dotIndsN(wN,:) = sub2ind(frameSz, tmpDotPositionsN(wN,1), tmpDotPositionsN(wN,2));
        wN = tmpDotPositionsN(:,2)>cStim(2)+sqOffset;
    end
    
    
    rotMatN = [cos(stim.bdirection(t)), -sin(stim.bdirection(t)); sin(stim.bdirection(t)), cos(stim.bdirection(t))];
    tmpDotPositionsN = round((tmpDotPositionsN-cStim)*rotMatN+cStim);
  
    if(~strcmp(stim.bmaskRad,'none'))
     idxN=((tmpDotPositionsN(:,2)-(cStim(2))).^2+(tmpDotPositionsN(:,1)-(cStim(1))).^2) < stim.bmaskRad.^2;
     clear tmpDotPositionsCN;
     tmpDotPositionsCN(:,1)=tmpDotPositionsN(idxN,1);
     tmpDotPositionsCN(:,2)=tmpDotPositionsN(idxN,2);
     tmpDotPositionsN=tmpDotPositionsCN;
    end
    
    sizesN=repmat(stim.dotSz*2,1,length(tmpDotPositionsN));
    
    colorN=repmat(stim.color'*2,1,length(tmpDotPositionsN));
   
    colorGN=repmat([0 stim.color(2)*2 0]',1,length(tmpDotPositionsN));
   
    colorRN=repmat([stim.color(1)*2 0 0]',1,length(tmpDotPositionsN));
   
    
if disp(2)>0
    clear  tmpDotPositionsNG
    clear  tmpDotPositionsNR
    tmpDotPositionsNG(:,2)=tmpDotPositionsN(:,2)-offset(2);
    tmpDotPositionsNR(:,2)=tmpDotPositionsN(:,2)+offset(2);
    Screen('DrawDots',display.windowPtr,[tmpDotPositionsNG(:,2)';tmpDotPositionsN(:,1)'], sizesN, colorGN, center,stim.dot_type);
    Screen('DrawDots',display.windowPtr,[tmpDotPositionsNR(:,2)';tmpDotPositionsN(:,1)'], sizesN, colorRN, center,stim.dot_type);
elseif disp(2)<0
    clear  tmpDotPositionsNG
    clear  tmpDotPositionsNR
    tmpDotPositionsNG(:,2)=tmpDotPositionsN(:,2)+offset(2);
    tmpDotPositionsNR(:,2)=tmpDotPositionsN(:,2)-offset(2);
    Screen('DrawDots',display.windowPtr,[tmpDotPositionsNG(:,2)';tmpDotPositionsN(:,1)'], sizesN, colorGN, center,stim.dot_type);
    Screen('DrawDots',display.windowPtr,[tmpDotPositionsNR(:,2)';tmpDotPositionsN(:,1)'], sizesN, colorRN, center,stim.dot_type);
else
   Screen('DrawDots',display.windowPtr,[tmpDotPositionsN(:,2)';tmpDotPositionsN(:,1)'], sizesN, colorN, center,stim.dot_type);
end



Screen('Flip',display.windowPtr);
stim.noise_level=stim.noise_level+stim.noise(t);
 
 end
  
stim.noise_level=stim.noise_level/stim.size(3);


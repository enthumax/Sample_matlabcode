% This code generates transparent motion stimuli- 2 random dot patches moving in two directions at two disparities. You can either
% use the output structure 's' or create a video to save the stimuli.


% Set monitor/display parameters
frameRate = 100;         %frames/sec
stim.duration=0.1;
 %20 px/degree, duration = stimSz(3)/frameRate ([301,301,100]) - Changed to 337 to fit new screen size (LCD, 45.39 x 26.53)
  %disparity (degrees) [signal, noise] with negative closer to fixation pt
stim.dotSz = 4;               %dot size px
ppdeg=21; %pixels per degree for the monitor 

%RDP 1
stim.sigradindeg=10;
stim.bradindeg=0;
stim.maskRad =stim.sigradindeg;

   %direction 0< dir <2*pi
stim.speed = 10.*ppdeg/frameRate;   %2           %speed px/frame = (deg/sec)*(px/degree)*frameRate %Originally set at 1.5 (4.58 for 5deg/s and 9.16 for 10deg/s)
stim.density = 0.0005;           % (dots/px) %0.0045== 2 dots per deg sq
stim.luminance = [1, 0];        %luminance (0-1) only change first element
stim.coherence = 1;             %coherence 0-1
stim.loc = [0,0];   %[Y,X] stimulus center coordinates


%RDP 2
%stim.bMaskRad = 'full';
%stim.bMaskRad =stim.bradindeg*ppdeg;
stim.bspeed = stim.speed;
stim.bdirection = 120;
stim.bdirection = ((stim.bdirection)/180)*pi;
stim.bdensity = 0;% stim.density;
stim.bcoherence = 1;
stim.bluminance = [1,0];
%stim.bMaskRad=(1+2.5)*ppdeg;
%stim.size = [stim.bMaskRad*2+1, stim.bMaskRad*2+1,stim.duration*frameRate]; 
stim.size = [449, 449,stim.duration*frameRate]; 
%stim.size = [200, 200,stim.duration*frameRate]; 

%stim.bLoc=[0,2.5*ppdeg];

stim.refdir=0;
stim.direction = ((stim.refdir)/180)*pi;
stim.disparity(1) = 0.1; %RDP 1 disparity
stim.disparity(2) = -0.1; %RDP 2 disparity

%stim.noisedisp=-0.1;
%stim.dirstep=[-30 -22.5 -15 -10 -5 -3 -1.5 1.5 3 5 10 15 22.5 30];
%stim.dirstep=[60 3 5 10 15 22.5 30]; %-30 -22.5 -15 -10 -5 -3 -60 0 
stim.dirstep=0;

% stim.dirstep= 30 ;
%c=['N' 'F'];


s = mkDispDotsRotMovie(stim); 
%s=mktwopatchDotsRotMovie(stim);
shMovie(s,sprintf('test'),frameRate);

 
%This is a PsychoToolbox experiment script that runs a fine direction discrimination task of a signal RDP at a disparity overlapping with a noise RDP at another disparity. 
% Subjects press "Leftarrow" for counterclockwise and "Rightarrow" for clockwise' direction reports.
clear all
addpath(genpath('P:\labFolderNew\PsychoTBox'))
subject='dummie';
date='12april2020';
display.dist = 61;  %cm
display.width = 52.5; %cm
%display.resolution=[angle2pix(display,52.5) angle2pix(display,29.5)];
display.skipChecks = 0; %avoid Screen's timing checks and verbosity
%display.resolution = [1024 768];
display.bkColor= [0 0 0];

%Specify stimulus parameters%
stim.dotSz =  3;     %dot size px
stim.lifetime=Inf; 


%Signal patch
stim.disparity(1) =0;    %disparity of signal patch (degrees)
stim.maskRad = 2.5;             %radius of signal patch in deg

stim.density = 0.015;           %signal patch density (dots/deg)
stim.luminance = [1, 0];        %luminance (0-1) only change first element
stim.coherence = 0.8;         %coherence 0-1
stim.loc = [0,0];               %[Y,X] stimulus center coordinates


stim.bmaskRad = 2.5;    %radius of noise patch
stim.bspeed = 0;

stim.bdirection = 0;
stim.bdensity =0.015;
stim.bcoherence = 0.4;
stim.bluminance = [1,0];
stat_stim=stim;
stat_stim.speed=0;
stat_stim.bcoherence=0;
%stat_stim.dir=[-180:45:0 45:45:135] ;
%stat_stim.dir=[-45 45 -135 135] ;
stat_stim.dir=-90; %reference direction


design.stepSize = 1;
%dir=[-60:-46 -44:-30]; %-180 is leftward 90 is downward 0 rightward -90 upward
stim.dir=[-10 -7 -6:-3 -1.5 0 1.5 3:6 7 10];
%stim.dir=[-4 4];

noise_disparity=[-0.1 0.1]; %N or F





design.nTrials = length(stat_stim.dir)*length(stim.dir)*length(noise_disparity);
design.nblocks=1;
results.RT=zeros(design.nTrials,design.nblocks);


%design.revnum=3;


% We'll use the same timing parameters:
design.stimDur = 0.15;  
design.segl = 0.5;  %duration of stimulus (sec)
design.segu = 0.5;
design.responseInterval = 1.2;  %window in time to record the subject's response
design.ITI = 0.5;             %inter-trial-interval (time between the start of each trial)


%keychars=['l' 'r' 'esc'];  %% l=leftarrow r=rightarrow
choice1='left';
choice2='right';
choice3='esc';

%Let's open the screen and get going!
try
    screens = Screen('Screens');
    
    display.screenNum= max(screens);
    display.windowPtr=display.screenNum;
    %display.screenNum=1;
    
    display.skipChecks=1;
    
    display = OpenWindow(display);
    
   
    nFrames = secs2frames(display,design.stimDur);
    
    stim.speed = 5/display.frameRate;              %speed of signal dots in deg/frame
    %stim.bspeed = stim.speed;
    stat_stim.speed=stim.speed;



 %stim.size = [display.resolution(2) display.resolution(1) nFrames];
 
    stim.size = [25 25 nFrames]; %frame size in deg; make sure this is bigger than either of the signal or noise patches


    %drawText(display,[0,6],'Press "u" for up and "d" for
    %down',[255,255,255]);e
    drawText(display,[0,6],'Press "Leftarrow" for counterclockwise and "Rightarrow" for clockwise',[255,255,255]);
    drawText(display,[0,5],'Press Any Key to Begin.',[255,255,255]);
display.fixation.shape='FillRect';
    display= drawFixation(display);
    while KbCheck; end
    KbWait;

    %record the clock time at the beginning
    startTime = GetSecs;
    
   
    jumbocomb=[];
    %loop through the trials
    for block_num=1:design.nblocks
        for i=1:length(stat_stim.dir)
        comb=vertcat([stim.dir' repmat(noise_disparity(1),length(stim.dir),1)],[stim.dir' repmat(noise_disparity(2),length(stim.dir),1)]);
        comb2=horzcat(comb,stat_stim.dir(i)*ones(size(comb,1),1));
        jumbocomb=vertcat(jumbocomb,comb2);
        end
        rand_comb=Shuffle(jumbocomb,2);
        trial_num=1;
    while trial_num<design.nTrials+1
    
    stat_stim.direction=(rand_comb(trial_num,3)/180)*pi;   %reference direction 
    stim.dirdiff=(rand_comb(trial_num,1)/180)*pi;
    stim.direction = stat_stim.direction+(rand_comb(trial_num,1)/180)*pi;
    stim.disparity(2)=rand_comb(trial_num,2);
    
    stat_stim.disparity(2)=stim.disparity(2);
    results.trialdir_ref(trial_num,block_num)=rand_comb(trial_num,3);
    results.trialdir_diff(trial_num,block_num) = rand_comb(trial_num,1);
    results.trialnoisedisp(trial_num,block_num) = stim.disparity(2);
    design.stimDur_stat(trial_num,block_num) = design.segl+(design.segu-design.segl)*rand(1);
    nFrames_stat=secs2frames(display,design.stimDur_stat(trial_num,block_num));
    stat_stim.size=[stim.size(1:2) nFrames_stat];
    
    
    
   % while rev<design.revnum
%         if rev>2
%            design.stepSize = 5; 
%         end
        
        


        %Show the stimulus
        %movingDots(display,dots,design.stimDur);
         
         display.fixation.shape='FillRect';

           startTime1 = GetSecs;  
%            fromH=angle2pix(display,stim.loc(2));
%            fromV=angle2pix(display,stim.loc(1));
%            toH=angle2pix(display,fromH+stim.maskRad*cos(stat_stim.direction));
%             toV=angle2pix(display,fromV+stim.maskRad*sin(stat_stim.direction));
%             penWidth=2;
%            Screen('DrawLine', display.windowPtr,80, fromH, fromV, toH, toV,penWidth);
%            waitTill(trial_num*design.ITI,startTime);
          %movingDispDots(stat_stim,display,design.stimDur_stat(trial_num,block_num));
          [final_stim,rt,keyCode]=movingDispDots(stim,display,design.stimDur);
         results.noise_level(trial_num,block_num)=final_stim.noise_level;
         
%     if isnan(round(keyCode))
%         keys=[];
%     else
% 
%         keys=KbName(keyCode);  % KbName('KeyNames') to see the whole list of keycode and keyname mapping
%          %[keys, rt] = waitTill(secs);
%     end       
         display.fixation.color = [255,255,0];
         display.fixation.shape='FillOval';
         drawFixation(display);  
         starttime= GetSecs;
%          while GetSecs-starttime<design.responseInterval
%           [keyIsDown, secs, keyCode] = KbCheck; 
%   
% if keyIsDown
% resp=find(keyCode);
% resp=resp(1);
% rt=secs-starttime;
% keys=KbName(resp); % KbName('KeyNames') to see the whole list of keycode and keyname mapping
% else
%     keys=[];
% end 
%end
         
          [keys, rt] = waitTill(design.responseInterval);
          keys=char(keys);
         if size(keys,1)>1
         keys=keys(1);
         end
               

        if isempty(keys)  
            drawFixation(display);
            results.response(trial_num,block_num) = 0;
            %display.fixation.color= [255,0,0];
            trial_num=trial_num+1;
        elseif isequal(keys,choice1)  %counterclockwise
                results.response(trial_num,block_num) = -1;
                results.RT(trial_num,block_num)=rt(1);
                trial_num=trial_num+1;
        elseif isequal(keys,choice2)   %clockwise
                results.response(trial_num,block_num) = 1;
                results.RT(trial_num,block_num)=rt(1);
                trial_num=trial_num+1;
        elseif isequal(keys,choice3)  %e key for esc
            Screen('CloseAll');
        else
            drawFixation(display);
            results.response(trial_num,block_num) = NaN;
            display.fixation.color= [255,0,0];
            trial_num=trial_num+1;
        end

        display.fixation.color = [255,255,255];
        drawFixation(display);

        %Now wait for the clock to catch up to the time for the next trial
        waitTill(trial_num*design.ITI,startTime);

       end

        
    end
    catch ME
    Screen('CloseAll');
    rethrow(ME)
end
Screen('CloseAll');
%cd P:\labFolderNew\PsychoTBox\results
save(sprintf('results%s%s',date,subject), 'results','design','stim','display')

% Save the results

% save(sprintf('resultsStaircase%s%s',date,subject), 'results','design')
% figure(1)
% clf
% stairs(results.greenlum);
% 
% % Let's get fancy and plot green and red symbols where trials had correct
% % and incorrect responses respectively:
% 
% upTrials = results.response==0;
% hold on
% plot(find(upTrials),results.greenlum(upTrials),'ko','MarkerFaceColor','g');
% 
% downTrials = results.response==1;
% hold on
% plot(find(downTrials),results.greenlum(downTrials),'ko','MarkerFaceColor','r');
% % upTrials = results.response==0;
% % plot(find(upTrials),results.greenlum(upTrials),'ko','MarkerFaceColor','k');
% %set(gca,'YTick',log(2.^[-4:0]))
% %logy2raw;
% 
% xlabel('Trial Number')
% ylabel('Green_LUT value');
% 
% 	
% 
% %Then we'll loop through these intensities calculating the proportion of
% %times that 'response' is equal to 1:
% 
% greenlums = unique(results.greenlum);
% nup = zeros(1,length(greenlums));
% nTrials = zeros(1,length(greenlums));
% 
% for i=1:length(greenlums)
%     id = results.greenlum ==greenlums(i) & isreal(results.response);
%     nTrials(i) = sum(id);
%     if results.response(id)==0
%     nup(i) = nup(i)+1;
%     end
% end
% 
% pup = nup./nTrials;
% 
% figure(2)
% clf
% hold on
%  plot(greenlums,100*pup,'o-','MarkerFaceColor','b');
%  %loop through each intensity so each data point can have it's own size.
% for i=1:length(greenlums)
%     sz = nTrials(i)+2;
%     plot(greenlums(i),100*pup(i),'ko-','MarkerFaceColor','b','MarkerSize',sz);
% end
% 
% set(gca,'XTick',greenlums);
% 
% set(gca,'YLim',[0,100]);
% xlabel('Greenlum');
% ylabel('Percent VA');

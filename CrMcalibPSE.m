function CrMcalibPSE(varargin)
% 200121 edited by jyh - coherence level (QUEST), data structure (coh0)
% 190619 edited by jyh - coherence level, visual stim time, auditory stim
% time, zittering time has been changed.
% 190612 written by jyh
% visual motion detection task (up or down)
% IV: stim Type (only new for now) X coherence (0/0.2/0.4/0.6) X congruency(cong/incong/noise)
% Each condition is repeated 20 times.
% All of audio files will be randomly sampled for 4 times for each subject.
% & each sound direction (up/down/noise).
% visual motion will be shown for 150 ms
% while playing sound for 3.5 s.
% Participant have to listen to melodies for at least 1 s.
% stimOnset time is between 1~1.2s.
% Response window is 2 s from stimOnset.
% zittering time is 200 ms at maximum. 
% SN, trial, motDir, cong, stimType, coher, sdDir, sdFile, vStartTime, dResp, dRT, dCorr

try
    main
catch me % error handling
    clearEnvironment;
    fprintf(2, '\n\n???:: %s\n\n', me.message); % main error message
    for k = 1:(length(me.stack) - 1)
        current = me.stack(k);
        fprintf('Error in ==> ');
        fprintf('<a href="matlab: opentoline(''%s'',%d,0)">%s at %d</a>\n\n',current.file, current.line, current.name, current.line);
    end
    fclose('all');
end
return


function main(varargin)
testMode = [0, 0, 0]; % set all 0 at real game.

fastMode = testMode(1); %1 fast
skipSync = testMode(2);
winSize  = testMode(3); % 0 for full size, 1 for 800x600 size

fclose('all');
format short;
ClockRandSeed;

% screen setup
global cx cy

% sound setup
InitializePsychSound;

%% --------------------------- DIRECTORIES INFO
curDir = pwd; %code directory
stimDir = fullfile(curDir, 'Stim');
subDir = fullfile(curDir, 'subcodes');
datDir = fullfile(pwd, 'dataCrMcalibPSE');
if ~exist(datDir,'dir'); mkdir(datDir); end

addpath(genpath(stimDir));
addpath(genpath(subDir));

%% --------------------------- SUBJECT & FILE INFO
prompt = {'enter subject number: ','Initials (e.g. HJY)',...
    'Mac OSX or Windows? (0-OSX, 1-Windows)', ... 
    'threshold for 75% acc.'}; %200121 edited by jyh 
defaults = {'99','','1', ''};
answer = inputdlg(prompt, 'experimental setup information', 1, defaults);
[SN, NM, wdw, t] = deal(answer{:});

SN = str2double(SN);   wdw = str2double(wdw); t = str2double(t);

% open files and specify file format
baseName = ['CrMcalibPSE' sprintf('%02d%s', SN, NM)];
fileName = ['datRaw' baseName '.csv'];	dataFile = fopen(fullfile(datDir, fileName), 'a');
fileName = ['datStm' baseName '.txt'];	stimFile = fopen(fullfile(datDir, fileName), 'a');
paramName =['datVar' baseName '.mat'];

%% --------------------------- KEYBOARD INFO
% [keyboardIndex, deviceName, allInfo] = GetKeyboardIndices
% ===== Experimenter
device(1).product = 'Apple Internal Keyboard / Trackpad';	%
device(1).vendorID= 1452;
% ===== iMac
device(2).product = 'Magic Keyboard';
device(2).vendorID= 76;
% ===== K375s keyboard
device(3).product = 'Keyboard K375s';
device(3).vendorID= 1133;
% ===== Participant
% Empirisoft Research Software DirectIN`
device(4).product = 'DirectIN Hardware';
device(4).vendorID= 6017;

if wdw %windows7
Participant = []; KbName('UnifyKeyNames');
elseif ~wdw %Mac 
Participant = IDKeyboards(device(1));  % See below.
end

SPC = KbName('1!'); % proceed to next
DUKey = KbName({'8*', '2@'}); %8 - down, 2 - up %edited by jyh 200221

%% --------------------------- SCREEN SETUP

% Removes the blue screen flash and minimize extraneous warnings.
oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
oldSkipSyncTests = Screen('Preference', 'SkipSyncTests', skipSync); % When you use CRT monitor, set 0.
Screen('Preference', 'TextRenderer', 0);   % When you use high quality Korean font, set 1.
% Screen('Preference', 'ConserveVRAM', 64);

whichScreen = max(Screen('Screens')); % determining a dual monitor setup.

if ~winSize
    [D.w, screenRect] = Screen('OpenWindow', whichScreen);
else
    [D.w, screenRect] = Screen('OpenWindow', whichScreen);
    sca; smallRect = screenRect*0.8; newRect = CenterRect(smallRect,screenRect);
    [D.w, screenRect] = Screen('OpenWindow', whichScreen, [], newRect);
end

%D.bcolor	 = BlackIndex(D.w);
D.bcolor     = GrayIndex(D.w);
Screen('BlendFunction', D.w, GL_ONE, GL_ONE); % make gabor patch smooth
Screen('BlendFunction', D.w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); % for using PNG's transparency
%Screen('HideCursorHelper', D.w);
Screen('FillRect', D.w, D.bcolor);
Screen('Flip', D.w);

[cx, cy]     = RectCenter(screenRect);
D.resolution = screenRect([3,4]);

%single door room
D.width = 40;
D.height = 30;
D.dist = 60; %for the double-door room, distance should be 90

%fixation parameters
D.fixation.size = .2; %fixation dot size(degree)
D.fixation.mask = .8;  %oval mask size(degree)
%for timing setup
blanking = Screen('GetFlipInterval', D.w);
D.frameRate = 1/blanking; %Hz

%% --------------------------- TIMING SETUP

durWAIT = (round(2.0/blanking)-0.5)*blanking; %before trial started
durAUD = 3.5; %(round(4.5/blanking)-0.5)*blanking; % audio playtime = stim time 
durVIS = 0.15; %(round(1.0/blanking)-0.5)*blanking; %dot motion
durRESP = 2; %(round(2.0/blanking)-0.5)*blanking; %response window
durBREAK = (round(2.0/blanking)-0.5)*blanking; %after break
durZIT = 0.2; %(round(0.1/blanking)-0.5)*blanking; %zittering range
durITI = (round(1.5/blanking)-0.5)*blanking; %interval between trials
%calculate frames
tFrames = secs2frames(D, durAUD); %total frame length
vFrames = secs2frames(D, durVIS); %visual frame length
rFrames = secs2frames(D, durRESP); %response window frame length 
ztFrames = secs2frames(D, durZIT); %zittering frames maximum (from middle frame)
mFrame = secs2frames(D, 1.1); %round(tFrames/2); %middle

ztFrameRange = [round(mFrame-ztFrames/2), round(mFrame+ztFrames/2)];

%% --------------------------- EXPERIMENTAL DESIGN
% 200121 edited by jyh
% Reduced number of trials for coherence 0: 40 to 20
% This is because there is no motion direction for coherence 0 cond.
% we build a data structure for coherence levels except 0 first.
% Then add a structure for 0 coherence (60 trials)

%congruency X stim type X coherence (3 X 2 X 5 )

% build data structure for coherence level [x, y, z]
MD = 2; %visual motion direction: down vs. up
CG = 3; %congruency: congruent vs. incongruent vs. noise
ST = 1; %stimulus type: new(glissando)
CH = 3; %coherence
Rep = 20;
data.nTRIAL = MD*CG*ST*CH*Rep;

xIndex=0:(data.nTRIAL-1);
data.cMD = 0:(data.nTRIAL-1);
data.cMD(mod(xIndex,MD)==0)=1; % down
data.cMD(mod(xIndex,MD)==1)=2; % up

data.cCG = 0:(data.nTRIAL-1);
data.cCG(mod(fix(xIndex/MD),CG)==0)=1; % congruent
data.cCG(mod(fix(xIndex/MD),CG)==1)=2; % incongruent
data.cCG(mod(fix(xIndex/MD),CG)==2)=3; % noise

data.cSD = repmat([1 2 2 1 3 3], [1, data.nTRIAL/(MD*CG)]); %sound direction according to cMD & cCG

data.cST = 0:(data.nTRIAL-1);
data.cST(mod(fix(xIndex/(MD*CG)),ST)==0)=1; % original
data.cST(mod(fix(xIndex/(MD*CG)),ST)==1)=2; % new

data.cCH = 0:(data.nTRIAL-1);
data.cCH = mod(fix(xIndex/(MD*CG*ST)), CH)+2;  %1 should be coherence 0 later

% data structure for 0 coherence - coh0.~ 
tmp.MD = 1; % no motion direction
tmp.SD = 3; % sound direction
tmp.CH = 1;
tmp.ST = 1;
tmp.Rep = 20;
coh0.nTRIAL = tmp.MD*tmp.SD*tmp.CH*tmp.Rep;

xIndex=0:(coh0.nTRIAL-1);
coh0.cMD = 0:(coh0.nTRIAL-1);
coh0.cMD(mod(xIndex,tmp.MD)==0)=1; % down
coh0.cMD(mod(xIndex,tmp.MD)==1)=2; % up

coh0.cSD = 0:(coh0.nTRIAL-1);
coh0.cSD(mod(fix(xIndex/tmp.MD),tmp.SD)==0)=1; % down
coh0.cSD(mod(fix(xIndex/tmp.MD),tmp.SD)==1)=2; % up
coh0.cSD(mod(fix(xIndex/tmp.MD),tmp.SD)==2)=3; % noise

coh0.cCG = zeros(1, coh0.nTRIAL); %congruency is marked as 0 

coh0.cST = 0:(coh0.nTRIAL-1);
coh0.cST(mod(fix(xIndex/(tmp.MD*tmp.SD)),tmp.ST)==0)=1; % original
coh0.cST(mod(fix(xIndex/(tmp.MD*tmp.SD)),tmp.ST)==1)=2; % new

coh0.cCH = 0:(coh0.nTRIAL-1);
coh0.cCH = mod(fix(xIndex/(tmp.MD*tmp.SD*tmp.ST)), tmp.CH)+1;  %coherence level


%merge two datastruct 

data.cMD = [data.cMD, coh0.cMD];
data.cCG = [data.cCG, coh0.cCG];
data.cSD = [data.cSD, coh0.cSD];
data.cST = [data.cST, coh0.cST];
data.cCH = [data.cCH, coh0.cCH];
data.nTRIAL = data.nTRIAL + coh0.nTRIAL; %Final nTRIAL
block = data.nTRIAL/(Rep*3); 
data.bTRIAL = data.nTRIAL/block; %trials in a block

%
data.dResp   = ones(1,data.nTRIAL)*7; % 2AFC, coherence yes vs. no
data.dRT     = ones(1,data.nTRIAL)*7; % response time
data.dCorr   = ones(1,data.nTRIAL)*7; % 1 for correct, 2 for incorrect.
data.xOrder  = Shuffle(1:data.nTRIAL);

nSF = 4; %sound file samples
data.cSF = 0:(data.nTRIAL-1);
data.cSF = mod(Shuffle(data.cSF), nSF)+1;  %list 1 2 3 4..

%% --------------------------- DOT PARAMETERS
tmp.cohLevel = linspace(0,t, 3);
dots.cohLevel = round([tmp.cohLevel tmp.cohLevel(3)+tmp.cohLevel(2)],3);

%%%
dots.nDots = 500;  % number of dots.
dots.color = [255 255 255];
dots.size = angle2pix(D, 0.12);  % size of dots in degrees, converted to pixels.

dots.center = [0,0];  % center of the field of dots (x,y)
dots.speed = 6;       %degrees/second.
dots.ovalSize = [13.5, 13.5];     % size of elliptical aperture [x,y] in degrees. 3.75 diameter was used in paper.
dots.lifetime = 18; % the maximum number of frames a dot can survive for. used to be 1800.
dots.life = ceil(rand(1,dots.nDots)*dots.lifetime);

nDots = sum([dots.nDots]);
colors = zeros(3,nDots);
sizes = zeros(1,nDots);

% Pulsation parameters
dots.pulsation = 0; %set dots to pulsate
dots.pulsationLength = 0.10; % Length of time the dots stay big for.
dots.pulsationPeriod = 0.34; % Time between pulsations
dots.pulsationFactor = 2; % Factor by which the dots' size will increase. Diameter.

%Calculate total number of temporal frames and convert seconds to frames
%for other relevant fields.
pulsationLengthFrames = secs2frames(D,dots.pulsationLength);
pulsationPeriodFrames = secs2frames(D,dots.pulsationPeriod);
pulsationNum = tFrames/pulsationPeriodFrames;
pulsationStart = (1:pulsationNum)*pulsationPeriodFrames;
pulsationEnd = pulsationStart + pulsationLengthFrames;


%% --------------------------- APERTURE CONFIGURATION
%Calculate the left, right top and bottom of each aperture (in degrees)

l = dots.center(1)-dots.ovalSize(1)/2;
r = dots.center(1)+dots.ovalSize(1)/2;
b = dots.center(2)-dots.ovalSize(2)/2;
t = dots.center(2)+dots.ovalSize(2)/2;

%% --------------------------- STIM

folderName = {'01Down', '02Up', '03Noise'};
audioExtension = {'*.wav', '*.wav'};
freq = 44100; nrchannels = 2; reqlatencyclass = 2; buffer = [];
MySongHandle = PsychPortAudio('Open', [], [], reqlatencyclass, freq, nrchannels);

buffer = []; subDir = {}; samples = {}; audiodataC = {}; 
subDir{1} = 'new'; %audacity
%subDir{2} = 'original'; %

soundDuration = round(durAUD);

for stim = 1:length(subDir)
    for fd = 1:length(folderName)
        audioDir = fullfile(stimDir, subDir{stim}, folderName{fd});
        cd(audioDir);
        audioList = dir(audioExtension{stim});
        soundName = audioList(1).name;
        for cc = 1:nSF
            [audiodata, infreq] = audioread(soundName);
            if infreq ~= freq
                audiodata = resample(audiodata, freq, infreq);
            end
            if fd ~= 3 %congruent, incongruent
                durationInSec = size(audiodata, 1)/freq;
                durMax = round(durationInSec) - soundDuration;
                durMin = 1;
                range = durMax/nSF;
                rdmStart = randi([round(durMax/nSF*(cc-1)+durMin), round(durMax/nSF*cc)]);
                samples{cc} = [rdmStart*freq, rdmStart*freq+soundDuration*freq];
                [audiodataC{cc}, freq] = audioread(soundName, samples{cc});
                [~, ninchannels] = size(audiodataC{cc});
                audiodataC{cc} = repmat(transpose(audiodataC{cc}), nrchannels / ninchannels, 1);
                buffer(stim, fd, cc) = PsychPortAudio('CreateBuffer', [], audiodataC{cc}); %a queue of buffer
            elseif fd == 3 %noise
                [~, ninchannels] = size(audiodata);
                audiodata = repmat(transpose(audiodata), nrchannels / ninchannels, 1);
                buffer(stim, fd, cc) = PsychPortAudio('CreateBuffer', [], audiodata); %a queue of buffer
                
            end
        end
    end
end
clear samples audiodata audiodataC; 
cd(curDir);
%PsychPortAudio('UseSchedule', MySongHandle, 1);

%% --------------------------- GET STARTED

% Logging
fprintf(		  '\n***** Task begins at %s\n', datestr(now, 0));
fprintf(stimFile, '\n***** Task begins at %s\n', datestr(now, 0));

fprintf(          'SN, trial, motDir, cong, stimType, coher, sdDir, sdFile, vStartTime, dResp, dRT, dCorr\n');
fprintf(stimFile, 'SN, trial, motDir, cong, stimType, coher, sdDir, sdFile, vStartTime, dResp, dRT, dCorr\n');

showInstruction('inst_start.txt', D);

keysOfInterest = zeros(1,256);
keysOfInterest(SPC) = 1;
KbQueueCreate(Participant, keysOfInterest);
KbQueueWait(Participant);

KbQueueRelease(Participant);
Screen('Flip', D.w);
WaitSecs(durWAIT);

keysOfInterest = zeros(1,256);
keysOfInterest(DUKey) = 1;
KbQueueCreate(Participant, keysOfInterest);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRIAL LOOP

for trial = 1:data.nTRIAL
    %% ------------------------ CURRENT CONDITION
    pressed = 0; showMotion = 0;
    it = data.xOrder(trial);
    xMD = data.cMD(it); % 1down vs. 2up
    xCG = data.cCG(it); % 1cong vs. 2incong. vs. 3noise vs. 0coherence 0
    xST = data.cST(it); % 1new
    xCH = data.cCH(it); % coherence [0 x y z]
    xSD = data.cSD(it); % 1down vs. 2up sound direction
    xSF = data.cSF(it); % sound file no. (1~4)
    vStartFrame = randi(ztFrameRange); %start visual frames
    
    % apply conditions
    dots.direction = mod(xMD,2)*180;% down 1 should be 180, 2 should be 0. 180903 edited
    dots.coherence = dots.cohLevel(xCH);
    dots.pulsation = xST-1; %set dots to pulsate. stim1: no pulsation
    
    nCoherent = ceil(dots.coherence*dots.nDots);  %num. of coherent dots
    
    % create response queue
    KbQueueStart(Participant);
    KbQueueFlush(Participant);
    
    %% ------------------------ DOTS CURRENT CONDITION
    order=  randperm(nDots);
    % Intitialize the dot positions and define some other initial parameters
    %Generate random starting positions
    dots.x = (rand(1,dots.nDots)-.5)*dots.ovalSize(1) + dots.center(1);
    dots.y = (rand(1,dots.nDots)-.5)*dots.ovalSize(2) + dots.center(2);
    %Create a direction vector for a given coherence level
    direction = rand(1,dots.nDots)*360;
    direction(1:nCoherent) = dots.direction;  %Set the 'coherent' directions
    
    %Calculate dx and dy vectors in real-world coordinates, for both
    %regular and oddball conditions
    dots.dx = dots.speed*sin(direction*pi/180)/D.frameRate;
    dots.dy = -dots.speed*cos(direction*pi/180)/D.frameRate;
    
    %Fill in the 'colors' and 'sizes' vectors for this field
    id = 1:nDots;  %index into the nDots length vector
    colors(:,order(id)) = repmat(dots.color(:),1,dots.nDots);
    sizes(order(id)) = repmat(dots.size,1,dots.nDots);
    sizes_orig = sizes;
    pulsationSizes = sizes*dots.pulsationFactor;
    
    dots.dx_orig = dots.dx;
    dots.dy_orig = dots.dy;
    
    initialDotPositions.x = dots.x;
    initialDotPositions.y = dots.y;
    %Zero out the screen position vectors and the 'goodDots' vector
    pixpos.x = zeros(1,nDots);
    pixpos.y = zeros(1,nDots);
    goodDots = zeros(1,nDots);
    
    %% ------------------------ PLAY
    %PsychPortAudio('AddToSchedule', MySongHandle, buffer(xST, xSD, xSF)); %, 1, 0.0, 1.0, 1
    PsychPortAudio('FillBuffer', MySongHandle, buffer(xST, xSD, xSF));
    PsychPortAudio('Start', MySongHandle, 1, 0, 1);
    %PsychPortAudio('DeleteBuffer', MySongHandle);
 
    for frameNum=pulsationStart(1):(tFrames+pulsationStart(1)-1) % was 1:nFrames  frameNum=pulsationStart(pulStart):(playLength+pulsationStart(pulStart)-1)
        
        %% ------------------------ SHOW
        if frameNum == pulsationStart(1)+vStartFrame | showMotion == 1
            showMotion = 1;
            
            %Update the dot position's real-world coordinates
            dots.x = dots.x + dots.dx;
            dots.y = dots.y + dots.dy;
            
            %Move the dots that are outside the aperture back one aperture width.
            dots.x(dots.x<l) = dots.x(dots.x<l) + dots.ovalSize(1);
            dots.x(dots.x>r) = dots.x(dots.x>r) - dots.ovalSize(1);
            dots.y(dots.y<b) = dots.y(dots.y<b) + dots.ovalSize(2);
            dots.y(dots.y>t) = dots.y(dots.y>t) - dots.ovalSize(2);
            
            %Increment the 'life' of each dot
            dots.life = dots.life+1;
            
            %Find the 'dead' dots
            deadDots = mod(dots.life,dots.lifetime)==0;
            
            %Replace the positions of the dead dots to random locations
            dots.x(deadDots) = (rand(1,sum(deadDots))-.5)*dots.ovalSize(1) + dots.center(1);
            dots.y(deadDots) = (rand(1,sum(deadDots))-.5)*dots.ovalSize(2) + dots.center(2);
            
            %Calculate the index for this field's dots into the whole list of
            %dots.  Using the vector 'order' means that, for example, the first
            %field is repscreenRectented not in the first n values, but rather is
            %distributed throughout the whole list.
            id = order(id);
            
            %Calculate the screen positions for this field from the real-world coordinates
            pixpos.x(id) = angle2pix(D,dots.x)+ D.resolution(1)/2;
            pixpos.y(id) = angle2pix(D,dots.y)+ D.resolution(2)/2;
            
            %Determine which of the dots in this field are outside this field's
            %elliptical aperture
            goodDots(id) = (dots.x-dots.center(1)).^2/(dots.ovalSize(1)/2)^2 +...
                (dots.y-dots.center(2)).^2/(dots.ovalSize(2)/2)^2 < 1;
            goodDots = logical(goodDots);
            
            if (find(frameNum==pulsationStart) & dots.pulsation == 1)
                sizes = pulsationSizes;
            elseif find(frameNum == pulsationEnd)
                sizes = sizes_orig;
            end
            
            Screen('DrawDots',D.w,[pixpos.x(goodDots);pixpos.y(goodDots)], sizes(goodDots), colors(:,goodDots),[0,0],1);
        end
        drawFixationDot2(D);
        %Screen('Close'); %delete back buffer
        
        
        %% ------------------------ CHECK THE LAST FRAME
        
        if frameNum == pulsationStart(1)+vStartFrame %start point
            KbQueueStart(Participant);
            KbQueueFlush(Participant);
            stimOnset = GetSecs;
        end
        
        if frameNum == pulsationStart(1)+vStartFrame+vFrames-1
            showMotion = 0;
        end
        
        %% ------------------------ COLLECT RESPONSE
        
        if (frameNum > pulsationStart(1)+vStartFrame) &&...
                (frameNum <= pulsationStart(1)+vStartFrame+rFrames -1) && ~pressed
            [pressed, firstPress] = KbQueueCheck(Participant);
            if pressed
                sT  = min(firstPress(firstPress ~= 0));	% shortest time
                xKey = find(firstPress == sT);		% key with the shortest time
                
                if length(xKey) > 1 % if more than one keys are pressed at the same time, get 6.
                    data.dCorr(it) = 6;
                else
                    data.dRT(it)	 = sT - stimOnset;
                    data.dResp(it)   = find(DUKey == xKey);
                    data.dCorr(it)   = (find(DUKey == xKey) == xMD); %should change xAlpha
                end
            end
        end
        
        
    end %end of frames
    
    %% ------------------------ STOP PLAYING
    
    PsychPortAudio('Stop', MySongHandle);
    %PsychPortAudio('UseSchedule', MySongHandle, 2);
    Screen('Flip', D.w);
    WaitSecs(durITI);
    
    %% ------------------------ LOGGING RESPONSE
    
    fprintf(           'sn: %02d, %d, con: %d, %d, %d, %1.3f, %d, %d, %1.2f, rsp: %d, %2.4f, %d\n', SN, trial, ...
        xMD, xCG, xST, dots.coherence, xSD, xSF, vStartFrame/D.frameRate, data.dResp(it), data.dRT(it), data.dCorr(it));
    fprintf(dataFile, '%02d, %d, %d, %d, %d, %1.3f, %d, %d, %1.2f, %d, %2.4f, %d\n', SN, trial, ...
        xMD, xCG, xST, dots.coherence, xSD, xSF, vStartFrame/D.frameRate, data.dResp(it), data.dRT(it), data.dCorr(it));
    fprintf(stimFile, '%02d, %d, %d, %d, %d, %1.3f, %d, %d, %1.2f, %d, %2.4f, %d\n', SN, trial, ...
        xMD, xCG, xST, dots.coherence, xSD, xSF, vStartFrame/D.frameRate, data.dResp(it), data.dRT(it), data.dCorr(it));
    
    %% ------------------------ CHECK BREAK
    if mod(trial, data.bTRIAL)==0 && (trial < data.nTRIAL)
        cntBlk = (data.nTRIAL-trial)/data.bTRIAL;
        
        Screen('TextColor', D.w, 0);
        Screen('TextSize',  D.w, 25);
        Screen('TextStyle', D.w, 1);
        
        rText = sprintf('Take a rest. Press 1 to continue...');
        if cntBlk==1
            bText = sprintf('%d block remains.', cntBlk);
        else
            bText = sprintf('%d blocks remain.', cntBlk);
        end
        rBound = Screen('TextBounds', D.w, rText)/2;
        bBound = Screen('TextBounds', D.w, bText)/2;
        Screen('DrawText', D.w, rText, cx-rBound(3), cy-rBound(4)-100);
        Screen('DrawText', D.w, bText, cx-bBound(3), cy-bBound(4)-200);
        Screen('Flip', D.w);
        
        bGO = 1;
        if fastMode, bGO = 0; end
        while bGO
            [keyIsDown, ~, keyCode] = KbCheck(Participant);
            if keyIsDown
                if keyCode(SPC)
                    bGO = 0;
                end
            end
        end
        Screen('Flip', D.w);
        WaitSecs(durBREAK);
    end
    
end

%% ------------------------ LETS WRAP UP
KbQueueRelease(); PsychPortAudio('Close', MySongHandle);
showInstruction('inst_ending.txt', D);
fprintf(           '\n***** Experiment ended at %s\n', datestr(now, 0));
fprintf(stimFile, '\n***** Experiment ended at %s\n', datestr(now, 0));
clear buffer;
save(fullfile(datDir, paramName));


eGO = 1;
if fastMode, eGO = 0; end
while eGO
    [keyIsDown, ~, keyCode] = KbCheck(Participant);
    if keyIsDown
        if keyCode(SPC)
            eGO = 0;
        end
    end
end

% clean up and go home
Screen('CloseAll');
fclose('all');
ListenChar(0);
return


%% SubFunctions

function deviceIndex = IDKeyboards (kbStruct)
devices	= PsychHID('Devices');
kbs		= find([devices(:).usageValue]==6); % value of keyboard

deviceIndex = [];
for mm=1:length(kbs)
    if strcmp(devices(kbs(mm)).product,kbStruct.product) && ...
            ismember(devices(kbs(mm)).vendorID, kbStruct.vendorID)
        deviceIndex = kbs(mm);
    end
end

if isempty(deviceIndex)
    error('No %s detected on the system',kbStruct.product);
end
return

function showInstruction(txtfile, window)
txtdir = fullfile(pwd, 'txtdir');

% Text files edited with a text editor such as vi or nano.
fid = fopen(fullfile(txtdir, txtfile), 'r', 'n','UTF-8');
xtext = fread(fid, '*char');
fclose(fid);
xtext = double(transpose(xtext));

%         Screen('TextFont', window.w, allFonts(idx).name);
Screen('TextColor', window.w, 254);
Screen('TextSize', window.w, 13);
Screen('TextStyle',window.w, 1);

Screen('FillRect', window.w, window.bcolor);
DrawFormattedText(window.w, xtext, 'center', 'center');
Screen('Flip', window.w);

return

%% EOF.

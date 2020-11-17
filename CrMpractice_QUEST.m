function CrMpractice_QUEST(varargin)
%200120 written by jyh
%This code is for estimating a coherence level for targeted acc.
%Targeted accuracy is 75%
%Visual motion detection task (up or down)
%all conditions (cong/incong/ambig) appears
%Total of 40 trials
%Initial parameters for QUEST:
%guess: 0.3; guessSd: 0.15;
% All of audio files will be randomly sampled for 4 times for each subject.
% & each sound direction (up/down/noise).
% visual motion will be shown for 150 ms
% while playing sound for 3.5 s.
% Participant have to listen to melodies for at least 1 s.
% stimOnset time is between 1~1.2s.
% Response window is 2 s from stimOnset.
% zittering time is 200 ms at maximum.
% SN, trial, stimIntensity, dCorr
% Final product: t (threshold), sd (standard deviation)

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
testMode = [0, 0, 1]; % set all 0 at real game.

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
datDir = fullfile(pwd, 'dataQUEST');
if ~exist(datDir,'dir'); mkdir(datDir); end

addpath(genpath(stimDir));
addpath(genpath(subDir));

%% --------------------------- SUBJECT & FILE INFO
prompt = {'enter subject number: ',...
    'Initials (e.g. HJY)',...
    'Mac OSX or Windows? (0-OSX, 1-Windows)',...
    'Estimate threshold (e.g. 0.2)',...
    'Estimate the standard deviation of your guess (e.g. 0.16)'};

defaults = {'99','','1','0.4','0.35'};
answer = inputdlg(prompt, 'experimental setup information', 1, defaults);
[SN, NM, wdw, tGuess, tGuessSd] = deal(answer{:});

SN = str2double(SN); wdw = str2double(wdw);
tGuess = str2double(tGuess); tGuessSd = str2double(tGuessSd);

% open files and specify file format
baseName = ['CrMpractice_QUEST' sprintf('%02d%s', SN, NM)];
fileName = ['datRaw' baseName '.csv'];	dataFile = fopen(fullfile(datDir, fileName), 'a');
fileName = ['datStm' baseName '.txt'];	stimFile = fopen(fullfile(datDir, fileName), 'a');
paramName =['datVar' baseName '.mat'];

% open files and specify file format
thresFile = fopen(fullfile(datDir, 'allSub_75threshold.csv'), 'a');

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
DUKey = KbName({'8*', '7&'}); %8 - down, 7 - up

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
    whichScreen = min(Screen('Screens')); % determining a dual monitor setup.
    
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

MD = 2; %visual motion direction: down vs. up
CG = 3; %congruency: congruent vs. incongruent vs. noise
ST = 1; %stimulus type: new(glissando)
Rep = 7;
block = 8;
data.nTRIAL = MD*CG*ST*Rep;
data.bTRIAL = data.nTRIAL;

xIndex=0:(data.nTRIAL-1);
data.cMD = 0:(data.nTRIAL-1);
data.cMD(mod(xIndex,MD)==0)=1; % down
data.cMD(mod(xIndex,MD)==1)=2; % up

data.cCG = 0:(data.nTRIAL-1);
data.cCG(mod(fix(xIndex/MD),CG)==0)=1; % congruent
data.cCG(mod(fix(xIndex/MD),CG)==1)=2; % incongruent
data.cCG(mod(fix(xIndex/MD),CG)==2)=3; % noise

data.cSD = repmat([1 2 2 1 3 3], [1, data.nTRIAL]); %sound direction according to cMD & cCG

data.dResp   = ones(1,data.nTRIAL)*7; % 2AFC, coherence yes vs. no
data.dRT     = ones(1,data.nTRIAL)*7; % response time
data.dCorr   = zeros(1,data.nTRIAL); % 
data.xOrder  = Shuffle(1:data.nTRIAL);

nSF = 4; %sound file samples
data.cSF = 0:(data.nTRIAL-1);
data.cSF = mod(Shuffle(xIndex), nSF)+1;  %list 1 2 3 4..

%% --------------------------- DOT PARAMETERS

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

%% --------------------------- QUEST PARAMS
wrongRight={'wrong', 'right'};

pThreshold=0.75; %target acc
beta=3; delta=0.01; gamma=0.5; %Weibull dist. params
grain = 0.01; %stepsize
%range = 0.3; % -range/2:grain:range/2

q=QuestCreate(tGuess, tGuessSd, pThreshold, beta, delta, gamma, grain);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRIAL LOOP

for trial = 1:data.nTRIAL
    %% ------------------------ CURRENT CONDITION
    pressed = 0; showMotion = 0;
    it = data.xOrder(trial);
    xMD = data.cMD(it); % 1down vs. 2up
    xCG = data.cCG(it); % 1cong vs. 2incong. vs. 3noise
    xST = 1; % 1new
    xCH = QuestQuantile(q); % coherence level
    xSD = data.cSD(it); % 1down vs. 2up sound direction
    xSF = data.cSF(it); % sound file no. (1~4)
    vStartFrame = randi(ztFrameRange); %start visual frames
    
    % apply conditions
    dots.direction = mod(xMD,2)*180;% down 1 should be 180, 2 should be 0. 180903 edited
    dots.coherence = xCH;
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
                    data.dCorr(it) = 0;
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
    fprintf('Trial %1d at %5.2f is %s\n',trial,xCH,char(wrongRight(data.dCorr(it)+1)));
    %fprintf(           'sn: %02d, %d, con: %d, %d, %d, %1.2f, %d, %d, %1.2f, rsp: %d, %2.4f, %d\n', SN, trial, ...
    %    xMD, xCG, xST, xCH, xSD, xSF, vStartFrame/D.frameRate, data.dResp(it), data.dRT(it), data.dCorr(it));
    fprintf(dataFile, '%02d, %d, %d, %d, %d, %1.2f, %d, %d, %1.2f, %d, %2.4f, %d\n', SN, trial, ...
        xMD, xCG, xST, xCH, xSD, xSF, vStartFrame/D.frameRate, data.dResp(it), data.dRT(it), data.dCorr(it));
    fprintf(stimFile, '%02d, %d, %d, %d, %d, %1.2f, %d, %d, %1.2f, %d, %2.4f, %d\n', SN, trial, ...
        xMD, xCG, xST, xCH, xSD, xSF, vStartFrame/D.frameRate, data.dResp(it), data.dRT(it), data.dCorr(it));
    
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
    
    %% ------------------------ UPDATE PDF
    q = QuestUpdate(q,xCH,data.dCorr(it));
    
end

%% ------------------------ LETS WRAP UP
KbQueueRelease(); PsychPortAudio('Close', MySongHandle);
showInstruction('inst_ending.txt', D);
fprintf(           '\n***** QUEST ended at %s\n', datestr(now, 0));
fprintf(stimFile, '\n***** QUEST ended at %s\n', datestr(now, 0));
clear buffer;

% Ask Quest for the final estimate of threshold.
t=QuestMean(q);		% Recommended by Pelli (1989) and King-Smith et al. (1994). Still our favorite.
sd=QuestSd(q);
fprintf('Final threshold estimate (mean+-sd) is %.2f +- %.2f\n',t,sd);
fprintf(thresFile, '%d, %.2f, %.2f', SN,t,sd);

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

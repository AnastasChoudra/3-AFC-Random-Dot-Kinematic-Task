% Instructions function added
% changed the fonts and sizes
% added more instructions between practice trials and main
% Added breaks after each block
% data saved in a Contingency table 
% Confidence function updated to a likert mouse hover scale
% Thank you message added at the end
% Demographics added to contingency table & Converted to .csv
% Results saved based on participant number
% Visual angle set to 5degrees for the apperture size!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NO EYETRACKING - PILOTING DOT TASK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Clear the workspace
close all;
clear;
% Collect participant information
participantNumber = input('Enter participant number: ');
age = input('Enter participant age: ');
sex = input('Enter participant sex (M/F): ', 's');

% Suppress PTB messages
oldLevel = Screen('Preference', 'Verbosity', 0);
% Initialize Psychtoolbox
Screen('Preference', 'SkipSyncTests', 1);
PsychDefaultSetup(2);
% Get the screen numbers
screens = Screen('Screens');
% Select the external screen if it is present, else revert to the native screen
screenNumber = max(screens);
% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
% Open an on-screen window
[window, windowRect] = Screen('OpenWindow', screenNumber, black);
% Get the size of the on-screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
% Get the center coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);
% Query the frame duration
ifi = Screen('GetFlipInterval', window);

clear display dots

% Set display parameters
display.dist = 60;  % Set at 60cm for all!
display.width = 34.5; % cm
display.skipChecks = 1; % avoid Screen's timing checks and verbosity
display.windowPtr = window;
display.frameRate = Screen('FrameRate', window);
display.resolution = [screenXpixels, screenYpixels];
display.bkColor = black;

% Define dot fields
dots.coherence = 0; % or any other valid initial value
dots.nDots = round(16.7 * pi * (2.5^2)); % number of dots in 5 degrees diameter aperture
dots.speed = 4; % degrees per second
dots.direction = 0; % will be set to left or right
dots.lifetime = 3;  % Updated lifetime to 3 frames
dots.apertureSize = [5, 5]; % diameter in degrees
dots.center = [0, 0]; % Center in degrees, not pixels
dots.color = [255, 255, 255]; % white
dots.size = 2; % pixels

% Define fixation point parameters
display.fixation.size = 0.25;  % default is .5 degrees
display.fixation.color = [255, 255, 255];  % default is white
display.fixation.flip = 0;  % do not flip yet

try
    HideCursor;
    % Display instructions
    displayInstructions(display);
    
    % Run practice trials
    runPracticeTrials(display, dots, ifi);

    % Indicate practice trials are finished
    DrawFormattedText(display.windowPtr, 'Practice trials finished. Press any key to begin the main experiment.', 'center', 'center', [255, 255, 255]);
    Screen('Flip', display.windowPtr);
    KbWait;  % Wait for any key press
    % Clear any previous key presses
    FlushEvents('keyDown');

    % Run main experiment
    DrawFormattedText(display.windowPtr, 'Main experiment is starting now. Press any key to continue.', 'center', 'center', [255, 255, 255]);
    Screen('Flip', display.windowPtr);
    KbWait;  % Wait for any key press

    runMainExperiment(display, dots, ifi, participantNumber, age, sex);
    
    % Display thank you message at the end
    DrawFormattedText(display.windowPtr, 'Thank you for participating. The experiment is now complete.', 'center', 'center', [255, 255, 255]);
    Screen('Flip', display.windowPtr);
    KbWait;  % Wait for any key press
catch ME
    Screen('CloseAll');
    rethrow(ME);
end
Screen('CloseAll');

% Restore PTB message level
Screen('Preference', 'Verbosity', oldLevel);

function [key, RT, correctDirection] = movingDots(display, dots)
    % Initialize dot positions and directions
    dots.x = (rand(1, dots.nDots) - 0.5) * dots.apertureSize(1) + dots.center(1);
    dots.y = (rand(1, dots.nDots) - 0.5) * dots.apertureSize(2) + dots.center(2);
    dots.life = randi([0, dots.lifetime], 1, dots.nDots);

    % Set correctDirection based on coherence
    if dots.coherence == 0
        correctDirection = 'n';
        coherentDirection = NaN; % No coherent direction for 0% coherence
    else
        % Determine the direction of coherent motion
        if dots.coherenceDirection == "right"
            coherentDirection = 0; % Rightward motion
            correctDirection = 'RightArrow';
        elseif dots.coherenceDirection == "left"
            coherentDirection = 180; % Leftward motion
            correctDirection = 'LeftArrow';
        else
            error('Invalid coherenceDirection. Choose "left" or "right".');
        end
    end

    % Calculate coherent and random directions
    nCoherentDots = round(dots.coherence * dots.nDots);
    nRandomDots = dots.nDots - nCoherentDots;

    % Coherent motion
    dots.dx = zeros(1, dots.nDots);
    dots.dy = zeros(1, dots.nDots);
    if ~isnan(coherentDirection)
        dots.dx(1:nCoherentDots) = dots.speed * cosd(coherentDirection) / display.frameRate;
        dots.dy(1:nCoherentDots) = -dots.speed * sind(coherentDirection) / display.frameRate;
    end

    % Random motion
    randomDirections = rand(1, nRandomDots) * 360;
    dots.dx(nCoherentDots+1:end) = dots.speed * cosd(randomDirections) / display.frameRate;
    dots.dy(nCoherentDots+1:end) = -dots.speed * sind(randomDirections) / display.frameRate;

    % Pre-render stimulus for 1 second
    for preFrame = 1:round(display.frameRate * 1)
        dots.x = dots.x + dots.dx;
        dots.y = dots.y + dots.dy;

        % Lifetime update
        dots.life = dots.life - 1;
        replife = dots.life <= 0;
        dots.life(replife) = dots.lifetime;
        dots.x(replife) = (rand(1, sum(replife)) - 0.5) * dots.apertureSize(1) + dots.center(1);
        dots.y(replife) = (rand(1, sum(replife)) - 0.5) * dots.apertureSize(2) + dots.center(2);
    end
    
    % Generate the circular Gabor visual stimulus with a specific orientation
    gaborSizePix = angle2pix(display, dots.apertureSize(1)); % Convert aperture size to pixels
    orientation = 45; % Set the desired orientation of the Gabor patch (adjust as needed)
    gradedImage = generateCircularGabor(gaborSizePix, gaborSizePix, orientation);

    % Create a texture from the graded image
    gradedImageTexture = Screen('MakeTexture', display.windowPtr, gradedImage);

    % Determine the destination rectangle for the graded image (centered at aperture)
    dstRect = CenterRectOnPointd([0 0 gaborSizePix gaborSizePix], display.resolution(1)/2, display.resolution(2)/2);

    % Initialize fade-out parameters for the graded image
    fadeOutTime = 1; % duration of fade-out in seconds
    fadeOutFrames = round(display.frameRate * fadeOutTime);
    alphaValues = linspace(1, 0, fadeOutFrames); % alpha values from 1 to 0 over fade-out frames

    frameCount = 0;  % Initialize frameCount
    key = '';
    startTime = GetSecs;

    while isempty(key)
        frameCount = frameCount + 1;

        % Clear screen
        Screen('FillRect', display.windowPtr, display.bkColor);

        % Update dot positions
        dots.x = dots.x + dots.dx;
        dots.y = dots.y + dots.dy;

        % Replot every third frame
        if mod(frameCount, 3) == 0
            dots.life = dots.life - 1;
            replife = dots.life <= 0;
            dots.life(replife) = dots.lifetime;
            dots.x(replife) = (rand(1, sum(replife)) - 0.5) * dots.apertureSize(1) + dots.center(1);
            dots.y(replife) = (rand(1, sum(replife)) - 0.5) * dots.apertureSize(2) + dots.center(2);
        end

        % Convert to pixel positions
        pixpos.x = angle2pix(display, dots.x) + display.resolution(1) / 2;
        pixpos.y = angle2pix(display, dots.y) + display.resolution(2) / 2;

        % Use the equation of an ellipse to determine which dots fall inside the circular aperture
        goodDots = (dots.x - dots.center(1)).^2 / (dots.apertureSize(1) / 2)^2 + ...
                   (dots.y - dots.center(2)).^2 / (dots.apertureSize(2) / 2)^2 < 1;

        % Draw only the 'good' dots
        Screen('DrawDots', display.windowPtr, [pixpos.x(goodDots); pixpos.y(goodDots)], dots.size, dots.color, [0, 0], 1);

        % Apply the graded image overlay with decreasing transparency
        if frameCount <= fadeOutFrames
            currentAlpha = alphaValues(frameCount);
            Screen('DrawTexture', display.windowPtr, gradedImageTexture, [], dstRect, [], [], currentAlpha);
        end

        % Flip screen
        Screen('Flip', display.windowPtr);

        % Check for key press
        [keyIsDown, timeSecs, keyCode] = KbCheck;
        if keyIsDown
            pressedKey = KbName(find(keyCode));
            if ischar(pressedKey) && any(strcmp(pressedKey, {'LeftArrow', 'RightArrow', 'n'}))
                key = pressedKey;
                RT = timeSecs - startTime;
            end
        end
    end

    % Close the texture
    Screen('Close', gradedImageTexture);
end

% Function to convert angle to pixels
function pix = angle2pix(display, angle)
    % ANGLE2PIX Converts visual angle in degrees to pixels
    % pix = angle2pix(display, angle)
    % Input:
    %   display - structure containing display properties
    %       display.dist - distance from screen (cm)
    %       display.width - width of screen (cm)
    %       display.resolution - resolution of screen (pixels)
    %   angle - visual angle (degrees)
    % Output:
    %   pix - corresponding number of pixels

    % Calculate pixels per degree
    ppd = pi * display.resolution(1) / atan(display.width / (2 * display.dist)) / 360;

    % Convert angle to pixels
    pix = angle * ppd;
end

function display = drawFixation(display)
    % Define default values for fixation if not provided
    if ~isfield(display, 'fixation')
        display.fixation.size = 0.5; % degrees
        display.fixation.color = [255, 255, 255]; % white
        display.fixation.flip = 1;
    else
        if ~isfield(display.fixation, 'size')
            display.fixation.size = 0.5; % degrees
        end
        if ~isfield(display.fixation, 'color')
            display.fixation.color = [255, 255, 255]; % white
        end
        if ~isfield(display.fixation, 'flip')
            display.fixation.flip = 1;
        end
    end

    % Calculate fixation size in pixels
    fixSize = angle2pix(display, display.fixation.size);
    % Center of the screen
    [xCenter, yCenter] = RectCenter([0 0 display.resolution(1) display.resolution(2)]);
    % Draw fixation cross
    Screen('DrawLine', display.windowPtr, display.fixation.color, xCenter - fixSize, yCenter, xCenter + fixSize, yCenter, 2);
    Screen('DrawLine', display.windowPtr, display.fixation.color, xCenter, yCenter - fixSize, xCenter, yCenter + fixSize, 2);
    % Flip screen if needed
    if display.fixation.flip
        Screen('Flip', display.windowPtr);
    end
end

function runPracticeTrials(display, dots, ifi)
    numPracticeTrials = 10; %Set to 10!
    for trial = 1:numPracticeTrials
        % Show fixation cross for 1,000 ms
        display = drawFixation(display);
        Screen('Flip', display.windowPtr);
        WaitSecs(1);

        % Randomly choose coherence direction
        if rand < 0.5
            dots.coherenceDirection = 'left';
        else
            dots.coherenceDirection = 'right';
        end

        % Run the moving dots stimulus
        [key, RT] = movingDots(display, dots); 

        % Display feedback (e.g., correct/incorrect, but omitted for brevity)
        fprintf('Practice Trial %d: You pressed the "%s" key after %5.2f seconds\n', trial, key, RT);

        % Collect confidence rating
        confidence = getConfidenceRating(display, ifi);
        fprintf('Confidence rating: %d\n', confidence);
    end
end

function runMainExperiment(display, dots, ifi, participantNumber, age, sex)
    coherenceLevels = [0.000, 0.007, 0.013, 0.027, 0.053, 0.107, 0.213, 0.427, 1.000]; 
    numBlocks = length(coherenceLevels);
    trialsPerBlock = 49; % Set to 49!!
    totalTrials = numBlocks * trialsPerBlock;

    coherenceAssignments = repelem(coherenceLevels, totalTrials / length(coherenceLevels));
    coherenceAssignments = coherenceAssignments(randperm(length(coherenceAssignments)));

    results = struct('participantNumber', [], 'trial', [], 'block', [], 'coherence', [], 'direction', [], 'response', [], 'correct', [], 'confidence', [], 'RT', [], 'age', [], 'sex', []);
    
    trialIndex = 1;
    for block = 1:numBlocks
        for trial = 1:trialsPerBlock
            display = drawFixation(display);
            Screen('Flip', display.windowPtr);
            WaitSecs(1);

            dots.coherence = coherenceAssignments(trialIndex);

            if dots.coherence == 0
                dots.coherenceDirection = 'none'; % No coherent direction
            elseif rand < 0.5
                dots.coherenceDirection = 'left';
            else
                dots.coherenceDirection = 'right';
            end

            [key, RT, correctDirection] = movingDots(display, dots);

            % Determine if the response was correct
            isCorrect = strcmpi(key, correctDirection);

            confidence = getConfidenceRating(display, ifi);

            results.participantNumber(trialIndex) = participantNumber;
            results.trial(trialIndex) = trialIndex;
            results.block(trialIndex) = block;
            results.coherence(trialIndex) = dots.coherence;
            results.direction{trialIndex} = correctDirection; 
            results.response{trialIndex} = key;
            results.correct(trialIndex) = isCorrect;
            results.confidence(trialIndex) = confidence;
            results.RT(trialIndex) = RT;
            results.age(trialIndex) = age;
            results.sex{trialIndex} = sex;

            fprintf('Block %d, Trial %d: You pressed the "%s" key after %5.2f seconds\n', block, trial, key, RT);
            fprintf('Confidence rating: %d\n', confidence);

            trialIndex = trialIndex + 1;
        end

        if block < numBlocks
            DrawFormattedText(display.windowPtr, 'You can take a short break. Press any key to continue to the next block.', 'center', 'center', [255, 255, 255]);
            Screen('Flip', display.windowPtr);
            KbStrokeWait;
        end
    end
    
    contingencyTables = constructContingencyTables(results, participantNumber, age, sex);
    save(['P', num2str(participantNumber), '.mat'], 'results', 'contingencyTables');
end

function contingencyTables = constructContingencyTables(results, participantNumber, age, sex)
    % Initialize contingency table structure
    contingencyTables = struct('participantNumber', [], 'trialNumber', [], 'block', [], 'coherence', [], 'Type2hits', [], 'Type2misses', [], 'Type2falseAlarms', [], 'Type2correctRejections', [], ...
        'confidence', [], 'RT', [], 'correct', [], 'response', [], 'direction', [], 'age', [], 'sex', []);

    trialNumbers = (1:length(results.trial))';  % Create a sequence of trial numbers

    % Iterate over each trial and construct the contingency table
    for i = 1:length(trialNumbers)
        idx = trialNumbers(i);

        Type2hits = sum(results.correct(idx) & results.confidence(idx) >= 4);
        Type2misses = sum(results.correct(idx) & results.confidence(idx) < 4);
        Type2falseAlarms = sum(~results.correct(idx) & results.confidence(idx) >= 4);
        Type2correctRejections = sum(~results.correct(idx) & results.confidence(idx) < 4);

        % Assign values to the contingency table structure
        contingencyTables(i).participantNumber = participantNumber;
        contingencyTables(i).trialNumber = idx;
        contingencyTables(i).block = results.block(idx);
        contingencyTables(i).coherence = results.coherence(idx);
        contingencyTables(i).Type2hits = Type2hits;
        contingencyTables(i).Type2misses = Type2misses;
        contingencyTables(i).Type2falseAlarms = Type2falseAlarms;
        contingencyTables(i).Type2correctRejections = Type2correctRejections;
        contingencyTables(i).confidence = results.confidence(idx);
        contingencyTables(i).RT = results.RT(idx);
        contingencyTables(i).correct = results.correct(idx); 
        contingencyTables(i).direction = results.direction{idx};
        contingencyTables(i).response = results.response{idx};
        contingencyTables(i).age = age;
        contingencyTables(i).sex = {sex};
    end

    % Convert the structure to a table
    contingencyTable = struct2table(contingencyTables);
    % Convert 'sex' column to a cell array of strings to ensure 'M' and 'F' are not converted to logical
    contingencyTable.sex = cellfun(@char, contingencyTable.sex, 'UniformOutput', false);
    % Save the table to a CSV file
    writetable(contingencyTable, fullfile('G:\Other computers\My Computer\MSc\Project; PYM0PP(prep) & PYM0EP\Data', ['P', num2str(participantNumber), '_Data.csv']));
end

function confidence = getConfidenceRating(display, ifi)
    % Set the desired font and size
    Screen('TextFont', display.windowPtr, 'Arial');
    Screen('TextSize', display.windowPtr, 44); % Set the font size

    % Likert scale settings
    scaleLengthPix = display.resolution(2) / 1.5;
    scaleHLengthPix = scaleLengthPix / 2;
    xCenter = display.resolution(1) / 2;
    yCenter = display.resolution(2) / 2;
    leftEnd = [xCenter - scaleHLengthPix yCenter];
    rightEnd = [xCenter + scaleHLengthPix yCenter];
    scaleLineCoords = [leftEnd' rightEnd'];
    scaleLineWidth = 10;
    dim = 16;
    hDim = dim / 2;
    numScalePoints = 6;
    xPosScalePoints = linspace(xCenter - scaleHLengthPix, xCenter + scaleHLengthPix, numScalePoints);
    yPosScalePoints = repmat(yCenter, 1, numScalePoints);
    xyScalePoints = [xPosScalePoints; yPosScalePoints];
    sliderLabels = {'Very Low', 'Very High'};

    % Get bounding boxes for the scale end label text
    textBoundsAll = nan(2, 4);
    for i = 1:2
        [~, ~, textBoundsAll(i, :)] = DrawFormattedText(display.windowPtr, sliderLabels{i}, 0, 0, [255, 255, 255]);
    end

    % Width and height of the scale end label text bounding boxs
    textWidths = textBoundsAll(:, 3)';
    halfTextWidths = textWidths / 2;
    textHeights = range([textBoundsAll(:, 2) textBoundsAll(:, 4)], 2)';
    halfTextHeights = textHeights / 2;

    % Position of the scale text
    textPixGap = 50;
    leftTextPosX = xCenter - scaleHLengthPix - hDim - textWidths(1) - textPixGap;
    rightTextPosX = xCenter + scaleHLengthPix + hDim + textPixGap;
    leftTextPosY = yCenter + halfTextHeights(1);
    rightTextPosY = yCenter + halfTextHeights(2);

    % Colors for the likert scale buttons when pressed (blue to red)
    br = linspace(0, 1, numScalePoints);
    bg = zeros(1, numScalePoints);
    bb = abs(1 - br);
    bRGB = [br; bg; bb];

    % Number of frames to wait before updating the screen
    waitframes = 1;

    % Sync us and get a time stamp. We blank the window first to remove the text that we drew to get the bounding boxes.
    Screen('FillRect', display.windowPtr, display.bkColor)
    vbl = Screen('Flip', display.windowPtr);

    % Loop the animation until a key is pressed
    confidence = 0;
    while confidence == 0
        % Get the current position of the mouse
        [mx, my, buttons] = GetMouse(display.windowPtr);

        % Check if the mouse is within any of the circles
        inCircles = sqrt((xPosScalePoints - mx).^2 + (yPosScalePoints - my).^2) < hDim;
        weInCircle = sum(inCircles) > 0;
        if weInCircle == 1
            [~, posCircle] = max(inCircles);
            coordsCircle = xyScalePoints(:, posCircle);
        end

        % Draw the scale line
        Screen('DrawLines', display.windowPtr, scaleLineCoords, scaleLineWidth, [128, 128, 128]);

        % Text for the ends of the slider
        DrawFormattedText(display.windowPtr, sliderLabels{1}, leftTextPosX, leftTextPosY, [0, 0, 255]);
        DrawFormattedText(display.windowPtr, sliderLabels{2}, rightTextPosX, rightTextPosY, [255, 0, 0]);

        % Draw the title for the slider
        DrawFormattedText(display.windowPtr, 'Rate your confidence (1-6):', 'center', display.resolution(2) * 0.25, [255, 255, 255]);

        % If we are in a circle, identify it with a frame
        if weInCircle == 1
            Screen('DrawDots', display.windowPtr, coordsCircle, dim * 1.2, [255, 255, 255], [], 2);
        end

        % Draw the likert scale points
        Screen('DrawDots', display.windowPtr, xyScalePoints, dim, [64, 64, 64], [], 2);

        % If we are clicking a circle, highlight the pressed button
        if weInCircle == 1 && sum(buttons) > 0
            Screen('DrawDots', display.windowPtr, coordsCircle, dim * 1.2, bRGB(:, posCircle), [], 2);
            confidence = posCircle;  % Set confidence to the selected scale point
        end

        % Draw the current rating
        Screen('DrawDots', display.windowPtr, [mx my], 10, [255, 255, 255], [], 2);

        % Flip to the screen
        vbl  = Screen('Flip', display.windowPtr, vbl + (waitframes - 0.5) * ifi);
    end

    % Clear screen
    Screen('FillRect', display.windowPtr, display.bkColor);
    Screen('Flip', display.windowPtr);
end

function displayInstructions(display)
    % Set the desired font and size
    Screen('TextFont', display.windowPtr, 'Arial');
    Screen('TextSize', display.windowPtr, 44); % Set the font size to 44

    instructions = [
        'In this experiment, you will see a field of moving dots.\n\n' ...
        'Your task is to determine whether the coherent motion is to the left or right.\n\n' ...
        'Press the Left Arrow key if you think the motion is to the left.\n\n' ...
        'Press the Right Arrow key if you think the motion is to the right.\n\n' ...
        'After each trial, you will rate your confidence on a scale from 1 to 6.\n\n' ...
        'Press any key to start the practice trials.'
    ];
    DrawFormattedText(display.windowPtr, instructions, 'center', 'center', [255, 255, 255]);
    Screen('Flip', display.windowPtr);
    KbStrokeWait;
end

function gradedImage = generateCircularGabor(width, height, orientation)
    % Create a meshgrid for spatial coordinates
    [x, y] = meshgrid(linspace(-1, 1, width), linspace(-1, 1, height));
    
    % Parameters for the Gabor-like stimulus
    frequency = 5;  % Spatial frequency (cycles per degree), adjust this value
    sigma = 0.3;  % Standard deviation of the Gaussian envelope, adjust this value
    phase = 0;  % Phase offset of the sinusoidal grating
    contrast = 1;  % Contrast of the grating
    
    % Convert orientation to radians
    orientationRad = orientation * pi / 180;
    
    % Generate the oriented Gabor-like stimulus
    xTheta = x * cos(orientationRad) + y * sin(orientationRad);
    yTheta = -x * sin(orientationRad) + y * cos(orientationRad);
    grating = contrast * sin(2 * pi * frequency * xTheta + phase);
    gaussianEnvelope = exp(-(x.^2 + y.^2) / (2 * sigma^2));
    gabor = grating .* gaussianEnvelope;
    
    % Create a circular mask to make the Gabor patch circular
    circularMask = sqrt(x.^2 + y.^2) <= 1; % Define a circular mask
    
    % Apply the circular mask to the Gabor patch
    gabor(~circularMask) = NaN; % Set the values outside the circular area to NaN
    
    % Convert NaNs to a background color (e.g., mid-gray)
    gabor(isnan(gabor)) = -1;
    
    % Convert to 8-bit image format with the circular mask applied
    gradedImage = uint8(255 * (gabor + 1) / 2);
end

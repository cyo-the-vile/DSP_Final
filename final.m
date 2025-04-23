
%%%%%%%%%%%%%% " the quick brown fox jumps over the lazy dog
%%%%%%%%%%%%%%%%%% so far we have 95% accuracy rating

%%%%%%%%% first load all audio

[x_gabe_1, fs] = audioread('gabe_1.wav');
[x_gabe_2, ~] = audioread('gabe_2.wav');
[x_gabe_3, ~] = audioread('gabe_3.wav');
[x_gabe_4, ~] = audioread('gabe_4.wav');
[x_gabe_5, ~] = audioread('gabe_5.wav');  

[x_gracelyn_1, ~] = audioread('gracelyn_1.wav');
[x_gracelyn_2, ~] = audioread('gracelyn_2.wav');
[x_gracelyn_3, ~] = audioread('gracelyn_3.wav');
[x_gracelyn_4, ~] = audioread('gracelyn_4.wav');
[x_gracelyn_5, ~] = audioread('gracelyn_5.wav');

[x_laura_1, ~] = audioread('laura_1.wav');
[x_laura_2, ~] = audioread('laura_2.wav');
[x_laura_3, ~] = audioread('laura_3.wav');
[x_laura_4, ~] = audioread('laura_4.wav');
[x_laura_5, ~] = audioread('laura_5.wav');

[x_lydia_1, ~] = audioread('lydia_1.wav');
[x_lydia_2, ~] = audioread('lydia_2.wav');
[x_lydia_3, ~] = audioread('lydia_3.wav');
[x_lydia_4, ~] = audioread('lydia_4.wav');
[x_lydia_5, ~] = audioread('lydia_5.wav');

[x_mark_1, ~] = audioread('mark_1.wav');
[x_mark_2, ~] = audioread('mark_2.wav');
[x_mark_3, ~] = audioread('mark_3.wav');
[x_mark_4, ~] = audioread('mark_4.wav');
[x_mark_5, ~] = audioread('mark_5.wav');

[x_micah_1, ~] = audioread('micah_1.wav');
[x_micah_2, ~] = audioread('micah_2.wav');
[x_micah_3, ~] = audioread('micah_3.wav');
[x_micah_4, ~] = audioread('micah_4.wav');
[x_micah_5, ~] = audioread('micah_5.wav');

[x_sana_1, ~] = audioread('sana_1.wav');
[x_sana_2, ~] = audioread('sana_2.wav');
[x_sana_3, ~] = audioread('sana_3.wav');
[x_sana_4, ~] = audioread('sana_4.wav');
[x_sana_5, ~] = audioread('sana_5.wav');

[x_ben_1, ~] = audioread('ben_1.wav');
[x_ben_2, ~] = audioread('ben_2.wav');
[x_ben_3, ~] = audioread('ben_3.wav');
[x_ben_4, ~] = audioread('ben_4.wav');
[x_ben_5, ~] = audioread('ben_5.wav');

[x_blessing_1, ~] = audioread('blessing_1.wav');
[x_blessing_2, ~] = audioread('blessing_2.wav');
[x_blessing_3, ~] = audioread('blessing_3.wav');
[x_blessing_4, ~] = audioread('blessing_4.wav');
[x_blessing_5, ~] = audioread('blessing_5.wav');

[x_dj_1, ~] = audioread('dj_1.wav');
[x_dj_2, ~] = audioread('dj_2.wav');
[x_dj_3, ~] = audioread('dj_3.wav');
[x_dj_4, ~] = audioread('dj_4.wav');
[x_dj_5, ~] = audioread('dj_5.wav');

%%%%%%%%%% intruder audio. we need to reject these.

[x_silvario_1, ~] = audioread('silvario_1.wav');
[x_silvario_2, ~] = audioread('silvario_2.wav');
[x_silvario_3, ~] = audioread('silvario_3.wav');

[x_crawley_1, ~] = audioread('crawley_1.wav');
[x_crawley_2, ~] = audioread('crawley_2.wav');
[x_crawley_3, ~] = audioread('crawley_3.wav');

[x_holden_1, ~] = audioread('holden_1.wav');
[x_holden_2, ~] = audioread('holden_2.wav');
[x_holden_3, ~] = audioread('holden_3.wav');

[x_tim_1, ~] = audioread('tim_1.wav');
[x_tim_2, ~] = audioread('tim_2.wav');
[x_tim_3, ~] = audioread('tim_3.wav');

[x_kris_1, ~] = audioread('kris_1.wav');
[x_kris_2, ~] = audioread('kris_2.wav');
[x_kris_3, ~] = audioread('kris_3.wav');

[noise_background, ~] = audioread('empty_room.wav');

%%%%%%%%% batch process the matrixes from stereo to mono.  all the matrices
%%%%%%%%% are in a stereo format and mono closely represents how we
%%%%%%%%% physically speak.  ie, we are not speaking from two mouths

vars = whos; 

for i = 1:length(vars)
    name = vars(i).name;

    % all initial matrixes basically will be grabbed. that's why we went
    % with the naming convention
    if startsWith(name, 'x_')
        data = eval(name);

        if ismatrix(data) && size(data, 2) == 2
            mono = mean(data, 2);  
            assignin('base', name, mono); 
        end
    end
end

%%%%%%%%%% background noise to mono for our first filter
if size(noise_background, 2) == 2
    noise_background = mean(noise_background, 2);
end
noise_background = noise_background / max(abs(noise_background));

%%%%%%%%%%% filter background noise.

N = length(noise_background);
Y = abs(fft(noise_background));
Y = Y(1:floor(N/2));
f = linspace(0, fs/2, length(Y));

[~, idx] = max(Y);      
f_noise = f(idx);           % we have 28.37hz


%%%%%%%%%%% attempting to run an fir filter with the matlab function
%%%%%%%%%%% fir1.  
%%%%%%%%%%%  https://www.mathworks.com/help/signal/ref/fir1.html

bandwidth = 3; %%%%%%%%%%% parameter for a tolerance 

f_low = (f_noise - bandwidth/2) / (fs/2); 
f_high = (f_noise + bandwidth/2) / (fs/2);

f_low = max(f_low, 0);
f_high = min(f_high, 1);
filter_order = 128;              %%%%%%%%%% amount of coefficients is actually filter_order+1. so 129

b = fir1(filter_order, [f_low f_high], 'stop');  

%%%%%%%%%% clean everything now with FIR using another loop similar to
%%%%%%%%%% earlier one.

vars = whos;

for i = 1:length(vars)
    name = vars(i).name;

    if startsWith(name, 'x_') && ~contains(name, '_clean')
        x = eval(name);

        if size(x, 2) == 2
            x = mean(x, 2);
        end

        x = x(:);
        x = x / max(abs(x));

        x_clean = filter(b, 1, x);
        x_clean = x_clean / max(abs(x_clean));

        clean_name = [name '_clean'];
        assignin('base', clean_name, x_clean);
    end
end


%%%%%%%%%%%%%%%  grabs magnitude.  we can't do voice recognition on
%%%%%%%%%%%%%%%  frequency alone unfortunately.  I spent a significant
%%%%%%%%%%%%%%%  amount of time concluding this.
features = struct(); 
vars = whos;

for i = 1:length(vars)
    name = vars(i).name;
    if startsWith(name, 'x_') && endsWith(name, '_clean') %%%%%% grabs the cleaned up ones and doesn't mess up golden copies
        x = eval(name);
        N = 2^nextpow2(length(x));

        X = abs(fft(x, N));
        X = X(1:floor(N/2));    
        X = X / max(X);          
        features.(name) = X;
    end
end

%%%%%%%%% so next we want to try averaging the cleaned recordings together.
%%%%%%%%%  This took a bit longer specifically on finding a way to traverse
%%%%%%%%%  everything.

%%%%%%%%%%%%%%%%%%%% these are the people we want to take.
speakers = {'gabe', 'gracelyn', 'laura', 'lydia', 'mark', 'micah', 'sana',... 
    'ben', 'blessing', 'dj'};

speaker_profiles = struct();

for i = 1:length(speakers)
    spk = speakers{i};
    fft_list = [];

    max_len = 0;
    temp_ffts = {};

    for j = 1:5
        key = sprintf('x_%s_%d_clean', spk, j);
        if isfield(features, key)
            temp = features.(key);
            max_len = max(max_len, length(temp));
            temp_ffts{end+1} = temp;
        end
    end

    %%%%%%%%%%% had to pad with zeros and the reason for this is that
    %%%%%%%%%%% all the audio lengths are different from one another
    for k = 1:length(temp_ffts)
        vec = temp_ffts{k}(:)';     %%%%%%%% matrices werent matching
        padded = [vec, zeros(1, max_len - length(vec))];
        fft_list(end+1, :) = padded;  % Now dimensions match
    end
    speaker_profiles.(spk) = mean(fft_list, 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%   master console input %%%%%%%%%%%%


fprintf('\n Password required! \n');

while true
    test_name = input('Enter cleaned test sample name: ', 's');

    if strcmpi(test_name, 'exit')
        disp('exiting');
        break;
    end

    if ~isfield(features, test_name)
        fprintf(' audio recording "%s" not found. \n\n', test_name);
        continue;
    end

    test_fft = features.(test_name);

    % Pad test_fft if needed with zeros
    test_fft = [test_fft(:)', zeros(1, max_len - length(test_fft))];

    
    min_dist = inf;
    best_match = '';
    speaker_names = fieldnames(speaker_profiles);

    for i = 1:length(speaker_names)
        spk = speaker_names{i};
        profile = speaker_profiles.(spk);

        % Padding with zeros again
        profile = [profile(:)', zeros(1, max_len - length(profile))];
%%%%%%%%%% this is very important we are using this value to calculate amplitude distance. Very important.
        dist = norm(test_fft - profile); 

        if dist < min_dist
            min_dist = dist;
            best_match = spk;
        end
    end

    % Threshold check
    threshold = 8.5;

    if min_dist > threshold
        fprintf(' Not Authorized \n\n');
    %    fprintf('DEBUG - Closest Match: %s\n', best_match);
    %    fprintf('DEBUG - Distance to Match: %.4f\n', min_dist);

    else
        fprintf(' Authorized. You may enter \n\n');
   %     fprintf('DEBUG - Closest Match: %s\n', best_match);
   %     fprintf('DEBUG - Distance to Match: %.4f\n', min_dist);

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% this was used for debugging and not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% required for functionality

fprintf('\n Cleaned Samples vs Speaker Profiles\n');

sample_names = fieldnames(features);
profile_names = fieldnames(speaker_profiles);
%max_len = max([structfun(@length, speaker_profiles), structfun(@length, features)]);
all_lengths = [structfun(@length, speaker_profiles); structfun(@length, features)];
max_len = max(all_lengths);

%max_len = max(structfun(@length, speaker_profiles)); 

threshold = 8.5;  %%%%%%%%%%%%%% this was tweaked until we have our confidence range.

for i = 1:length(sample_names)
    sample = sample_names{i};
    test_fft = features.(sample);
    test_fft = [test_fft(:)', zeros(1, max_len - length(test_fft))];

    min_dist = inf;
    best_match = '';

    for j = 1:length(profile_names)
        profile = speaker_profiles.(profile_names{j});
        profile = [profile(:)', zeros(1, max_len - length(profile))];

        dist = norm(test_fft - profile);

        if dist < min_dist
            min_dist = dist;
            best_match = profile_names{j};
        end
    end
    if min_dist > threshold
        result = 'INTRUDER';
    else
        result = 'AUTHORIZED';
    end

    fprintf('%-25s Closest Match: %-10s amplitude: %.4f %s\n', sample, best_match, min_dist, result);
end

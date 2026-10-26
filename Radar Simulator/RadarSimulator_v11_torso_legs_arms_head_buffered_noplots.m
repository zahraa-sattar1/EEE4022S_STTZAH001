%% Walking radar simulator
% Zahraa Sattar | EEE4022S
%
% Generate 24GHz CW radar returns from the moving torso, legs, arms and head.
% Joint motion follows the walking model, with stationary buffers and adjustable noise.

clear all;
close all;

%% Optional settings
if ~exist('AmplitudeMode','var')
    AmplitudeMode = "A";
end

if ~exist('SimulatorR0_m','var')
    SimulatorR0_m = 10;
end

if ~exist('SimulatorVerbose','var')
    SimulatorVerbose = true;
end

if ~exist('SimulatorMakePlots','var')
    SimulatorMakePlots = true;
end

AmplitudeMode = string(AmplitudeMode);

assert(isscalar(AmplitudeMode) && ...
    any(AmplitudeMode == ["A","B","Old"]), ...
    'AmplitudeMode must be "A", "B", or "Old".');

validateattributes(SimulatorR0_m, {'numeric'}, ...
    {'real','finite','scalar','positive'});

%% Physical Constants
c = 3e8;                        % speed of light, m/s

%% Display Settings
maxSpeed_m_s = 10; % (m/s) maximum speed to display
minSpeed_m_s = -2;

%% Radar Parameters
fc = 24e9;                       % carrier frequency

RadarHeight_m = 0.75;              % radar height, m
RadarCord_ZX = [RadarHeight_m, 0]; % Coordinate order: [Z, X]

%% Human: base walking geometry
R0 = SimulatorR0_m;

if ~exist('v0','var')
    v0 = 1.1;
end

validateattributes(v0, {'numeric'}, ...
    {'real','finite','scalar','positive'});

A  = 1;                          % return signal amplitude (arbitrary units, same for every scatterer)

%% Human - Boulic Walking Model Parameters (Appendix A)
Height_m = 1.7;                 % subject height, m
H_Spine = 0.58*Height_m;
H_t = 0.53*Height_m;            % hip height, m
L_torso = 0.288*Height_m;       % torso length
L_upperLeg = 0.245*Height_m;              % Length of Upper Leg
L_lowerLeg = 0.246*Height_m;              % Length of Lower Leg
L_legHeight = L_upperLeg + L_lowerLeg;    % Length of Leg
L_foot = 0.143*Height_m;                  % Length of foot, ankle to toe
L_upperArm = 0.188*Height_m;            % Length of upper arm, shoulder to elbow
L_forearm = 0.145*Height_m;             % Length of forearm, elbow to wrist
L_head = 0.13*Height_m;         % length of head segment
Vr = v0/L_legHeight;            % Relative velocity normalised by the leg height
RLc = 1.346*sqrt(Vr);           % relative cycle length, Inman et al. (1981), Eq. 1
Tc  = RLc/Vr;                   % gait-cycle duration, s
Ts  = 0.752*Tc - 0.143;         % stance duration, s, Murray et al. (1964)
Ts_pct = Ts/Tc;                 % relative stance duration

%% Simulation Settings
fs = 44.1e3;                     % sample rate, Hz
SimulationTime_s = 4;            % total length of the simulated signal, s

if v0 < 1
    SimulationTime_s = 8;              % extend the walk at slow speeds to include more cycles
end

bufferTime_s = 0.5;              % stationary time before and after walking, s
N = round((SimulationTime_s + 2*bufferTime_s)*fs);  % samples including both buffers
if ~exist('SNR_dB','var')
    SNR_dB = 40;
end

%% Processing parameters

if exist('CPI_s','var')
    CPI = CPI_s;
else
    CPI = 0.025;
end

CPI_requested_s_Sim = CPI;

if v0 < 1.0                      % retain the low-speed simulator CPI override
    CPI = 0.075;
end

NumSamplesPerFrame = 2^(nextpow2(round(CPI*fs)));       % actual window length after rounding up to a power of 2

CPI_actual_s_Sim = NumSamplesPerFrame/fs;

assert(NumSamplesPerFrame <= N, ...
    'CPI exceeds the recording duration.');

%% Time Vector and Normalised Gait Time
lamda = c/fc;       % radar wavelength, m
t = (0:1:(N-1))*1/fs;           % real time vector, seconds
t_walk = min(max(t - bufferTime_s, 0), SimulationTime_s);   % hold the walking clock at each buffer boundary
t_pct = mod(t_walk/Tc,1);            % Relative gait phase in the interval 0 <= tR < 1
%% Spine Movement Calculations
Av = 0.015*Vr;                  % vertical bob amplitude, units of Ht
Spine_Z = H_Spine + ( -Av + Av*sin(2*pi*(2*t_pct - 0.35)) );    % Vertical translation of spine
if Vr > 0.5
    Aa = -0.021;                 % fore/aft amplitude, units of Ht
else
    Aa = -0.084*Vr^2 + 0.084*Vr;
end
phi_a  = 0.625 - Ts_pct;         % phase term, from the paper's forward/backward translation expression
x_offset = Aa*sin(2*pi*(2*t_pct + 2*phi_a));
Spine_X = R0 - v0*t_walk + x_offset;     % Spine's forward position (changed to + x_offset)
if Vr > 0.5
    Ar  = 2;
else
    Ar = -8*Vr^2 + 8*Vr;
end
theta_r = -Ar + Ar*sin(2*pi*(2*t_pct - 0.1));   % pelvis forward/backward rotation angle, in degrees
Torso_X = Spine_X + (L_torso/2)*sin(theta_r*pi/180);
Torso_Z = Spine_Z + (L_torso/2)*cos(theta_r*pi/180);

%% Torso Radar Return Signal
Torso_R = sqrt((RadarCord_ZX(2) - Torso_X).^2 + (RadarCord_ZX(1) - Torso_Z).^2);    % range from radar to torso point

%% Flexing at the Hip (Chen Appendix 4A / Boulic Appendix B.1)
Hip_Z = H_t + ( -Av + Av*sin(2*pi*(2*t_pct - 0.35)) );   % vertical position of the hip joint
Hip_X = R0 - v0*t_walk + x_offset;                            % fore/aft position of the hip joint
if Vr < 0.5
    gaitCharacteristic = 'a';
elseif Vr < 1.3
    gaitCharacteristic = 'b';
else
    gaitCharacteristic = 'c';
end
theta_deg = computeHipFlexionDeg(Vr,gaitCharacteristic,t_pct);   % hip flexion angle, degrees, one value per t_pct sample
theta_rad = theta_deg*pi/180;                                     % convert to radians for sin/cos below
L_thigh_midpoint = L_upperLeg/2;                 % lever arm, hip to thigh midpoint
Thigh_X = Hip_X - L_thigh_midpoint*sin(theta_rad);
Thigh_Z = Hip_Z - L_thigh_midpoint*cos(theta_rad);

%% Thigh Radar Return Signal
Thigh_R = sqrt((RadarCord_ZX(2) - Thigh_X).^2 + (RadarCord_ZX(1) - Thigh_Z).^2);    % range from radar to thigh midpoint, m

%% Knee Movement Calculations (Chen Appendix 4A / Boulic Appendix B.2)
Knee_X = Hip_X - L_upperLeg*sin(theta_rad);
Knee_Z = Hip_Z - L_upperLeg*cos(theta_rad);
theta_knee_deg = computeKneeFlexionDeg(Vr,gaitCharacteristic,t_pct);   % how much the shank bends back relative to the thigh, degrees
theta_knee_rad = theta_knee_deg*pi/180;
theta_shank_rad = theta_rad - theta_knee_rad;       % shank angle relative to the thigh
L_shank_midpoint = L_lowerLeg/2;
Shank_X = Knee_X - L_shank_midpoint*sin(theta_shank_rad);
Shank_Z = Knee_Z - L_shank_midpoint*cos(theta_shank_rad);

%% Shank Radar Return Signal
Shank_R = sqrt((RadarCord_ZX(2) - Shank_X).^2 + (RadarCord_ZX(1) - Shank_Z).^2);

%% Ankle Movement Calculations (Chen Appendix 4A / Boulic Appendix B.3)
Ankle_X = Knee_X - L_lowerLeg*sin(theta_shank_rad);
Ankle_Z = Knee_Z - L_lowerLeg*cos(theta_shank_rad);
theta_ankle_deg = computeAnkleFlexionDeg(Vr,gaitCharacteristic,Ts_pct,t_pct);   % ankle flexion angle, degrees
theta_ankle_rad = theta_ankle_deg*pi/180;
theta_foot_rad = theta_rad - theta_knee_rad + theta_ankle_rad;  % Relative foot angle

Toe_X = Ankle_X - L_foot*cos(theta_foot_rad);
Toe_Z = Ankle_Z + L_foot*sin(theta_foot_rad);

Foot_X = (Ankle_X + Toe_X)/2;
Foot_Z = (Ankle_Z + Toe_Z)/2;

%% Foot Radar Return Signal
Foot_R = sqrt((RadarCord_ZX(2) - Foot_X).^2 + (RadarCord_ZX(1) - Foot_Z).^2);

%% Right Leg (Chen Appendix 4A / Boulic Appendix B.1-B.3)
t_pct_right = mod(t_pct + 0.5, 1);   % right-leg phase, half a cycle from the left
RHip_X = Hip_X;
RHip_Z = Hip_Z;
RHipFlex_deg = computeHipFlexionDeg(Vr,gaitCharacteristic,t_pct_right);     % right hip flexion angle, degrees
RHipFlex_rad = RHipFlex_deg*pi/180;
RKneeFlex_deg = computeKneeFlexionDeg(Vr,gaitCharacteristic,t_pct_right);   % right knee flexion angle, degrees
RKneeFlex_rad = RKneeFlex_deg*pi/180;
RAnkleFlex_deg = computeAnkleFlexionDeg(Vr,gaitCharacteristic,Ts_pct,t_pct_right);  % right ankle flexion angle, degrees
RAnkleFlex_rad = RAnkleFlex_deg*pi/180;
theta_Rthigh_rad = RHipFlex_rad;                                  % right thigh angle
theta_Rshank_rad = RHipFlex_rad - RKneeFlex_rad;                  % right shank angle
theta_Rfoot_rad  = RHipFlex_rad - RKneeFlex_rad + RAnkleFlex_rad; % right foot angle
RKnee_X = RHip_X - L_upperLeg*sin(theta_Rthigh_rad);
RKnee_Z = RHip_Z - L_upperLeg*cos(theta_Rthigh_rad);
RAnkle_X = RKnee_X - L_lowerLeg*sin(theta_Rshank_rad);
RAnkle_Z = RKnee_Z - L_lowerLeg*cos(theta_Rshank_rad);
RToe_X = RAnkle_X - L_foot*cos(theta_Rfoot_rad);
RToe_Z = RAnkle_Z + L_foot*sin(theta_Rfoot_rad);
RThigh_X = (RHip_X + RKnee_X)/2;
RThigh_Z = (RHip_Z + RKnee_Z)/2;
RShank_X = (RKnee_X + RAnkle_X)/2;
RShank_Z = (RKnee_Z + RAnkle_Z)/2;
RFoot_X  = (RAnkle_X + RToe_X)/2;
RFoot_Z  = (RAnkle_Z + RToe_Z)/2;

%% Right Leg Radar Return Signals
RThigh_R = sqrt((RadarCord_ZX(2) - RThigh_X).^2 + (RadarCord_ZX(1) - RThigh_Z).^2);
RShank_R = sqrt((RadarCord_ZX(2) - RShank_X).^2 + (RadarCord_ZX(1) - RShank_Z).^2);
RFoot_R  = sqrt((RadarCord_ZX(2) - RFoot_X).^2  + (RadarCord_ZX(1) - RFoot_Z).^2);

%% Left Arm (Boulic Appendix D.2/D.3)
Shoulder_X = Spine_X + L_torso*sin(theta_r*pi/180);   % shoulder position at the top of the torso segment
Shoulder_Z = Spine_Z + L_torso*cos(theta_r*pi/180);
As = 9.88*Vr;   % shoulder swing amplitude shared by both arms
LShoulderFlex_deg = 3  - As*cos(2*pi*t_pct_right);   % left arm, opposite of right foot
LShoulderFlex_rad = LShoulderFlex_deg*pi/180;
LUpperArmAngle_rad = LShoulderFlex_rad - theta_r*pi/180;   % torso's own lean subtracted out
LElbowFlex_deg = computeElbowFlexionDeg(Vr,gaitCharacteristic,t_pct_right);   % left elbow, same phase as left shoulder
LElbowFlex_rad = LElbowFlex_deg*pi/180;
LForeArmAngle_rad = LUpperArmAngle_rad + LElbowFlex_rad;
LElbow_X = Shoulder_X - L_upperArm*sin(LUpperArmAngle_rad);
LElbow_Z = Shoulder_Z - L_upperArm*cos(LUpperArmAngle_rad);
LWrist_X = LElbow_X - L_forearm*sin(LForeArmAngle_rad);
LWrist_Z = LElbow_Z - L_forearm*cos(LForeArmAngle_rad);
LUpperArm_X = (Shoulder_X + LElbow_X)/2;
LUpperArm_Z = (Shoulder_Z + LElbow_Z)/2;
LForeArm_X  = (LElbow_X + LWrist_X)/2;
LForeArm_Z  = (LElbow_Z + LWrist_Z)/2;

%% Left Arm Radar Return Signals
LUpperArm_R = sqrt((RadarCord_ZX(2) - LUpperArm_X).^2 + (RadarCord_ZX(1) - LUpperArm_Z).^2);
LForeArm_R  = sqrt((RadarCord_ZX(2) - LForeArm_X).^2  + (RadarCord_ZX(1) - LForeArm_Z).^2);

%% Right Arm (Boulic Appendix D.2/D.3)
RShoulderFlex_deg = 3  - As*cos(2*pi*t_pct);   % right arm, opposite of left foot
RShoulderFlex_rad = RShoulderFlex_deg*pi/180;
RUpperArmAngle_rad = RShoulderFlex_rad - theta_r*pi/180;   % torso's own lean subtracted out
RElbowFlex_deg = computeElbowFlexionDeg(Vr,gaitCharacteristic,t_pct);   % right elbow, using the right-shoulder phase
RElbowFlex_rad = RElbowFlex_deg*pi/180;
RForeArmAngle_rad = RUpperArmAngle_rad + RElbowFlex_rad;
RElbow_X = Shoulder_X - L_upperArm*sin(RUpperArmAngle_rad);
RElbow_Z = Shoulder_Z - L_upperArm*cos(RUpperArmAngle_rad);
RWrist_X = RElbow_X - L_forearm*sin(RForeArmAngle_rad);
RWrist_Z = RElbow_Z - L_forearm*cos(RForeArmAngle_rad);
RUpperArm_X = (Shoulder_X + RElbow_X)/2;
RUpperArm_Z = (Shoulder_Z + RElbow_Z)/2;
RForeArm_X  = (RElbow_X + RWrist_X)/2;
RForeArm_Z  = (RElbow_Z + RWrist_Z)/2;

%% Right Arm Radar Return Signals
RUpperArm_R = sqrt((RadarCord_ZX(2) - RUpperArm_X).^2 + (RadarCord_ZX(1) - RUpperArm_Z).^2);
RForeArm_R  = sqrt((RadarCord_ZX(2) - RForeArm_X).^2  + (RadarCord_ZX(1) - RForeArm_Z).^2);

%% HEAD (Boulic Appendix D)
Head_X = Spine_X + (L_torso + L_head)*sin(theta_r*pi/180);
Head_Z = Spine_Z + (L_torso + L_head)*cos(theta_r*pi/180);

%% Head Radar Return Signal % NEW
Head_R = sqrt((RadarCord_ZX(2) - Head_X).^2 + (RadarCord_ZX(1) - Head_Z).^2);

%% Toe Radar Return Signal (the actual tip of the toe, not the ankle-toe midpoint)
Toe_R = sqrt((RadarCord_ZX(2) - Toe_X).^2 + (RadarCord_ZX(1) - Toe_Z).^2);
RToe_R = sqrt((RadarCord_ZX(2) - RToe_X).^2 + (RadarCord_ZX(1) - RToe_Z).^2);

%% Calculate RCS from the moving segment geometry

AmplitudeNames = ...
    ["torso","head","thigh","shank","foot","uarm","farm"];

SegmentNames = ...
    ["torso","head","Lthigh","Lshank","Lfoot", ...
     "Luarm","Lfarm","Rthigh","Rshank","Rfoot","Ruarm","Rfarm"];

SegmentRadii_m = [0.15,0.10,0.07,0.06,0.05,0.06,0.05];

SegmentType = [1,2,3,4,5,6,7,3,4,5,6,7];

SegmentEnds = {
    Spine_X,    Spine_Z,    Shoulder_X, Shoulder_Z
    Shoulder_X, Shoulder_Z, Head_X,     Head_Z
    Hip_X,      Hip_Z,      Knee_X,     Knee_Z
    Knee_X,     Knee_Z,     Ankle_X,    Ankle_Z
    Ankle_X,    Ankle_Z,    Toe_X,      Toe_Z
    Shoulder_X, Shoulder_Z, LElbow_X,   LElbow_Z
    LElbow_X,   LElbow_Z,   LWrist_X,   LWrist_Z
    RHip_X,     RHip_Z,     RKnee_X,    RKnee_Z
    RKnee_X,    RKnee_Z,    RAnkle_X,   RAnkle_Z
    RAnkle_X,   RAnkle_Z,   RToe_X,     RToe_Z
    Shoulder_X, Shoulder_Z, RElbow_X,   RElbow_Z
    RElbow_X,   RElbow_Z,   RWrist_X,   RWrist_Z
    };

WalkMask = t > bufferTime_s & ...
    t < bufferTime_s + SimulationTime_s;

assert(any(WalkMask), 'No walking samples found.');

SegmentRCS_m2 = cell(12,1);
SegmentSqrtRCS = cell(12,1);

MeanSegmentRCS_m2 = zeros(12,1);

for segIdx = 1:12

    SegmentRCS_m2{segIdx} = v11EllipsoidRCS( ...
        SegmentEnds{segIdx,1}, ...
        SegmentEnds{segIdx,2}, ...
        SegmentEnds{segIdx,3}, ...
        SegmentEnds{segIdx,4}, ...
        SegmentRadii_m(SegmentType(segIdx)), ...
        RadarCord_ZX);

    SegmentSqrtRCS{segIdx} = sqrt(SegmentRCS_m2{segIdx});

    MeanSegmentRCS_m2(segIdx) = ...
        mean(SegmentRCS_m2{segIdx}(WalkMask));
end

TorsoMeanSqrtRCS = mean(SegmentSqrtRCS{1}(WalkMask));

assert(isfinite(TorsoMeanSqrtRCS) && TorsoMeanSqrtRCS > 0, ...
    'Invalid torso amplitude normalisation.');

RCSAmplitudeWeights = zeros(12,1);

for segIdx = 1:12
    RCSAmplitudeWeights(segIdx) = ...
        mean(SegmentSqrtRCS{segIdx}(WalkMask))/TorsoMeanSqrtRCS;
end

AmplitudeWeights = RCSAmplitudeWeights(1:7).';

OldAmplitudeWeights = ...
    [1,0.5,0.45,0.35,0.25,0.3,0.2, ...
     0.45,0.35,0.25,0.3,0.2];

AppliedAmplitude = cell(12,1);
AppliedAmplitudeMean = zeros(12,1);

for segIdx = 1:12

    switch AmplitudeMode

        case "A"
            AppliedAmplitude{segIdx} = ...
                AmplitudeWeights(SegmentType(segIdx));

            AppliedAmplitudeMean(segIdx) = ...
                AppliedAmplitude{segIdx};

        case "B"
            AppliedAmplitude{segIdx} = ...
                SegmentSqrtRCS{segIdx}/TorsoMeanSqrtRCS;

            AppliedAmplitudeMean(segIdx) = ...
                RCSAmplitudeWeights(segIdx);

        case "Old"
            AppliedAmplitude{segIdx} = OldAmplitudeWeights(segIdx);

            AppliedAmplitudeMean(segIdx) = ...
                OldAmplitudeWeights(segIdx);
    end
end

AmplitudeSummary = table( ...
    SegmentNames(:), ...
    MeanSegmentRCS_m2, ...
    RCSAmplitudeWeights, ...
    AppliedAmplitudeMean, ...
    'VariableNames', ...
    {'Segment','MeanRCS_m2','RCSWeightMean','AppliedWeightMean'});

%% Amplitude variables for compatibility
A = AppliedAmplitude{1};
A_head = AppliedAmplitude{2};

A_thigh = AppliedAmplitude{3};
A_shank = AppliedAmplitude{4};
A_foot = AppliedAmplitude{5};

A_upperarm = AppliedAmplitude{6};
A_forearm = AppliedAmplitude{7};

A_Rthigh = AppliedAmplitude{8};
A_Rshank = AppliedAmplitude{9};
A_Rfoot = AppliedAmplitude{10};

A_Rupperarm = AppliedAmplitude{11};
A_Rforearm = AppliedAmplitude{12};

%% Noise
win = hamming(NumSamplesPerFrame); 
processingGain = sum(win)^2 / sum(win.^2); % retained window-gain diagnostic, not used in noisePower

noisePower = 10^(-(SNR_dB - 10*log10(NumSamplesPerFrame))/10); % - compute noise power
noise = sqrt(noisePower/2)*(randn(1,N) + 1i*randn(1,N));       % independent Gaussian I/Q noise samples

%% Combined radar return signals

%% Clean radar returns
RdrRxSignal_Torso = A .* exp(-1i*4*pi*fc*Torso_R/c);

RdrRxSignal_Head = A_head .* exp(-1i*4*pi*fc*Head_R/c);

RdrRxSignal_Thigh = A_thigh .* exp(-1i*4*pi*fc*Thigh_R/c);

RdrRxSignal_Shank = A_shank .* exp(-1i*4*pi*fc*Shank_R/c);

RdrRxSignal_Foot = A_foot .* exp(-1i*4*pi*fc*Foot_R/c);

RdrRxSignal_RThigh = A_Rthigh .* exp(-1i*4*pi*fc*RThigh_R/c);

RdrRxSignal_RShank = A_Rshank .* exp(-1i*4*pi*fc*RShank_R/c);

RdrRxSignal_RFoot = A_Rfoot .* exp(-1i*4*pi*fc*RFoot_R/c);

RdrRxSignal_LUpperArm = A_upperarm .* exp(-1i*4*pi*fc*LUpperArm_R/c);

RdrRxSignal_LForeArm = A_forearm .* exp(-1i*4*pi*fc*LForeArm_R/c);

RdrRxSignal_RUpperArm = A_Rupperarm .* exp(-1i*4*pi*fc*RUpperArm_R/c);

RdrRxSignal_RForeArm = A_Rforearm .* exp(-1i*4*pi*fc*RForeArm_R/c);

RdrRxSignal_Toe = A_foot .* exp(-1i*4*pi*fc*Toe_R/c);

RdrRxSignal_RToe = A_Rfoot .* exp(-1i*4*pi*fc*RToe_R/c);

%% Amplitude-independent toe references
RdrRxSignal_Toe_GT = 1*exp(-1i*4*pi*fc*Toe_R/c);

RdrRxSignal_RToe_GT = 1*exp(-1i*4*pi*fc*RToe_R/c);

RdrRxSignal_clean = ...
    RdrRxSignal_Torso + ...
    RdrRxSignal_Thigh + ...
    RdrRxSignal_Shank + ...
    RdrRxSignal_Foot + ...
    RdrRxSignal_RThigh + ...
    RdrRxSignal_RShank + ...
    RdrRxSignal_RFoot + ...
    RdrRxSignal_LUpperArm + ...
    RdrRxSignal_LForeArm + ...
    RdrRxSignal_RUpperArm + ...
    RdrRxSignal_RForeArm + ...
    RdrRxSignal_Head;

RdrRxSignal = RdrRxSignal_clean + noise;

MeasuredSampleSNR_dB = 10*log10(mean(abs(RdrRxSignal_clean).^2)/mean(abs(noise).^2));

%% Spectrogram Settings
OverlapFactor = 0.8;
noverlap = floor(OverlapFactor*NumSamplesPerFrame);
nfft = NumSamplesPerFrame;
[S, FrequencyAxis, TimeAxis] = Zahraa_Spectrogram(RdrRxSignal, win, noverlap, nfft, fs);  % shared spectrogram helper
VelocityAxis = FrequencyAxis*lamda/2;

%% Crop To Speed Range Of Interest
speed_m_per_sec_Idx = find((VelocityAxis <= maxSpeed_m_s) & (VelocityAxis >= minSpeed_m_s));
SpeedVectorOfInterest = VelocityAxis(speed_m_per_sec_Idx);
S_OfInterest = S(speed_m_per_sec_Idx, :);
S_OfInterestToPlot = abs(S_OfInterest)/max(max(abs(S_OfInterest)));

%% Plot: Spectrogram

dV = H_Spine*abs(Aa)*(4*pi/Tc);

%% Theoretical Torso Speed Prediction
fprintf('Cycle duration Dc = %.3f s, torso oscillation freq = %.3f Hz\n', Tc, 2/Tc);
fprintf('Predicted Vr range: [%.3f, %.3f] m/s  (v0 = %.3f, dV = %.3f)\n', v0-dV, v0+dV, v0, dV);

%% Plot: X-Z Position History
if SimulatorMakePlots
    figure;
    plot(Torso_X, Torso_Z, '.');
    hold on;
    plot(Thigh_X, Thigh_Z, 'g.');
    plot(Shank_X, Shank_Z, 'm.');
    plot(Foot_X, Foot_Z, 'k.');
    plot(RThigh_X, RThigh_Z, 'm.');   % right thigh
    plot(RShank_X, RShank_Z, 'g.');   % right shank
    plot(RFoot_X, RFoot_Z, 'r.');     % right foot
    plot(LUpperArm_X, LUpperArm_Z, 'c.');
    plot(LForeArm_X, LForeArm_Z, 'b.');
    plot(RUpperArm_X, RUpperArm_Z, 'w.');
    plot(RForeArm_X, RForeArm_Z, 'g.');
    plot(Head_X, Head_Z, 'y.');     % head
    plot(RadarCord_ZX(2), RadarCord_ZX(1), 'r^', 'MarkerFaceColor', 'r');
    axis xy;
    ylim([0 5]);
    xlabel('X (m)');
    ylabel('Z (m)');
    title('CW spectrogram of Full Body (Torso + Legs + Arms + Head)');
    legend('Torso centre','L Thigh','L Shank','L Foot','R Thigh','R Shank','R Foot','L UpperArm','L ForeArm','R UpperArm','R ForeArm','Head','radar');
    grid on;
end

%% Ground Truth Verifications

%% Ground Truth Verifications
CPI_groundTruth_s = 0.025;   % fixed reference window, independent of the swept detector CPI
NumSamplesPerFrame_gt = 2^(nextpow2(round(CPI_groundTruth_s*fs)));
win_gt = hamming(NumSamplesPerFrame_gt);
noverlap_gt = floor(OverlapFactor*NumSamplesPerFrame_gt);
nfft_gt = NumSamplesPerFrame_gt;

RdrRxSignal_LegsOnly_clean = RdrRxSignal_Toe_GT + RdrRxSignal_RToe_GT;
[S_legs, freq_axis_legs, time_axis_legs] = Zahraa_Spectrogram(RdrRxSignal_LegsOnly_clean, win_gt, noverlap_gt, nfft_gt, fs);

P_legs = abs(S_legs).^2;
legFreq_clean_Hz = zeros(1, size(P_legs,2));
for tIdx = 1:size(P_legs,2)
    frameMax = max(P_legs(:,tIdx));
    idx = find(P_legs(:,tIdx) > 0.01*frameMax, 1, 'last');   % highest bin above 1% of this frame maximum
    if ~isempty(idx)
        legFreq_clean_Hz(tIdx) = freq_axis_legs(idx);
    end
end

minDistance_s = 0.2;   % minimum separation to suppress numerical peak jitter, s
minProminence_Hz = 0.25*(max(legFreq_clean_Hz) - min(legFreq_clean_Hz));   % prominence relative to the track range
[peakVals, peakIdx] = findpeaks(legFreq_clean_Hz, 'MinPeakDistance', round(minDistance_s/(time_axis_legs(2)-time_axis_legs(1))), ...
    'MinPeakProminence', minProminence_Hz);

trueStepTimes_s_corrected = time_axis_legs(peakIdx);
edgeMargin_s = 0.1;   % exclude buffer-transition peaks within this margin, s
trueStepTimes_s_corrected = trueStepTimes_s_corrected(trueStepTimes_s_corrected > bufferTime_s + edgeMargin_s & trueStepTimes_s_corrected < bufferTime_s + SimulationTime_s - edgeMargin_s);
fprintf('Corrected ground-truth step count: %d\n', numel(trueStepTimes_s_corrected));
disp(trueStepTimes_s_corrected);

%% Optional console output
if SimulatorVerbose

    fprintf( ...
        '\nV11 mode %s | speed %.2f m/s | SNR setting %.1f dB\n', ...
        char(AmplitudeMode), v0, SNR_dB);

    fprintf( ...
        'Simulator CPI: requested %.6f s, actual %.6f s\n', ...
        CPI_requested_s_Sim, CPI_actual_s_Sim);

    fprintf( ...
        'Gait cycle %.6f s | reference steps %d\n', ...
        Tc, numel(trueStepTimes_s_corrected));

    disp(AmplitudeSummary);
end

%% Optional spectrogram plot
if SimulatorMakePlots

    figure;

    imagesc( ...
        TimeAxis, SpeedVectorOfInterest, ...
        20*log10(max(S_OfInterestToPlot,realmin)), [-50 0]);

    axis xy;

    xlabel('Time (s)');
    ylabel('Velocity (m/s)');

    colorbar;
    colormap('jet');

    title(sprintf( ...
        'Radar simulator V11: mode %s', char(AmplitudeMode)));
end

%% Local functions: Chen Appendix 4A joint-angle trajectories
function hipFlexion_deg = computeHipFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,tR)
if gaitCharacteristic == 'a'
    hip_x1 = -0.1; hip_x2 = 0.5; hip_x3 = 0.9;
    hip_y1 = 50*relativeWalkingVelocity;
    hip_y2 = -30*relativeWalkingVelocity;
    hip_y3 = 50*relativeWalkingVelocity;
elseif gaitCharacteristic == 'b'
    hip_x1 = -0.1; hip_x2 = 0.5; hip_x3 = 0.9;
    hip_y1 = 25; hip_y2 = -15; hip_y3 = 25;
else
    hip_x1 = 0.2*(relativeWalkingVelocity-1.3)/1.7 - 0.1;
    hip_x2 = 0.5;
    hip_x3 = 0.9;
    hip_y1 = 5*(relativeWalkingVelocity-1.3)/1.7 + 25;
    hip_y2 = -15;
    hip_y3 = 6*(relativeWalkingVelocity-1.3)/1.7 + 25;
end
if abs((hip_x1 + 1) - hip_x3) < 1e-12
    hipPhasePoints = [hip_x1-1,hip_x2-1,hip_x1,hip_x2,hip_x3,hip_x2+1,hip_x3+1];
    hipAnglePoints_deg = [hip_y1,hip_y2,hip_y1,hip_y2,hip_y3,hip_y2,hip_y3];
else
    hipPhasePoints = [hip_x1-1,hip_x2-1,hip_x3-1,hip_x1,hip_x2,hip_x3, ...
                      hip_x1+1,hip_x2+1,hip_x3+1];
    hipAnglePoints_deg = [hip_y1,hip_y2,hip_y3,hip_y1,hip_y2,hip_y3, ...
                          hip_y1,hip_y2,hip_y3];
end
hipFlexion_deg = pchip(hipPhasePoints,hipAnglePoints_deg,tR);
end
function kneeFlexion_deg = computeKneeFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,tR)
if gaitCharacteristic == 'a'
    knee_x1 = 0.17; knee_x2 = 0.4; knee_x3 = 0.75; knee_x4 = 1;
    knee_y1 = 3;
    knee_y2 = 3;
    knee_y3 = 140*relativeWalkingVelocity;
    knee_y4 = 3;
elseif gaitCharacteristic == 'b'
    knee_x1 = 0.17; knee_x2 = 0.4; knee_x3 = 0.75; knee_x4 = 1;
    knee_y1 = 3; knee_y2 = 3; knee_y3 = 70; knee_y4 = 3;
else
    knee_x1 = -0.05*(relativeWalkingVelocity-1.3)/1.7 + 0.17;
    knee_x2 = 0.4;
    knee_x3 = -0.05*(relativeWalkingVelocity-1.3)/1.7 + 0.75;
    knee_x4 = -0.03*(relativeWalkingVelocity-1.3)/1.7 + 1;
    knee_y1 = 22*(relativeWalkingVelocity-1.3)/1.7 + 3;
    knee_y2 = 3;
    knee_y3 = -5*(relativeWalkingVelocity-1.3)/1.7 + 70;
    knee_y4 = 3*(relativeWalkingVelocity-1.3)/1.7 + 3;
end
kneePhasePoints = [knee_x1-1,knee_x2-1,knee_x3-1,knee_x4-1, ...
                   knee_x1,knee_x2,knee_x3,knee_x4, ...
                   knee_x1+1,knee_x2+1,knee_x3+1,knee_x4+1];
kneeAnglePoints_deg = [knee_y1,knee_y2,knee_y3,knee_y4, ...
                       knee_y1,knee_y2,knee_y3,knee_y4, ...
                       knee_y1,knee_y2,knee_y3,knee_y4];
kneeFlexion_deg = pchip(kneePhasePoints,kneeAnglePoints_deg,tR);
end
function ankleFlexion_deg = computeAnkleFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,relativeSupportDuration,tR)
if gaitCharacteristic == 'a'
    ankle_x1 = 0;
    ankle_x2 = 0.08;
    ankle_x3 = 0.5;
    ankle_x4 = relativeSupportDuration;
    ankle_x5 = 0.85;
    ankle_y1 = -3;
    ankle_y2 = -30*relativeWalkingVelocity - 3;
    ankle_y3 = 22*relativeWalkingVelocity - 3;
    ankle_y4 = -34*relativeWalkingVelocity - 3;
    ankle_y5 = -3;
elseif gaitCharacteristic == 'b'
    ankle_x1 = 0;
    ankle_x2 = 0.08;
    ankle_x3 = 0.5;
    ankle_x4 = relativeSupportDuration;
    ankle_x5 = 0.85;
    ankle_y1 = -3;
    ankle_y2 = -18;
    ankle_y3 = 8;
    ankle_y4 = -20;
    ankle_y5 = -3;
else
    ankle_x1 = 0;
    ankle_x2 = 0.08;
    ankle_x3 = -0.1*(relativeWalkingVelocity-1.3)/1.7 + 0.5;
    ankle_x4 = relativeSupportDuration;
    ankle_x5 = 0.85;
    ankle_y1 = 5*(relativeWalkingVelocity-1.3)/1.7 - 3;
    ankle_y2 = 4*(relativeWalkingVelocity-1.3)/1.7 - 18;
    ankle_y3 = -3*(relativeWalkingVelocity-1.3)/1.7 + 8;
    ankle_y4 = -8*(relativeWalkingVelocity-1.3)/1.7 - 20;
    ankle_y5 = 5*(relativeWalkingVelocity-1.3)/1.7 - 3;
end
anklePhasePoints = [ankle_x1-1,ankle_x2-1,ankle_x3-1,ankle_x4-1,ankle_x5-1, ...
                    ankle_x1,ankle_x2,ankle_x3,ankle_x4,ankle_x5, ...
                    ankle_x1+1,ankle_x2+1,ankle_x3+1,ankle_x4+1,ankle_x5+1];
ankleAnglePoints_deg = [ankle_y1,ankle_y2,ankle_y3,ankle_y4,ankle_y5, ...
                        ankle_y1,ankle_y2,ankle_y3,ankle_y4,ankle_y5, ...
                        ankle_y1,ankle_y2,ankle_y3,ankle_y4,ankle_y5];
ankleFlexion_deg = pchip(anklePhasePoints,ankleAnglePoints_deg,tR);
end
function elbowFlexion_deg = computeElbowFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,tR)
if gaitCharacteristic == 'a'
    elbow_x1 = 0.05; elbow_x2 = 0.5; elbow_x3 = 0.9;
    elbow_y1 = 6*relativeWalkingVelocity + 3;
    elbow_y2 = 34*relativeWalkingVelocity + 3;
    elbow_y3 = 10*relativeWalkingVelocity + 3;
elseif gaitCharacteristic == 'b'
    elbow_x1 = 0.05;
    elbow_x2 = 0.01*(relativeWalkingVelocity-0.5)/0.8 + 0.5;
    elbow_x3 = 0.9;
    elbow_y1 = 8*(relativeWalkingVelocity-0.5)/0.8 + 6;
    elbow_y2 = 24*(relativeWalkingVelocity-0.5)/0.8 + 20;
    elbow_y3 = 9*(relativeWalkingVelocity-0.5)/0.8 + 8;
else
    elbow_x1 = 0.05;
    elbow_x2 = 0.04*(relativeWalkingVelocity-1.3)/1.7 + 0.51;
    elbow_x3 = -0.1*(relativeWalkingVelocity-1.3)/1.7 + 0.9;
    elbow_y1 = -6*(relativeWalkingVelocity-1.3)/1.7 + 14;
    elbow_y2 = 26*(relativeWalkingVelocity-1.3)/1.7 + 44;
    elbow_y3 = -6*(relativeWalkingVelocity-1.3)/1.7 + 17;

end
elbowPhasePoints = [elbow_x1-1,elbow_x2-1,elbow_x3-1, ...
                    elbow_x1,elbow_x2,elbow_x3, ...
                    elbow_x1+1,elbow_x2+1,elbow_x3+1];
elbowAnglePoints_deg = [elbow_y1,elbow_y2,elbow_y3, ...
                        elbow_y1,elbow_y2,elbow_y3, ...
                        elbow_y1,elbow_y2,elbow_y3];
elbowFlexion_deg = pchip(elbowPhasePoints,elbowAnglePoints_deg,tR);
end

function sigma_m2 = v11EllipsoidRCS( ...
    x1, z1, x2, z2, a_m, radarZX)

ax = x2 - x1;
az = z2 - z1;

segmentLength_m = hypot(ax, az);
cSemi_m = segmentLength_m/2;

mx = (x1 + x2)/2;
mz = (z1 + z2)/2;

lx = radarZX(2) - mx;
lz = radarZX(1) - mz;

losLength_m = hypot(lx, lz);

assert( ...
    all(isfinite(segmentLength_m)) && ...
    all(segmentLength_m > 0) && ...
    all(isfinite(losLength_m)) && ...
    all(losLength_m > 0), ...
    'Invalid ellipsoid geometry.');

cosTheta = (ax.*lx + az.*lz)./ ...
    (segmentLength_m.*losLength_m);

cosTheta = max(-1, min(1, cosTheta));

sinSquared = max(0, 1 - cosTheta.^2);

sigma_m2 = pi*a_m^4.*cSemi_m.^2 ./ ...
    (a_m^2.*sinSquared + cSemi_m.^2.*cosTheta.^2).^2;

assert( ...
    all(isfinite(sigma_m2)) && all(sigma_m2 > 0), ...
    'Invalid ellipsoid RCS.');

end

function [S, freq_axis, time_axis] = Zahraa_Spectrogram(x, win, noverlap, nfft, fs)

    x   = x(:).';       % row vector without conjugation
    win = win(:).';     % row window without conjugation
    L = length(win);        % frame length = window length
    N = length(x);
    hop = L - noverlap;     % step between successive frame starts

    M = floor((N - L) / hop) + 1;    % number of frames
    S = zeros(nfft, M);

    frame_starts = zeros(1, M);      % frame start indices for the time axis

    for m = 1:M
        start_idx = (m - 1) * hop + 1;
        stop_idx  = start_idx + L - 1;
        frame_starts(m) = start_idx;
        segment = x(start_idx:stop_idx) .* win;     % window the current segment
        S(:, m) = fft(segment, nfft).';
    end

    time_axis = (frame_starts - 1 + L/2) / fs;  % frame-centre time, s

    if isreal(x)
        freq_axis = (0:nfft-1) * (fs / nfft);
    else
        S = fftshift(S, 1);
        freq_axis = (-nfft/2 : nfft/2 - 1) * (fs / nfft);
    end
end

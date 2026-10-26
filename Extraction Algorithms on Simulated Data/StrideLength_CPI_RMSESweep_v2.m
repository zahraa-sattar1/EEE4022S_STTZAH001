%% Stride length CPI sweep
% Zahraa Sattar | EEE4022S
%
% Compare CPI settings across repeated simulated walking trials.
% Report stride length RMSE and percentage error against the simulator reference.

CPIValues = [0.025, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];    % same set used for the other CPI sweeps
SNR_dB_fixed = 30;                          % SNR held constant while CPI is the thing being varied
numTrials_cpi = 15;
rmse_m = zeros(size(CPIValues));            % one RMSE result per CPI value
meanPercentError = zeros(size(CPIValues));  % one percent-error result per CPI value, for context alongside the RMSE

stdError_m = zeros(size(CPIValues));        % trial-to-trial spread of the error at each CPI value, for the error bar

CPI_actual_SL_s = zeros(size(CPIValues));   % will hold the CPI actually achieved for each requested CPI_s value

%% Run the sweep:

set(0, 'DefaultFigureVisible', 'off');  % suppress figures

for cpiIdx = 1:length(CPIValues)
    CPI_s = CPIValues(cpiIdx);                          % set this CPI value before the trials below use it
    SNR_dB = SNR_dB_fixed;                              % keep SNR fixed for every trial in this sweep

    trialRMSE_m = zeros(1, numTrials_cpi);

    allStrideErrors_m = [];          % pooled matched per-stride errors across all trials, at this CPI
    allGroundTruthMatched_m = [];

for trial = 1:numTrials_cpi     % repeat numTrials_cpi times, since the noise is different every run

        % Finding_StrideLength_v2_noplots_torsovelocity
        %% Step 1: Step recognition: get detected step-event timestamps
        % Finding_StepRecognition_v4_noplots
        % RadarSimulator_v11_torso_legs_arms_head_buffered_noplots
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

        RadarHeight_m = 0.75;              % Replace with actual measured height
        RadarCord_ZX = [RadarHeight_m, 0]; % Coordinate order: [Z, X]

        %% Human: base walking geometry
        R0 = SimulatorR0_m;

        if ~exist('v0','var')
            v0 = 1.2;
        end

        validateattributes(v0, {'numeric'}, ...
            {'real','finite','scalar','positive'});

        A  = 1;                          % return signal amplitude (arbitrary units, same for every scatterer)

        %% Human - Boulic Walking Model Parameters (Appendix A)
        Height_m = 1.7;                 % subject height, m
        H_Spine = 0.58*Height_m;
        H_t = 0.53*Height_m;            % Height of hip joint
        L_torso = 0.288*Height_m;       % torso length
        L_upperLeg = 0.245*Height_m;              % Length of Upper Leg
        L_lowerLeg = 0.246*Height_m;              % Length of Lower Leg
        L_legHeight = L_upperLeg + L_lowerLeg;    % Length of Leg
        L_foot = 0.143*Height_m;                  % Length of foot, ankle to toe
        L_upperArm = 0.188*Height_m;            % Length of upper arm, shoulder to elbow
        L_forearm = 0.145*Height_m;             % Length of forearm, elbow to wrist
        L_head = 0.13*Height_m;         % length of head segment
        Vr = v0/L_legHeight;            % Relative velocity normalised by the leg height
        RLc = 1.346*sqrt(Vr);
        Tc  = RLc/Vr;                   % Gait cycle duration, seconds
        Ts  = 0.752*Tc - 0.143;
        Ts_pct = Ts/Tc;                 % relative stance duration

        %% Simulation Settings
        fs = 44.1e3;                     % sample rate, Hz
        SimulationTime_s = 4;            % total length of the simulated signal, s

        if v0 < 1
            SimulationTime_s = 8;
        end

        bufferTime_s = 0.5;
        N = round((SimulationTime_s + 2*bufferTime_s)*fs);  % total number of samples, now includes both buffers
        if ~exist('SNR_dB','var')
            SNR_dB = 40;
        end

        %% Processing parameters

        if exist('CPI_s','var')
            CPI = CPI_s;
        else
            CPI = 0.03;
        end

        CPI_requested_s_Sim = CPI;

        if v0 < 1.0                      % slower speeds need a longer CPI to spectrogram clearly
            CPI = 0.075;
        end

        NumSamplesPerFrame = 2^(nextpow2(round(CPI*fs)));       % actual window length after rounding up to a power of 2

        CPI_actual_s_Sim = NumSamplesPerFrame/fs;

        assert(NumSamplesPerFrame <= N, ...
            'CPI exceeds the recording duration.');

        %% Time Vector and Normalised Gait Time
        lamda = c/fc;       % radar wavelength, m
        t = (0:1:(N-1))*1/fs;           % real time vector, seconds
        t_walk = min(max(t - bufferTime_s, 0), SimulationTime_s);
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
        theta_deg = d3_computeHipFlexionDeg(Vr,gaitCharacteristic,t_pct);   % hip flexion angle, degrees, one value per t_pct sample
        theta_rad = theta_deg*pi/180;                                     % convert to radians for sin/cos below
        L_thigh_midpoint = L_upperLeg/2;                 % lever arm, hip to thigh midpoint
        Thigh_X = Hip_X - L_thigh_midpoint*sin(theta_rad);
        Thigh_Z = Hip_Z - L_thigh_midpoint*cos(theta_rad);

        %% Thigh Radar Return Signal
        Thigh_R = sqrt((RadarCord_ZX(2) - Thigh_X).^2 + (RadarCord_ZX(1) - Thigh_Z).^2);    % range from radar to torso point

        %% Knee Movement Calculations (Chen Appendix 4A / Boulic Appendix B.2)
        Knee_X = Hip_X - L_upperLeg*sin(theta_rad);
        Knee_Z = Hip_Z - L_upperLeg*cos(theta_rad);
        theta_knee_deg = d3_computeKneeFlexionDeg(Vr,gaitCharacteristic,t_pct);   % how much the shank bends back relative to the thigh, degrees
        theta_knee_rad = theta_knee_deg*pi/180;
        theta_shank_rad = theta_rad - theta_knee_rad;       % The shanks angles is the knee angle relative to the hip angle
        L_shank_midpoint = L_lowerLeg/2;
        Shank_X = Knee_X - L_shank_midpoint*sin(theta_shank_rad);
        Shank_Z = Knee_Z - L_shank_midpoint*cos(theta_shank_rad);

        %% Shank Radar Return Signal
        Shank_R = sqrt((RadarCord_ZX(2) - Shank_X).^2 + (RadarCord_ZX(1) - Shank_Z).^2);

        %% Ankle Movement Calculations (Chen Appendix 4A / Boulic Appendix B.3)
        Ankle_X = Knee_X - L_lowerLeg*sin(theta_shank_rad);
        Ankle_Z = Knee_Z - L_lowerLeg*cos(theta_shank_rad);
        theta_ankle_deg = d3_computeAnkleFlexionDeg(Vr,gaitCharacteristic,Ts_pct,t_pct);   % ankle flexion angle, degrees
        theta_ankle_rad = theta_ankle_deg*pi/180;
        theta_foot_rad = theta_rad - theta_knee_rad + theta_ankle_rad;  % Relative foot angle

        Toe_X = Ankle_X - L_foot*cos(theta_foot_rad);
        Toe_Z = Ankle_Z + L_foot*sin(theta_foot_rad);

        Foot_X = (Ankle_X + Toe_X)/2;
        Foot_Z = (Ankle_Z + Toe_Z)/2;

        %% Foot Radar Return Signal
        Foot_R = sqrt((RadarCord_ZX(2) - Foot_X).^2 + (RadarCord_ZX(1) - Foot_Z).^2);

        %% Right Leg (Chen Appendix 4A / Boulic Appendix B.1-B.3)
        t_pct_right = mod(t_pct + 0.5, 1);   % right leg's own gait phase
        RHip_X = Hip_X;
        RHip_Z = Hip_Z;
        RHipFlex_deg = d3_computeHipFlexionDeg(Vr,gaitCharacteristic,t_pct_right);     % right hip flexion angle, degrees
        RHipFlex_rad = RHipFlex_deg*pi/180;
        RKneeFlex_deg = d3_computeKneeFlexionDeg(Vr,gaitCharacteristic,t_pct_right);   % right knee flexion angle, degrees
        RKneeFlex_rad = RKneeFlex_deg*pi/180;
        RAnkleFlex_deg = d3_computeAnkleFlexionDeg(Vr,gaitCharacteristic,Ts_pct,t_pct_right);  % right ankle flexion angle, degrees
        RAnkleFlex_rad = RAnkleFlex_deg*pi/180;
        theta_Rthigh_rad = RHipFlex_rad;                                  % right thigh's angle in space = right hip flexion angle directly
        theta_Rshank_rad = RHipFlex_rad - RKneeFlex_rad;                  % right shank angle = thigh angle, bent back by the knee
        theta_Rfoot_rad  = RHipFlex_rad - RKneeFlex_rad + RAnkleFlex_rad; % right foot angle = shank angle, bent again by the ankle
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
        Shoulder_X = Spine_X + L_torso*sin(theta_r*pi/180);   % shoulder sits at the top of the torso segment, full L_torso
        Shoulder_Z = Spine_Z + L_torso*cos(theta_r*pi/180);
        As = 9.88*Vr;   % shoulder swing amplitude, shared by both arms, only phase differs
        LShoulderFlex_deg = 3  - As*cos(2*pi*t_pct_right);   % left arm, opposite of right foot
        LShoulderFlex_rad = LShoulderFlex_deg*pi/180;
        LUpperArmAngle_rad = LShoulderFlex_rad - theta_r*pi/180;   % torso's own lean subtracted out
        LElbowFlex_deg = d3_computeElbowFlexionDeg(Vr,gaitCharacteristic,t_pct_right);   % left elbow, same phase as left shoulder
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
        RElbowFlex_deg = d3_computeElbowFlexionDeg(Vr,gaitCharacteristic,t_pct);   % left elbow, same phase as left shoulder
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

            SegmentRCS_m2{segIdx} = d3_v11EllipsoidRCS( ...
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
        processingGain = sum(win)^2 / sum(win.^2);

        noisePower = 10^(-(SNR_dB - 10*log10(NumSamplesPerFrame))/10); % - compute noise power
        noise = sqrt(noisePower/2)*(randn(1,N) + 1i*randn(1,N));       % - generate complex noise samples. I, Q from Guassian distribution

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
        [S, FrequencyAxis, TimeAxis] = Zahraa_Spectrogram(RdrRxSignal, win, noverlap, nfft, fs);
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
            plot(RThigh_X, RThigh_Z, 'm.');   % left thigh
            plot(RShank_X, RShank_Z, 'g.');   % left shank
            plot(RFoot_X, RFoot_Z, 'r.');     % left foot
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
        CPI_groundTruth_s = 0.025;
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
            idx = find(P_legs(:,tIdx) > 0.01*frameMax, 1, 'last');   % 1% of this frame's own peak, not an absolute value
            if ~isempty(idx)
                legFreq_clean_Hz(tIdx) = freq_axis_legs(idx);
            end
        end

        minDistance_s = 0.2;   % keeps numerical jitter from being counted as separate peaks
        minProminence_Hz = 0.25*(max(legFreq_clean_Hz) - min(legFreq_clean_Hz));
        [peakVals, peakIdx] = findpeaks(legFreq_clean_Hz, 'MinPeakDistance', round(minDistance_s/(time_axis_legs(2)-time_axis_legs(1))), ...
            'MinPeakProminence', minProminence_Hz);

        trueStepTimes_s_corrected = time_axis_legs(peakIdx);
        edgeMargin_s = 0.1;
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

        %% Step 0 - define input parameters
        freq_high_Yardibi = 350*24/5.8;

        c = 299792458;
        fc = 24e9;
        lambda = c/fc;

        %% Step 1: average filter on raw radar data
        marginFactor = 1.5;
        avgSamples = floor(fs / (marginFactor * freq_high_Yardibi));
        xAvg = d2_AverageFilterRadarData(RdrRxSignal, avgSamples);
        sampleWindow = 1:500;   % zoom into a short segment to actually see the smoothing
        fprintf('Std of raw signal (this segment): %.4f\n', std(real(RdrRxSignal(sampleWindow))));
        fprintf('Std of averaged signal (this segment): %.4f\n', std(real(xAvg(sampleWindow))));

        %% Step 2: compute spectrogram
        if ~exist('CPI_s','var')
            CPI_s = 0.1;   % Window-length parameter, matching the walking-speed algorithm
        end
        [S_step, freq_axis_step_full, time_axis_step, fsProc_step] = d2_ComputeGaitSpectrogram(xAvg, fs, CPI_s); % change freq_axis_step to freq_axis_step_full

        CPI_actual_s = length(freq_axis_step_full) / fsProc_step;     % CPI_actual_s = N*(1/fs)
        fprintf('Requested CPI %.4f:\n', CPI_s);
        fprintf('Actual CPI %.4f:\n', CPI_actual_s);

        %% Step 3: background noise subtraction
        nonzeroIndex = find( (freq_axis_step_full > 0) & (freq_axis_step_full <= freq_high_Yardibi));
        freq_axis_step = freq_axis_step_full(nonzeroIndex);
        S_step_leg = S_step(nonzeroIndex,:);
        P_step2 = abs(S_step_leg).^2;
        Pm = d2_SubtractBackgroundNoise(S_step_leg);
        fprintf('Check: Pm never negative -> %d\n', all(Pm(:) >= 0));
        fprintf('Check: Pm never exceeds original power -> %d\n', all(Pm(:) <= P_step2(:) + eps));

        %% Step 4: median filter on spectral data
        frame_period_s = time_axis_step(2) - time_axis_step(1);   % actual spacing between spectrogram frames
        medianWindow_s = Tc/8;
        medianFrames = max(1, round(medianWindow_s / frame_period_s));
        Pmedian = d2_MedianFilterSpectrogram(Pm, medianFrames);
        noiseBinIdx = 1;
        fprintf('Std at edge frequency bin BEFORE median filter: %.4f\n', std(Pm(noiseBinIdx,:)));
        fprintf('Std at edge frequency bin AFTER  median filter: %.4f\n', std(Pmedian(noiseBinIdx,:)));
        fprintf('Pmedian range (dB): [%.1f, %.1f]\n', ...
            10*log10(min(Pmedian(Pmedian>0)) + eps), 10*log10(max(Pmedian(:)) + eps));

        %% Step 5: Select Leg Frequency

        if exist('PthOverride','var')
            Pth = PthOverride;
        elseif exist('SNR_dB','var')
            Pth = 10 - SNR_dB;
            Pth = max(Pth, -40);
        else
            Pth = -30;
        end

        [legFreq_Hz, isDetected] = d2_SelectLegFrequency(Pmedian, freq_axis_step, Pth);
        fprintf('Frames with a leg-frequency detection: %d / %d\n', sum(isDetected), numel(isDetected));

        %% Step 6
        [legFreq_refined, isOutlier] = d2_RefineLegFrequencyOutliers(legFreq_Hz, isDetected);
        fprintf('Frames flagged as outliers: %d\n', sum(isOutlier));

        figure('Name','Step 6 verification - remove outliers');
        plot(time_axis_step, legFreq_Hz, 'b.');
        hold on;
        plot(time_axis_step, legFreq_refined, 'r-', 'LineWidth', 1.2);
        xlabel('Time (s)');
        ylabel('Leg frequency (Hz)');
        legend('Raw Step 5 detections', 'Step 6 refined (outliers removed)');
        title('Steps 5-6: leg frequency selection and outlier rejection');
        grid on;

        %% Step 7
        smoothWindow_s = Tc/8;
        avgFrames7 = max(1, round(smoothWindow_s / frame_period_s));
        F_smooth = d2_SmoothLegFrequency(legFreq_refined, avgFrames7);

        %% Step 8
        frequencyNormalization_Hz = max(F_smooth,[],'omitnan');
        normalizedLegFrequency =  F_smooth/frequencyNormalization_Hz;   % normalise leg frequency vector so maximum value is 1

        DopplerRes_Hz = 1/CPI_actual_s;

        cPeak =  DopplerRes_Hz/freq_high_Yardibi;
        stepIdx8 = find( (normalizedLegFrequency(2:end-1) >= normalizedLegFrequency(1:end-2)) & ...
                          (normalizedLegFrequency(2:end-1) >= normalizedLegFrequency(3:end) + cPeak) ) + 1;
        fprintf('Step 8 raw peaks: %d\n', numel(stepIdx8));

        %% Step 9
        minInterval_s = 0.35;   % Wang et al., 2014
        stepIdx9 = d2_RemoveCloseSteps(stepIdx8, time_axis_step, F_smooth, minInterval_s);
        detectedStepTimes_s = time_axis_step(stepIdx9);
        fprintf('Step 9 final step count: %d\n', numel(stepIdx9));

        %% Step 10: Calculate step and stride parameters
        stepSpeeds_m_s =  F_smooth(stepIdx9)*lambda/2;
        stepDurations_s = diff(detectedStepTimes_s);
        if numel(detectedStepTimes_s) >= 3
            strideDurations_s = detectedStepTimes_s(3:end) - detectedStepTimes_s(1:end-2);
        else
            strideDurations_s = [];
        end


        %% FUNCTIONS

        %% Step 1: Average Filter

        %% Step 2: Torso velocity-vs-time

        torsoVelocityMethod = 'walkingspeed';

        if strcmp(torsoVelocityMethod, 'walkingspeed')   % swap in the more refined walking-speed pipeline instead
            
            % Finding_WalkingSpeed_v3_noplots
            %% 1. Initial passband (Step 1's K-scaling, using the real fc)
            fc_paper = 5.8e9;              % Yardibi et al. (2011)'s radar carrier, Hz

            if ~exist('CPI_s','var')
                CPI_s = 0.1;                % Introduce CPI input variable
            end

            N_actual = 2^(nextpow2(round(CPI_s*fs)));
            CPI_actual_s_WS = N_actual / fs;        % CPI_actual = N * (1/fs)
            fprintf('Requested CPI %.4f:\n', CPI_s);
            fprintf('Actual CPI %.4f:\n', CPI_actual_s_WS);

            K = fc / fc_paper;              % fc comes from the simulator workspace
            f_low0  = 5 * K;    % multiply low and high frequencies by K scaler
            f_high_min = 50*K;
            maxExpectedSpeed_mps = 4.0;             % covers brisk walking/light jogging
            f_high0 = maxExpectedSpeed_mps * 2 * fc / c;
            fprintf('K = %.6f, initial band = [%.4f, %.4f] Hz\n', K, f_low0, f_high0);

            %% 2. Steps 1-6
            numIterations = 3;
            resultsIter = d4_RunGaitBandRefinement(RdrRxSignal, fs, f_low0, f_high_min, f_high0, numIterations, CPI_s); % f_high_min is also an input
            for iter = 1:numIterations
                fprintf('Iteration %d: band = [%.2f, %.2f] Hz, centre = %.4f Hz\n', ...
                    iter, resultsIter(iter).freq_low, resultsIter(iter).freq_high, resultsIter(iter).centreFreq_Hz);
            end
            final = resultsIter(end);   % most refined pass

            %% 3. Step 7: suppress no-motion frames
            finalValidFrequencyIndex = find( abs(final.freq_axis) >= final.freq_low &  abs(final.freq_axis) <= final.freq_high);

            energyPerFrame = sum(final.Ssmooth(finalValidFrequencyIndex, :), 1);
            energyBaseline = median(energyPerFrame,'omitnan');
            energySpread = median(abs(energyPerFrame-energyBaseline),'omitnan');
            energyThreshold = max(0,energyBaseline-2*energySpread);
            isNoMotion = energyPerFrame < energyThreshold;                          % identify frames with no motion (ie below threshold)
            torsoFreq_suppressed = final.torsoFreq_refined_Hz;
            torsoFreq_suppressed(isNoMotion) = 0;

            %% 4. Step 8: compute walking speed
            c = 3e8;   % speed of light, m/s
            velocityTrack_mps = torsoFreq_suppressed * c / (2 * fc);
            velocityTrack_mps(isNoMotion) = NaN;                              % instantaneous velocity
            movingFrames = torsoFreq_suppressed ~= 0;
            meanSpeed_mps = mean(velocityTrack_mps(movingFrames), 'omitnan');   % averages speed that are not equal to 0
            fprintf('\nExtracted mean walking speed: %.4f m/s\n', meanSpeed_mps);


            %% Ground-truth check
            if exist('v0','var')
                percentError = 100 * (meanSpeed_mps - v0) / v0;
                fprintf('Ground-truth mean speed v0        : %.3f m/s\n', v0);
                fprintf('Extracted mean speed              : %.3f m/s\n', meanSpeed_mps);
                fprintf('Error vs. ground truth            : %+.1f %%\n', percentError);
            end
            if all(torsoFreq_suppressed == 0)
                warning('Every frame was suppressed as no-motion, energyThreshold is set too high.');
            end

            %% FUNCTIONS

            %% Step 1: Ideal Bandpass Filter Function
            torsoVelocity_m_s = velocityTrack_mps;
            torsoVelocity_m_s = fillmissing(torsoVelocity_m_s, 'linear', 'EndValues', 'nearest');
            torsoVelocityTime_s = final.time_axis;              % matching time axis from that script's output

        else   

            [S_torso, freq_axis_torso, time_axis_torso] = Zahraa_Spectrogram(RdrRxSignal, win, noverlap, nfft, fs);

            torsoBand_Hz = [5, 350*fc/5.8e9];   % widened band at 350 Hz, restricts search to torso return bandwidth

            bandIdx = find(freq_axis_torso >= torsoBand_Hz(1) & freq_axis_torso <= torsoBand_Hz(2));    % frequency-axis indices that fall inside the torso band

            torsoFreq_Hz = zeros(1, size(S_torso,2));                   % one torso frequency value per time frame, filled in below
            
            for tIdx = 1:size(S_torso,2)                                % loop over every time frame in the spectrogram
                P_col = abs(S_torso(bandIdx,tIdx)).^2;                  % power in the torso band, this frame only
                [~, relIdx] = max(P_col);
                torsoFreq_Hz(tIdx) = freq_axis_torso(bandIdx(relIdx));
            end

            torsoVelocity_m_s = torsoFreq_Hz*lamda/2;   % v = f*lambda/2 (standard Doppler equation)
            torsoVelocityTime_s = time_axis_torso;      % keeps the matching time axis alongside the velocity track

        end

        %% Step 3: Ground truth (heel-strike based) - correct reference for stride LENGTH

        trueStepTimes_s = bufferTime_s : (Tc/2) : (bufferTime_s + SimulationTime_s);

        %% Step 4: Ground-truth stride length, from the simulator's own torso trajectory
        [groundTruth_stepLengths_m, groundTruth_strideLengths_m] = d1_ComputeGroundTruthStrideLength(t, Torso_R, trueStepTimes_s);

        %% Step 5: Estimated stride length, from integrated Doppler velocity + detected steps
        [estimated_stepLengths_m, estimated_strideLengths_m, displacement_m] = d1_EstimateStrideLength(torsoVelocity_m_s, torsoVelocityTime_s, detectedStepTimes_s);

        %% Step 6: Compare
        fprintf('--- Stride length comparison ---\n');
        fprintf('Ground truth: mean %.4f m, std %.4f m (n=%d)\n', mean(groundTruth_strideLengths_m), std(groundTruth_strideLengths_m), numel(groundTruth_strideLengths_m));     % summary stats for the true stride lengths
        fprintf('Estimated:    mean %.4f m, std %.4f m (n=%d)\n', mean(estimated_strideLengths_m), std(estimated_strideLengths_m), numel(estimated_strideLengths_m));           % same summary stats for the estimated stride lengths
        fprintf('Constant speed check, v0*Tc = %.4f m\n', v0*Tc);

        %% Local functions

        estimatedStrideStartTimes_s = detectedStepTimes_s(1:end-2);
        groundTruthStrideStartTimes_s = trueStepTimes_s(1:end-2);
        matchTolerance_s = Tc/4;
        [estimatedMatched_m, groundTruthMatched_m] = MatchStrides(estimated_strideLengths_m, estimatedStrideStartTimes_s, groundTruth_strideLengths_m, groundTruthStrideStartTimes_s, matchTolerance_s);
        strideErrors_m = estimatedMatched_m - groundTruthMatched_m;     % per-stride error, one value per matched pair

        trialRMSE_m(trial) = sqrt(mean(strideErrors_m.^2, 'omitnan'));   % this trial's RMSE across its own matched strides

        allStrideErrors_m = [allStrideErrors_m, strideErrors_m];        % pooled across all trials at this CPI
        allGroundTruthMatched_m = [allGroundTruthMatched_m, groundTruthMatched_m];
end

    rmse_m(cpiIdx) = mean(trialRMSE_m, 'omitnan');
    stdError_m(cpiIdx) = std(trialRMSE_m, 'omitnan');

    meanPercentError(cpiIdx) = 100 * mean(allStrideErrors_m, 'omitnan') / mean(allGroundTruthMatched_m, 'omitnan');  % average signed percent error, for context

    CPI_actual_SL_s(cpiIdx) = 2^(nextpow2(round(CPI_s*fs))) / fs;   % (CPI = N*(1/fs))
    fprintf('CPI = %.3f s -> RMSE = %.4f m (%.1f%% mean error)\n', CPI_s, rmse_m(cpiIdx), meanPercentError(cpiIdx));    % progress readout while the sweep runs
end

set(0, 'DefaultFigureVisible', 'on');

close all;

%% Tidy up the arrays and catch a partial run before it causes a confusing error

CPIValues = CPIValues(:).';
CPI_actual_SL_s = CPI_actual_SL_s(:).';
rmse_m = rmse_m(:).';

if length(CPIValues) ~= length(rmse_m)  
    error('CPIValues has %d entries but rmse_m has %d - the script was probably only partly run. Re-run the WHOLE file from the top, not a highlighted section.', ...
        length(CPIValues), length(rmse_m));
end

%% Plot the results: RMSE against the actual CPI

figure('Name','CPI vs stride-length RMSE');
errorbar(CPI_actual_SL_s, rmse_m, stdError_m, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor','b', 'CapSize', 8);
xlabel('Actual CPI (s)');
ylabel('RMSE (m)');
title(sprintf('Stride-length RMSE vs. CPI (SNR fixed at %d dB)', SNR_dB_fixed));
grid on;
[bestRMSE, bestIdx] = min(rmse_m);      % finds the lowest RMSE and which CPI gave it
bestCPI = CPI_actual_SL_s(bestIdx);     % the actual CPI value at that lowest point
hold on;
plot(bestCPI, bestRMSE, 'rp', 'MarkerSize', 14, 'MarkerFaceColor','r');     % marks the best point on the plot

%% Display final results

fprintf('\n--- Results: CPI sweep (SNR fixed at %d dB) ---\n', SNR_dB_fixed);

for i = 1:length(CPIValues)
    fprintf('CPI_s requested = %.3f s (actual = %.3f s) -> RMSE = %.4f m\n', CPIValues(i), CPI_actual_SL_s(i), rmse_m(i));
end
fprintf('Best CPI = %.3f s (actual), with RMSE = %.4f m\n', bestCPI, bestRMSE);     % prints the answer directly

%% Local functions

function hipFlexion_deg = d3_computeHipFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,tR)
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
function kneeFlexion_deg = d3_computeKneeFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,tR)
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
function ankleFlexion_deg = d3_computeAnkleFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,relativeSupportDuration,tR)
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
function elbowFlexion_deg = d3_computeElbowFlexionDeg(relativeWalkingVelocity,gaitCharacteristic,tR)
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

function sigma_m2 = d3_v11EllipsoidRCS( ...
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

function xAvg = d2_AverageFilterRadarData(x, avgSamples)
    isRowInput = isrow(x);
    x = x(:).';
    xAvg = movmean(x, avgSamples);
    if ~isRowInput
        xAvg = xAvg.';
    end
end

%% Step 2: Compute Spectogram
function [S, freq_axis, time_axis, fsProc] = d2_ComputeGaitSpectrogram(x, fs, CPI_s)  % CPI request, s
    Seifert_approach = 0;
    if (Seifert_approach == 1)
        processingFs = 2560;   % Hz
        if fs > processingFs
            [p, q]  = rat(processingFs / fs, 1e-6);
            xProc   = resample(x, p, q);
            fsProc  = fs * p / q;
        else
            xProc  = x;
            fsProc = fs;
        end
        win      = hamming(256);   % Seifert et al. (2020): M = 256
        noverlap = 255;              % Seifert et al. (2020): M-1 samples
        nfft     = 2048;             % zero-padded
    else
        xProc  = x;
        fsProc = fs;
        N = 2^(nextpow2(round(CPI_s*fsProc)));
        win      = hamming(N);
        noverlap = N*3/4;                           % Wang uses 75% overlap
        nfft     = N;
    end
    [S, freq_axis, time_axis] = d2_Zahraa_Spectrogram(xProc, win, noverlap, nfft, fsProc);
end

%% Step 3: Subtract background noise
function Pm = d2_SubtractBackgroundNoise(S)
    P = abs(S).^2;                  % power, not magnitude (Wang's P(t,f) is a power quantity)
    Pback = zeros(size(P,1),1);                                 % intialise variable containing background power per Doppler bin
    numberBackgroundSamples = max(1,round(0.20*size(P,2)));
    for frequencyIndex = 1:size(P,1)
        sortedPower = sort(P(frequencyIndex,:),'ascend');                        % sort values in ascending order
        Pback(frequencyIndex) =  mean(sortedPower(1:numberBackgroundSamples));
    end
    backgroundMatrix =  repmat(Pback,1,size(P,2));
    Pm = P - backgroundMatrix;                                                   % removes background noise from signal to get Pm
    Pm(Pm < 0) = eps;
end

%% Step 4: Median Filter Spectogram
function Pmedian = d2_MedianFilterSpectrogram(Pm, medianFrames)
    Pmedian = movmedian(Pm, medianFrames, 1);  % median applied along Doppler dimension
end

%% Step 5: Select highest frequency for Leg
function [legFreq_Hz, isDetected] = d2_SelectLegFrequency(Pmedian, freq_axis, Pth)
    PmedianNorm =  Pmedian/(max(Pmedian(:)) + eps);    % normalise Spectrogram
    PmedianDb =  10*log10(PmedianNorm + eps);          % convert to dB, where maximum value is 0 dB
    aboveThresh = PmedianDb > Pth;         % logical mask, same size as Pmedian
    legFreq_Hz = NaN(1, size(Pmedian,2));    % initialise to NaN
    isDetected = false(1, size(Pmedian,2));
    for t = 1:size(Pmedian,2)
        idx = find(aboveThresh(:,t), 1, 'last');   % freq_axis is ascending, so 'last' = highest frequency
        if ~isempty(idx)
            legFreq_Hz(t) = freq_axis(idx);
            isDetected(t) = true;
        end
    end
end

%% Step 6: Finding frequency Outliers
function [legFreq_interp, isOutlier] = d2_RefineLegFrequencyOutliers(legFreq_Hz, isDetected)   % returns function legFreq_interp and not legFreq_refined
    mu    = mean(legFreq_Hz(isDetected));   % Finding average of leg frequencies where Step 5 found actual detections
    sigma = std(legFreq_Hz(isDetected));    % Finding the standard deviation
    isOutlier = false(size(legFreq_Hz));    % Array size of legFreq_Hz set to false
    isOutlier(isDetected) = abs(legFreq_Hz(isDetected) - mu) > 3*sigma; % Among detected frames, checks how far value is from mean
    legFreq_refined = legFreq_Hz;
    legFreq_refined(isOutlier) = NaN;   % Overwrites frequencies that are outliers
    legFreq_interp = fillmissing(legFreq_refined, 'linear', 'EndValues','nearest');  % fills the missing values using linear interpolation
end

%% Step 7: Smooth F(t) with an average filter
function F_smooth = d2_SmoothLegFrequency(F, avgFrames)
    F_smooth = movmean(F, avgFrames, 'omitnan');
end

%% Step 8: Obtain step instances through time domain peak selection for F (t),
function stepIdx = d2_DetectStepPeaks(F, c)
    N = length(F);                                  % Number of time frames in signal
    isPeak = false(1, N);                           % logical flag for each frame, all start as false
    for t = 2:N-1                                   % skip the first/last frame, both need a neighbour on each side
        if F(t) >= F(t-1) && F(t) >= F(t+1) + c     % rising/flat before t, and drops by at least c after t
            isPeak(t) = true;                       % if t satisfies both conditions, mark it a step
        end
    end
    stepIdx = find(isPeak);                         % convert the logical flags into a list of frame indices
end

%% Step 9
function stepIdx_refined = d2_RemoveCloseSteps(stepIdx, time_axis, F_smooth, minInterval_s)
    stepTimes = time_axis(stepIdx);              % look up the actual time of each candidate step
    keep = true(size(stepIdx));                  % logical flag for each candidate, all start as true
    lastKeptIdx = 1;
    for k = 2:length(stepIdx)                    % start from second candidate
        if (stepTimes(k) - stepTimes(lastKeptIdx)) <= minInterval_s    % too close to the last kept step
            if F_smooth(stepIdx(k)) > F_smooth(stepIdx(lastKeptIdx))
                 keep(lastKeptIdx) = false;
                 lastKeptIdx = k;
            else
                 keep(k) = false;
            end
        else
             lastKeptIdx = k;
        end
    end
    stepIdx_refined = stepIdx(keep);
end
function [S, freq_axis, time_axis] = d2_Zahraa_Spectrogram(x, win, noverlap, nfft, fs)
    x   = x(:).';
    win = win(:).';
    L = length(win);        % frame length = window length
    N = length(x);
    hop = L - noverlap;     % step between successive frame starts
    M = floor((N - L) / hop) + 1;    % number of frames
    S = zeros(nfft, M);
    frame_starts = zeros(1, M);      % keep track of where each frame starts (for time axis later)
    for m = 1:M
        start_idx = (m - 1) * hop + 1;
        stop_idx  = start_idx + L - 1;
        frame_starts(m) = start_idx;
        segment = x(start_idx:stop_idx) .* win;     % grab the chunk and window it
        S(:, m) = fft(segment, nfft).';
    end
    time_axis = (frame_starts - 1) / fs;
    if isreal(x)
        freq_axis = (0:nfft-1) * (fs / nfft);
    else
        S = fftshift(S, 1);
        freq_axis = (-nfft/2 : nfft/2 - 1) * (fs / nfft);
    end
end

function [y, f, mask] = d4_IdealBandpassFilter(x, fs, f_low, f_high)
    isRowInput = isrow(x);
    x = x(:).';
    N = length(x);
    X  = fft(x);
    Xs = fftshift(X);
    f = (-floor(N/2):ceil(N/2)-1) * (fs / N);
    mask = (abs(f) >= f_low) & (abs(f) <= f_high);
    Xf = Xs .* mask;
    Xu = ifftshift(Xf);
    y  = ifft(Xu);
    if ~isRowInput
        y = y.';        % restores the caller's original orientation
    end
end

%% Step 2: Compute Spectogram
function [S, freq_axis, time_axis, fsProc] = d4_ComputeGaitSpectrogram(x, fs, CPI_s, freq_high)  % included CPI_s
Seifert_approach = 0;  % Choose Seifert's approach of computing spectrogram
if (Seifert_approach == 1)
    processingFs = 2560;   % Hz
    if fs > processingFs
        [p, q]  = rat(processingFs / fs, 1e-6);
        xProc   = resample(x, p, q);
        fsProc  = fs * p / q;
    else
        xProc  = x;
        fsProc = fs;
    end
    N = 2^(nextpow2(round(CPI_s*fsProc)));
    win      = hamming(N);
    noverlap = N-1;
    nfft     = N*8;
else
        xProc  = x;
        fsProc = fs;
        N = 2^(nextpow2(round(CPI_s*fsProc)));
        win      = hamming(N);
        noverlap = N*3/4;                           % Wang uses 75% overlap
        nfft     = N;
end
    [S_temp, freq_axis_temp, time_axis] = d4_Zahraa_Spectrogram(xProc, win, noverlap, nfft, fsProc);
    validFrequencyIndex = find( (abs(freq_axis_temp) <= freq_high*3) );
    S = S_temp(validFrequencyIndex,:);                  % extract subset of spectogram
    freq_axis = freq_axis_temp(validFrequencyIndex);
end

%% Step 3: Reduce spectral noise with average filter
function Ssmooth = d4_SmoothSpectrogramNoise(S, avgFrames)
      Ssmooth = movmean(abs(S), avgFrames, 1);          % Apply moving average across Doppler bins (i.e dimension 1)
end

%% Step 4: torso Doppler frequency = frequency of peak
function [torsoFreq_Hz, peakVal] = d4_TrackPeakFrequency(Ssmooth, freq_axis, freq_low, freq_high)
    validFrequencyIndex = find( (abs(freq_axis) >= freq_low) & (abs(freq_axis) <= freq_high) );
    Subset_Ssmooth = Ssmooth(validFrequencyIndex,:);
    [peakVal, peakIdx] = max(Subset_Ssmooth, [], 1);
    torsoFreq_Hz = freq_axis(validFrequencyIndex(peakIdx));
end

%% Step 5: Outlier rejection
function [torsoFreq_refined, isOutlier] = d4_RefineOutliers(torsoFreq_Hz)
    mu    = mean(torsoFreq_Hz, 'omitnan');
    sigma = std(torsoFreq_Hz, 'omitnan');
    isOutlier = abs(torsoFreq_Hz - mu) > 3*sigma;
    torsoFreq_refined = torsoFreq_Hz;
    torsoFreq_refined(isOutlier) = NaN;
    torsoFreq_refined = fillmissing(torsoFreq_refined,'linear','EndValues','nearest'); % Fill in NaNs with 'linear' interpolated values
end

%% Step 6: Iterative band refinement (repeats Steps 1-5)
function results = d4_RunGaitBandRefinement(x, fs, freq_low0, f_high_min, freq_high0, numIterations, CPI_s)
    freq_low  = freq_low0;
    freq_high = freq_high0;
    results = struct('freq_low', {}, 'freq_high', {}, 'torsoFreq_Hz', {}, ...
        'torsoFreq_refined_Hz', {}, 'centreFreq_Hz', {});
    for iter = 1:numIterations
        xFilt = d4_IdealBandpassFilter(x, fs, freq_low, freq_high);
        [S, freq_axis, time_axis, fsProc] = d4_ComputeGaitSpectrogram(xFilt, fs, CPI_s, freq_high);
        Power_S = abs(S).^2;                                    % Compute Power of Spectrogram (linear units)
        Ssmooth = d4_SmoothSpectrogramNoise(Power_S, 21);          % Smooth the Power of Spetrogram (linear units)
        torsoFreq_Hz = d4_TrackPeakFrequency(Ssmooth, freq_axis, freq_low, freq_high);
        torsoFreq_refined_Hz = d4_RefineOutliers(torsoFreq_Hz);
        centreFreq_Hz = mean(torsoFreq_refined_Hz, 'omitnan');
        results(iter).freq_low             = freq_low;
        results(iter).freq_high            = freq_high;
        results(iter).torsoFreq_Hz         = torsoFreq_Hz;
        results(iter).torsoFreq_refined_Hz = torsoFreq_refined_Hz;
        results(iter).centreFreq_Hz        = centreFreq_Hz;
        results(iter).Ssmooth              = Ssmooth;
        results(iter).freq_axis            = freq_axis;
        results(iter).time_axis            = time_axis;
        results(iter).fsProc               = fsProc;
        if iter < numIterations
            freq_low = max(0.5*centreFreq_Hz,freq_low0);       % 0.5*MeanGaitFrequency (see Yardibi paper)
            freq_high = max(2.0*centreFreq_Hz,f_high_min);  % 2*MeanGaitFrequency (see Yardibi paper)
        end
    end
end

%% Step 7: Suppress no-motion frames
function torsoFreq_suppressed = d4_SuppressNoMotion(torsoFreq_Hz, Ssmooth, energyThreshold)
    energyPerFrame = sum(Ssmooth, 1);      % total energy across all bins, per frame
    isNoMotion = energyPerFrame < energyThreshold;
    torsoFreq_suppressed = torsoFreq_Hz;
    torsoFreq_suppressed(isNoMotion) = 0;   % wherever isNoMotion is true, that position is set to 0
end

%% Step 8: Compute walking speed
function [velocityTrack_mps, meanSpeed_mps] = d4_ComputeWalkingSpeed(torsoFreq_Hz, ft_Hz, c)
    velocityTrack_mps = torsoFreq_Hz * c / (2 * ft_Hz);
    movingFrames = torsoFreq_Hz ~= 0;
    meanSpeed_mps = mean(velocityTrack_mps(movingFrames), 'omitnan');   % averages speed that are not equal to 0
end

%% Zahraa_Spectrogram   Manually computed spectrogram (no built-in spectrogram).
function [S, freq_axis, time_axis] = d4_Zahraa_Spectrogram(x, win, noverlap, nfft, fs)
    x   = x(:).';
    win = win(:).';
    L = length(win);        % frame length = window length
    N = length(x);
    hop = L - noverlap;     % step between successive frame starts
    M = floor((N - L) / hop) + 1;    % number of frames
    S = zeros(nfft, M);
    frame_starts = zeros(1, M);      % keep track of where each frame starts (for time axis later)
    for m = 1:M
        start_idx = (m - 1) * hop + 1;
        stop_idx  = start_idx + L - 1;
        frame_starts(m) = start_idx;
        segment = x(start_idx:stop_idx) .* win;     % grab the chunk and window it
        S(:, m) = fft(segment, nfft).';
    end
    time_axis = (frame_starts - 1) / fs;
    if isreal(x)
        freq_axis = (0:nfft-1) * (fs / nfft);
    else
        S = fftshift(S, 1);
        freq_axis = (-nfft/2 : nfft/2 - 1) * (fs / nfft);
    end
end

function [stepLengths_m, strideLengths_m] = d1_ComputeGroundTruthStrideLength(t, Torso_R, stepTimes_s)

    Torso_R_atSteps = interp1(t, Torso_R, stepTimes_s, 'linear');               % reads the simulator's real torso position at each step time
    stepLengths_m = abs(diff(Torso_R_atSteps));                                 % distance moved between each pair of consecutive steps
    strideLengths_m = abs(Torso_R_atSteps(3:end) - Torso_R_atSteps(1:end-2));   % distance moved between every other step (same foot, one full stride)
end

function [stepLengths_m, strideLengths_m, displacement_m] = d1_EstimateStrideLength(velocity_m_s, velocityTime_s, stepTimes_s)

    displacement_m = cumtrapz(velocityTime_s, velocity_m_s);                                % integrates velocity over time to reconstruct a position track
    displacement_atSteps = interp1(velocityTime_s, displacement_m, stepTimes_s, 'linear');  % reads that reconstructed position at each detected step time
    stepLengths_m = abs(diff(displacement_atSteps));                                        % distance moved between each pair of consecutive detected steps
    strideLengths_m = abs(displacement_atSteps(3:end) - displacement_atSteps(1:end-2));     % distance moved between every other detected step
end

function [estimatedMatched_m, groundTruthMatched_m] = MatchStrides(estimatedStrideLengths_m, estimatedStrideStartTimes_s, groundTruthStrideLengths_m, groundTruthStrideStartTimes_s, matchTolerance_s)
    estimatedMatched_m = [];
    groundTruthMatched_m = [];
for i = 1:numel(estimatedStrideLengths_m)
        timeDiffs_s = abs(groundTruthStrideStartTimes_s - estimatedStrideStartTimes_s(i));
        [minDiff_s, matchIdx] = min(timeDiffs_s);         % closest ground truth stride in time
if minDiff_s <= matchTolerance_s
            estimatedMatched_m(end+1) = estimatedStrideLengths_m(i);
            groundTruthMatched_m(end+1) = groundTruthStrideLengths_m(matchIdx);
end
end
end

function [S, freq_axis, time_axis] = Zahraa_Spectrogram(x, win, noverlap, nfft, fs)

    x   = x(:).';
    win = win(:).';
    L = length(win);        % frame length = window length
    N = length(x);
    hop = L - noverlap;     % step between successive frame starts

    M = floor((N - L) / hop) + 1;    % number of frames
    S = zeros(nfft, M);

    frame_starts = zeros(1, M);      % keep track of where each frame starts (for time axis later)

    for m = 1:M
        start_idx = (m - 1) * hop + 1;
        stop_idx  = start_idx + L - 1;
        frame_starts(m) = start_idx;
        segment = x(start_idx:stop_idx) .* win;     % grab the chunk and window it
        S(:, m) = fft(segment, nfft).';
    end

    time_axis = (frame_starts - 1 + L/2) / fs;  % labelling the time_axis as the center of the window

    if isreal(x)
        freq_axis = (0:nfft-1) * (fs / nfft);
    else
        S = fftshift(S, 1);
        freq_axis = (-nfft/2 : nfft/2 - 1) * (fs / nfft);
    end
end

%Synergy Comparison
%---Code developed by Federico Tessari, PhD 
%---Massachusetts Institute of Technology
%---Latest Version, Nov-18-2024

%The present code serves to reproduce the results presented in the paper:
%"A Geometric Approach for the Comparison of Kinematic Synergy Postures" by
%Tessari, West and Hogan

%It also serves to familiarize with the use of the Geometric Configuration
%Similarity

% Initialization of Variables
clear, clc, close all

set(0, 'DefaultLineLineWidth', 2);
set(groot,'defaultAxesFontSize',12);
set(0,'defaultfigurecolor',[1 1 1]); % white figure background
set(groot,'defaultAxesBox','on'); % box on
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

%Loading Data
load('T_WH_piano_data.mat')

%Add Path to Functions and Models
addpath('functions')
addpath('functions\HandModel')
addpath('functions\HandModel\HumanHand')

%Specify Simulation Path and CAD parts Path
simulation_path = 'functions\HandModel\hand_simulator_v4.slx';
CAD_path = 'functions\HandModel\HumanHand\';

%Max and Min ROM ANGLES per Each Finger [deg] 
% For MCP, DIP, PIP fingers excluding thumb (Bain et al. 2015) 
% For CMC, MCP, IP fo Thum (Gracia-Ibanez et al. 2016)
% For abductions/adduction
% (https://www.physiotutors.com/wiki/wrist-hand-active-range-of-motion/)
%MCP DIP PIP ABD
IDX_min = (pi/180)*[-13 -5.3 -5.3 -12.5];
IDX_max = (pi/180)*[86 102 76.8 12.5];
MID_min = (pi/180)*[-17.5 -6.5 -5.7 -12.5];
MID_max = (pi/180)*[91.5 99.6 85 12.5];
RIG_min = (pi/180)*[-19.5 -7.3 -4.5 -12.5];
RIG_max = (pi/180)*[91.9 104.1 85.8 12.5];
LIT_min = (pi/180)*[-24.8 -6.5 -5.7 -12.5];
LIT_max = (pi/180)*[89.8 98.8 89.0 12.5];
THB_min = (pi/180)*[-26.2 -21.0 -12.4 -30];
THB_max = (pi/180)*[42.1 26.1 102.1 65];
%Matrix Form
ROM_min=[THB_min(1);THB_min(2);THB_min(3);THB_min(4);...
        IDX_min(1);IDX_min(2);IDX_min(3);IDX_min(4);...
        MID_min(1);MID_min(2);MID_min(3);MID_min(4);...
        RIG_min(1);RIG_min(2);RIG_min(3);RIG_min(4);...
        LIT_min(1);LIT_min(2);LIT_min(3);LIT_min(4)];
ROM_max=[THB_max(1);THB_max(2);THB_max(3);THB_max(4);...
        IDX_max(1);IDX_max(2);IDX_max(3);IDX_max(4);...
        MID_max(1);MID_max(2);MID_max(3);MID_max(4);...
        RIG_max(1);RIG_max(2);RIG_max(3);RIG_max(4);...
        LIT_max(1);LIT_max(2);LIT_max(3);LIT_max(4)];

%Load Simulink Model
handle = load_system(simulation_path);

exp_condition = input(['Select: [1] - Numerical Results, \n'...
                         '[2] - Experimental Results, \n' ...
                         '[3] - Whole Trajectory Simulation \n']);
switch exp_condition
    case 1
        %Comparing Numerical Synergies
        clc, close all

        avg_post = zeros(20,1);
        %Choosing Max and Min ROM
        S1 = ROM_max;
        S2 = ROM_min;

        %Cosine Similarity
        cosSim = abs(S1'*S2)/(norm(S1)*norm(S2));

        t = [0 0.01 0.02 0.03 0.04];
        trj1 = [avg_post';S1';S2';avg_post';avg_post'];
        trj2 = [avg_post';S2';S1';avg_post';avg_post'];
        %Joint Trajectories
        joint_trj1 = [t',trj1];
        joint_trj2 = [t',trj2];

        tf = 0.04;
        ts = 0.01;
        conf = 0;

        [GCS_mean_avgp,GCS_vec] = geomSim(tf,ts,conf,simulation_path);  %conf = 1 returns nice graphical outputs
        GCS_avgp = GCS_mean_avgp(1);
        GCS_maxp = GCS_mean_avgp(2);
        GCS_minp = GCS_mean_avgp(3);
    case 2
        % Comparing Experimental Synergies
        clc, close all

        example = input(['Select: [1] - Same subject 1st vs 2nd synergy, \n'...
            '[2] - Different subjects, same task, 1st synergies. \n']);
        switch example
            case 1
                %Subject 2 - Bach Prelude 1st Syn vs 2nd Syn
                ii = 97;
                jj = 97;
                S1 = T_WH_piano_data.(6){ii,1};
                S2 = T_WH_piano_data.(6){jj,1};
                %Select the Synergy Number from S1 (i) and S2 (j)
                i = 1;
                j = 2;
            case 2
                %Subject 2 - Bach Prelude 1st Syn vs 2nd Syn
                ii = 97;
                jj = 128;
                S1 = T_WH_piano_data.(6){ii,1};
                S2 = T_WH_piano_data.(6){jj,1};
                %Select the Synergy Number from S1 (i) and S2 (j)
                i = 1;
                j = 1;
        end

        %Cosine Similarity
        figure()
        cosineSimilarity = getCosineSimilarity(S1,S2,'Plot','on');


        % Max and Min Activations and Average Posture
        avgP1 = T_WH_piano_data.(10){ii,1};
        avgP2 = T_WH_piano_data.(10){jj,1};

        maxP1 = T_WH_piano_data.(11){ii,1};
        minP1 = T_WH_piano_data.(12){ii,1};

        maxP2 = T_WH_piano_data.(11){jj,1};
        minP2 = T_WH_piano_data.(12){jj,1};

        % Geometric Approach
        conf = input('Select the Simulation Method: [0] - Discrete, [1] - Continuous \n');
        %Simulation Configuration
        switch conf
            case 0
                %Simulation Time Configuration
                tf = 0.04;
                ts = 0.01;
            case 1
                %Simulation Time Configuration
                tf = 5;
                ts = 0.01;
        end

        %Temporal Evolution Configuration
        n = input('Select Temporal Evolution Coeff: [1] - max/min(U*S), [2] - max/min(ROM) \n');

        epsilon = 0.05;
        switch n
            case 1
                %Maximum and Minimum Postures extracted by the max and min of (U*S)
                if (S1(:,i)'*S2(:,j)) >= epsilon
                    joint_trj1 = joint_traj_gen(S1(:,i),avgP1,tf,ts,minP1(i),maxP1(i),conf);
                    joint_trj2 = joint_traj_gen(S2(:,j),avgP2,tf,ts,minP2(j),maxP2(j),conf);
                elseif  (S1(:,i)'*S2(:,j)) <= -epsilon
                    joint_trj1 = joint_traj_gen(S1(:,i),avgP1,tf,ts,minP1(i),maxP1(i),conf);
                    joint_trj2 = joint_traj_gen(-S2(:,j),avgP2,tf,ts,-maxP2(j),-minP2(j),conf);
                elseif abs((S1(:,i)'*S2(:,j))) < epsilon
                    joint_trj1 = joint_traj_gen(S1(:,i),avgP1,tf,ts,minP1(i),maxP1(i),conf);
                    joint_trj2 = joint_traj_gen(S2(:,j),avgP2,tf,ts,minP2(j),maxP2(j),conf);
                end
            case 2
                %Maximum and Minimum Postures to maximize ROM of largest element in
                %the Synergy Vector
                [v_S1max,i_S1max] = max(abs(S1(:,i)));
                [v_S2max,i_S2max] = max(abs(S2(:,j)));
                joint_trj1 = joint_traj_gen(S1(:,i),avgP1,tf,ts,ROM_min(i_S1max)/S1(i_S1max,i),ROM_max(i_S1max)/S1(i_S1max,i),conf);
                joint_trj2 = joint_traj_gen(S2(:,j),avgP2,tf,ts,ROM_min(i_S2max)/S2(i_S2max,j),ROM_max(i_S2max)/S2(i_S2max,j),conf);
        end

        [GCS_mean_avgp,GCS_vec] = geomSim(tf,ts,conf,simulation_path);  %conf = 1 returns nice graphical outputs
        GCS_avgp = GCS_mean_avgp(1);
        GCS_maxp = GCS_mean_avgp(2);
        GCS_minp = GCS_mean_avgp(3);

    case 3
        % Simulating Raw Data
        clc, close all

        conf = input('Select the Simulation Method: [0] - Discrete, [1] - Continuous \n');
        %Simulation Configuration
        switch conf
            case 0
                %Simulation Time Configuration
                tf = 0.04;
                ts = 0.01;
            case 1
                %Simulation Time Configuration
                tf = 5;
                ts = 0.001;
        end

        subj_sel = input('Select the Subject: [0] - Subj 2, [1] - Subj 6 \n');
        %Simulation Configuration
        switch subj_sel
            case 0
                ii = 97; %Subject 2 - Bach Prelude
            case 1
                ii = 128; %Subject 6 - Bach Prelude
        end


        X = T_WH_piano_data.Data{ii,1};
        tt = linspace(0,tf,length(X))';

        %Data Pre-Processing
        % NaN Removal and Wrapping to Pi
        for i = length(X):-1:1
            if sum(double(isnan(X(i,:)))) > 0
                tt(i) = [];
                X(i,:) = [];
            end
            for j = 1:20
                if abs(X(i,j)) > pi
                    X(i,j) = X(i,j) - 2*sign(X(i,j))*pi;
                end
            end
        end

        %Joint Trajectory Simulation
        joint_trj1 = [tt,X];
        joint_trj2 = [tt,X];

        %Geometric Similarity Calculation
        [GCS_mean_avgp,GCS_vec] = geomSim(tf,ts,conf,simulation_path); %conf = 1 returns nice graphical outputs
        GCS_avgp(i,j) = GCS_mean_avgp(1);
        GCS_maxp(i,j) = GCS_mean_avgp(2);
        GCS_minp(i,j) = GCS_mean_avgp(3);
end



%Hand Model Calibration
%---Code developed by Federico Tessari, PhD 
%---Massachusetts Institute of Technology
%---Latest Version, Nov-18-2024
clear, clc, close all

%This code serves to compute the maximum distance between the key-points of
%each finger (PIP,DIP,TIP)

%The results of this calibration are already available in the file:
%"mean_d_max_15_final.mat"

%If you change the hand model (dimensions or range of motions) in some way,
% you need to run this calibration and update the file.

%Add Path to Functions and Models
addpath('functions')
addpath('functions\HandModel')
addpath('functions\HandModel\HumanHand')

%Specify Simulation Path and CAD parts Path
simulation_path = 'functions\HandModel\hand_simulator_v4.slx';
CAD_path = 'functions\HandModel\HumanHand\';

%Load Simulink Model
handle = load_system(simulation_path);

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


% Simulation to Explore the whole task space

%Simulation Sampling Frequency
ts = 0.01;
%Simulation Half-Time
tf = 100;
%Simulation Total-Time
tT = 200;

t = 0:ts:tf;

joint_trj1_p1 = [t',zeros(length(t),20)];
joint_trj2_p1 = [t',zeros(length(t),20)];

%Joint Explorations from Min to Max on both hands
idx = 2;
for i = 1:5
    if i == 1
        for j = 1:4 %THUMB
            if j<=3
            p1 = polyfit([0 tf/4 tf/2 3*tf/4 tf],[0 THB_max(j) 0 THB_min(j) 0],4);
            p2 = polyfit([0 tf/4 tf/2 3*tf/4 tf],[0 THB_min(j) 0 THB_max(j) 0],4);
            joint_trj1_p1(:,idx) = p1(1)*t.^4+p1(2)*t.^3+p1(3)*t.^2+p1(4)*t+p1(5);
            joint_trj2_p1(:,idx) = p2(1)*t.^4+p2(2)*t.^3+p2(3)*t.^2+p2(4)*t+p2(5);
            else
            joint_trj1_p1(:,idx) = ((THB_max(j)-THB_min(j))/2)*sin(10*t)+((THB_max(j)+THB_min(j))/2);
            joint_trj2_p1(:,idx) = -((THB_max(j)-THB_min(j))/2)*sin(10*t)+((THB_max(j)+THB_min(j))/2);
            end
            idx=idx+1;
        end
    elseif i == 2 %INDEX
        for j = 1:4
            if j<= 3
            p1 = polyfit([0 tf/4 tf/2 3*tf/4 tf],[0 IDX_max(j) 0 IDX_min(j) 0],4);
            p2 = polyfit([0 tf/4 tf/2 3*tf/4 tf],[0 IDX_min(j) 0 IDX_max(j) 0],4);
            joint_trj1_p1(:,idx) = p1(1)*t.^4+p1(2)*t.^3+p1(3)*t.^2+p1(4)*t+p1(5);
            joint_trj2_p1(:,idx) = p2(1)*t.^4+p2(2)*t.^3+p2(3)*t.^2+p2(4)*t+p2(5);
            else
            joint_trj1_p1(:,idx) = IDX_max(j)*cos(10*t);
            joint_trj2_p1(:,idx) = -IDX_max(j)*cos(10*t);
            end
            idx=idx+1;
        end
    elseif i == 3 %MIDDLE
        for j = 1:4
            if j<=3
            p1 = polyfit([0 tf/4 tf/2 3*tf/4 tf],[0 MID_max(j) 0 MID_min(j) 0],4);
            p2 = polyfit([0 tf/4 tf/2 3*tf/4 tf],[0 MID_min(j) 0 MID_max(j) 0],4);
            joint_trj1_p1(:,idx) = p1(1)*t.^4+p1(2)*t.^3+p1(3)*t.^2+p1(4)*t+p1(5);
            joint_trj2_p1(:,idx) = p2(1)*t.^4+p2(2)*t.^3+p2(3)*t.^2+p2(4)*t+p2(5);
            else
            joint_trj1_p1(:,idx) = MID_max(j)*cos(10*t);
            joint_trj2_p1(:,idx) = -MID_max(j)*cos(10*t);
            end
            idx=idx+1;
        end
    elseif i == 4 %RING
        for j = 1:4
            if j<=3
            p1 = polyfit([0 tf/4 tf/2 3*tf/4 tf],[0 RIG_max(j) 0 RIG_min(j) 0],4);
            p2 = polyfit([0 tf/4 tf/2 3*tf/4 tf],[0 RIG_min(j) 0 RIG_max(j) 0],4);
            joint_trj1_p1(:,idx) = p1(1)*t.^4+p1(2)*t.^3+p1(3)*t.^2+p1(4)*t+p1(5);
            joint_trj2_p1(:,idx) = p2(1)*t.^4+p2(2)*t.^3+p2(3)*t.^2+p2(4)*t+p2(5);
            else
            joint_trj1_p1(:,idx) = RIG_max(j)*cos(10*t);
            joint_trj2_p1(:,idx) = -RIG_max(j)*cos(10*t);
            end
            idx=idx+1;
        end
    elseif i == 5 %LITTLE
        for j = 1:4
            if j<=3
            p1 = polyfit([0 tf/4 tf/2 3*tf/4 tf],[0 LIT_max(j) 0 LIT_min(j) 0],4);
            p2 = polyfit([0 tf/4 tf/2 3*tf/4 tf],[0 LIT_min(j) 0 LIT_max(j) 0],4);
            joint_trj1_p1(:,idx) = p1(1)*t.^4+p1(2)*t.^3+p1(3)*t.^2+p1(4)*t+p1(5);
            joint_trj2_p1(:,idx) = p2(1)*t.^4+p2(2)*t.^3+p2(3)*t.^2+p2(4)*t+p2(5);
            else
            joint_trj1_p1(:,idx) = LIT_max(j)*cos(10*t);
            joint_trj2_p1(:,idx) = -LIT_max(j)*cos(10*t);
            end
            idx=idx+1;
        end
    end
end

%Random Joint Exploration
t2 = (tf+ts):ts:tT;

joint_trj1_p2 = [t2',zeros(length(t2),20)];
joint_trj2_p2 = [t2',zeros(length(t2),20)];

joint_trj1_p2(:,1) = t2';
joint_trj2_p2(:,1) = t2';

%Repeat 10 times for improving the results
for k = 1:10
    % Joint Explorations Randomly from Min to Max
    idx = 2;
    for i = 1:5
        if i == 1
            for j = 1:4 %THUMB
                joint_trj1_p2(:,idx) = rand([length(t2),1])*(THB_max(j)-THB_min(j))+THB_min(j);
                joint_trj2_p2(:,idx) = rand([length(t2),1])*(THB_max(j)-THB_min(j))+THB_min(j);
                idx=idx+1;
            end
        elseif i == 2 %INDEX
            for j = 1:4
                joint_trj1_p2(:,idx) = rand([length(t2),1])*(IDX_max(j)-IDX_min(j))+IDX_min(j);
                joint_trj2_p2(:,idx) = rand([length(t2),1])*(IDX_max(j)-IDX_min(j))+IDX_min(j);
                idx=idx+1;
            end
        elseif i == 3 %MIDDLE
            for j = 1:4
                joint_trj1_p2(:,idx) = rand([length(t2),1])*(MID_max(j)-MID_min(j))+MID_min(j);
                joint_trj2_p2(:,idx) = rand([length(t2),1])*(MID_max(j)-MID_min(j))+MID_min(j);
                idx=idx+1;
            end
        elseif i == 4 %RING
            for j = 1:4
                joint_trj1_p2(:,idx) = rand([length(t2),1])*(RIG_max(j)-RIG_min(j))+RIG_min(j);
                joint_trj2_p2(:,idx) = rand([length(t2),1])*(RIG_max(j)-RIG_min(j))+RIG_min(j);
                idx=idx+1;
            end
        elseif i == 5 %LITTLE
            for j = 1:4
                joint_trj1_p2(:,idx) = rand([length(t2),1])*(LIT_max(j)-LIT_min(j))+LIT_min(j);
                joint_trj2_p2(:,idx) = rand([length(t2),1])*(LIT_max(j)-LIT_min(j))+LIT_min(j);
                idx=idx+1;
            end
        end
    end

    joint_trj1 = [joint_trj1_p1;joint_trj1_p2];
    joint_trj2 = [joint_trj2_p1;joint_trj2_p2];

    %Simulation
    out = sim(simulation_path,'StartTime','0','StopTime',num2str(tT),'FixedStep',num2str(ts));
    
    % Extraction of Maximum Distance
    max_distances(:,k) = [max(out.d_thb);max(out.d_idx);max(out.d_mid);max(out.d_rig);max(out.d_lit);
                          max(out.d_thb_dip);max(out.d_idx_dip);max(out.d_mid_dip);max(out.d_rig_dip);max(out.d_lit_dip);
                          max(out.d_thb_pip);max(out.d_idx_pip);max(out.d_mid_pip);max(out.d_rig_pip);max(out.d_lit_pip)];
end

% Max Average Finger Distance and Standard Deviation
mean_d_max = mean(max_distances,2);
std_d_max = std(max_distances');




%% EE 523 Assignment 01 - Sajjad Uddin Mahmud - Spring 2023 - WSU

%% Basic Initialization
clc;
clear all;
close all;

%% Adding Folder Path
Data_Path = [fileparts(pwd),'/Data'];
addpath(Data_Path);

Function_Path = [fileparts(pwd),'/Functions'];
addpath(Function_Path);

%% Setting Up The Input Data As Per Assignment
Problem = 'A';
if Problem == 'A'
    Excel_Worksheet = 'Problem_A';
elseif Problem == 'B'
    Excel_Worksheet = 'Problem_B';
elseif Problem == 'C'
    Excel_Worksheet = 'Problem_C';
end

%% Reading From Bus Data
%% Bus Number
All_Bus_Number = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'A4:A14'); % Reading All Bus ID Data
Total_Bus = length(All_Bus_Number); % Calculating Total Bus Number

%% Bus Type
All_Bus_Type = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'G4:G14'); % Reading All Bus Type Data
PQ_Bus_Type = 0;
PQ_Bus_Type_1 = 1;
PV_Bus_Type = 2; 
Slack_Bus_Type = 3;
Bus = All_Bus_Number(find(All_Bus_Type ~= Slack_Bus_Type)); % Bus Type Data Except the Slack Bus
Bus_Type = All_Bus_Type(Bus); % Bus Type Data Except the Slack Bus

%% Bus Information
% Slack_Bus_Number = 1

Base_MVA = 100; 
V_Desired = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'O4:O14'); % Given Desired Voltage
Delta_in_Rad = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'I4:I14'); % Given Voltage Angle
P_Load = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'J4:J14')/Base_MVA; % Load MW pu
Q_Load = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'K4:K14')/Base_MVA; % Load MVAR pu
P_Gen = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'L4:L14')/Base_MVA; % Generator MW pu
Q_Gen = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'M4:M14')/Base_MVA; % Generator MVAR pu

%% Reading from Branch Data
%% Branch Number
From_Bus = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'A18:A27');
To_Bus = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'B18:B27');

%% Bus Shunt Conductance and Shunt Susceptance
G_Shunt_Bus = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'R4:R14');
B_Shunt_Bus = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'S4:S14');

%% Calculating Bus Shunt Admittance
Y_Shunt_Bus = G_Shunt_Bus + j.*B_Shunt_Bus;

%% Branch Resistance Per Unit
R_Branch = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'G18:G27');

%% Branch Reactance Per Unit
X_Branch = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'H18:H27');

%% Line Charging B Per Unit
B_Branch = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'I18:I27');

%% Transformer Turns Ratio
XFR_TurnRatio = xlsread('Kundur_11Bus_System.xlsx',Excel_Worksheet,'O18:O27');

%% Calculating Branch Impedence and Admittance
for i=1:length(From_Bus)
    Z_Branch(i) = R_Branch(i) + j * X_Branch(i); % Per Unit Impedance
    Y_Branch(i) = 1 / Z_Branch(i); % Per Unit Admittance
end 

%% Tap Consideration
Tap_Consideration = 1; % 0 = Without Taps, 1 = With Taps
if (Tap_Consideration == 0)
    for i = 1:length(XFR_TurnRatio)
        XFR_TurnRatio(i) = 0; % If We Do Not Consider Tap, All the Turn Ratio of Transformer are 0
    end
end

%% Calculating Y Bus Matrix: 

% Initialization
Y_Bus = zeros(Total_Bus,Total_Bus);

% LOOP: Computing Off-Diagonal Elements
for i=1:length(Y_Branch)
    if (XFR_TurnRatio(i)==0)
        Y_Bus(From_Bus(i),To_Bus(i)) = - Y_Branch(i);
        Y_Bus(To_Bus(i),From_Bus(i)) = - Y_Branch(i);
        Y_Bus_Diag(From_Bus(i),To_Bus(i)) = - Y_Branch(i);
        Y_Bus_Diag(To_Bus(i),From_Bus(i)) = - Y_Branch(i);
    else
        T = (1/(XFR_TurnRatio(i)));
        Y_Bus(From_Bus(i),To_Bus(i)) = - Y_Branch(i) * (T);
        Y_Bus(To_Bus(i),From_Bus(i)) = - Y_Branch(i) * (T);
        Y_Bus_Diag(From_Bus(i),To_Bus(i)) = - Y_Branch(i);
        Y_Bus_Diag(To_Bus(i),From_Bus(i)) = - Y_Branch(i) * (T^2);
    end
end
 
% LOOP: Computing Diagonal Elements
Y_Bus_Sum = sum(Y_Bus_Diag);
for i=1:Total_Bus
Y_Bus(i,i) = -Y_Bus_Sum(i) + Y_Shunt_Bus(i); % Adding Shunt Capacitance
end

% LOOP: Adding Line Charaging Capacitance
for i=1:length(From_Bus)
    Y_Bus(From_Bus(i),From_Bus(i)) = Y_Bus(From_Bus(i),From_Bus(i)) + j * (B_Branch(i) / 2);
    Y_Bus(To_Bus(i),To_Bus(i)) = Y_Bus(To_Bus(i),To_Bus(i)) + j * (B_Branch(i) / 2);
end

% Converting Y Bus Data into Polar Form
Rho = abs(Y_Bus); % Magnitude of Y Bus Entries
Theta = angle(Y_Bus); % Angle of Y Bus Entries in radian
B = imag(Y_Bus); % Imaginary Part of Y Bus Entries
G = real(Y_Bus); % Real Part of Y Bus Entries

% End of Y Bus Formation. Y Bus is Ready

%% Power Flow

%% Method
    Task = 1; % Newton-Raphson = 1; Fast Decoupled = 2
    
    %% Power Flow - Newton-Raphson Method
    if Task == 1
        [V, Delta_in_Rad, Iteration] = Newton_Raphson_Function(Y_Bus, V_Desired, Delta_in_Rad, P_Gen, P_Load, Q_Gen, Q_Load, All_Bus_Number, All_Bus_Type);
    

    %%  Power Flow - Fast Decoupled
    elseif Task == 2
        [V_FD, Delta_in_Rad_FD, Iteration_FD] = Newton_Raphson_Function_1(Y_Bus, V_Desired, Delta_in_Rad, P_Gen, P_Load, Q_Gen, Q_Load, All_Bus_Number, All_Bus_Type); % Putting Values after 4 Iterations of NR as Input for Fast Decoupled
        
        V_FD = transpose(V_FD);
        Delta_in_Rad_FD = transpose(Delta_in_Rad_FD);

        [V, Delta_in_Rad, Iteration] = Fast_Decoupled_Function(Y_Bus, V_FD, Delta_in_Rad_FD, P_Gen, P_Load, Q_Gen, Q_Load, All_Bus_Number, All_Bus_Type);

    end
      

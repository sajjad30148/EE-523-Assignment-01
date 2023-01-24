%% EE 523 Assignment 01 - Sajjad Uddin Mahmud - Spring 2023 - WSU

%% Basic Initialization
clc;
clear all;
close all;


%% Setting Up The Input Data As Per Assignment
Problem = 'B';
if Problem == 'A'
    Excel_Worksheet = 'Problem_A';
end
if Problem == 'B'
    Excel_Worksheet = 'Problem_B';
end
if Problem == 'C'
    Excel_Worksheet = 'Problem_C';
end

%% Reading From Bus Data
%% Bus Number
All_Bus_Number = xlsread('Kundur_Two_Area_System.xlsx',Excel_Worksheet,'A4:A14'); % Reading All Bus ID Data
Total_Bus = length(All_Bus_Number); % Calculating Total Bus Number

%% Bus Type
All_Bus_Type = xlsread('Kundur_Two_Area_System.xlsx',Excel_Worksheet,'G4:G14'); % Reading All Bus Type Data
PQ_Bus_Type = 0;
PQ_Bus_Type_1 = 1;
PV_Bus_Type = 2; 
Slack_Bus_Type = 3;
Bus = All_Bus_Number(find(All_Bus_Type ~= Slack_Bus_Type)); % Bus Type Data Except the Slack Bus
Bus_Type = All_Bus_Type(Bus); % Bus Type Data Except the Slack Bus

%% Bus Information
% Slack_Bus_Number = 1

Base_MVA = 100; 
V_Desired = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','O4:O14'); % Given Desired Voltage
Delta_in_Rad = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','I4:I14'); % Given Voltage Angle
P_Load = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','J4:J14')/Base_MVA; % Load MW pu
Q_Load = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','K4:K14')/Base_MVA; % Load MVAR pu
P_Gen = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','L4:L14')/Base_MVA; % Generator MW pu
Q_Gen = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','M4:M14')/Base_MVA; % Generator MVAR pu

%% Reading from Branch Data
%% Branch Number
From_Bus = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','A18:A27');
To_Bus = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','B18:B27');

%% Bus Shunt Conductance and Shunt Susceptance
G_Shunt_Bus = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','R4:R14');
B_Shunt_Bus = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','S4:S14');

%% Calculating Bus Shunt Admittance
Y_Shunt_Bus = G_Shunt_Bus + j.*B_Shunt_Bus;

%% Branch Resistance Per Unit
R_Branch = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','G18:G27');

%% Branch Reactance Per Unit
X_Branch = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','H18:H27');

%% Line Charging B Per Unit
B_Branch = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','I18:I27');

%% Transformer Turns Ratio
XFR_TurnRatio = xlsread('Kundur_Two_Area_System.xlsx','Problem_A','O18:O27');

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

% Schedule Real and Reactive Power
P_Scheduled = transpose(P_Gen - P_Load);
Q_Scheduled = transpose(Q_Gen - Q_Load);

% Initial Voltage Magnitude
V = ones(1,length(All_Bus_Number));
V(1,find(V_Desired)) = V_Desired(find(V_Desired),1);

% Initialization
Iteration = 0;
Tolerance = 0.01;
while 1

    %% Calculating Real Power 
    
    % Initialization
    P_Calculated = zeros(1,Total_Bus);
    
    % LOOP: Computing Real Power
    for i=1:Total_Bus
        for n=1:Total_Bus
            P_Calculated(i) = P_Calculated(i) + (abs(abs(Y_Bus(i,n)) * V(i) * V(n))) * (cos(angle(Y_Bus(i,n)) + Delta_in_Rad(n) - Delta_in_Rad(i)));
        end
    end
    
    %% Calculating Reactive Power
    
    % Initialization
    Q_Calculated = zeros(1,Total_Bus);
    
    % LOOP: Computing Reactive Power 
    for i=1:Total_Bus
        for n=1:Total_Bus
            Q_Calculated(i) = Q_Calculated(i) + (abs(abs(Y_Bus(i,n)) * V(i) * V(n))) * (sin(angle(Y_Bus(i,n)) + Delta_in_Rad(n) - Delta_in_Rad(i)));
        end
        Q_Calculated(i) = - Q_Calculated(i);
    end
    
    %% Calculating Mismatch
    Delta_P = P_Scheduled - P_Calculated;
    Delta_Q = Q_Scheduled - Q_Calculated;
    
    %% Method
    Task = 2; % Newton-Raphson = 1; Fast Decoupled = 2
    
    %% Power Flow by Newton-Raphson Method
    if Task == 1
    
        % Initializating Jacobian Matrix
        J11 = zeros(length(Bus));
        J12 = zeros(length(Bus));
        J21 = zeros(length(Bus));
        J22 = zeros(length(Bus));
        
        % LOOP: Computing Jacobian Matrix for All the Buses Except Slack Bus
        for i=1:length(Bus)
            for j=1:length(Bus)
                if (i==j)
                    J11(i,j) = - Q_Calculated(Bus(i)) - ((V(Bus(i)))^2) * (imag(Y_Bus(Bus(i),Bus(i))));
                    J21(i,j) = P_Calculated(Bus(i)) - ((V(Bus(i)))^2) * (real(Y_Bus(Bus(i),Bus(i))));
                    J12(i,j) = P_Calculated(Bus(i)) + ((V(Bus(i)))^2) * (real(Y_Bus(Bus(i),Bus(i))));
                    J22(i,j) = Q_Calculated(Bus(i)) - ((V(Bus(i)))^2) * (imag(Y_Bus(Bus(i),Bus(i))));
                else
                    J11(i,j) = - abs(V(Bus(i)) * V(Bus(j)) * abs(Y_Bus(Bus(i),Bus(j)))) * sin(angle(Y_Bus(Bus(i),Bus(j))) + Delta_in_Rad(Bus(j)) - Delta_in_Rad(Bus(i)));
                    J21(i,j) = - abs(V(Bus(i)) * V(Bus(j)) * abs(Y_Bus(Bus(i),Bus(j)))) * cos(angle(Y_Bus(Bus(i),Bus(j))) + Delta_in_Rad(Bus(j)) - Delta_in_Rad(Bus(i)));
                    J12(i,j) = - J21(i,j);
                    J22(i,j) = J11(i,j);
                end
            end
        end
        
        % Removing Rows and Columns from Jacobian for PV Bus
        PV = find(Bus_Type==PV_Bus_Type);
        J12(:,PV) = [];
        J21(PV,:) = [];
        J22(:,PV) = [];
        J22(PV,:) = [];
        J = [J11 J12; J21 J22];
    
        % Delta
        Delta_J = Delta_in_Rad(find(All_Bus_Type ~= Slack_Bus_Type));
        V_J = V(find((All_Bus_Type == PQ_Bus_Type)));
        Delta_P_J = Delta_P(find(All_Bus_Type ~= Slack_Bus_Type));
        Delta_Q_J = Delta_Q(find((All_Bus_Type == PQ_Bus_Type)));
        Delta_P_Q = [transpose(Delta_P_J);transpose(Delta_Q_J)];
        
        %% Updating V and Delta through LU Factorization
    
        % Function Calling: LU Factorization Using Dolittle's Method
        [V_Delta_Corrected] = LU_Factorization_Dolittle_Function(J,Delta_P_Q);
    
        % LOOP: Sorting the Voltages and Angles after LU Factorization
        for i=1:length(V_Delta_Corrected)
            if (i <= length(Delta_P_J))
                Delta_Corrected(i) = V_Delta_Corrected(i);
            else
                V_Corrected(i-length(Delta_P_J)) = V_Delta_Corrected(i);
            end
        end
        
        % Updating Voltages and Angles
        Delta_Updated = Delta_J + Delta_Corrected;
        V_Updated = V_J .* (1 + V_Corrected);
        
        % Preparing for Next Iteration
        V_i = (find((All_Bus_Type == PQ_Bus_Type)));
        Delta_i = find(All_Bus_Type ~= Slack_Bus_Type);
        
        for i=1:length(Delta_i)
            Delta_New(Delta_i(i)) = Delta_Updated(i);
        end
        
        for i=1:length(V_i)
            V_Desired(V_i(i)) = V_Updated(i);
        end
        
        V = transpose(V_Desired);
        Delta_in_Rad = Delta_New;
        Delta_in_Degree = (180 / pi) * Delta_in_Rad;
        
        Iteration = Iteration + 1;
    
        %% Output
    
        %fprintf("YBus: \n")
        %Y_Bus
        
        fprintf("Number of Iteration: \n");
        Iteration
    
        fprintf("Voltage Magnitude: \n")
        V
    
        fprintf("Voltage Angles in Degree: \n")
        Delta_in_Degree
    
        fprintf("Real Power in MW: \n")
        P_Calculated * Base_MVA
    
        fprintf("Reactive Power in MVAR: \n")
        Q_Calculated * Base_MVA
        
        % Checking Tolerance Limit
        if (max(abs(Delta_P_J)) < Tolerance & max(abs(Delta_Q_J)) < Tolerance)
            break;
        end
        
    %%  Power Flow by Newton-Raphson Method - Fast Decoupled
    elseif Task == 2
    
        % Initialization of Jacobian in Fast Decoupled Method; J12=J21=0
        J11 = zeros(length(Bus));
        J22 = zeros(length(Bus));
        
        % LOOP: Computing Jacobian Matrix for All the Buses Except Slack Bus 
        for i=1:length(Bus)
            for j=1:length(Bus)
                    J11(i,j) = - (imag(Y_Bus(Bus(i),Bus(j))));
                    J22(i,j) = - (imag(Y_Bus(Bus(i),Bus(j))));
            end
        end
        
        % Removing Rows and Columns from Jacobian for PV Bus
        PV = find(Bus_Type == PV_Bus_Type);
        J22(:,PV) = [];
        J22(PV,:) = [];
       
        % Delta
        Delta_J = Delta_in_Rad(find(All_Bus_Type ~= Slack_Bus_Type));
        V_J = V(find((All_Bus_Type == PQ_Bus_Type)));
        Delta_P_J = Delta_P(find(All_Bus_Type ~= Slack_Bus_Type));
        Delta_Q_J = Delta_Q(find((All_Bus_Type == PQ_Bus_Type)));
    
        %% Updating V and Delta through LU Factorization
        
        % Function Calling: LU Factorization Using Dolittle's Method
        [Delta_Corrected] = LU_Factorization_Dolittle_Function(J11,Delta_P_J);
        [V_Corrected] = LU_Factorization_Dolittle_Function(J22,Delta_Q_J);
    
        % Updating Voltages and Angles
        Delta_Updated = Delta_J + transpose(Delta_Corrected);
        V_Updated = V_J.*(1  + transpose(V_Corrected));
        
        % Preparing for Next Iteration
        V_i = (find((All_Bus_Type == PQ_Bus_Type)));
        Delta_i = find(All_Bus_Type ~= Slack_Bus_Type);
    
        for i=1:length(Delta_i)
            Delta_New(Delta_i(i)) = Delta_Updated(i);
        end
        
        for i=1:length(V_i)
            V_Desired(V_i(i)) = V_Updated(i);
        end
        
        V = transpose(V_Desired);
        Delta_in_Rad = Delta_New;
        Delta_in_Degree = (180 / pi) * Delta_in_Rad;
        
        Iteration = Iteration + 1
    
        %% Output
    
        fprintf("YBus: \n")
        Y_Bus
        
        fprintf("Number of Iteration: \n")
        Iteration
    
        fprintf("Voltage Magnitude: \n")
        V
    
        fprintf("Voltage Angles in Degree: \n")
        Delta_in_Degree
    
        fprintf("Real Power in MW: \n")
        P_Calculated * Base_MVA
    
        fprintf("Reactive Power in MVAR: \n")
        Q_Calculated * Base_MVA
        
        % Checking Tolerance Limit
        if (max(abs(Delta_P_J)) < Tolerance & max(abs(Delta_Q_J)) < Tolerance)
            break;
        end
    
    end

end



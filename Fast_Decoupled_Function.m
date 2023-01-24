%% Power Flow Function: Fast Decoupled

function [V, Delta_in_Rad, Iteration] = Fast_Decoupled_Function(Y_Bus, V_Desired, Delta_in_Rad, P_Gen, P_Load, Q_Gen, Q_Load, All_Bus_Number, All_Bus_Type)

%% Basic Initialization
Total_Bus = length(All_Bus_Number);
PQ_Bus_Type = 0;
PQ_Bus_Type_1 = 1;
PV_Bus_Type = 2; 
Slack_Bus_Type = 3;
Bus = All_Bus_Number(find(All_Bus_Type ~= Slack_Bus_Type)); % Bus Type Data Except the Slack Bus
Bus_Type = All_Bus_Type(Bus); % Bus Type Data Except the Slack Bus
Base_MVA = 100; 

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
        
        Iteration = Iteration + 1;
    
        %% Output
    
        %fprintf("YBus: \n")
        %Y_Bus
        
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
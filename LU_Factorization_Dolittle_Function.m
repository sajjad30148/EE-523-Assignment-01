%% LU Factorization Function: Dolittle's Algorithm

function [ X_Matrix ] = LU_Factorization_Dolittle_Function(A_Matrix,B_Matrix)

% Getting the Size of Input Matrix
Length_A = length(A_Matrix);

% Initializing The Lower and Upper Triangular Matrices
Lower_Triangular_Matrix = zeros(Length_A,Length_A);
Upper_Triangular_Matrix = zeros(Length_A,Length_A);

% LOOP: Assigning 1 into All Diagonal Elements of Lower Traingular Matrix
for j = 1:Length_A
    Lower_Triangular_Matrix(j,j) = 1;
end

% Computing 1st Row of Upper Traingular Matrix
Upper_Triangular_Matrix(1,:) = A_Matrix(1,:);

% Computing 1st Column of Lower Traingular Matrix
Lower_Triangular_Matrix(:,1) = A_Matrix(:,1)/Upper_Triangular_Matrix(1,1);

% LOOP: Computing All Other Rows and Column of Upper and Lower Traingular Matrix
for j = 2:Length_A
    for k = j:Length_A
        Upper_Triangular_Matrix(j,k) = A_Matrix(j,k) - Lower_Triangular_Matrix(j,1:j-1) * Upper_Triangular_Matrix(1:j-1,k);
    end

    for l = j+1:Length_A
        Lower_Triangular_Matrix(l,j) = (A_Matrix(l,j) - Lower_Triangular_Matrix(l,1:j-1) * Upper_Triangular_Matrix(1:j-1,j)) / Upper_Triangular_Matrix(j,j);
    end
    
end

% Output
% A_Matrix
% Lower_Triangular_Matrix
% Upper_Triangular_Matrix

% Verification
% A_Matrix - (Lower_Triangular_Matrix * Upper_Triangular_Matrix)


%% Forward Substitution

% Initialization of Y Matrix
Y_Matrix = zeros(Length_A,1);

% Computing First Value of Y Matrix
Y_Matrix(1) = B_Matrix(1) / Lower_Triangular_Matrix(1,1);

% LOOP: Computing Rest of the Entries of Y Matrix
for j = 2:Length_A
    Y_Matrix(j) = (B_Matrix(j) - Lower_Triangular_Matrix(j,1:j-1) * Y_Matrix(1:j-1)) / Lower_Triangular_Matrix(j,j);
end

% Output
% Y_Matrix


%% Backward Substitution

% Initialization of X Matrix
X_Matrix = zeros(Length_A,1);

% Computing Last Value of X Matrix
X_Matrix(Length_A) = Y_Matrix(Length_A) / Upper_Triangular_Matrix(Length_A,Length_A);

% LOOP: Computing Rest of the Entries of X Matrix
for j = Length_A-1:-1:1 
    X_Matrix(j) = (Y_Matrix(j) - Upper_Triangular_Matrix(j,j+1:Length_A) * X_Matrix(j+1:Length_A)) / Upper_Triangular_Matrix(j,j);
end

% Output
% X_Matrix

end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paper Source
% "Step-by-Step HHL Algorithm Walkthrough to Enhance the Understanding of
% Critical Quantum Computing Concepts"
% https://arxiv.org/abs/2108.09004

%To cite

%Hector Morrell and Hiu Yung Wong, "Step-by-Step HHL Algorithm Walkthrough to Enhance the Understanding of Critical Quantum Computing Concepts", arXiv preprint arXiv:2108.09004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear all

display('This section is to calculate the neccessary components of the circuit');
% Classical Solution

%%%%%%%% Intiating the Matrix and the Vector %%%%%%%%
A = [1, -1/3; -1/3, 1] 

b = [0; 1]

%%%%%%%% Inverting the Matrix to find x %%%%%%%%
inv_A = inv(A)

x = inv_A*b


%%%%%%%% Find the expected ratio of measurement %%%%%%%%

norm_factor = norm(x)

x_norm = x/norm_factor

Prob_0 = x_norm(1,1)^2
Prob_1 = x_norm(2,1)^2

%%%%%%%% Initiazing t and finding the eigenvalues and vector of the Matrix %%%%%%%%
t = 3*pi/4

[v,e] = eig(A)

V_tran = transpose(v)

Md = v*A*V_tran

%% 
%%%%%%%% Diaginalizing the vectors needed for QPE rotation %%%%%%%%
y_diagonal = exp((i*Md*t));

y_diagonal(1,2) = 0;
y_diagonal(2,1) = 0;

% display(y_diagonal);

y_non_diagonal = V_tran*y_diagonal*v

y_diagonal_2 = exp((i*Md*t)*2);

y_diagonal_2(1,2) = 0;
y_diagonal_2(2,1) = 0;

% display(y_diagonal_2);

y_non_diagonal_2 = v*y_diagonal_2*V_tran

display('This section is to calculate qubit evolution step-by-step');

%% 
%%%%%%%%%%%%%%%% Quantum Gates 

I = [1,0;0,1];

H = (1/sqrt(2)).*[1,1;1,-1];

Not = [0,1;1,0];

%% 

%Initialization

Psi_0 = [1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]

%State Preparation

Psi_1 = kron(Not,kron(I,kron(I,I)))*Psi_0 

%% 
% Superposition of clock qubits

Psi_2 = kron(I, kron(kron(H,H),I))*Psi_1

%%

%%%%%%%%%%%%%%%%%%% Initiating the Control Gates

U2 = [0,-1;-1,0];

U1 = .5*[(-1+i),(1+i);...
        (1+i),(-1+i)];
Bra_0 = [1;0];
Bra_1 = [0;1];

TBra_0 = transpose(Bra_0);
TBra_1 = transpose(Bra_1);

outer_00 = Bra_0*TBra_0;
outer_11 = Bra_1*TBra_1;

%%%%%%%%%%%% Passing through first control

CU1 = kron(kron(kron(I,I),outer_00)+ kron(kron(U1,I),outer_11),I);

Psi_31 = CU1*Psi_2

%%%%%%%%%% Passing through Second control

CU2 = kron(kron(kron(I,outer_00),I)+ kron(kron(U2,outer_11),I),I);

Psi_3 = CU2*Psi_31

%%%%%%%%%%%%%% Inverse QFT
Swap = [1,0,0,0;...
         0,0,1,0;...
         0,1,0,0;...
         0,0,0,1];
Ru = [1,0,0,0;...
         0,1,0,0;...
         0,0,1,0;...
         0,0,0,-i];

% QFT Verification    
y = kron(H,I);    
y2 = kron(I,H);     
Gate_QFT = Swap*y2*Ru*y;
Gate_Final = kron(I,Gate_QFT);

% QFT Matrix Multiplication
inter1 = kron(kron(I,kron(H,I)),I);
inter2 = kron(kron(I,Ru),I);
inter3 = kron(kron(I,kron(I,H)),I);
inter4 = kron(kron(I, Swap),I);
Psi_4=inter4*inter3*inter2*inter1*Psi_3


%% 

%%%%%%%%%%%%%%%%%%%% Ancilla Rotation C = 1

%%{

%%%%%%%% Ry rotation 01
theta1 = 2*asin(1/1);
RY01 = [cos(theta1/2),-sin(theta1/2);...
        sin(theta1/2),cos(theta1/2)];

CRY01 = kron(I,kron(I,kron(outer_00,I))+kron(I,kron(outer_11,RY01)));

%%%%%%% Ry rotation 10
theta2 = 2*asin(1/2);

%%%%%% FIx this error Here
RY10 = [cos(theta2/2),-sin(theta2/2);...
        sin(theta2/2),cos(theta2/2)];

CRY10 = kron(I,kron(outer_00,kron(I,I))+ kron(outer_11,kron(I,RY10)));


Psi_5 = CRY10*CRY01*Psi_4
%}

%Measurement is not performed in this notebook
Psi_6=Psi_5

% IQFT

%%%%%%%%%%%%%% QFT
Ru_inv = [1,0,0,0;...
         0,1,0,0;...
         0,0,1,0;...
         0,0,0,i];
   
Gate_QFT = y*Ru_inv*y2*Swap;
% Inverse QFT Matrix Multiplication
inter4i = kron(kron(I,kron(H,I)),I);
inter3i = kron(kron(I,Ru_inv),I);
inter2i = kron(kron(I,kron(I,H)),I);
inter1i = kron(kron(I, Swap),I);
Psi_7=inter4i*inter3i*inter2i*inter1i*Psi_6

%%%%%%%%%%%%%%%%%%% Initiating the Inverse Control Gates

inv_U2 = [0,-1;-1,0];

inv_U1 = .5*[(-1-i),(1-i);...
        (1-i),(-1-i)];
    
%%%%%%%%%% Passing through Square control

inv_CU2 = kron(kron(I,outer_00),I)+ kron(kron(inv_U2,outer_11),I);

Gate_inv_CU2 = kron(inv_CU2,I);

%%%%%%%%%%%% Passing through Power 1 control

inv_CU1 = kron(kron(I,I),outer_00)+ kron(kron(inv_U1,I),outer_11);

Gate_inv_CU1 = kron(inv_CU1,I);

Psi_8 = Gate_inv_CU1*Gate_inv_CU2*Psi_7

% Apply Hadamard Gates

Psi_9=kron(I,kron(kron(H,H),I))*Psi_8


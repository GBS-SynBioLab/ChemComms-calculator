% recursive solver for the N rook problem
% it calculates and saves all possible chess board configurations
% the code is based on: https://stackoverflow.com/a/10289790
% ZAT 2018, Imperial College London

clear variables
clc

% global variables to keep track of the solutions
global boards
global solution_idx

% size of the problem
number_of_rooks = 8;
rows            = 8;
columns         = 8;
% the chessboard where the rooks are placed
board           = zeros(rows,columns);
% calculating the number of solutions to allocate memory
num_of_sol = nchoosek(rows,number_of_rooks)*nchoosek(columns,number_of_rooks)*factorial(number_of_rooks);
boards = zeros(rows,columns,num_of_sol);
solution_idx = 1;

%% start the recursion 
tic;
recurse_rooks(1,number_of_rooks,board);
total = toc;

assert(num_of_sol == size(boards,3))

%% generate report
fprintf('board size: %dx%d, rooks: %d, number of solutions: %d in %g sec\n',rows,columns,number_of_rooks,size(boards,3),total)
% This function generates all possible combinations of devices and AHLs for
% the given number of channels
% output: 3D matrix, each column is a possible device and AHL pair, the third slice contains the all possible combinations 
%
% ZAT 2018 Imperial College London
function all_possible_channel_comb = get_possible_channels(num_of_channels,num_of_devices,num_of_AHLs)

% global variables to keep track of the solutions
global boards
global solution_idx

% size of the problem
number_of_rooks = num_of_channels;
rows            = num_of_devices;
columns         = num_of_AHLs;
% the chessboard where the rooks are placed
board           = zeros(rows,columns);
% calculating the number of solutions to allocate memory
num_of_sol   = nchoosek(rows,number_of_rooks)*nchoosek(columns,number_of_rooks)*factorial(number_of_rooks);
boards       = zeros(rows,columns,num_of_sol);
solution_idx = 1;

%% start the recursion 
tic;
recurse_rooks(1,number_of_rooks,board);
total = toc;

assert(num_of_sol == size(boards,3))

%% process results
[r,c,~] = ind2sub(size(boards),find(boards));

all_possible_channel_comb = reshape([r,c]',2,num_of_channels,num_of_sol);

end
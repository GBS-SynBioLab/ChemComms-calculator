% ZAT 2018 Imperial College London
function recurse_rooks(row, numberOfRooks,board)
global boards
global solution_idx
if (numberOfRooks == 0)
    boards(:,:,solution_idx) = board;
    solution_idx = solution_idx + 1;
    return;
end
for k=row:size(board,1)
    for l = 1:size(board,2)
        if sum(board(k,:)) == 0 && sum(board(:,l)) == 0
            board(k,l) = 1;
            recurse_rooks(k+1, numberOfRooks-1,board)
            board(k,l) = 0;
        end
    end
end
end
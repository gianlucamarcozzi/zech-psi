%% Linear/non-linear fitting model
function winham = windowhamming(L, iStart)
    % Apply Hamming window function to a vector long L only from the index iStart onwards.
    winham = zeros(1, L);
    winham(1:iStart - 1) = ones(1, iStart - 1);
    tempham = hamming(2*(L - iStart));
    winham(iStart:end) = tempham(end/2:end);
end
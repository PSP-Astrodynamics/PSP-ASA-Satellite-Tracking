function val = nonuniqueAngle(array)
    % Round values to 4 decimal places to handle numerical errors
    array = mod(real(array), 2*pi);
    roundedArray = round(array, 8);
    
    % Find unique values and their counts
    [uniqueVals, ~, indices] = unique(roundedArray);
    counts = histc(indices, 1:numel(uniqueVals));
    
    % Identify values that are repeated
    val = uniqueVals(counts > 1);
end

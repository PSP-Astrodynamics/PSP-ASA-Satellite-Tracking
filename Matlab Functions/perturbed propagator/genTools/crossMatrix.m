function f_bar = crossMatrix(f_vec)
    f_bar = [0, -f_vec(3), f_vec(2);...
    f_vec(3), 0, -f_vec(1);...
    -f_vec(2), f_vec(1), 0];
end
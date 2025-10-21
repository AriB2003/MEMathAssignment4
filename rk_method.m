function BT_struct = rk_method(name)
    BT_struct = struct();
    switch name
        case "midpoint"
            BT_struct.A = [0, 0; 0.5, 0];
            BT_struct.B = [0; 1];
            BT_struct.C = [0; 0.5];
        case "kutta3rd"
        BT_struct.A = [0, 0, 0;
                       0.5, 0, 0;
                       -1, 2, 0];
        BT_struct.B = [1/6; 2/3; 1/6];
        BT_struct.C = [0; 0.5; 1];
        case "nystrom5th"
            BT_struct.A = [0, 0, 0, 0, 0, 0; 
                           1/4, 0, 0, 0, 0, 0;
                           4/25, 6/25, 0, 0, 0, 0;
                           1/4, -3, 15/4, 0, 0, 0;
                           2/27, 10/9, -50/81, 8/81, 0, 0;
                           2/25, 12/25, 2/15, 8/75, 0, 0];
            BT_struct.B = [23/192; 0; 125/192; 0; -27/64; 125/192];
            BT_struct.C = [0; 1/3; 2/5; 1; 2/3; 4/5];
    end
end
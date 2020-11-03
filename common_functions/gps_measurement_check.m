function new_data = gps_measurement_check(k, k_gps)
if k <= k_gps
    new_data = false;
elseif mod(k, k_gps) == 0
    new_data = true;
else
    new_data = false;
end
end
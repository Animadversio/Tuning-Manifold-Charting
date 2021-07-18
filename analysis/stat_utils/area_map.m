function area = area_map(chan)
if chan<=48 & chan>=33
area = "V1";
elseif chan>48
area = "V4";
elseif chan<33
area = "IT";
end
end
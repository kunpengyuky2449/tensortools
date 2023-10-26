all_b = [];
for i = 1:width(rez.dshift)
    all_b = cat(1, all_b, rez.dshift(:, i));
end
plot(rez.dshift)
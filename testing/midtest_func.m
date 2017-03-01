function result = midtest_func(n)
% at first i will put one thing then switch it mid-way

pause(5);
calc = sum(sum(magic(n)));
result = sprintf('the second thing i put and result is %d', calc);

end
function out=test_rand1(seed)

rs = RandStream('mt19937ar');
RandStream.setGlobalStream(rs);
reset(rs, seed);

for iter=1:5
    out(iter,:)=randperm(10);
end

end
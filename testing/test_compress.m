function out=test_compress(in)

in_attribs=whos('in');

out=dzip(in);

out_attribs=whos('out');

ratio=out_attribs.bytes/in_attribs.bytes;

fprintf('The compressed file is %0.2f the size of the input\n',ratio);

out_unzipped=dunzip(out);

if isequal(in,out_unzipped)
    disp('The zipping algorithm was truly lossless')
else
    disp('The zipping algorithm was NOT!!! truly lossless')
end

end
function refreshdisp(str,prevstr,iteration)

if ~exist('iteration','var')
    iteration=2;
end

if iteration==1
    fprintf(str);
else
    fprintf(repmat('\b',1,length(prevstr)));
    fprintf(str);
end


function q = GetCartesianVariable(fid,length_block_metadata,length_block,ndims,sof)

global string_length;

npts = fread(fid,ndims,'int32');

if sof == 4
    floatstring = 'float32';
elseif sof == 8
    floatstring = 'float64';
end
    
stagger = fread(fid,ndims,floatstring);
extents = fread(fid,2,floatstring);

q.('meshname') = char(fread(fid,string_length,'uchar'))';
q.('meshclass') = char(fread(fid,string_length,'uchar'))';
fseek(fid,length_block_metadata-4*ndims-ndims*sof-2*sof-2*string_length,'cof');

q.('data') = fread(fid,prod(npts),floatstring);

if ndims > 1
    q.('data') = reshape(q.('data'),npts');
end


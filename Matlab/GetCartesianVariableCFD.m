function q = GetCartesianVariable(fid, length_block_metadata, ...
        length_block, ndims, sof);

global string_length;

npts = fread(fid, ndims, 'int32');

if sof == 4
    typestring = 'float32';
elseif sof == 8
    typestring = 'float64';
end

stagger = fread(fid, ndims, typestring);
extents = fread(fid, 2, typestring);

q.meshname = char(fread(fid, string_length, 'uchar'))';
q.meshclass = char(fread(fid, string_length, 'uchar'))';
fseek(fid, length_block_metadata - (4 + sof) * ndims - 2 * sof ...
    - 2 * string_length, 'cof');

q.data = fread(fid, prod(npts), typestring);

if ndims > 1
    q.data = reshape(q.data, npts');
end

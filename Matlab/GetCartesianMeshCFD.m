function q = GetCartesianMeshCFD(fid, length_block_metadata, length_block, ...
        ndims, sof);

npts = fread(fid, ndims, 'int32');

if sof == 4
    typestring = 'float32';
elseif sof == 8
    typestring = 'float64';
end

extents = fread(fid, 2 * ndims, typestring);
q.x = fread(fid, npts(1), typestring);
if (ndims >= 2)
    q.y = fread(fid, npts(2), typestring);
end
if (ndims >= 3)
    q.z = fread(fid, npts(3), typestring);
end

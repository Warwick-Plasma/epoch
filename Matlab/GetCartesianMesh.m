function q = GetCartesianMesh(fid,length_block_metadata,length_block,ndims,sof)

npts = fread(fid,ndims,'int32');

if sof == 4
    floatstring = 'float32';
elseif sof == 8
    floatstring = 'float64';
end

extents = fread(fid,ndims*2,floatstring);
q.('x') = fread(fid,npts(1),floatstring);
if (ndims >= 2)
    q.('y') = fread(fid,npts(2),floatstring);
end
if (ndims >=3)
        q.('z') = fread(fid,npts(3),floatstring);
end

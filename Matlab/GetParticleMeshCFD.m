function q = GetParticleMesh(fid, length_block_metadata, length_block, ...
        ndims, sof);

PARTICLE_CARTESIAN = 0;
PARTICLE_POLAR = 1;
PARTICLE_CYLINDRICAL = 2;

coord_type = fread(fid, 1, 'int32');
npart = fread(fid, 1, 'int64');

if sof == 4
    typestring = 'float32';
elseif sof == 8
    typestring = 'float64';
end

extents = fread(fid, 2 * ndims, typestring);
q.data = fread(fid, [npart ndims], typestring);

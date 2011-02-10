function q = GetMesh(fid, length_block_metadata, length_block);

MESH_CARTESIAN = 0;
MESH_PARTICLE = 1;

type = fread(fid, 1, 'int32');
ndims = fread(fid, 1, 'int32');
sof = fread(fid, 1, 'int32');

if type == MESH_CARTESIAN
    q = GetCartesianMesh(fid, length_block_metadata - 12, ...
        length_block - 12, ndims, sof);
elseif type == MESH_PARTICLE
    q = GetParticleMesh(fid, length_block_metadata - 12, ...
        length_block - 12, ndims, sof);
else
    fseek(fid, length_block - 12, 'cof');
    q = 0;
end

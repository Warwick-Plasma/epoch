function q = GetMeshVariable(fid, length_block_metadata, length_block);

VAR_CARTESIAN = 0;
VAR_PARTICLE = 1;

type = fread(fid, 1, 'int32');
ndims = fread(fid, 1, 'int32');
sof = fread(fid, 1, 'int32');

if type == VAR_CARTESIAN
    q = GetCartesianVariable(fid, length_block_metadata - 12, ...
        length_block - 12, ndims, sof);
elseif type == VAR_PARTICLE
    q = GetParticleVariable(fid, length_block_metadata - 12, ...
        length_block - 12, ndims, sof);
else
    fseek(fid, length_block - 12, 'cof');
    q = 0;
end

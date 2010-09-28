function q = GetParticleVariable(fid,length_block_metadata,length_block,ndims,sof)

global string_length;

npart = fread(fid,1,'int64');

if sof == 4
    floatstring = 'float32';
elseif sof == 8
    floatstring = 'float64';
end
    
extents = fread(fid,2,floatstring);

q.('meshname') = char(fread(fid,string_length,'uchar'))';
q.('meshclass') = char(fread(fid,string_length,'uchar'))';
fseek(fid,length_block_metadata-8-2*sof-2*string_length,'cof');

q.('data') = fread(fid,npart,floatstring);

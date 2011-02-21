function q = GetDataCFD(filename)

%%%%%%%%%%%%%%%%%

global string_length;

TYPE_SCRIBBLE = -1;
TYPE_ADDITIONAL = 0;
TYPE_MESH = 1;
TYPE_MESH_VARIABLE = 2;
TYPE_SNAPSHOT = 3;
TYPE_STITCHED_VECTOR = 4;
TYPE_STITCHED_MAGNITUDE = 5;
TYPE_CONSTANT = 6;
TYPE_ARB_DB = 7;

MESH_CARTESIAN = 0;
MESH_PARTICLE = 1;

PARTICLE_CARTESIAN = 0;
PARTICLE_POLAR = 1;
PARTICLE_CYLINDRICAL = 2;

DIMENSION_IRRELEVANT = 0;
DIMENSION_1D = 1;
DIMENSION_2D = 2;
DIMENSION_3D = 3;

%%%%%%%%%%%%%%%%

fid = fopen(filename);

if fid == -1; disp('bad filename'); q = 'fail'; return; end

% File header
cfd_marker = char(fread(fid, 3, 'uchar'))';
length_file_header = fread(fid, 1, 'int32');
length_block_header = fread(fid, 1, 'int32');
version = fread(fid, 2, 'int32');
string_length = fread(fid, 1, 'int32');
num_blocks = fread(fid, 1, 'int32');

if ~(version(1) == 1 && version(2) == 0)
    endianness = fread(fid, 1, 'int32');
    start_sec = fread(fid, 1, 'int32');
    start_millisec = fread(fid, 1, 'int32');
    step = fread(fid, 1, 'int32');
    time = fread(fid, 1, 'float64');
    q.step = step;
    q.time = time;
end

% Now skipping rest of file header

fseek(fid, length_file_header, 'bof');

q.version = version;

for c = 1:num_blocks
    name = char(fread(fid, string_length, 'uchar'))';
    class = char(fread(fid, string_length, 'uchar'))';
    block_type = fread(fid, 1, 'int32');
    length_block_metadata = fread(fid, 1, 'int64');
    length_block = fread(fid, 1, 'int64');

    tag = deblank(name);
    tag = regexprep(tag, '\^', 'pow');
    tag = regexprep(tag, ' ', '_');
    tag = regexprep(tag, '/', '_');
    tag = regexprep(tag, '\W', '');

    blockknown = 0;
    if block_type == TYPE_SNAPSHOT
        if (version(1) == 1 && version(2) == 0)
            q.(tag) = GetSnapshot(fid, length_block_metadata, length_block);
            blockknown = 1;
        end
    elseif block_type == TYPE_CONSTANT
        q.(tag) = GetConstant(fid, length_block_metadata, length_block);
        blockknown = 1;
    elseif block_type == TYPE_MESH_VARIABLE
        q.(tag) = GetMeshVariable(fid, length_block_metadata, length_block);
        blockknown = 1;
    elseif block_type == TYPE_MESH
        q.(tag) = GetMesh(fid, length_block_metadata, length_block);
        blockknown = 1;
    end

    if blockknown == 0 % look for next block
        fseek(fid, length_block, 'cof');
    end
end

fclose(fid);

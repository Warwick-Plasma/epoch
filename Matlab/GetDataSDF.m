function q = GetDataSDF(filename)

global block;

%%%%%%%%%%%%%%%%%

h.ID_LENGTH = 32;
h.ENDIANNESS = 16911887;
h.VERSION = 1;
h.REVISION = 2;
h.MAGIC = 'SDF1';

h.BLOCKTYPE.SCRUBBED = -1;
h.BLOCKTYPE.NULL = 0;
h.BLOCKTYPE.PLAIN_MESH = 1;
h.BLOCKTYPE.POINT_MESH = 2;
h.BLOCKTYPE.PLAIN_VARIABLE = 3;
h.BLOCKTYPE.POINT_VARIABLE = 4;
h.BLOCKTYPE.CONSTANT = 5;
h.BLOCKTYPE.ARRAY = 6;
h.BLOCKTYPE.RUN_INFO = 7;
h.BLOCKTYPE.SOURCE = 8;
h.BLOCKTYPE.STITCHED_TENSOR = 9;
h.BLOCKTYPE.STITCHED_MATERIAL = 10;
h.BLOCKTYPE.STITCHED_MATVAR = 11;
h.BLOCKTYPE.STITCHED_SPECIES = 12;
h.BLOCKTYPE.FAMILY = 13;

h.BLOCKTYPE_NAME = { 'Invalid block'; 'Plain mesh'; 'Point mesh'; ...
    'Plain variable'; 'Point variable'; 'Constant'; 'Simple array'; ...
    'Run information'; 'Source code'; 'Stitched tensor'; ...
    'Stitched material'; 'Stitched material variable'; 'Stitched species'; ...
    'Particle family' };

h.DATATYPE.NULL = 0;
h.DATATYPE.INTEGER4 = 1;
h.DATATYPE.INTEGER8 = 2;
h.DATATYPE.REAL4 = 3;
h.DATATYPE.REAL8 = 4;
h.DATATYPE.REAL16 = 5;
h.DATATYPE.CHARACTER = 6;
h.DATATYPE.LOGICAL = 7;
h.DATATYPE.OTHER = 8;

%%%%%%%%%%%%%%%%

h.filename = filename;
h.fid = fopen(filename);

if h.fid == -1; disp('bad filename'); q = 'fail'; return; end

% File header
sdf_magic = char(fread(h.fid, 4, 'uchar'))';
endianness = fread(h.fid, 1, 'int32');
version = fread(h.fid, 1, 'int32');
revision = fread(h.fid, 1, 'int32');
code_name = deblank(strtrim(char(fread(h.fid, h.ID_LENGTH, 'uchar'))'));
first_block_location = fread(h.fid, 1, 'int64');
summary_location = fread(h.fid, 1, 'int64');
summary_size = fread(h.fid, 1, 'int32');
nblocks = fread(h.fid, 1, 'int32');
h.block_header_length = fread(h.fid, 1, 'int32');
q.step = fread(h.fid, 1, 'int32');
q.time = fread(h.fid, 1, 'float64');
jobid1 = fread(h.fid, 1, 'int32');
jobid2 = fread(h.fid, 1, 'int32');
string_length = fread(h.fid, 1, 'int32');
code_io_version = fread(h.fid, 1, 'int32');

% Now seek to first block
b.block_start = first_block_location;

for n = 1:nblocks
    fseek(h.fid, b.block_start, 'bof');
    b.next_block_location = fread(h.fid, 1, 'uint64');
    b.data_location = fread(h.fid, 1, 'uint64');
    b.id = deblank(strtrim(char(fread(h.fid, h.ID_LENGTH, 'uchar'))'));
    b.data_length = fread(h.fid, 1, 'uint64');
    b.blocktype = fread(h.fid, 1, 'uint32');
    b.datatype = fread(h.fid, 1, 'uint32');
    b.ndims = fread(h.fid, 1, 'uint32');
    b.name = deblank(strtrim(char(fread(h.fid, string_length, 'uchar'))'));
    b.mesh_id = '';
    b.var = 0;
    b.map = 0;

    block = b;

    switch b.blocktype
      case h.BLOCKTYPE.PLAIN_MESH
        block.var = GetPlainMeshSDF(h);
      case h.BLOCKTYPE.POINT_MESH
        block.var = GetPointMeshSDF(h);
      case h.BLOCKTYPE.PLAIN_VARIABLE
        block.var = GetPlainVariableSDF(h);
      case h.BLOCKTYPE.POINT_VARIABLE
        block.var = GetPointVariableSDF(h);
      case h.BLOCKTYPE.CONSTANT
        block.var = GetConstantSDF(h);
    end

    blocklist(n) = block;
    b.block_start = b.next_block_location;
end

for n = 1:nblocks
    b = blocklist(n);

    count = 0;
    rem = regexprep(b.name, ' ', '_');
    while size(rem) > 0
        [str, rem] = strtok(rem, '/');
        count = count + 1;
        s{count} = regexprep(str, '\W', '');
    end

    hasgrid = 0;
    add = 0;
    switch b.blocktype
      case h.BLOCKTYPE.PLAIN_MESH
          add = 1;
      case h.BLOCKTYPE.POINT_MESH
          add = 1;
      case h.BLOCKTYPE.PLAIN_VARIABLE
          add = 1;
          hasgrid = h.BLOCKTYPE.PLAIN_MESH;
      case h.BLOCKTYPE.POINT_VARIABLE
          add = 1;
          hasgrid = h.BLOCKTYPE.POINT_MESH;
      case h.BLOCKTYPE.CONSTANT
          add = 1;
    end

    got = 0;
    if hasgrid
        for gn=1:nblocks
            g = blocklist(gn);
            if g.blocktype == hasgrid && strcmp(b.mesh_id, g.id)
                grid = g.var;
                gname = regexprep(g.name, ' ', '_');
                gname = regexprep(gname, '/', '.');
                got = 1;
                break
            end
        end
    end

    if add
      switch count
        case 1
          q.(s{1}) = b.var;
          if got
              q.(s{1}).grid = grid;
              q.(s{1}).grid.name = gname;
          end
        case 2
          q.(s{1}).(s{2}) = b.var;
          if got
              q.(s{1}).(s{2}).grid = grid;
              q.(s{1}).(s{2}).grid.name = gname;
          end
        case 3
          q.(s{1}).(s{2}).(s{3}) = b.var;
          if got
              q.(s{1}).(s{2}).(s{3}).grid = grid;
              q.(s{1}).(s{2}).(s{3}).grid.name = gname;
          end
        case 4
          q.(s{1}).(s{2}).(s{3}).(s{4}) = b.var;
          if got
              q.(s{1}).(s{2}).(s{3}).(s{4}).grid = grid;
              q.(s{1}).(s{2}).(s{3}).(s{4}).grid.name = gname;
          end
        case 5
          q.(s{1}).(s{2}).(s{3}).(s{4}).(s{5}) = b.var;
          if got
              q.(s{1}).(s{2}).(s{3}).(s{4}).(s{5}).grid = grid;
              q.(s{1}).(s{2}).(s{3}).(s{4}).(s{5}).grid.name = gname;
          end
        case 6
          q.(s{1}).(s{2}).(s{3}).(s{4}).(s{5}).(s{6}) = b.var;
          if got
              q.(s{1}).(s{2}).(s{3}).(s{4}).(s{5}).(s{6}).grid = grid;
              q.(s{1}).(s{2}).(s{3}).(s{4}).(s{5}).(s{6}).grid.name = gname;
          end
      end
    end
end

fclose(h.fid);

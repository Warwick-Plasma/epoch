function q = GetPlainVariableSDF(h);

global block;

fseek(h.fid, block.block_start + h.block_header_length, 'bof');

mult = fread(h.fid, 1, 'float64');
units = deblank(strtrim(char(fread(h.fid, h.ID_LENGTH, 'uchar'))'));
block.mesh_id = deblank(strtrim(char(fread(h.fid, h.ID_LENGTH, 'uchar'))'));
npts = fread(h.fid, block.ndims, 'int32');
stagger = fread(h.fid, 1, 'int32');

if block.datatype == h.DATATYPE.REAL4
    typestring = 'single';
elseif block.datatype == h.DATATYPE.REAL8
    typestring = 'double';
end

offset = block.data_location;

tagname = 'data';
block.map = memmapfile(h.filename, 'Format', ...
        {typestring npts' tagname}, 'Offset', offset, ...
        'Repeat', 1, 'Writable', false);
q.(tagname) = block.map.data.(tagname);

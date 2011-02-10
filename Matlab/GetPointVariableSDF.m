function q = GetPointVariableSDF(h);

global block;

fseek(h.fid, block.block_start + h.block_header_length, 'bof');

mult = fread(h.fid, 1, 'float64');
units = deblank(strtrim(char(fread(h.fid, h.ID_LENGTH, 'uchar'))'));
block.mesh_id = deblank(strtrim(char(fread(h.fid, h.ID_LENGTH, 'uchar'))'));
npart = fread(h.fid, 1, 'int64');

if block.datatype == h.DATATYPE.REAL4
    typestring = 'single';
elseif block.datatype == h.DATATYPE.REAL8
    typestring = 'double';
end

offset = block.data_location;

tagname = 'data';
block.map = memmapfile(h.filename, 'Format', ...
        {typestring npart tagname}, 'Offset', offset, ...
        'Repeat', 1, 'Writable', false);
q.(tagname) = block.map.data.(tagname);

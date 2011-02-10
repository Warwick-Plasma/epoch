function constant = GetConstant(fid, length_block_metadata, length_block);

constant = fread(fid, 1, 'float64');

fseek(fid, length_block_metadata - 8, 'cof');

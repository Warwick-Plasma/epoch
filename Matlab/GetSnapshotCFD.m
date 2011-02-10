function q = GetSnapshotCFD(fid, length_block_metadata, length_block);

q.Snapshot = fread(fid, 1, 'int32');
q.Time = fread(fid, 1, 'float64');

fseek(fid, length_block_metadata - 12, 'cof');

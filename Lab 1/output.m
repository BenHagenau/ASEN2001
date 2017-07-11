function [] = output(outfile_path, data_cell)
  % fprintf(fileID,formatSpec,A1,...,An) writes to a file nicely


  fileID = fopen(outfile_path, 'w'); % discard previous contents of file

  if(fileID == -1) % coudln't open file
    fileID = 1;    % redirect to STDOUT
  end

  num_lines = size(data_cell, 1);

  for i = 1:num_lines
    force_or_moment = data_cell{i, 1};
    magnitude       = data_cell{i, 2};
    fprintf(fileID, '%s | %d\n', force_or_moment, magnitude);
  end

  fclose(fileID);
end

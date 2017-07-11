% Copyright Jacob Killelea 2016
% MIT License
% returns the next line that doesn't begin with a #
% if this function is called repeatedly, it will return a new line each time, until EOF
function line = next_non_comment_line(fileID)
  line = fgets(fileID);
  while (line(1) == '#')
    line = fgets(fileID);
  end
end

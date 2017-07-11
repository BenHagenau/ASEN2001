% Copyright Jacob Killelea, Ben Hagenau
% MIT License
% NOTE:: str2num('2 1') => [2, 1]. Very useful.
% function [num_forces, force_application_coords, force_vector_coords, num_moments, moment_application_coords, moment_vector_coords, NUM_SUPPORTS, support_coords, support_reaction_data] = get_file_input(filepath) % I'm kinda sorry. No other way to return all that data without lumping it together
function force_application_coords = get_file_input(filepath) % I'm kinda sorry. No other way to return all that data without lumping it together
  NUM_SUPPORTS = 6; %const? can't be determined from input data

  fileID = fopen(filepath); % open file

  % first data line is the number of forces and the number of moments. Cast them as ints
  line = str2num(next_non_comment_line(fileID));
  num_forces  = line(1);
  num_moments = line(2);

  % allocate matrices for storing data about the forces
  force_application_coords = cell(num_forces, 3);
  force_vector_coords      = cell(num_forces, 4);

  % allocate matrices for storing data about the moments
  moment_application_coords = cell(num_moments, 3);
  moment_vector_coords      = cell(num_moments, 4);
  
   % allocate matrices for storing datat about the supports
  support_coords        = cell(NUM_SUPPORTS, 3);
  support_reaction_data = cell(NUM_SUPPORTS, 4);

% FORCE DATA
  for i = 1:num_forces % get the coordinates of where the forces are applied
    force_application_coords(i, :) = num2cell(double(str2num(next_non_comment_line(fileID))));
  end
  for i = 1:num_forces % get the vector data (magnitude, direction)
    force_vector_coords(i, :) = num2cell(double(str2num(next_non_comment_line(fileID))));
  end

% MOMENT DATA
  for i = 1:num_moments % get the coordinates of where the moments are applied
    moment_application_coords(i, :) = num2cell(double(str2num(next_non_comment_line(fileID))));
  end
  for i = 1:num_moments % get their data (magnitude, direction)
    moment_vector_coords(i, :) = num2cell(double(str2num(next_non_comment_line(fileID))));
  end

% SUPPORT DATA
  for i = 1:NUM_SUPPORTS % get the coordinates of the supports
    support_coords(i, :) = num2cell(double(str2num(next_non_comment_line(fileID))));
  end

  for i = 1:NUM_SUPPORTS % get the data about the reactions at each support (force or moment, vector coords)
    tmp = str2mat(next_non_comment_line(fileID)); % get the line as a character array (since it's mixed characters and numbers)
    tmp = strsplit(tmp);    % split into cell, which eliminates whitespace nicely (ex: 'F' '1.0' '6.0' '-7.0' '')
                            %                                                                                 ^ extra column for some reason
    for j = 1:length(tmp)-1 % ignore the last column in tmp, since it doesn't contain any data
      support_reaction_data{i, j} = tmp{j}; % enter the data into each row of the support_reaction_data cell-array
    end
    for j = 2:length(tmp)-1
      support_reaction_data{i,j} = str2double(support_reaction_data{i,j}); % cast everything except the first column as a double, since that's the form of data we want.
    end
  end
  % disp(support_reaction_data);

end
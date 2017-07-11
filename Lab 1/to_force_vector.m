function vector = to_force_vector(force)
    % Returns a vector based on the direction and magnitude of the input vector, [f, dx, dy, dz]
    % [x, y, z] = to_force_vector([f, dx, dy, dz])
    % the magnitudes of dx, dy, and dz do not matter. Only their porportions relative to each other are important.
    % The returned vector is of magnitude f.
    % this function is input-agnostic, and could be called with parameter {f, dx, dy, dz}
    if isa(force, 'cell') % if force is cell, change to vector
      force = cell2mat(force);
    end

    magnitude_force = force(1); % get value inside vector
    vector = magnitude_force .* to_unit_vector(force(2:4));
end

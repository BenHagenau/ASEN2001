%create the b matrix of known moments and forces
function [force_mat,moment_mat,b] = force_vector(force_vector_coords, num_forces, force_application_coords, moment_application_coords, num_moments,moment_vector_coords)
%cell to matrix
force_vector_coords = cell2mat(force_vector_coords);
force_application_coords = cell2mat(force_application_coords);
moment_vector_coords = cell2mat(moment_vector_coords);
moment_application_coords = cell2mat(moment_application_coords);
%preallocate force vector matrix
force_mat = zeros(num_forces,3);
position_f = zeros(num_forces,3);
moment_mat = zeros(num_forces+num_moments,3);
%determine unit vector of forces
for i=1:num_forces
      dx = force_vector_coords(i,2);
      dy = force_vector_coords(i,3);
      dz = force_vector_coords(i,4);
%determine unit vector components and save in matrix
    vector_magnitude = magnitude(dx, dy, dz);
      x = dx/vector_magnitude;
      y = dy/vector_magnitude;
      z = dz/vector_magnitude;
%multiply the magnitude of the forces by the force unit vector components
     force_mat(i,:) = [x, y, z]*force_vector_coords(i,1);
%compute position vectors of forces
position_f(i,:) = force_application_coords(i,:);
%compute moment: r x F
moment_mat(i,:)=cross(position_f(i,:),force_mat(i,:));
end

%determine moments given in file
%preallocate moment vector matrix
given_moment_mat = zeros(num_forces,3);
position_m = zeros(num_forces,3);
%determine unit vector of forces
for i=1:num_moments
      dx = moment_vector_coords(i,2);
      dy = moment_vector_coords(i,3);
      dz = moment_vector_coords(i,4);
%determine unit vector components and save in matrix
    vector_magnitude = magnitude(dx, dy, dz);
      x = dx/vector_magnitude;
      y = dy/vector_magnitude;
      z = dz/vector_magnitude;
%multiply the magnitude of the forces by the force unit vector components
     moment_mat(num_forces+i,:) = [x, y, z]*moment_vector_coords(i,1);
end





%sum the components of the the momemnts and the forces to create matrix x
%preallocate x
b = zeros(6,1);
%sum the force components and place them in matrix x [Fx;Fy;Fz]
for i=1:3
    b(i) = sum(force_mat(:,i));
end
%sum the moment components and place them in matrix x [Mx;My;Mz]
for i=1:3
    b(i+3) = sum(moment_mat(:,i));
end

end

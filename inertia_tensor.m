%% helper functions (don't use/edit - only used to help initialize I)
function tensor = inertia_tensor(num)

n = num2str(num);

tensor = [sym(['Ixx' n]) sym(['Ixy' n]) sym(['Ixz' n]);
          sym(['Iyx' n]) sym(['Iyy' n]) sym(['Iyz' n]);
          sym(['Izx' n]) sym(['Izy' n]) sym(['Izz' n])];

assume(tensor, 'real');

end
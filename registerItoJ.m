function J = registerItoJ(ops, I)

dsnew  = registration_offsets(I, ops, 0);
J  = register_movie(I, ops, dsnew);
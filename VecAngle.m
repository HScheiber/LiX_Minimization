% Function that outputs the angle in degrees between two vectors
function Angle = VecAngle(a_vec,b_vec)
    Angle = atan2d(norm(cross(a_vec,b_vec)),dot(a_vec,b_vec));
end

function clear_persistents()
% Clear the persistent variables used in nonlinear optimization

f_clear_persistent = true;

deriv_patterns([],[],[],[], f_clear_persistent);
interpolate_image_stacks([],[],[], f_clear_persistent);
brightness_error_vec([],[],[], f_clear_persistent);
smoothness_error_vec([],[],[], f_clear_persistent);

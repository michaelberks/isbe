#include "ncm_apt_null_controller.h"

#include <vcl_cassert.h>
#include <vcl_iostream.h>

//
// Public member functions
//

ncm_apt_null_controller::ncm_apt_null_controller()
{
}

//
//: Return ID of the controller.
long ncm_apt_null_controller::id() const
{
  do_nothing();
  return -1L;
}

//
//: Make LED on controller flash
void ncm_apt_null_controller::identify() const
{
  do_nothing();
}

////
////: Set offset of home position from endstop
//void ncm_apt_null_controller::set_home_offset(float home_offset)
//{
//  do_nothing();
//}

//: Reverse the direction of motion.
void ncm_apt_null_controller::set_reverse(bool reverse /* = true */)
{
  do_nothing();
}

//
//: Go to the zero position (with respect to the reverse limit switch).
void ncm_apt_null_controller::move_zero(bool return_immediately /* = true */)
{
  do_nothing();
}
 
//: Set offset of home position from endstop.
bool ncm_apt_null_controller::set_home_position(float home_offset)
{
  return true;
}

//
//: Go to the home position. The offset from the limit switch can be set by
//  giving home_offset a positive value.
void ncm_apt_null_controller::move_home(float home_offset,
                                        bool return_immediately /* = true */)
{
  do_nothing();
}

//
//: Go to the home position. The offset from the limit switch can be set by
//  giving home_offset a positive value.
void ncm_apt_null_controller::move_home(bool return_immediately /* = true */)
{
  do_nothing();
}

//
//: Move stage to an absolute position with respect to home position.
int ncm_apt_null_controller::move_to(float absolute_position,
                                 bool return_immediately /* = true */)
{
  do_nothing();
  return -1;
}

//
//: Move stage by a relative distance with respect to current position.
int ncm_apt_null_controller::move_by(float relative_distance, 
                                 bool return_immediately /* = true */)
{
  do_nothing();
  return -1;
}

//
//: Move stage at a 
int ncm_apt_null_controller::move_at_velocity(float velocity)
{
  do_nothing();
  return -1;
}

//: Move the stage at a specified multiple of the maximum velocity.
void ncm_apt_null_controller::move_at_relative_velocity(float velocity_scale)
{
  do_nothing();
}
 
//: Buffer the current position.
void ncm_apt_null_controller::buffer_position()
{
  do_nothing();
}

//
//: Return position (in mm) of the stage with respect to its home position.
float ncm_apt_null_controller::position() const
{
  do_nothing();
  return 0.0f;
}

//
//: Return position in absolute coordinates (i.e. relative to endstop)
float ncm_apt_null_controller::absolute_position() const
{
  do_nothing();
  return 12.5f;
}

//
//: Return maximum position allowed
float ncm_apt_null_controller::min_position() const
{
  do_nothing();
  return 0.0f;
}

//
//: Return midpoint of travel
float ncm_apt_null_controller::mid_position() const
{
  do_nothing();
  return 0.0f;
}

//
//: Return minimum position allowed
float ncm_apt_null_controller::max_position() const
{
  do_nothing();
  return 25.0f;
}

//
//: Return distance to nearest endstop.
float ncm_apt_null_controller::distance_to_endstop() const
{
  do_nothing();
  return 0.0f;
}

//
//: Return estimated instantaneous velocity
float ncm_apt_null_controller::velocity()
{
  do_nothing();
  return 0.0f;
}

//
//: Set the maximum velocity for moves
void ncm_apt_null_controller::set_velocity(float velocity)
{
  do_nothing();
}

//
//: Set the maximum acceleration for moves
void ncm_apt_null_controller::set_acceleration(float acceleration)
{
  do_nothing();
}


//
//: Stop gracefully.
int ncm_apt_null_controller::stop() const
{
  do_nothing();
  return -1;
}

//
//: Return if homed.
bool ncm_apt_null_controller::homed() const
{
  do_nothing();
  return true;
}

//
//: Return if reversed.
bool ncm_apt_null_controller::reversed() const
{
  do_nothing();
  return true;
}

//
// Private member functions
//

void ncm_apt_null_controller::do_nothing() const
{
  // Do nothing! (For now.)
}
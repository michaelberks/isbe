#include "ncm_apt_controller.h"

#include <vcl_cassert.h>
#include <vcl_iostream.h>
#include <vcl_cmath.h>

#include <nailfold/ncm_encrypt.h>

// Header files needed for the APT API
#include <windows.h>
#include <APTAPI.h>

//
// Public member functions
//

ncm_apt_controller::ncm_apt_controller(long serial_number)
: id_(serial_number),
  reversed_(false),
  velocity_(1.0f),
  acceleration_(1.0f),
  position_(0.0f),
  home_offset_(0.0f), // read from device immediately
  home_velocity_(0.0f), // read from device immediately
  max_velocity_(0.0),
  max_acceleration_(0.0),
  min_position_(-1.0), // read from device immediately
  max_position_(-1.0), // read from device immediately
  // Constants
  min_position_abs_(0.0f),
  max_position_abs_(25.0f),
  min_velocity_(0.0),
  min_acceleration_(0.0),
	homed_(false)
{
  assert(id_is_valid());

  // It's questionable how useful this lot is. Most of the parameters are reset
  // to predefined values during APTInit() anyway.

  InitHWDevice(id_);

  // Get limits on acceleration and velocity.
  // These are then fixed for the class and should not be changed.
  MOT_GetVelParamLimits(id_, &max_acceleration_, &max_velocity_);

  // Get home move parameters
  long home_direction = 0; // not used
  long home_lim_switch = 0; // not used
  MOT_GetHomeParams(id_, &home_direction, &home_lim_switch, 
                         &home_velocity_, &home_offset_);

	//Set backlash to zero
	float backlash_dist = 0;
	MOT_SetBLashDist(id_, backlash_dist);
	
  // Set home velocity.
  home_velocity_ = 2.0f;
  MOT_SetHomeParams(id_, HOME_REV, HOMELIMSW_REV,
                         home_velocity_, home_offset_);

  // Get limits on position (and others of less interest)
  // It's questionable how useful these limits are, since they are redefined
  // when you set the home position.
  long units;
  float pitch;
  MOT_GetStageAxisInfo(id_, &min_position_, &max_position_, 
                            &units, &pitch);

  //MOT_SetVelParams(serialNumber, /* min velocity = */ 0.0f, 
  //                               /* acceleration = */ maxAcceleration, 
  //                               /* max velocity = */ 0.5f);

  // Check if the stage needs homing (a lot of commands won't work if the 
  // stage isn't homed, though there's no good reason why they shouldn't).
  // If the position is *exactly* zero, it is likely that the controller
  // had been disconnected and needs rehoming. 
  MOT_GetPosition(id_, &position_);
  if (vcl_abs(position_) == 0)
    move_home();
}

//
//: Return ID of controller
long ncm_apt_controller::id() const
{
  return id_;
}

//
//: Make LED on controller flash
void ncm_apt_controller::identify() const
{
  assert(id_is_valid());

  MOT_Identify(id_);
}

//
//: Reverse the direction of motion for move_by() and move_at_velocity().
//  This is useful if the platform is at 180 degrees to the camera coordinate
//  frame and you want to align the two.
void ncm_apt_controller::set_reverse(bool reverse /* = true */)
{
  reversed_ = reverse;
}
 
//
//: Go to the zero position (with respect to the reverse limit switch).
void ncm_apt_controller::move_zero(bool return_immediately /* = true */)
{
  assert(id_is_valid());

  // Store copy of current home position
  float home_offset = 0.0f;
  long direction = 0L;
  long switchtype = 0L;
  MOT_GetHomeParams(id_, &direction, &switchtype,
                         &home_velocity_, &home_offset);

  // Set home position to zero and move there
  const float zero_offset = 0.0f;
  MOT_SetHomeParams(id_, HOME_REV, HOMELIMSW_REV,
                         home_velocity_, zero_offset);

  const bool bWait = !return_immediately;
  long return_code = MOT_MoveHome(id_, bWait);

  // Reinstate previous home position
  MOT_SetHomeParams(id_, HOME_REV, HOMELIMSW_REV,
                         home_velocity_, home_offset);
}
 
//: Set the home position, and update the limits of travel (from the home
//  position) accordingly.
bool ncm_apt_controller::set_home_position(float home_offset)
{
  if ((min_position_abs_ <= home_offset) && 
                           (home_offset <= max_position_abs_))
  {
    home_offset_ = home_offset;

    // Get units and pitch
    long units;
    float pitch;
    MOT_GetStageAxisInfo(id_, &min_position_, &max_position_, 
                              &units, &pitch);

    // Redefine limits relative to home position
    min_position_ = min_position_abs_ - home_offset_;
    max_position_ = max_position_abs_ - home_offset_;

    MOT_SetStageAxisInfo(id_, min_position_, max_position_,
                              units, pitch);

    // Set new home position in controller.
    MOT_SetHomeParams(id_, HOME_REV, HOMELIMSW_REV,
                           home_velocity_, home_offset_);

    return true;
  }
  else
  {
    vcl_cout << "Home position not in valid range ["
             << min_position_ << ", " << max_position_ 
             << "]" << vcl_endl;

    return false;
  }
}
 
//: Go to the existing home position.
void ncm_apt_controller::move_home(bool return_immediately /* = true */)
{
  assert(id_is_valid());

  vcl_cout << "Homing " << id_ << vcl_endl;

  // Move to zero position first.
  const bool bWait = !return_immediately;
  long return_code = MOT_MoveHome(id_, bWait);

	//Set homed flag
	homed_ = true;
}
 
//: Go to a given home position.
void ncm_apt_controller::move_home(float home_offset,
                                   bool return_immediately /* = true */)
{
  assert(id_is_valid());

  if ( set_home_position(home_offset) )
    move_home(return_immediately);
}
 
//: Move stage to an absolute position with respect to home position.
int ncm_apt_controller::move_to(float position,
                                bool return_immediately /* = true */)
{
  assert(id_is_valid());

  const bool bWait = !return_immediately;
  int return_code = MOT_MoveAbsoluteEx(id_, position, bWait);

  return return_code;
}

//
//: Move stage by a relative distance with respect to current position.
int ncm_apt_controller::move_by(float relative_distance, 
                                bool return_immediately /* = true */)
{
  assert(id_is_valid());
  
  if (reversed_)
    relative_distance = -relative_distance;

	float new_position = position_ + relative_distance;
	if (new_position > min_position_ && new_position < max_position_)
	{
		const bool bWait = !return_immediately;
		int return_code = MOT_MoveRelativeEx(id_, relative_distance, bWait);

		return return_code;
	}
	else return -1;
}

//
//: Move stage at a given velocity
int ncm_apt_controller::move_at_velocity(float velocity)
{
  assert(id_is_valid());

  int return_code = -1;

  // I don't know why the direction of movement is opposite to what you'd
  // expect for the MT25 stage.
  if (velocity >= 0.0f)
  {
    set_velocity(velocity);

    if (reversed_ && position_ > min_position_)
      return_code = MOT_MoveVelocity(id_, MOVE_FWD);
    else if (position_ < max_position_)
      return_code = MOT_MoveVelocity(id_, MOVE_REV);
  }
  else //if (velocity < 0.0f)
  {
    set_velocity(-velocity);

    if (reversed_ && position_ < max_position_)
      return_code = MOT_MoveVelocity(id_, MOVE_REV);
    else if (position_ > min_position_)
      return_code = MOT_MoveVelocity(id_, MOVE_FWD);
  }

  return return_code;
}

void ncm_apt_controller::move_at_relative_velocity(float velocity_scale)
{
  // move_at_velocity() will apply appropriate limits.
  move_at_velocity( velocity_scale * max_velocity_ );
}
 
//: Read the position from the motors and store it for future use.
//  This is surprisingly expensive, and should be run within the APT thread.
void ncm_apt_controller::buffer_position()
{
  assert(id_is_valid());

  MOT_GetPosition(id_, &position_);
}

//
//: Return position (in mm) of the stage with respect to its home position.
float ncm_apt_controller::position() const
{
  return position_;
}

//
//: Return position in absolute coordinates (i.e. relative to endstop)
float ncm_apt_controller::absolute_position() const
{
  return position_ - min_position_;
}

//
//: Return minimum position allowed
float ncm_apt_controller::min_position() const
{
  return min_position_;
}

//
//: Return midpoint of travel
float ncm_apt_controller::mid_position() const
{
  return 0.5f * (min_position_ + max_position_);
}

//
//: Return maximum position allowed
float ncm_apt_controller::max_position() const
{
  return max_position_;
}

//
//: Return distance to nearest endstop.
float ncm_apt_controller::distance_to_endstop() const
{
  const float d_rev = position_ - min_position_;
  const float d_fwd = max_position_ - position_;

  if (d_rev < d_fwd)
    return d_rev;
  else
    return d_fwd;
}

//
//: Return estimated instantaneous velocity
float ncm_apt_controller::velocity()
{
  return 0.0f;
}

//
//: Set the maximum velocity for moves
void ncm_apt_controller::set_velocity(float velocity)
{
  assert(id_is_valid());
  
  const bool too_high = (velocity > max_velocity_);
  const bool too_low = (velocity < min_velocity_);

  if (!too_high && !too_low)
    velocity_ = velocity;
  else if (too_high)
    velocity_ = max_velocity_;
  else
    velocity_ = min_velocity_;

  MOT_SetVelParams(id_, min_velocity_, acceleration_, velocity_);
}

//
//: Set the maximum acceleration for moves
void ncm_apt_controller::set_acceleration(float acceleration)
{
  assert(id_is_valid());
  
  const bool too_high = (acceleration > max_acceleration_);
  const bool too_low = (acceleration < min_acceleration_);

  if (!too_high && !too_low)
    acceleration_ = acceleration;
  else if (too_high)
    acceleration_ = max_acceleration_;
  else
    acceleration_ = min_acceleration_;

  MOT_SetVelParams(id_, min_velocity_, acceleration_, velocity_);
}

//
//: Stop gracefully.
int ncm_apt_controller::stop() const
{
  assert(id_is_valid());
  
  int return_code = MOT_StopProfiled(id_);

  return return_code; 
}

//
//: Return if homed.
bool ncm_apt_controller::homed() const
{
  return homed_;
}

//
//: Return if reversed
bool ncm_apt_controller::reversed() const
{
	return reversed_;
}

//
// Private member functions
//

//
//: Prevent construction of a controller without an ID by making the default
//  constructor private.
ncm_apt_controller::ncm_apt_controller()
: // Constants
  min_position_abs_(0.0f),
  max_position_abs_(25.0f),
  min_velocity_(0.0),
  min_acceleration_(0.0)
{
}

//
//: Return true only if ID is within a valid range.
bool ncm_apt_controller::id_is_valid() const
{
  return ((83000000 <= id_) && (id_ <=83999999));
}


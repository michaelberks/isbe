#ifndef ncm_apt_controller_base_h_
#define ncm_apt_controller_base_h_

#include <vcl_string.h>

class ncm_apt_controller_base
{
// INTERFACE

public:
  // No member variables here, please

  //: Default constructor
  ncm_apt_controller_base();

  //: Return ID of the controller.
  virtual long id() const = 0;

  //: Make controller LED flash.
  virtual void identify() const = 0;

  //: Reverse the direction of motion.
  virtual void set_reverse(bool reverse = true) = 0;

  //: Move stage to the zero position (with respect to the reverse limit switch).
  virtual void move_zero(bool return_immediately = true) = 0;

  //: Set offset of home position from endstop.
  virtual bool set_home_position(float home_offset) = 0;

  //: Move stage to the existing home position.
  virtual void move_home(bool return_immediately = true) = 0;

  //: Move stage to the given home position.
  virtual void move_home(float home_offset,
                         bool return_immediately = true) = 0;

  //: Move stage to given position (relative to home position).
  virtual int move_to(float position,
                      bool return_immediately = true) = 0;

  //: Move stage by a relative distance.
  virtual int move_by(float relative_distance,
                      bool return_immediately = true) = 0;

  //: Move the stage at a specified velocity until stopped.
  virtual int move_at_velocity(float velocity) = 0;

  //: Move the stage at a specified multiple of the maximum velocity.
  virtual void move_at_relative_velocity(float velocity_scale) = 0;

  //: Buffer the current position.
  virtual void buffer_position() = 0;

  //: Return position relative to home position.
  virtual float position() const = 0;

  //: Return position in absolute coordinates (i.e. relative to endstop).
  virtual float absolute_position() const = 0;

  virtual float min_position() const = 0;
  virtual float mid_position() const = 0;
  virtual float max_position() const = 0;

  //: Return distance to nearest endstop.
  virtual float distance_to_endstop() const = 0;

  //: Return estimated instantaneous velocity.
  virtual float velocity() = 0;

  //: Set the maximum velocity for moves.
  virtual void set_velocity(float velocity) = 0;

  //: Set the acceleration for moves.
  virtual void set_acceleration(float acceleration) = 0;

  //: Stop gracefully.
  virtual int stop() const = 0;

	//: Return if homed.
  virtual bool homed() const = 0;

	//: Return if homed.
  virtual bool reversed() const = 0;


// IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes

private:
};

#endif // ncm_apt_controller_base_h_
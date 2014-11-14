#ifndef ncm_apt_controller_h_
#define ncm_apt_controller_h_

#include <vcl_string.h>

#include <nailfold/ncm_apt_controller_base.h>

class ncm_apt_controller : public ncm_apt_controller_base
{
// INTERFACE

public:
  // No member variables here, please

  //: Default ctor
  ncm_apt_controller(long serial_number);

  //: Return ID of the controller.
  virtual long id() const;

  //: Make controller LED flash.
  virtual void identify() const;

  //: Reverse the direction of motion.
  virtual void set_reverse(bool reverse = true);

  //: Move stage to the zero position (with respect to the reverse limit switch).
  virtual void move_zero(bool return_immediately = true);

  //: Set offset of home position from endstop.
  virtual bool set_home_position(float home_offset);

  //: Move stage to the existing home position.
  virtual void move_home(bool return_immediately = true);

  //: Move stage to a given home position.
  virtual void move_home(float home_offset,
                         bool return_immediately = true);
  
  //: Move stage to given position (relative to home position).
  virtual int move_to(float position,
                      bool return_immediately = true);

  //: Move stage by a relative distance.
  virtual int move_by(float relative_distance,
                      bool return_immediately = true);

  //: Move the stage at a specified velocity until a stopped.
  virtual int move_at_velocity(float velocity);

  //: Move the stage at a specified multiple of the maximum velocity.
  virtual void move_at_relative_velocity(float velocity_scale);

  //: Buffer the current position.
  virtual void buffer_position();

  //: Return position relative to home position.
  virtual float position() const;

  //: Return position in absolute coordinates (i.e. relative to endstop).
  virtual float absolute_position() const;

  virtual float min_position() const;
  virtual float mid_position() const;
  virtual float max_position() const;

  //: Return distance to nearest endstop.
  virtual float distance_to_endstop() const;

  //: Return estimated instantaneous velocity.
  virtual float velocity();

  //: Set the maximum velocity for moves.
  virtual void set_velocity(float velocity);

  //: Set the acceleration for moves.
  virtual void set_acceleration(float acceleration);

  //: Stop gracefully.
  virtual int stop() const;

	//: Return if homed.
  virtual bool homed() const;

	//: Return if homed.
  virtual bool reversed() const;


// IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes

private:
  // Members and functions visible only to objects of this class

  //: Prevent construction without a serial number by making the default
  //  constructor private.
  ncm_apt_controller();

  //: Check if ID is within a valid range
  virtual bool id_is_valid() const;


  //: Serial number of the controller
  long id_;

  //: Direction of movement in velocity move commands.
  bool reversed_;

  //: Offset of 'home' position from the REV limit switch
  float home_offset_;

  //: Maximum velocity at which to move to the home position
  float home_velocity_;

  //: Minimum position allowed (determined by polling the stage)
  float min_position_;

  //: Maximum position allowed (determined by polling the stage)
  float max_position_;

  //: Minimum position allowed (determined by polling the stage)
  const float min_position_abs_;

  //: Maximum position allowed (determined by polling the stage)
  const float max_position_abs_;

  //: Buffered position.
  float position_;

  //: Minimum velocity allowed by the stage/controller (fixed at zero).
  const float min_velocity_;

  //: Maximum velocity allowed by the stage/controller (determined by polling 
  //  the stage/controller).
  float max_velocity_;

  //: Maximum velocity reached in moves
  float velocity_;

  //: Minimum acceleration allowed by the stage/controller.
  const float min_acceleration_;

  //: Maximum acceleration allowed by the stage/controller (determined by 
  //  polling the stage/controller).
  float max_acceleration_;

  //: Acceleration to use in trapezoid profile.
  float acceleration_;

	//: Flag to set if homed (by this software)
	bool homed_;
};

#endif // ncm_apt_controller_h_
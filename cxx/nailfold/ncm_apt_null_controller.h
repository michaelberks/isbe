#ifndef ncm_apt_null_controller_h_
#define ncm_apt_null_controller_h_

#include <vcl_string.h>

#include <nailfold/ncm_apt_controller_base.h>

class ncm_apt_null_controller : public ncm_apt_controller_base
{
// INTERFACE

public:
  // No member variables here, please

  //: Default ctor
  ncm_apt_null_controller();

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
  
  //: Move stage to the home position.
  virtual void move_home(bool return_immediately = true);

  //: Move stage to the home position.
  virtual void move_home(float home_offset,
                         bool return_immediately = true);

  //: Move stage to given position (relative to home position).
  virtual int move_to(float absolute_position,
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

	//: Return if reversed.
  virtual bool reversed() const;


// IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes

private:
  // Members and functions visible only to objects of this class

  //: Do nothing special (though it can be modified to do 
  //  something - such as throw an error - at a later date).
  void do_nothing() const;
};

#endif // ncm_apt_null_controller_h_
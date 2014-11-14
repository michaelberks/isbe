#ifndef ncm_apt_server_h_
#define ncm_apt_server_h_

#include <vcl_string.h>
#include <vcl_vector.h>

// Forward declarations
class ncm_apt_controller_base;
class ncm_apt_controller;
class ncm_apt_null_controller;

class ncm_apt_server
{
// INTERFACE

public:
  //: Default constructor
  ncm_apt_server();

  //: Destructor
  ~ncm_apt_server();

  //: Return number of controllers in list.
  unsigned n_controllers() const;

  //: Return vector of controller IDs available.
  vcl_vector<long> ids();

  //: Assign a controller to the X axis
  void set_X_id(long id);

  //: Assign a controller to the Y axis
  void set_Y_id(long id);
  
  //: Assign a controller to the Z axis
  void set_Z_id(long id);

  //: Get pointer to controller assigned to X axis
  ncm_apt_controller_base* X() const;

  //: Get pointer to controller assigned to Y axis
  ncm_apt_controller_base* Y() const;

  //: Get pointer to controller assigned to Z axis
  ncm_apt_controller_base* Z() const;

  //: Move stages to (x,y,z) location
  void move_to_xyz(float x, float y, float z);

  //: Move all stages to zero position
  void move_all_zero(bool return_immediately = true);

  //: Move all stages to home position
  void move_all_home(bool return_immediately = true);

  //: Enable/Disable event dialog box.
  void set_event_dialogs(bool dialogs_enabled = true) const;

	//: Populate controllers_ with new instances of ncm_apt_controller.
  void create_controllers();

	//: Clear controllers_.
  void destroy_controllers();

// IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes

	

private:
  // Members and functions visible only to objects of this class

  //: Find a controller with a given ID.
  ncm_apt_controller* controller_with_id(long id);


  //: Vector of pointers to controllers found.
  vcl_vector<ncm_apt_controller*> controllers_;

  //: Null controller - returned from X(), Y() or Z() if no controller is
  //  assigned.
  ncm_apt_null_controller* null_controller_;

  //: Pointer to controller assigned to X axis
  ncm_apt_controller* X_;

  //: Pointer to controller assigned to Y axis
  ncm_apt_controller* Y_;

  //: Pointer to controller assigned to Z axis
  ncm_apt_controller* Z_;
};

#endif // ncm_apt_server_h_
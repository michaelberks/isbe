#include "ncm_apt_server.h"
#include <vcl_iostream.h>
#include <nailfold/ncm_encrypt.h>

#include <nailfold/ncm_apt_controller.h>
#include <nailfold/ncm_apt_null_controller.h>

// Header files needed for the APT API
#include <windows.h>
#include <APTAPI.h>

//
// Public member functions
//

//
//: Default constructor
ncm_apt_server::ncm_apt_server()
: controllers_(0),
  null_controller_(new ncm_apt_null_controller()),
  X_(NULL),
  Y_(NULL),
  Z_(NULL)
{
  
  //create_controllers(); No longer do this on object creation
	//Wait for apt_server_ to be asked to create controllers
}

//
//: Destructor
ncm_apt_server::~ncm_apt_server()
{
  destroy_controllers();
	delete null_controller_;
  APTCleanUp();
}

//
//: Enable/Disable event dialog box.
void ncm_apt_server::set_event_dialogs(bool dialogs_enabled /* = true */) const
{
  EnableEventDlg(dialogs_enabled);
}

//
//: Return vector of controller IDs available.
vcl_vector<long> ncm_apt_server::ids()
{
  vcl_vector<long> id_list(n_controllers());
  for (unsigned i = 0; i < n_controllers(); ++i)
  {
    id_list[i] = controllers_[i]->id();
  }

  return id_list;
}

//
//: Assign a controller to the X axis
void ncm_apt_server::set_X_id(long id)
{
  // Unassign controller if ID is -1
  if (id == -1)
  {
    X_ = NULL;
    return;
  }

  ncm_apt_controller* apt = controller_with_id(id);

  // Do not allow the same controller to be assigned to two different axes.
  // Otherwise, very different positions could be requested of the same
  // controller almost instantaneously.
  if ((apt != NULL) && (Y() != apt) && (Z() != apt))
    X_ = apt;
  else
    X_ = NULL;
}

//
//: Assign a controller to the Y axis
void ncm_apt_server::set_Y_id(long id)
{
  // Unassign controller if ID is -1
  if (id == -1)
  {
    Y_ = NULL;
    return;
  }

  ncm_apt_controller* apt = controller_with_id(id);

  // Do not allow the same controller to be assigned to two different axes.
  // Otherwise, very different positions could be requested of the same
  // controller almost instantaneously.
  if ((apt != NULL) && (X() != apt) && (Z() != apt))
    Y_ = apt;
  else
    Y_ = NULL;
}

//
//: Assign a controller to the Z axis
void ncm_apt_server::set_Z_id(long id)
{
  // Unassign controller if ID is -1
  if (id == -1)
  {
    Z_ = NULL;
    return;
  }

  ncm_apt_controller* apt = controller_with_id(id);

  // Do not allow the same controller to be assigned to two different axes.
  // Otherwise, very different positions could be requested of the same
  // controller almost instantaneously.
  if ((apt != NULL) && (X() != apt) && (Y() != apt))
    Z_ = apt;
  else
    Z_ = NULL;
}

//
//: Get pointer to controller assigned to X axis
ncm_apt_controller_base* ncm_apt_server::X() const
{
  if (X_ != NULL)
    return X_;
  else
    return null_controller_;
}

//
//: Get pointer to controller assigned to Y axis
ncm_apt_controller_base* ncm_apt_server::Y() const
{
  if (Y_ != NULL)
    return Y_;
  else
    return null_controller_;

}

//
//: Get pointer to controller assigned to Z axis
ncm_apt_controller_base* ncm_apt_server::Z() const
{
  if (Z_ != NULL)
    return Z_;
  else
    return null_controller_;

}

//
//: Move stages to (x,y,z) location
void ncm_apt_server::move_to_xyz(float x, float y, float z)
{
  // TODO
}

//
//: Move stages to zero position (relative to reverse limit switch)
void ncm_apt_server::move_all_zero(bool return_immediately /* = true */)
{
  for (unsigned i = 0; i < n_controllers(); ++i)
    controllers_[i]->move_zero(return_immediately);
}

//
//: Move stages to home position
void ncm_apt_server::move_all_home(bool return_immediately /* = true */)
{
  for (unsigned i = 0; i < n_controllers(); ++i)
    controllers_[i]->move_home(return_immediately);
}

//
// Private member functions
//

//
//: Populate controllers_ with pointers to new ncm_apt_controller instances.
void ncm_apt_server::create_controllers()
{
	vcl_cout << "Initialising APT motors = " << APTInit() << vcl_endl;
	set_event_dialogs(false);

  long n_controllers = 0;
  GetNumHWUnitsEx(HWTYPE_TDC001, &n_controllers);
  
  for (int i = 0; i < n_controllers; ++i)
  {
    long serial_number = -1;
    GetHWSerialNumEx(HWTYPE_TDC001, i, &serial_number);

    ncm_apt_controller* new_controller = 
        new ncm_apt_controller(serial_number);

    controllers_.push_back(new_controller);
  }
}

//
//: Delete pointers to all ncm_apt_controller and clear controllers_.
void ncm_apt_server::destroy_controllers()
{
  while (!controllers_.empty())
  {
    controllers_.pop_back();
  }
}

//
//: Return number of controllers.
unsigned ncm_apt_server::n_controllers() const
{
  return controllers_.size();
}

//
//: Find a controller with a given ID.
ncm_apt_controller* ncm_apt_server::controller_with_id(long id)
{
  for (vcl_vector<ncm_apt_controller*>::iterator it = controllers_.begin();
       it != controllers_.end();
       ++it)
  {
    if ((*it)->id() == id)
      return *it;
  }

  // Return null controller if no matching controller found.
  return NULL;
}
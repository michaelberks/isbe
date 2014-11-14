//:
// \file
// \brief 
//
// \author Phil Tresadern

#include "ncm_qapt_server.h"

#include <nailfold/ncm_apt_controller_base.h>
 
ncm_qapt_server::ncm_qapt_server()
{
  const float reference_position = 0.0f;
  reference_position_[AptX] = reference_position;
  reference_position_[AptY] = reference_position;
  reference_position_[AptZ] = reference_position;

  const float max_displacement = 1e9;
  set_max_displacement(AptX, max_displacement);
  set_max_displacement(AptY, max_displacement);
  set_max_displacement(AptZ, max_displacement);
}
 
//: Move specified axis to the zero position (with respect to the reverse 
//  limit switch).
void ncm_qapt_server::move_zero(apt_axis axis, 
                                bool return_immediately /* = true */)
{
  get_axis(axis)->move_zero(return_immediately);

  if (!return_immediately)
    emit move_complete(axis);
}
 
//: Move specified axis to the given home position.
void ncm_qapt_server::move_home(apt_axis axis, 
                                float home_position,
                                bool return_immediately /* = true */)
{
  get_axis(axis)->move_home(home_position, return_immediately);

  if (!return_immediately)
    emit move_complete(axis);
}
 
//: Move specified axis to the existing home position.
void ncm_qapt_server::move_home(apt_axis axis, 
                                bool return_immediately /* = true */)
{
  get_axis(axis)->move_home(return_immediately);

  if (!return_immediately)
    emit move_complete(axis);
}
 
//: Move specified axis to given position (relative to home position).
void ncm_qapt_server::move_to(ncm_qapt_server::apt_axis axis, 
                              float position,
                              bool return_immediately /* = true */)
{
  int return_code = get_axis(axis)->move_to(position, return_immediately);

  if (!return_immediately)
    emit move_complete(axis);
}
 
//: Move specified axis by a given distance.
void ncm_qapt_server::move_by(ncm_qapt_server::apt_axis axis,
                              float distance,
                              bool return_immediately /* = true */)
{
  int return_code = get_axis(axis)->move_by(distance, return_immediately);

  if (!return_immediately)
    emit move_complete(axis);
}

void ncm_qapt_server::set_reference_position(ncm_qapt_server::apt_axis axis)
{
  if ((First <= axis) && (axis <= Last))
    reference_position_[axis] = get_axis(axis)->absolute_position();
}
 
void ncm_qapt_server::set_max_displacement(ncm_qapt_server::apt_axis axis,
                                           float max_displacement /* = 1e9 */)
{
  if ((First <= axis) && (axis <= Last))
    max_displacement_[axis] = max_displacement;
}

void ncm_qapt_server::create_controllers()
{
	ncm_apt_server::create_controllers();
	QObject::connect(&apt_timer_, SIGNAL(timeout()),
                   this, SLOT(on_apt_timer()) );

  const int clicksPerSecond = 20;
  apt_timer_.start(/* interval = */ 1000.0 / clicksPerSecond);
}

//: Clear controllers_.
void ncm_qapt_server::destroy_controllers()
{
	QObject::disconnect(&apt_timer_, SIGNAL(timeout()),
                   this, SLOT(on_apt_timer()) );
  apt_timer_.stop();
	ncm_apt_server::destroy_controllers();
}

//
// Private slots
//

void ncm_qapt_server::on_apt_timer()
{
  // Read and store the motor positions from the units.
  X()->buffer_position();
  Y()->buffer_position();
  Z()->buffer_position();

  //for (apt_axis axis = First; axis < Last; ++axis)
  apt_axis axis = AptZ;
  {
    // Compute absolute displacement from reference position.
    float displacement = get_axis(axis)->absolute_position() - 
                         reference_position_[axis];

    if (displacement < 0)
      displacement = -displacement;

    if (displacement > max_displacement_[axis])
    {
      vcl_cout << "Displacement = " << displacement << vcl_endl;
      emit move_complete(axis);
    }
  }
}

//
// Private methods
//

//
//: Return the motor controller corresponding to the specified axis.
ncm_apt_controller_base*
ncm_qapt_server::get_axis(ncm_qapt_server::apt_axis axis)
{
  switch (axis)
  {
    case AptX:
      return X();
    case AptY:
      return Y();
    case AptZ:
      return Z();
    default:
      // This will cause problems when dereferenced. 
      // Maybe find a way to return the null_controller instead.
      return NULL;
  }
}


#include "ncm_camera.h"

//: Define static class members here
// ncm_camera::member_ = 0;

//: Define static const class members here
// ncm_camera::const_member_;

//: Define class member functions


ncm_camera::ncm_camera()
{
}

void ncm_camera::set_queue(ncm_video_frame_queue *_queue)
{
	queue_ = _queue;
}

////: Assignment operator
//ncm_camera&
//ncm_camera::operator=(ncm_camera const& rhs)
//{
//  // Check for self-assignment only if necessary (i.e. if self-assignment would
//  // cause the program to behave badly)
//  if (this == &rhs)
//    return *this;
//
//  // Call base class assignment if this class is derived
//  Base::operator=(rhs);
//
//  // Copy the bits needed from rhs
//  //*this = rhs;
//
//  return *this;
//}


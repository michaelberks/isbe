#include <nailfold/qtools/ncm_video_frame.h>

//: Define static class members here
// ncm_lumenera_camera::member_ = 0;

//: Define static const class members here
// ncm_lumenera_camera::const_member_;

//: Define class member functions


//: Default constructor.
ncm_video_frame::ncm_video_frame(QTime frame_time /* = QTime() */)
: frame_(0),
  header_(frame_time)
{
}

ncm_video_frame::ncm_video_frame(unsigned ni, unsigned nj, 
                                 QTime frame_time /* = QTime() */)
: frame_(new vil_image_view<vxl_byte>(ni, nj)),
  header_(frame_time)
{
}

vil_image_view<vxl_byte>* ncm_video_frame::frame()
{
	return &frame_;
}
 
const vil_image_view<vxl_byte>* ncm_video_frame::frame() const
{
	return &frame_;
}
//Set/Get metadata fo the frame
ncm_video_frame_header ncm_video_frame::get_frame_header() const
{
	return header_;
}

//: Set the filename for the frame.
void ncm_video_frame::set_frame_name( vcl_string frame_name )
{
	header_.frame_name_ = frame_name;
}
vcl_string ncm_video_frame::frame_name() const
{
	return header_.frame_name_;
}

//: Set/get the frame number
void ncm_video_frame::set_frame_num( int frame_num )
{
	header_.frame_num_ = frame_num;
}
	
int ncm_video_frame::frame_num() const
{
	return header_.frame_num_;
}


//: Set the timestamp for the frame.
void ncm_video_frame::set_frame_time( QTime frame_time )
{
	header_.frame_time_ = frame_time;
}
 
//: Return the frame's timestamp.
QTime ncm_video_frame::frame_time() const
{
	return header_.frame_time_;
}
 
//: Set the motor positions.
void ncm_video_frame::set_motor_xyz(float x, float y, float z)
{
  header_.motor_position_[0] = x;
  header_.motor_position_[1] = y;
  header_.motor_position_[2] = z;
}
 
//: Put motor positions into x, y and z.
void ncm_video_frame::get_motor_xyz(float& x, float& y, float& z) const
{
  x = header_.motor_position_[0];
  y = header_.motor_position_[1];
  z = header_.motor_position_[2];
}
 
//: Set the sharpness of the image.
void ncm_video_frame::set_sharpness(double sharpness)
{
  header_.sharpness_ = sharpness;
}
 
//: Return sharpness of the image.
double ncm_video_frame::sharpness() const
{
  return header_.sharpness_;
}
 
//: Indicate that this is the last frame in the sequence.
void ncm_video_frame::set_last_frame(bool is_last_frame /* = true */)
{
  header_.is_last_frame_ = is_last_frame;
}
 
//: Return true if this is the last frame in the sequence.
bool ncm_video_frame::is_last_frame() const
{
  return header_.is_last_frame_;
}

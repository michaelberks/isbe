#ifndef NCM_VIDEO_FRAME_H
#define NCM_VIDEO_FRAME_H

#include <vcl_vector.h>
#include <vil/vil_image_view.h>
#include <QTime>
#include <QMetaType>

class ncm_video_frame_header
	{
	public:
		ncm_video_frame_header(QTime frame_time = QTime())
			: frame_time_(frame_time),
				frame_name_(""),
				frame_num_(0),
			  motor_position_(3,0),
				sharpness_(-1.0),
				is_last_frame_(false)
		{
		}
		
		ncm_video_frame_header(const ncm_video_frame_header &other)
		{
			frame_num_ = other.frame_num_;
			frame_name_ = other.frame_name_;
			frame_time_ = other.frame_time_;
			motor_position_ = other.motor_position_;
			sharpness_ = other.sharpness_;
			is_last_frame_ = other.is_last_frame_;
		}

		~ncm_video_frame_header()
		{
		}

		//Number of frame in image sequence
		int frame_num_;

		vcl_string frame_name_; //Empty unless tagged for saving

		//: Timestamp for the frame.
		QTime frame_time_;

		//: Motor (x,y,z) positions at the point of capture.
		vcl_vector<float> motor_position_;

		//: Sharpness of the image.
		double sharpness_;

		//: Flag indicating this is the last frame to be saved.
		bool is_last_frame_;

};

Q_DECLARE_METATYPE(ncm_video_frame_header);

class ncm_video_frame
{
// INTERFACE

public:

  //: Constructor (default).
	ncm_video_frame(QTime frame_time = QTime());

  //: Constructor.
	ncm_video_frame(unsigned ni, unsigned nj, QTime frame_time = QTime());

  //: Return a pointer to the frame.
	vil_image_view<vxl_byte>* frame();
	const vil_image_view<vxl_byte>* frame() const;

  // Set/Get metadata
	ncm_video_frame_header get_frame_header() const;

	//: Set the filename for the frame.
	void set_frame_name( vcl_string frame_name );
	vcl_string frame_name() const;

	//: Set the filename for the frame.
	void set_frame_num( int frame_num );
	int frame_num() const;

  //: Set the timestamp for the frame.
	void set_frame_time( QTime frame_time );

  //: Return the frame's timestamp.
	QTime frame_time() const;

  //: Set the motor positions.
  void set_motor_xyz(float x, float y, float z);

  //: Put motor positions into x, y and z.
  void get_motor_xyz(float& x, float& y, float& z) const;

  //: Set the sharpness of the image.
  void set_sharpness(double sharpness);

  //: Return sharpness of the image.
  double sharpness() const;

  //: Indicate that this is the last frame in the sequence.
  void set_last_frame(bool is_last_frame = true);

  //: Return true if this is the last frame in the sequence.
  bool is_last_frame() const;


// IMPLEMENTATION

protected:

  // Members and functions visible to objects of this class and derived classes

private: // Variables

  //: The frame itself.
	vil_image_view<vxl_byte> frame_;

  //  Metadata
	ncm_video_frame_header header_;

};

#endif // NCM_VIDEO_FRAME_H
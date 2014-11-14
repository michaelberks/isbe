#include "ncm_qdata_stream.h"

 
//  Public Methods
 
//: Default constructor.
ncm_qdata_stream::ncm_qdata_stream()
{
  camera_.set_queue(&main_queue_);
  
  processor_.set_input_queue(&main_queue_);
  processor_.set_output_queue(&save_queue_);

  saver_.set_queue(&save_queue_);

  // FIXME
  //QObject::connect( &camera_, SIGNAL(frame_ready()),
  //                  &processor_, SLOT(process_frame()) );
}
 
//: Return pointer to camera.
ncm_dmk_camera& ncm_qdata_stream::camera()
{
  return camera_;
}
const ncm_dmk_camera& ncm_qdata_stream::camera() const
{
  return camera_;
}
 
//: Return pointer to the main queue (frame buffer).
ncm_video_frame_queue& ncm_qdata_stream::main_queue()
{
	return main_queue_;
}
const ncm_video_frame_queue& ncm_qdata_stream::main_queue() const
{
	return main_queue_;
}
 
//: Return pointer to processor.
ncm_qdata_processor& ncm_qdata_stream::processor()
{
  return processor_;
}
const ncm_qdata_processor& ncm_qdata_stream::processor() const
{
  return processor_;
}
 
//: Return pointer to the save queue (frame buffer).
ncm_video_frame_queue& ncm_qdata_stream::save_queue()
{
	return save_queue_;
}
const ncm_video_frame_queue& ncm_qdata_stream::save_queue() const
{
	return save_queue_;
}
 
//: Return pointer to saver.
ncm_qdata_saver& ncm_qdata_stream::saver()
{
  return saver_;
}
const ncm_qdata_saver& ncm_qdata_stream::saver() const
{
  return saver_;
}
 
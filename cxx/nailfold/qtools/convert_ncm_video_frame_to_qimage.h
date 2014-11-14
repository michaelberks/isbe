#include <QImage>
#include <QSharedPointer>

#include <vcl_cstring.h> // for vcl_memcpy

#include <nailfold/qtools/ncm_video_frame.h>

// Declarations

void convert_ncm_video_frame_to_qimage(
  QImage &qimage, 
  QSharedPointer<ncm_video_frame> vid_frame);


// Definitions

inline 
void convert_ncm_video_frame_to_qimage(
  QImage &qimage, 
  QSharedPointer<ncm_video_frame> vid_frame)
{
	int src_cols = vid_frame->frame()->ni();
	int src_rows = vid_frame->frame()->nj();

	if ( qimage.isNull() || 
      (qimage.width() != src_cols) || 
      (qimage.height() != src_rows) ) 
	{
		// Initialise qimage
		qimage = QImage( src_cols, src_rows, QImage::Format_Indexed8 );
		
		// Set colour values
		const unsigned n_colours = 256;
		qimage.setNumColors( n_colours );

		for (unsigned i = 0; i < n_colours; ++i)
		{
			qimage.setColor(i, qRgb(i, i, i));
		}

		vcl_cout << "Initialised qimage: " 
             << qimage.width() << " x " << qimage.height() << vcl_endl;
	}

  // Copy pixel data from frame to qimage
  vcl_memcpy( qimage.bits(), 
              vid_frame->frame()->top_left_ptr(), 
              src_cols*src_rows );
}
